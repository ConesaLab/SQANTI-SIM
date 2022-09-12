#!/usr/bin/env python3
"""
classify_gtf.py

Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account itself in the reference.
The classification step is adapted from the original SQANTI3 QC v5.0.
(sqanti3_qc.py -> https://github.com/ConesaLab/SQANTI3)

Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)
"""

import os
import sys
import distutils.spawn
import bisect
import subprocess
import itertools
import multiprocessing as mp
from collections import defaultdict, namedtuple
from tqdm import tqdm

try:
    from bx.intervals import Interval, IntervalTree
except ImportError:
    print("Unable to import bx-python! Please make sure bx-python is installed.", file=sys.stderr)
    sys.exit(-1)

try:
    from cupcake.tofu.compare_junctions import compare_junctions
except ImportError:
    print("Unable to import cupcake.tofu! Please make sure you install cupcake.", file=sys.stderr)
    sys.exit(-1)


utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "SQANTI3/utilities")
GTF2GENEPRED_PROG = os.path.join(utilitiesPath, "gtfToGenePred")

if distutils.spawn.find_executable(GTF2GENEPRED_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GTF2GENEPRED_PROG), file=sys.stderr)
    sys.exit(-1)


#####################################
#                                   #
#          DEFINE CLASSES           #
#                                   #
#####################################

class genePredReader(object):
    """Gets gtfToGenePred output and builds genePredRecord objects"""

    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)


class genePredRecord(object):
    """Saves the features of each transcript read by gtfTogenePred"""

    def __init__(self, id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene=None):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart  # 1-based start
        self.txEnd = txEnd  # 1-based end
        self.cdsStart = cdsStart  # 1-based start
        self.cdsEnd = cdsEnd  # 1-based end
        self.exonCount = exonCount
        self.exonStarts = exonStarts  # 0-based starts
        self.exonEnds = exonEnds  # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s, e in zip(exonStarts, exonEnds):
            self.length += e - s
            self.exons.append(Interval(s, e))

        # junctions are stored (1-based last base of prev exon, 1-based first base of next exon)
        self.junctions = [
            (self.exonEnds[i], self.exonStarts[i + 1])
            for i in range(self.exonCount - 1)
        ]

    @property
    def segments(self):
        return self.exons

    @classmethod
    def from_line(cls, line):
        raw = line.strip().split("\t")
        return cls(
            id=raw[0],
            chrom=raw[1],
            strand=raw[2],
            txStart=int(raw[3]),
            txEnd=int(raw[4]),
            cdsStart=int(raw[5]),
            cdsEnd=int(raw[6]),
            exonCount=int(raw[7]),
            exonStarts=[
                int(x) for x in raw[8][:-1].split(",")
            ],  # exonStarts string has extra , at end
            exonEnds=[
                int(x) for x in raw[9][:-1].split(",")
            ],  # exonEnds string has extra , at end
            gene=raw[11] if len(raw) >= 12 else None,
        )


class myQueryTranscripts:
    """Features of the query transcript and its associated reference"""

    def __init__(self, id, gene_id, tss_diff, tts_diff, num_exons, length,
                 str_class, subtype=None, genes=None, transcripts=None,
                 chrom=None, strand=None, refLen="NA",refExons="NA",
                 refStart="NA", refEnd="NA", q_splicesite_hit=0,
                 q_exon_overlap=0,junctions=None, tss=None, tts=None,
                 intergenic_assoc=None):

        self.id = id
        self.gene_id = gene_id  # gene where this transcript cames from
        self.tss_diff = tss_diff  # distance to TSS of best matching ref
        self.tts_diff = tts_diff  # distance to TTS of best matching ref
        self.tss_gene_diff = "NA"  # min distance to TSS of genes matching ref
        self.tts_gene_diff = "NA"  # min distance to TTS of genes matching ref
        self.genes = genes if genes is not None else []
        self.AS_genes = set()  # ref genes that are hit on the opposite strand
        self.transcripts = transcripts if transcripts is not None else []
        self.num_exons = num_exons
        self.length = length
        self.str_class = str_class  # structural classification of the isoform
        self.chrom = chrom
        self.strand = strand
        self.subtype = subtype
        self.refLen = refLen
        self.refExons = refExons
        self.refStart = refStart
        self.refEnd = refEnd
        self.q_splicesite_hit = q_splicesite_hit
        self.q_exon_overlap = q_exon_overlap
        self.junctions = junctions if len(junctions) > 0 else ["NA", "NA"]
        self.tss = tss if self.strand == "+" else tts
        self.tts = tts if self.strand == "+" else tss
        self.intergenic_assoc = intergenic_assoc

    def get_total_diff(self):
        return abs(self.tss_diff) + abs(self.tts_diff)

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons):
        self.transcripts = [ref_transcript]
        self.genes = [ref_gene]
        self.tss_diff = tss_diff
        self.tts_diff = tts_diff
        self.refLen = refLen
        self.refExons = refExons


#####################################
#                                   #
#         DEFINE FUNCTIONS          #
#                                   #
#####################################

def gtf_parser(gtf_name: str) -> defaultdict:
    """Parse input isoforms from GTF to dict grouped by chromosome regions

    Args:
        gtf_name (str) the GTF file location to be parsed

    Returns:
        defaultdict: a dict of genePredRecord sorted by chr regions
    """

    global queryFile
    queryFile = os.path.splitext(gtf_name)[0] + ".genePred"

    # run gtf to genePred
    cmd = (
        GTF2GENEPRED_PROG
        + " {0} {1} -genePredExt -allErrors -ignoreGroupsWithoutExons".format(
            gtf_name, queryFile
        )
    )
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("ERROR running cmd: {0}".format(cmd), file=sys.stderr)
        sys.exit(-1)

    isoforms_list = defaultdict(lambda: [])

    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)

    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    isoform_list_by_reg = defaultdict(lambda: {})
    for chrom in isoforms_list:
        c = 0
        for t in isoforms_list[chrom]:
            for k, r in isoform_list_by_reg[chrom].items():
                if t.txStart <= r[1] and r[0] <= t.txEnd:
                    isoform_list_by_reg[chrom][k][0] = min(r[0], t.txStart)
                    isoform_list_by_reg[chrom][k][1] = max(r[1], t.txEnd)
                    isoform_list_by_reg[chrom][k][2].append(t)
                    break
            else:
                isoform_list_by_reg[chrom][chrom + str(c)] = [t.txStart, t.txEnd, [t]]
                c += 1

    isoforms_list = defaultdict(lambda: [])
    for chrom in isoform_list_by_reg:
        for k, reg in isoform_list_by_reg[chrom].items():
            isoforms_list[chrom].append(reg[2])

    return isoforms_list


def transcript_classification(trans_by_region: list) -> dict:
    """Classify transcripts in SQANTI3 structural categories

    Given a list of transcripts from the same chromosomic region classify each
    transcript into the most suitable SQANTI3 structural category comparing it
    with the rest transcripts of the region

    Args:
        trans_by_region (list): transcripts from the same region of the chr

    Returns:
        defaultdict: transcripts classified with its associated reference
    """

    res = defaultdict(lambda: [])
    for trans in trans_by_region:
        # dict of junctions by region in the chr
        junctions_by_chr = defaultdict(
            lambda: {"donors": set(), "acceptors": set(), "da_pairs": set()}
        )
        # dict of junctions by gene
        junctions_by_gene = defaultdict(lambda: set())
        # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
        known_5_3_by_gene = defaultdict(lambda: {"begin": set(), "end": set()})

        for r in trans_by_region:
            if trans.id == r.id or r.length < min_ref_len:
                continue
            known_5_3_by_gene[r.gene]["begin"].add(r.txStart)
            known_5_3_by_gene[r.gene]["end"].add(r.txEnd)
            if r.exonCount >= 2:
                for d, a in r.junctions:
                    junctions_by_chr[r.chrom]["donors"].add(d)
                    junctions_by_chr[r.chrom]["acceptors"].add(a)
                    junctions_by_chr[r.chrom]["da_pairs"].add((d, a))
                    junctions_by_gene[r.gene].add((d, a))

        for k in junctions_by_chr:
            junctions_by_chr[k]["donors"] = list(junctions_by_chr[k]["donors"])
            junctions_by_chr[k]["donors"].sort()
            junctions_by_chr[k]["acceptors"] = list(
                junctions_by_chr[k]["acceptors"]
            )
            junctions_by_chr[k]["acceptors"].sort()
            junctions_by_chr[k]["da_pairs"] = list(
                junctions_by_chr[k]["da_pairs"]
            )
            junctions_by_chr[k]["da_pairs"].sort()

        start_ends_by_gene = dict(known_5_3_by_gene)
        junctions_by_chr = dict(junctions_by_chr)
        junctions_by_gene = dict(junctions_by_gene)

        # Find best reference hit
        isoform_hit = transcriptsKnownSpliceSites(
            trans, trans_by_region, start_ends_by_gene, min_ref_len
        )

        if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
            # not FSM or ISM --> see if it is NIC, NNC, or fusion
            isoform_hit = novelIsoformsKnownGenes(
                isoform_hit, trans, junctions_by_chr, junctions_by_gene
            )
        elif isoform_hit.str_class in ("", "geneOverlap"):
            # possibly NNC, genic, genic intron, anti-sense, or intergenic
            isoform_hit = associationOverlapping(
                isoform_hit, trans, junctions_by_chr, trans_by_region
            )

        # Save trans classification
        res[isoform_hit.chrom].append(isoform_hit)

    return dict(res)


def transcriptsKnownSpliceSites(trec: genePredRecord, ref_chr: list, start_ends_by_gene: dict, min_ref_len: int) -> myQueryTranscripts:
    """Find best reference hit for the query transcript

    Checks for full-splice-match, incomplete-splice-match, anyKnownJunction,
    anyKnownSpliceSite or geneOverlap.

    Args:
        trec (genePredRecord): query transcript to be classified
        ref_chr (list): list of reference transcript from the same region
        start_ends_by_gene (dict): begins and ends from genes
        min_ref_len (int): minimum length of a trancript to consider it as ref

    Returns:
        myQueryTranscripts: best reference hit(s) for the query transcript
    """

    def calc_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
        """Calculates the number of nucleotides that overlap"""
        if s1 == "NA" or s2 == "NA":
            return 0
        if s1 > s2:
            s1, e1, s2, e2 = s2, e2, s1, e1
        return max(0, min(e1, e2) - max(s1, s2))

    def calc_splicesite_agreement(query_exons: list, ref_exons: list) -> int:
        """Returns the number of matching splicesites"""
        q_sites = {}
        for e in query_exons:
            q_sites[e.start] = 0
            q_sites[e.end] = 0
        for e in ref_exons:
            if e.start in q_sites:
                q_sites[e.start] = 1
            if e.end in q_sites:
                q_sites[e.end] = 1
        return sum(q_sites.values())

    def calc_exon_overlap(query_exons: list, ref_exons: list) -> int:
        """The number of nucleotides in exons that overlap"""
        q_bases = {}
        for e in query_exons:
            for b in range(e.start, e.end):
                q_bases[b] = 0

        for e in ref_exons:
            for b in range(e.start, e.end):
                if b in q_bases:
                    q_bases[b] = 1
        return sum(q_bases.values())

    def get_diff_tss_tts(trec: genePredRecord, ref: genePredRecord) -> tuple:
        """Calculates the absolute difference in the TSS and TTS of 2 trans"""
        if trec.strand == "+":
            diff_tss = trec.txStart - ref.txStart
            diff_tts = ref.txEnd - trec.txEnd
        else:
            diff_tts = trec.txStart - ref.txStart
            diff_tss = ref.txEnd - trec.txEnd
        return diff_tss, diff_tts

    def get_gene_diff_tss_tts(isoform_hit: myQueryTranscripts):
        """Computes the difference to gene start and end

        Now that we know the reference (isoform) it hits,
        add the nearest start/end site for that gene (all isoforms of the gene)
        """

        nearest_start_diff, nearest_end_diff = float("inf"), float("inf")
        for ref_gene in isoform_hit.genes:
            for x in start_ends_by_gene[ref_gene]["begin"]:
                d = trec.txStart - x
                if abs(d) < abs(nearest_start_diff):
                    nearest_start_diff = d
            for x in start_ends_by_gene[ref_gene]["end"]:
                d = trec.txEnd - x
                if abs(d) < abs(nearest_end_diff):
                    nearest_end_diff = d

        if trec.strand == "+":
            isoform_hit.tss_gene_diff = (
                nearest_start_diff
                if nearest_start_diff != float("inf")
                else "NA"
            )
            isoform_hit.tts_gene_diff = (
                nearest_end_diff if nearest_end_diff != float("inf") else "NA"
            )
        else:
            isoform_hit.tss_gene_diff = (
                -nearest_end_diff
                if nearest_start_diff != float("inf")
                else "NA"
            )
            isoform_hit.tts_gene_diff = (
                -nearest_start_diff
                if nearest_end_diff != float("inf")
                else "NA"
            )

    def categorize_incomplete_matches(trec: genePredRecord, ref: genePredRecord)-> str:
        """Returns the subcategory of incomplete splice matches

        intron_retention --- at least one trec exon covers at least two adjacent ref exons
        complete --- all junctions agree and is not IR
        5prime_fragment --- all junctions agree but trec has less 5' exons. The isoform is a 5' fragment of the reference transcript
        3prime_fragment --- all junctions agree but trec has less 3' exons. The isoform is a 3' fragment of the reference transcript
        internal_fragment --- all junctions agree but trec has less 5' and 3' exons
        """
        # check intron retention
        ref_exon_tree = IntervalTree()
        for i, e in enumerate(ref.exons):
            ref_exon_tree.insert(e.start, e.end, i)
        for e in trec.exons:
            if (
                len(ref_exon_tree.find(e.start, e.end)) > 1
            ):  # multiple ref exons covered
                return "intron_retention"

        agree_front = trec.junctions[0] == ref.junctions[0]
        agree_end = trec.junctions[-1] == ref.junctions[-1]
        if agree_front:
            if agree_end:
                return "complete"
            else:  # front agrees, end does not
                return (
                    "5prime_fragment"
                    if trec.strand == "+"
                    else "3prime_fragment"
                )
        else:
            if agree_end:  # front does not agree, end agrees
                return (
                    "3prime_fragment"
                    if trec.strand == "+"
                    else "5prime_fragment"
                )
            else:
                return "internal_fragment"

    isoform_hit = myQueryTranscripts(
        id=trec.id,
        gene_id=trec.gene,
        tts_diff="NA",
        tss_diff="NA",
        num_exons=trec.exonCount,
        length=trec.length,
        str_class="",
        chrom=trec.chrom,
        strand=trec.strand,
        subtype="no_subcategory",
        junctions=trec.junctions,
        tss=trec.txStart,
        tts=trec.txEnd,
    )

    cat_ranking = {
        "full-splice_match": 5,
        "incomplete-splice_match": 4,
        "anyKnownJunction": 3,
        "anyKnownSpliceSite": 2,
        "geneOverlap": 1,
        "": 0,
    }

    # ----------------------------------#
    #       SPLICED TRANSCRIPT         #
    # ----------------------------------#
    if trec.exonCount >= 2:
        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        for ref in ref_chr:
            if trec.id == ref.id or ref.length < min_ref_len:
                continue  # cannot use itself as reference

            if trec.txStart < ref.txEnd and ref.txStart < trec.txEnd: # It is < and not <= because start is 0-based and end 1-based
                hits_by_gene[ref.gene].append(ref)

        if len(hits_by_gene) == 0:
            return isoform_hit

        for ref_gene in sorted(hits_by_gene):
            isoform_hit = myQueryTranscripts(
                id=trec.id,
                gene_id=trec.gene,
                tts_diff="NA",
                tss_diff="NA",
                num_exons=trec.exonCount,
                length=trec.length,
                str_class="",
                chrom=trec.chrom,
                strand=trec.strand,
                subtype="no_subcategory",
                junctions=trec.junctions,
                tss=trec.txStart,
                tts=trec.txEnd,
            )

            for ref in hits_by_gene[ref_gene]:
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue

                # --MONO-EXONIC REFERENCE--#
                if ref.exonCount == 1:
                    if (
                        ref.exonCount == 1
                    ):  # mono-exonic reference, handle specially here
                        if (
                            calc_exon_overlap(trec.exons, ref.exons) > 0
                            and cat_ranking[isoform_hit.str_class]
                            < cat_ranking["geneOverlap"]
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                trec.gene,
                                "NA",
                                "NA",
                                trec.exonCount,
                                trec.length,
                                "geneOverlap",
                                subtype="mono-exon",
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=[ref.id],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=0,
                                q_exon_overlap=calc_exon_overlap(
                                    trec.exons, ref.exons
                                ),
                                junctions=trec.junctions,
                                tss=trec.txStart,
                                tts=trec.txEnd,
                            )
                # --MULTI-EXONIC REFERENCE--#
                else:
                    match_type = compare_junctions(
                        trec,
                        ref,
                        internal_fuzzy_max_dist=0,
                        max_5_diff=999999,
                        max_3_diff=999999,
                    )

                    if match_type not in (
                        "exact",
                        "subset",
                        "partial",
                        "concordant",
                        "super",
                        "nomatch",
                    ):
                        raise Exception(
                            "Unknown match category {0}!".format(match_type)
                        )

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                    # #############################
                    # SQANTI's full-splice_match
                    # #############################
                    if match_type == "exact":
                        subtype = "multi-exon"
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["full-splice_match"]
                            or abs(diff_tss) + abs(diff_tts)
                            < isoform_hit.get_total_diff()
                        ):
                            # subcategory for matching 5' and matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) <= 50:
                                subtype = "reference_match"
                            # subcategory for matching 5' and non-matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) > 50:
                                subtype = "alternative_3end"
                            # subcategory for matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) <= 50:
                                subtype = "alternative_5end"
                            # subcategory for non-matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) > 50:
                                subtype = "alternative_3end5end"
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                trec.gene,
                                diff_tss,
                                diff_tts,
                                trec.exonCount,
                                trec.length,
                                str_class="full-splice_match",
                                subtype=subtype,
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=[ref.id],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(
                                    trec.exons, ref.exons
                                ),
                                junctions=trec.junctions,
                                tss=trec.txStart,
                                tts=trec.txEnd,
                            )

                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subtype = categorize_incomplete_matches(trec, ref)
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking[
                            "incomplete-splice_match"
                        ] or (
                            isoform_hit.str_class == "incomplete-splice_match"
                            and abs(diff_tss) + abs(diff_tts)
                            < isoform_hit.get_total_diff()
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                trec.gene,
                                diff_tss,
                                diff_tts,
                                trec.exonCount,
                                trec.length,
                                str_class="incomplete-splice_match",
                                subtype=subtype,
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=[ref.id],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(
                                    trec.exons, ref.exons
                                ),
                                junctions=trec.junctions,
                                tss=trec.txStart,
                                tts=trec.txEnd,
                            )

                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ("partial", "concordant", "super"):
                        q_sp_hit = calc_splicesite_agreement(
                            trec.exons, ref.exons
                        )
                        q_ex_overlap = calc_exon_overlap(trec.exons, ref.exons)
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["anyKnownJunction"]
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit > isoform_hit.q_splicesite_hit
                            )
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit == isoform_hit.q_splicesite_hit
                                and q_ex_overlap > isoform_hit.q_exon_overlap
                            )
                            or (
                                isoform_hit.str_class == "anyKnownJunction"
                                and q_sp_hit == isoform_hit.q_splicesite_hit
                                and q_exon_d
                                < abs(trec.exonCount - isoform_hit.refExons)
                            )
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                trec.gene,
                                "NA",
                                "NA",
                                trec.exonCount,
                                trec.length,
                                str_class="anyKnownJunction",
                                subtype="no_subcategory",
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=["novel"],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(
                                    trec.exons, ref.exons
                                ),
                                junctions=trec.junctions,
                                tss=trec.txStart,
                                tts=trec.txEnd,
                            )

                    else:  # must be nomatch
                        assert match_type == "nomatch"
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if (
                            cat_ranking[isoform_hit.str_class]
                            < cat_ranking["anyKnownSpliceSite"]
                            and calc_splicesite_agreement(
                                trec.exons, ref.exons
                            )
                            > 0
                        ):
                            isoform_hit = myQueryTranscripts(
                                trec.id,
                                trec.gene,
                                "NA",
                                "NA",
                                trec.exonCount,
                                trec.length,
                                str_class="anyKnownSpliceSite",
                                subtype="no_subcategory",
                                chrom=trec.chrom,
                                strand=trec.strand,
                                genes=[ref.gene],
                                transcripts=["novel"],
                                refLen=ref.length,
                                refExons=ref.exonCount,
                                refStart=ref.txStart,
                                refEnd=ref.txEnd,
                                q_splicesite_hit=calc_splicesite_agreement(
                                    trec.exons, ref.exons
                                ),
                                q_exon_overlap=calc_exon_overlap(
                                    trec.exons, ref.exons
                                ),
                                junctions=trec.junctions,
                                tss=trec.txStart,
                                tts=trec.txEnd,
                            )

                        if isoform_hit.str_class == "":  # still not hit yet, check exonic overlap
                            if (
                                cat_ranking[isoform_hit.str_class]
                                < cat_ranking["geneOverlap"]
                                and calc_exon_overlap(trec.exons, ref.exons)
                                > 0
                            ):
                                isoform_hit = myQueryTranscripts(
                                    trec.id,
                                    trec.gene,
                                    "NA",
                                    "NA",
                                    trec.exonCount,
                                    trec.length,
                                    str_class="geneOverlap",
                                    subtype="no_subcategory",
                                    chrom=trec.chrom,
                                    strand=trec.strand,
                                    genes=[ref.gene],
                                    transcripts=["novel"],
                                    refLen=ref.length,
                                    refExons=ref.exonCount,
                                    refStart=ref.txStart,
                                    refEnd=ref.txEnd,
                                    q_splicesite_hit=calc_splicesite_agreement(
                                        trec.exons, ref.exons
                                    ),
                                    q_exon_overlap=calc_exon_overlap(
                                        trec.exons, ref.exons
                                    ),
                                    junctions=trec.junctions,
                                    tss=trec.txStart,
                                    tts=trec.txEnd,
                                )
            best_by_gene[ref_gene] = isoform_hit

        # now we have best_by_gene:
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
        geneHitTuple = namedtuple(
            "geneHitTuple", ["score", "rStart", "rEnd", "rGene", "iso_hit"]
        )
        best_by_gene = [
            geneHitTuple(
                cat_ranking[iso_hit.str_class],
                iso_hit.refStart,
                iso_hit.refEnd,
                ref_gene,
                iso_hit,
            )
            for ref_gene, iso_hit in best_by_gene.items()
        ]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene))
        if len(best_by_gene) == 0:  # no hit
            try:
                isoform_hit.intergenic_assoc = ref_gene
            except:
                isoform_hit.intergenic_assoc = "novel"
            return isoform_hit

        best_by_gene.sort(
            key=lambda x: (
                x.score,
                x.iso_hit.q_splicesite_hit
                + (x.iso_hit.q_exon_overlap)
                * 1.0
                / sum(e.end - e.start for e in trec.exons)
                + calc_overlap(x.rStart, x.rEnd, trec.txStart, trec.txEnd)
                * 1.0
                / (x.rEnd - x.rStart)
                - abs(trec.exonCount - x.iso_hit.refExons),
            ),
            reverse=True,
        )  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end = best_by_gene[0].rStart, best_by_gene[0].rEnd
        for t in best_by_gene[1:]:
            if t.score == 0:
                break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0:
                isoform_hit.genes.append(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(
                    cur_end, t.rEnd
                )

    # ----------------------------------#
    #       UNSPLICED TRANSCRIPT       #
    # ----------------------------------#
    else:
        for ref in ref_chr:
            if trec.id == ref.id or ref.length < min_ref_len:
                continue  # cannot use itself as reference
            # if hits_exon(trec, ref) and ref.exonCount == 1:
            if (
                trec.txStart < ref.txEnd
                and ref.txStart < trec.txEnd
                and ref.exonCount == 1
            ):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue

                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                if isoform_hit.str_class == "":  # no match so far
                    isoform_hit = myQueryTranscripts(
                        trec.id,
                        trec.gene,
                        diff_tss,
                        diff_tts,
                        trec.exonCount,
                        trec.length,
                        "full-splice_match",
                        subtype="mono-exon",
                        chrom=trec.chrom,
                        strand=trec.strand,
                        genes=[ref.gene],
                        transcripts=[ref.id],
                        refLen=ref.length,
                        refExons=ref.exonCount,
                        junctions=trec.junctions,
                        tss=trec.txStart,
                        tts=trec.txEnd,
                    )
                elif abs(diff_tss) + abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(
                        ref.id,
                        ref.gene,
                        diff_tss,
                        diff_tts,
                        ref.length,
                        ref.exonCount,
                    )

        if isoform_hit.str_class == "":
            for ref in ref_chr:
                if trec.id == ref.id or ref.length < min_ref_len:
                    continue  # to not match with itself
                # if hits_exon(trec, ref) and ref.exonCount >= 2:
                if (
                    trec.txStart < ref.txEnd
                    and ref.txStart < trec.txEnd
                    and ref.exonCount >= 2
                ):
                    if (
                        calc_exon_overlap(trec.exons, ref.exons) == 0
                    ):  # no exonic overlap, skip!
                        continue

                    if ref.strand != trec.strand:
                        # opposite strand, just record it in AS_genes
                        isoform_hit.AS_genes.add(ref.gene)
                        continue

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                    for e in ref.exons:
                        if e.start <= trec.txStart < trec.txEnd <= e.end:
                            isoform_hit.str_class = "incomplete-splice_match"
                            isoform_hit.subtype = "mono-exon"
                            isoform_hit.modify(
                                ref.id,
                                ref.gene,
                                diff_tss,
                                diff_tts,
                                ref.length,
                                ref.exonCount,
                            )
                            # this is as good a match as it gets, we can stop the search here
                            get_gene_diff_tss_tts(isoform_hit)
                            return isoform_hit

                    # if we haven't exited here, then ISM hit is not found yet
                    # instead check if it's NIC by intron retention
                    # but we don't exit here since the next gene could be a ISM hit
                    if ref.txStart <= trec.txStart < trec.txEnd <= ref.txEnd:
                        isoform_hit.str_class = "novel_in_catalog"
                        isoform_hit.subtype = "mono-exon"
                        # check for intron retention
                        if len(ref.junctions) > 0:
                            for (d, a) in ref.junctions:
                                if trec.txStart < d < a < trec.txEnd:
                                    isoform_hit.subtype = (
                                        "mono-exon_by_intron_retention"
                                    )
                                    break
                        isoform_hit.modify(
                            "novel",
                            ref.gene,
                            "NA",
                            "NA",
                            ref.length,
                            ref.exonCount,
                        )
                        get_gene_diff_tss_tts(isoform_hit)
                        return isoform_hit

                    # if we get to here, means neither ISM nor NIC, so just add a ref gene and categorize further later
                    isoform_hit.genes.append(ref.gene)

    get_gene_diff_tss_tts(isoform_hit)
    isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]["begin"])
    return isoform_hit


def novelIsoformsKnownGenes(
    isoforms_hit: myQueryTranscripts,
    trec: genePredRecord,
    junctions_by_chr: dict,
    junctions_by_gene: dict,
) -> myQueryTranscripts:
    """Check for NIC, NNC or fusion
    At this point definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    Args:
        isoforms_hit (myQueryTranscript): best isoform hit at the moment
        trec (genePredRecord): query transcript
        junctions_by_chr (dict): dictionary with all the junctions in the region
        junction_by_gene (dict): dictionary with all the junctions by gene of the current region

    Returns:
        isoforms_hit (myQueryTranscript): updated isoforms hit (myQueryTranscripts object)
    """

    def has_intron_retention():
        for e in trec.exons:
            m = bisect.bisect_left(
                junctions_by_chr[trec.chrom]["da_pairs"], (e.start, e.end)
            )
            if (
                m < len(junctions_by_chr[trec.chrom]["da_pairs"])
                and e.start
                <= junctions_by_chr[trec.chrom]["da_pairs"][m][0]
                < junctions_by_chr[trec.chrom]["da_pairs"][m][1]
                < e.end
            ):
                return True
        return False

    ref_genes = list(set(isoforms_hit.genes))

    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to 'NA' for non-FSM/ISM matches)
    isoforms_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d, a in trec.junctions:
            all_junctions_known = (
                all_junctions_known
                and (d in junctions_by_chr[trec.chrom]["donors"])
                and (a in junctions_by_chr[trec.chrom]["acceptors"])
            )
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and (
                (d, a) in ref_gene_junctions
            )
        if all_junctions_known:
            isoforms_hit.str_class = "novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "combination_of_known_splicesites"
                # For NIC the splice site can come from different ref genes
                # Add those to ref
                st_other_gene = []
                other_ref_genes = set()
                gene_st = [st for sj in junctions_by_gene[ref_genes[0]] for st in sj]
                for d, a in trec.junctions:
                    if d not in gene_st:
                        st_other_gene.append(d)
                    if a not in gene_st:
                        st_other_gene.append(a)
                
                for g in junctions_by_gene:
                    gene_st = [st for sj in junctions_by_gene[g] for st in sj]
                    for st in gene_st:
                        if st in st_other_gene:
                            other_ref_genes.add(g)
                            st_other_gene.remove(st)
                    if len(st_other_gene) == 0:
                        break
                
                for g in other_ref_genes:
                    isoforms_hit.genes.append(g)
          
        else:
            isoforms_hit.str_class = "novel_not_in_catalog"
            isoforms_hit.subtype = "at_least_one_novel_splicesite"
    else:  # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(
            itertools.chain(
                junctions_by_gene[ref_gene]
                for ref_gene in ref_genes
                if ref_gene in junctions_by_gene
            )
        )

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict(
            (i, all_ref_junctions.count(junc))
            for i, junc in enumerate(trec.junctions)
        )

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = (
                "mono-exon" if trec.exonCount == 1 else "multi-exon"
            )

    if has_intron_retention():
        isoforms_hit.subtype = "intron_retention"

    return isoforms_hit


def associationOverlapping(
    isoforms_hit: myQueryTranscripts,
    trec: genePredRecord,
    junctions_by_chr: dict,
    trans_by_region: list,
) -> myQueryTranscripts:
    """Check for antisense, genic, genic-intron or intergenic

    At this point definitely not FSM or ISM or NIC or NNC
    possibly (in order of preference assignment):
        - antisense  (on opp strand of a known gene)
        - genic (overlaps a combination of exons and introns), ignore strand
        - genic_intron  (completely within an intron), ignore strand
        - intergenic (does not overlap a gene at all), ignore strand

    Args:
        isoforms_hit (myQueryTranscript): best isoform hit at the moment
        trec (genePredRecord): query transcript
        junctions_by_chr (dict): dictionary with all the junctions in the region
        trans_by_region (list): list of all transcripts in the region

    Returns:
        isoforms_hit (myQueryTranscript): updated isoforms hit (myQueryTranscripts object)
    """

    isoforms_hit.str_class = "intergenic"
    isoforms_hit.transcripts = ["novel"]
    isoforms_hit.subtype = "mono-exon" if trec.exonCount == 1 else "multi-exon"

    if len(isoforms_hit.genes) == 0:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if len(isoforms_hit.AS_genes) == 0:
            if trec.chrom in junctions_by_chr:
                # no hit even on opp strand
                # see if it is completely contained within a junction
                da_pairs = junctions_by_chr[trec.chrom]["da_pairs"]
                i = bisect.bisect_left(da_pairs, (trec.txStart, trec.txEnd))
                while i < len(da_pairs) and da_pairs[i][0] <= trec.txStart:
                    if (
                        da_pairs[i][0]
                        <= trec.txStart
                        <= trec.txStart
                        <= da_pairs[i][1]
                    ):
                        isoforms_hit.str_class = "genic_intron"
                        for ref in trans_by_region:
                            if (
                                da_pairs[i][0] in ref.exonEnds
                                and da_pairs[i][1] in ref.exonStarts
                            ):
                                if ref.exonEnds.index(da_pairs[i][0]) == (
                                    ref.exonStarts.index(da_pairs[i][1]) - 1
                                ):
                                    isoforms_hit.genes = [ref.gene]
                        break
                    i += 1
            else:
                pass  # remain intergenic
        else:
            # hits one or more genes on the opposite strand
            isoforms_hit.str_class = "antisense"
            isoforms_hit.genes = isoforms_hit.AS_genes
    else:
        isoforms_hit.str_class = "genic"

    return isoforms_hit


def trans_overlap(r1: genePredRecord, r2: genePredRecord) -> bool:
    """Check if 2 trans overlap

    Args:
        r1 (genePredRecord): transcript 1
        r2 (genePredRecord): transcript 2

    Returns:
        bool: True if they overlap
              False if they do not
    """

    if r1.txStart <= r2.txEnd and r2.txStart <= r1.txEnd:
        return True
    else:
        return False


def hits_exon(r1: genePredRecord, r2: genePredRecord) -> bool:
    """Check if any exon of r2 is overlaped by r1 start and end

    IntervalTree.find(end, start): Return a sorted list of all intervals overlapping [start,end)
    Overlap, not within that's the point

    Args:
        r1 (genePredRecord): transcript 1
        r2 (genePredRecord): transcript 2

    Returns:
        bool: True if it hits an exon
              False if it does not
    """

    for e in r2.exons:
        if r1.txStart <= e.end and e.start <= r1.txEnd:
            return True
    return False


def summary_table_cat(data: dict):
    """Prints a summary table of the classification

    Args:
        data (dict): dictionary with all the transcripts by region
    """

    counts = defaultdict(
        lambda: 0,
        {
            "full-splice_match": 0,
            "incomplete-splice_match": 0,
            "novel_in_catalog": 0,
            "novel_not_in_catalog": 0,
            "fusion": 0,
            "antisense": 0,
            "genic_intron": 0,
            "genic": 0,
            "intergenic": 0,
        },
    )
    for chrom in data.values():
        for trans in chrom:
            counts[trans.str_class] += 1

    print("\033[94m_" * 79 + "\033[0m")
    print("\033[92mS Q A N T I - S I M\033[0m \U0001F4CA")
    print()
    print("Classification summary Table \U0001F50E")
    print("\033[94m_" * 79 + "\033[0m")
    for k, v in counts.items():
        print("\033[92m|\033[0m " + k + ": " + str(v))


def write_category_file(data: dict, out_name: str):
    """
    Writes the file with the structural category of each transcript and its reference

    Args:
        data (dict) dictionary with tanscripts classified
        out_name (str) out file name
    """

    f_out = open(out_name, "w")
    f_out.write("transcript_id\tgene_id\tstructural_category\tassociated_gene\tassociated_trans\tchrom\tstrand\texons\tdonors\tacceptors\tTSS_genomic_coord\tTTS_genomic_coord\tlength\n")

    for chrom in data.values():
        for trans in chrom:
            donors = []
            acceptors = []
            for d, a in trans.junctions:
                donors.append(d)
                acceptors.append(a)
            if isinstance(donors[0], str):
                donors = ["NA"]
                acceptors = ["NA"]
            else:
                donors = [str(d) for d in donors]
                acceptors = [str(a + 1) for a in acceptors]  # Change to 1-based exon start
            trans.tss += 1  # Change to 1-based transcript start

            if trans.str_class == "intergenic":
                if not trans.intergenic_assoc:
                    trans.intergenic_assoc = "novel"
                f_out.write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (
                        trans.id,
                        trans.gene_id,
                        trans.str_class,
                        trans.intergenic_assoc,
                        "_".join(trans.transcripts),
                        trans.chrom,
                        trans.strand,
                        trans.num_exons,
                        ",".join(donors),
                        ",".join(acceptors),
                        trans.tss,
                        trans.tts,
                        trans.length,
                    )
                )
            else:
                f_out.write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                    % (
                        trans.id,
                        trans.gene_id,
                        trans.str_class,
                        "_".join(trans.genes),
                        "_".join(trans.transcripts),
                        trans.chrom,
                        trans.strand,
                        trans.num_exons,
                        ",".join(donors),
                        ",".join(acceptors),
                        trans.tss,
                        trans.tts,
                        trans.length,
                    )
                )

    f_out.close()


def classify_gtf(args):
    """Classifies all transcripts from a GTF annotation

    Given a GTF classifies all transcripts in its potential structural category

    Returns:
        trans_info (dict): all transcripts classified in SQANTI3 structural categories
    """

    def initializer():
        global min_ref_len
        min_ref_len = 0

    # parsing transcripts from GTF
    print("[SQANTI-SIM] Parsing transcripts from GTF reference annotation file")
    trans_by_chr = gtf_parser(args.gtf)

    # classify transcripts
    print("[SQANTI-SIM] Classifying transcripts according to its SQANTI3 structural category")
    trans_info = defaultdict(lambda: [])

    if args.cores <= 1:
        initializer()
        for chrom in trans_by_chr:
            print(chrom)
            for record in tqdm(range(len(trans_by_chr[chrom]))):
                trans_by_region = trans_by_chr[chrom][record]
                tmp = transcript_classification(trans_by_region)
                for k in tmp:
                    trans_info[k].extend(tmp[k])

    else:  # multiprocessing
        all_regions = []
        for chrom in trans_by_chr:
            all_regions.extend(trans_by_chr[chrom])

        pool = mp.Pool(args.cores, initializer, ())
        tmp = pool.map(transcript_classification, all_regions)
        for x in tmp:
            for k, v in x.items():
                trans_info[k].extend(v)

    # Write category file
    print("[SQANTI-SIM] Writting structural category file")
    cat_out = os.path.join(args.dir, (args.output + "_index.tsv"))
    write_category_file(trans_info, cat_out)

    return trans_info
