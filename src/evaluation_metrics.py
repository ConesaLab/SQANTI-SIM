#!/usr/bin/env python3
"""
evaluation_metrics.py

Generate SQANTI-SIM report and metrics

Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)
"""

import os
import pandas
import subprocess
import sys
from collections import defaultdict
from src.SQANTI3.utilities.short_reads import get_TSS_bed, get_ratio_TSS, get_bam_header
from src.SQANTI3.sqanti3_qc import CAGEPeak, STARcov_parser
from time import strftime


def sqanti3_stats(args):
    """Runs SQANTI3 and generates SQANTI-SIM report

    Given the reconstructed transcripts in GTF format it runs the SQANTI3
    pipeline and computes the SQANTI-SIM metrics to evaluate how well it
    recovered the novel transcripts

    Args:
        args (list): arguments to parse
    """

    def write_whithin_cage(row):
        return within_cage_dict[row["transcript_id"]]

    def write_dist_cage(row):
        return dist_cage_dict[row["transcript_id"]]

    def write_ratio_TSS(row):
        if row["transcript_id"] in ratio_TSS_dict:
            return ratio_TSS_dict[row["transcript_id"]]["max_ratio_TSS"]
        else:
            return 1

    def write_SJ_cov(row):
        min_cov = "NA"
        if row["exons"] == 1:
            return min_cov
        d = row["donors"].split(",")
        a = row["acceptors"].split(",")
        for i in range(int(row["exons"]) - 1):
            sample_cov = SJcovInfo[row["chrom"], row["strand"]][(int(d[i]), (int(a[i])-1))] # make exon starts (SJ acceptors) 0 based
            total_coverage_unique = (
                sum(
                    [cov_uniq for (cov_uniq, cov_multi) in sample_cov.values()]
                )
                if SJcovInfo is not None
                else "NA"
            )
            if min_cov == "NA" or min_cov > total_coverage_unique:
                min_cov = total_coverage_unique
        return min_cov
    
    def write_TSS(row):
        return trans_start_end[row["isoform"]][0]
    def write_TTS(row):
        return trans_start_end[row["isoform"]][1]

    print("[SQANTI-SIM][%s] Running SQANTI3" %(strftime("%d-%m-%Y %H:%M:%S")))
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, "SQANTI3/sqanti3_qc.py")

    MIN_REF_LEN = 0
    cmd = [
        sqanti3,
        args.transcriptome,
        args.gtf,
        args.genome,
        "-o",
        args.output,
        "-d",
        os.path.join(args.dir, "sqanti3"),
        "--cpus",
        str(args.cores),
        "--min_ref_len",
        str(MIN_REF_LEN),
        "--report",
        "skip",
        "--force_id_ignore",
        "--skipORF"
    ]

    if args.CAGE_peak:
        cmd.append("--CAGE_peak")
        cmd.append(args.CAGE_peak)

    if args.short_reads:
        cmd.append("--short_reads")
        cmd.append(args.short_reads)

    if args.fasta:
        cmd.append("--aligner_choice")
        cmd.append(args.aligner_choice)
        cmd.append("--fasta")

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running SQANTI3: {0}".format(cmd), file=sys.stderr)
        #sys.exit(1)

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    if args.CAGE_peak:
        print("[SQANTI-SIM][%s] Parsing CAGE Peak data" %(strftime("%d-%m-%Y %H:%M:%S")))
        cage_peak_data = CAGEPeak(args.CAGE_peak)

        within_cage_dict = defaultdict(lambda: False)
        dist_cage_dict = defaultdict(lambda: False)
        with open(args.trans_index, "r") as index_file:
            header_names = index_file.readline()
            header_names = header_names.split()
            id_pos = header_names.index("transcript_id")
            chrom_pos = header_names.index("chrom")
            strand_pos = header_names.index("strand")
            start_pos = header_names.index("TSS_genomic_coord") # No need to swap start and end coordinates -> already swapped in the index file for negative strand
            for line in index_file:
                line = line.split()
                within_cage, dist_cage = cage_peak_data.find(
                    line[chrom_pos], line[strand_pos], (int(line[start_pos])-1)
                ) # 0 based TSS
                within_cage_dict[line[id_pos]] = within_cage
                dist_cage_dict[line[id_pos]] = dist_cage
        index_file.close()

        trans_index["dist_to_CAGE_peak"] = trans_index.apply(
            write_dist_cage, axis=1
        )
        trans_index["within_CAGE_peak"] = trans_index.apply(
            write_whithin_cage, axis=1
        )

    # Short Read Coverage
    if args.coverage:
        print("[SQANTI-SIM][%s] Parsing Coverage data" %(strftime("%d-%m-%Y %H:%M:%S")))
        
        SJcovNames, SJcovInfo = STARcov_parser(args.coverage)
        trans_index["min_cov"] = trans_index.apply(write_SJ_cov, axis=1)

    elif args.short_reads:
        print("[SQANTI-SIM][%s] Parsing Short Read data" %(strftime("%d-%m-%Y %H:%M:%S")))
        star_out = os.path.join(args.dir, "sqanti3/STAR_mapping/")
        star_index = os.path.join(args.dir, "sqanti3/STAR_index/")

        SJcovNames, SJcovInfo = STARcov_parser(star_out)
        trans_index["min_cov"] = trans_index.apply(write_SJ_cov, axis=1)
    
    # Short reads ratio TSS - Code adapted from SQANTI3
    if args.SR_bam:
        print("[SQANTI-SIM][%s] Calculating ratio TSS from provided BAM" %(strftime("%d-%m-%Y %H:%M:%S")))

        if os.path.isdir(args.SR_bam):
            bams = []
            for files in os.listdir(args.SR_bam):
                if files.endswith('.bam'):
                    bams.append(args.SR_bam + '/' + files)
        else:
            b = open(args.SR_bam , "r")
            bams = []
            for file in b:
                bams.append(file.rstrip())
        chr_order = get_bam_header(bams[0])
        inside_bed, outside_bed = get_TSS_bed(args.gtf, chr_order)
        ratio_TSS_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order)
        trans_index["ratio_TSS"] = trans_index.apply(write_ratio_TSS, axis=1)

    elif args.short_reads:
        chr_order = os.path.join(star_index, "chrNameLength.txt")
        inside_bed, outside_bed = get_TSS_bed(args.gtf, chr_order)
        bams = []
        for filename in os.listdir(star_out):
            if filename.endswith(".bam"):
                bams.append(star_out + "/" + filename)
        ratio_TSS_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order)
        trans_index["ratio_TSS"] = trans_index.apply(write_ratio_TSS, axis=1)

    trans_index.to_csv(
        args.trans_index, sep="\t", na_rep="NA", header=True, index=False
    )

    print("[SQANTI-SIM][%s] Generating SQANTI-SIM report" %(strftime("%d-%m-%Y %H:%M:%S")))
    src_dir = os.path.dirname(os.path.realpath(__file__))
    classification_file = os.path.join(args.dir, "sqanti3/",(args.output + "_classification.txt"))
    junctions_file = os.path.join(args.dir, "sqanti3/", (args.output + "_junctions.txt"))
    corrected_genePred = os.path.join(args.dir, "sqanti3/", (args.output + "_corrected.genePred"))

    # Add TSS and TTS genomic coords to classification file
    # GenePred format -> https://genome.ucsc.edu/FAQ/FAQformat.html#format9
    trans_start_end = defaultdict(lambda: [None, None])
    with open(corrected_genePred, "r") as gp_file:
        for line in gp_file:
            line = line.split()
            name = line[0]
            strand = line[2]
            if strand == "+":
                txStart = int(line[3]) + 1 # Turn genePred 0-based start to 1-based 
                txEnd = line[4]
            else:
                txStart = line[4]
                txEnd = int(line[3]) + 1

            trans_start_end[name] = [txStart, txEnd]

    classif_f = pandas.read_csv(classification_file, sep="\t", header=0, dtype={"chrom":str})
    classif_f["TSS_genomic_coord"] = classif_f.apply(write_TSS, axis=1)
    classif_f["TTS_genomic_coord"] = classif_f.apply(write_TTS, axis=1)
    classif_f.to_csv(classification_file, sep="\t", header=True, index=False, na_rep="NA")

    # Generate SQANTI-SIM report
    cmd = [
        "Rscript",
        os.path.join(src_dir, "SQANTI-SIM_report.R"),
        classification_file,
        junctions_file,
        args.trans_index,
        str(args.min_support),
        src_dir,
        args.expression
    ]

    cmd = " ".join(cmd)
    if subprocess.check_call(cmd, shell=True) != 0:
        print(
            "[SQANTI-SIM] ERROR running SQANTI-SIM report generation: {0}".format(cmd),
            file=sys.stderr,
        )
        sys.exit(1)
