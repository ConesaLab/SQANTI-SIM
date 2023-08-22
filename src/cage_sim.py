#!/usr/bin/env python3
"""
cage_sim.py

SQANTI-SIM CAGE bed file simulator.

Author: Jorge Mestre Tomas (jorge.mestre.tomas@csic.es)
"""

import argparse
import csv
import os
import pandas
import pickle
import pybedtools
import numpy
import random
import re
import subprocess
import sys
from BCBio import GFF as BCBio_GFF
from bx.intervals import IntervalTree
from collections import defaultdict
from sklearn.linear_model import LogisticRegression


class CAGEPeak:
    """
    This class is adapted from the original SQANTI3 QC.
    (sqanti3_qc.py -> https://github.com/ConesaLab/SQANTI3)
    """
    def __init__(self, cage_bed_filename):
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks
        self.cage_length = defaultdict(lambda: 0)

        self.read_bed()

    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            cageid = str(raw[3])
            strand = raw[5]
            tss0 = int((start0+end1)/2)
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (cageid, tss0, start0, end1))
            self.cage_length[cageid] = (abs(end1-start0))-1

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within a cage peak>, <nearest dist to TSS>
        dist to TSS is 0 if right on spot
        dist to TSS is + if downstream, - if upstream (watch for strand!!!)
        """
        matched_CAGE = set()
        hit, within_peak, dist_peak = "NA", False, "NA"
        for (cageid, tss0,start0,end1) in self.cage_peaks[(chrom,strand)].find(query-search_window, query+search_window):
            # Skip those cage peaks that are downstream the detected TSS because degradation just make the transcript shorter
            if strand=="+" and start0>int(query) and end1>int(query):
                continue
            if strand=="-" and start0<int(query) and end1<int(query):
                continue

            within_out = (start0<=query<end1)
            if within_out:
                w = True
                matched_CAGE.add(cageid)
            else:
                w = False

            if not within_peak==True:
                hit, within_peak, dist_peak = cageid, w, (query - tss0) * (-1 if strand=="-" else +1)
            else:
                d = (query - tss0) * (-1 if strand=="-" else +1)
                if abs(d) < abs(dist_peak):
                   hit, within_peak, dist_peak = cageid, w, d
        return within_peak, dist_peak, hit, matched_CAGE

def get_bam_header(bam):
    # This code is taken from the original SQANTI3 QC
    o_dir=os.path.dirname(bam)
    out=o_dir + "/chr_order.txt"
    os.system("samtools view -H {b} | grep '^@SQ' | sed 's/@SQ\tSN:\|LN://g'  > {o}".format(b=bam, o=out))
    return(out)


def get_TSS_bed(corrected_gtf, chr_order):
    # This code is taken from the original SQANTI3 QC
    limit_info = dict(gff_type=["transcript"])
    out_directory=os.path.dirname(corrected_gtf)
    tmp_in= out_directory + "/coverage_inside_TSS.bed_tmp"
    tmp_out = out_directory + "/coverage_outside_TSS.bed_tmp"
    in_handle=open(corrected_gtf)
    with open(tmp_in,"w") as inside:
        with open(tmp_out,"w") as outside:
            for rec in BCBio_GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
                chr=rec.id
                iso_id=rec.features[0].qualifiers["transcript_id"][0]
                loc=str(rec.features[0].location)
                loc=re.split("[\(\)\[\]\:]",loc)
                loc=list(filter(None,loc))
                strand=str(loc[2])
                if strand=="+":
                    start_in=int(loc[0])
                    end_in=int(loc[0])+100
                    start_out=int(loc[0])-101
                    end_out=int(loc[0])-1
                    if start_out<0:
                        start_out=0
                else:
                    start_in=int(loc[1])-100
                    end_in=int(loc[1])
                    start_out=int(loc[1])+1
                    end_out=int(loc[1])+101
                if end_out<=0 or start_in<=0: 
                    print('{iso} will not be included in TSS ratio calculation since its TSS is located at the very beginning of the chromosome'.format(iso=iso_id))
                else:
                    inside.write(chr + "\t" + str(start_in) + "\t" + str(end_in) + "\t" + iso_id + "\t0\t" + strand + "\n")
                    outside.write(chr + "\t" + str(start_out) + "\t" + str(end_out) + "\t" + iso_id + "\t0\t" + strand + "\n")
    in_handle.close()
    i = pybedtools.BedTool(tmp_in)
    o = pybedtools.BedTool(tmp_out)
    inside_sorted = out_directory + "/inside_TSS.bed" 
    outside_sorted = out_directory + "/outside_TSS.bed"
    i.sort(g=chr_order, output=inside_sorted)
    o.sort(g=chr_order, output=outside_sorted) 
    os.system("rm {i} {o}".format(i=tmp_in , o=tmp_out))
    return(inside_sorted, outside_sorted)


def get_ratio_TSS(inside_bed, outside_bed, replicates, chr_order):
    # This code is taken from the original SQANTI3 QC 
    ## the idea would be to first calculate the average coverage per sample for in and out beds. Calculate each ratio
    # ## get the maximum the ratios across replicates and return it as a dictionary
    out_TSS_file = os.path.dirname(inside_bed) + "/ratio_TSS.csv"
    in_bed = pybedtools.BedTool(inside_bed)
    out_bed = pybedtools.BedTool(outside_bed)
    for b in [*range(0,len(replicates))]:
        bam_file=replicates[b] 
        in_cov = in_bed.coverage(bam_file, sorted=True, g=chr_order)
        out_cov = out_bed.coverage(bam_file, sorted=True, g=chr_order)
        inside_df = pandas.DataFrame(columns=["id","inside"])
        for entry in in_cov:
            new_entry = pandas.DataFrame({"id" : [entry.name] , "inside" : [float(entry[6])]})
            inside_df = pandas.concat([inside_df,new_entry], ignore_index=True)
        outside_df = pandas.DataFrame(columns=["id","outside"])
        for entry in out_cov:
            new_entry = pandas.DataFrame({"id" : [entry.name] , "outside" : [float(entry[6])]})
            outside_df = pandas.concat([outside_df, new_entry], ignore_index=True)
        merged = pandas.merge(inside_df, outside_df, on="id")
        merged["ratio_TSS"] = (merged["inside"]+0.01)/(merged["outside"]+0.01)
        merged["ratio_TSS"] = pandas.to_numeric(merged["ratio_TSS"])
        if b == 0 :
            ratio_rep_df = merged["id"]
        ratio_rep_df = pandas.merge(ratio_rep_df, merged[["id","ratio_TSS"]], on="id")
    #ratio_rep_df.to_csv(out_TSS_file, index=False)
    ratio_rep_df["max_ratio_TSS"]=ratio_rep_df.max(axis=1, numeric_only=True)
    ratio_rep_df = ratio_rep_df[["id","max_ratio_TSS"]]
    ratio_rep_dict = ratio_rep_df.set_index("id").T.to_dict()
    os.system('rm {i} {o}'.format(i=inside_bed, o=outside_bed))
    print('Temp files removed.\n')
    return(ratio_rep_dict)

def get_TSS_cov(gtf, bam_files, chr_order):
    """
    Calculate TSS coverage for a set of BAM files using a corrected GTF file.

    Args:
        gtf (str): Path to the GTF file.
        bam_files (list): List of BAM file paths.
        chr_order (str): Path to the chromosome order file.

    Returns:
        dict: Dictionary mapping isoform IDs to TSS coverage values.
    """

    # Define the limit information for GFF parsing and open GTF
    out_directory=os.path.dirname(gtf)
    tmp_bed= out_directory + "/TSS_cov.bed_tmp"
    limit_info = dict(gff_type=["transcript"])

    in_handle=open(gtf)
    with open(tmp_bed,"w") as TSS_cov:
        for rec in BCBio_GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
            chr=rec.id
            iso_id=rec.features[0].qualifiers["transcript_id"][0]
            loc=str(rec.features[0].location)
            loc=re.split("[\(\)\[\]\:]",loc)
            loc=list(filter(None,loc))
            strand=str(loc[2])
            if strand=="+":
                start=int(loc[0])
                end=int(loc[0]) + 20 # Coverage 20 bp downstream
            else:
                start=int(loc[1]) - 20
                end=int(loc[1])
            if start < 0:
                start = 0
            else:
                # Write the TSS information to the temporary BED file
                TSS_cov.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + iso_id + "\t0\t" + strand + "\n")
    in_handle.close()

    # Create and sort the BED file based on the chromosome order
    bed = pybedtools.BedTool(tmp_bed)
    bed_sorted = bed.sort(g=chr_order)
    
    # Calculate coverage of the sorted BED file using the BAM files
    bed_cov = bed_sorted.coverage(bam_files, sorted=True, counts=True, split=True, g=chr_order)
    TSS_cov_df = pandas.DataFrame(columns=["id","TSS_cov"])
    for entry in bed_cov:
        new_entry = pandas.DataFrame({"id" : [entry.name] , "TSS_cov" : [float(entry[6])]})
        TSS_cov_df = pandas.concat([TSS_cov_df,new_entry], ignore_index=True)
    TSS_cov_df["TSS_cov"] = pandas.to_numeric(TSS_cov_df["TSS_cov"])
    TSS_cov_dict = TSS_cov_df.set_index("id").T.to_dict()
        
    os.system("rm {b}".format(b=tmp_bed))

    return(TSS_cov_dict)


def star(genome, SR_fofn, output_dir, cores):
    # This code is taken from the original SQANTI3 QC 
    fasta_genome = genome #Fasta Format already checked
    index_dir = output_dir + "/STAR_index/"
    index_dir_tmp = index_dir + "/_STARtmp/"
    index_dir_o = index_dir + "SAindex" 
    mapping_dir = output_dir + "/STAR_mapping/"
    print('[SQANTI-SIM] START running STAR...')
    if not os.path.exists(mapping_dir):
        os.makedirs(mapping_dir)
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)
        if not os.path.exists(index_dir_o):
            print('Running indexing...')
            subprocess.call(["STAR", "--runThreadN", str(cores), "--runMode", "genomeGenerate", "--genomeDir", index_dir, "--genomeFastaFiles", fasta_genome, "--outTmpDir", index_dir_tmp])
            print('Indexing done.')
    else:
        print('Index identified. Proceeding to mapping.')
    star_mapping(index_dir, SR_fofn, output_dir, cores)
    return(mapping_dir, index_dir)


def star_mapping(index_dir, SR_fofn, output_dir, cores):
    # This code is taken from the original SQANTI3 QC 
    mapping_dir = output_dir + "/STAR_mapping"
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = [x.strip() for x in line.split(' ')]
            compressed = False
            if files[0][-3:] == ".gz":
                compressed = True
            if compressed :
                sample_name = os.path.splitext(files[0])[-2].split("/")[-1]
                sample_name = sample_name.split(".")[:-1][0]
            else:
                sample_name = os.path.splitext(files[0])[-2].split("/")[-1]
            sample_prefix = mapping_dir + "/" + sample_name
            if not os.path.exists(sample_prefix + "Log.final.out"):
                print('[SQANTI-SIM] Mapping for ', sample_name, ': in progress...')
                if not compressed:
                    if len(files) == 1:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(["STAR", "--runThreadN",  str(cores), "--genomeDir", index_dir, "--readFilesIn", files[0], "--outFileNamePrefix", sample_prefix,"--alignSJoverhangMin", "8", "--alignSJDBoverhangMin", "1", "--outFilterType", "BySJout", "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "20", "--outFilterMismatchNoverLmax", "0.04", "--outFilterMismatchNmax", "999", "--alignIntronMin", "20", "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000", "--sjdbScore", "1", "--genomeLoad", "NoSharedMemory", "--outSAMtype", "BAM", "SortedByCoordinate", "--twopassMode", "Basic"])
                    else:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(["STAR", "--runThreadN", str(cores),  "--genomeDir", index_dir, "--readFilesIn", files[0], files[1], "--outFileNamePrefix", sample_prefix,"--alignSJoverhangMin", "8", "--alignSJDBoverhangMin", "1", "--outFilterType", "BySJout", "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "20", "--outFilterMismatchNoverLmax", "0.04", "--outFilterMismatchNmax", "999", "--alignIntronMin", "20", "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000", "--sjdbScore", "1", "--genomeLoad", "NoSharedMemory", "--outSAMtype", "BAM", "SortedByCoordinate", "--twopassMode", "Basic"])
                else:
                    if len(files) == 1:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(["STAR", "--runThreadN",  str(cores), "--genomeDir", index_dir, "--readFilesIn", files[0], "--outFileNamePrefix", sample_prefix,"--alignSJoverhangMin", "8", "--alignSJDBoverhangMin", "1", "--outFilterType", "BySJout", "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "20", "--outFilterMismatchNoverLmax", "0.04", "--outFilterMismatchNmax", "999", "--alignIntronMin", "20", "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000", "--sjdbScore", "1", "--genomeLoad", "NoSharedMemory", "--outSAMtype", "BAM", "SortedByCoordinate", "--readFilesCommand", "zcat", "--twopassMode", "Basic"])
                    else:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(["STAR", "--runThreadN", str(cores), "--genomeDir", index_dir, "--readFilesIn", files[0], files[1], "--outFileNamePrefix", sample_prefix,"--alignSJoverhangMin", "8", "--alignSJDBoverhangMin", "1", "--outFilterType", "BySJout", "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "20", "--outFilterMismatchNoverLmax", "0.04", "--outFilterMismatchNmax", "999", "--alignIntronMin", "20", "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000", "--sjdbScore", "1", "--genomeLoad", "NoSharedMemory", "--outSAMtype", "BAM", "SortedByCoordinate", "--readFilesCommand", "zcat", "--twopassMode", "Basic"])


def sample_position(chrom_TSS_pos: dict, dist2TSS: list, cage_length: list) -> tuple:
    """
    Samples positions based on the distance to transcription start sites (TSS) in each chromosome.

    Args:
        chrom_TSS_pos (dict): Dictionary containing chromosome-wise TSS positions.
            Keys: Tuple of (chromosome, strand).
            Values: List of TSS positions for the corresponding chromosome and strand.
        dist2TSS (list): List of distances to TSS. Specifies the desired distance from the TSS for each position.
        cage_length (list): List of CAGE peak lengths. Specifies the peak length for each sampled position.

    Returns:
        Tuple: A tuple containing lists of sampled CAGE peak chromosome, strand, coordinate, distance to TSS, and peak length.
            - cage_peak_chrom (list): List of chromosomes for sampled CAGE peaks.
            - cage_peak_strand (list): List of strands for sampled CAGE peaks.
            - cage_peak_coord (list): List of coordinates for sampled CAGE peaks.
            - cage_peak_distance (list): List of distances to TSS for sampled CAGE peaks.
            - cage_peak_length (list): List of peak lengths for sampled CAGE peaks.
    """
     
    chrom_intervals = {}
    cage_peak_chrom = []
    cage_peak_strand = []
    cage_peak_coord = []
    cage_peak_distance = []
    cage_peak_length = []

    # Calculate chromosome intervals
    for chrom_strand, TSS_pos in chrom_TSS_pos.items():
        unique_sorted_pos = numpy.sort(numpy.unique(TSS_pos))
        chrom_TSS_pos[chrom_strand] = unique_sorted_pos
        chrom_intervals[chrom_strand] = numpy.diff(unique_sorted_pos)

    # Sample coordinates that are at least k bp away from existing TSS positions
    for i in range(len(dist2TSS)):
        k = dist2TSS[i]
        # Find potential positions based on distance
        potential_pos = {chrom_strand: numpy.where(intervals > abs(2*k))[0] for chrom_strand, intervals in chrom_intervals.items() if numpy.any(intervals > abs(2*k))}
        result_dict = {}
        if k > 0:
            for chrom_strand, positions in potential_pos.items():
                if chrom_strand[1] == "+":
                    result_dict[chrom_strand] = numpy.asarray(chrom_TSS_pos[chrom_strand])[positions+1].tolist()
                else:
                    result_dict[chrom_strand] = numpy.asarray(chrom_TSS_pos[chrom_strand])[positions].tolist()
        else:
            for chrom_strand, positions in potential_pos.items():
                if chrom_strand[1] == "+":
                    result_dict[chrom_strand] = numpy.asarray(chrom_TSS_pos[chrom_strand])[positions].tolist()
                else:
                    result_dict[chrom_strand] = numpy.asarray(chrom_TSS_pos[chrom_strand])[positions+1].tolist()
        potential_chrom_coord = [(chrom, strand, coord) for (chrom, strand), coords in result_dict.items() for coord in coords]
        
        if len(potential_chrom_coord) == 0:
            continue
        
        # Randomly sample a position
        index = random.randint(0, len(potential_chrom_coord) - 1)
        sampled_chrom, sampled_strand, sampled_coord = list(potential_chrom_coord)[index]
        cage_peak_chrom.append(sampled_chrom)
        cage_peak_strand.append(sampled_strand)
        cage_peak_coord.append(sampled_coord)
        cage_peak_distance.append(k)
        cage_peak_length.append(cage_length[i])

    return cage_peak_chrom, cage_peak_strand, cage_peak_coord, cage_peak_distance, cage_peak_length

def CAGEsim_bed(sim_trans, pa_model, true_peaks_df, false_peaks_df, prop_false_peaks: float = 0.2):
    """
    Simulates CAGE peaks based on the provided simulation parameters and data.

    Args:
        sim_trans (pd.DataFrame): DataFrame containing simulated transcript data.
        pa_model: The predictive model for determining if a CAGE peak is present.
        true_peaks_df (pd.DataFrame): DataFrame containing the length of true CAGE peaks.
        false_peaks_df (pd.DataFrame): DataFrame containing the length of false CAGE peaks.
        prop_false_peaks (float): Proportion of false CAGE peaks to simulate (default: 0.2).

    Returns:
        Tuple: A tuple containing the following:
            - DataFrame: DataFrame containing the simulated CAGE peaks in BED format.
            - List of tuples: Each tuple contains a CAGE peak ID and its corresponding transcript ID.
    """

    # Dict with isoform TSS coord (value) in each (chrom, strand) (key)
    chrom_TSS_pos = sim_trans.groupby(["chrom", "strand"])["TSS_genomic_coord"].apply(list).to_dict()

    # Predict if there is a CAGE peak supporting each simulated isoform
    sim_trans["ratio_TSS"] = numpy.log(sim_trans["ratio_TSS"].astype('float'))
    sim_trans["TSS_cov"] = numpy.log((sim_trans["TSS_cov"].astype('float') + 1) / numpy.sum(sim_trans["TSS_cov"].astype('float') + 1))
    sim_trans["within_CAGE_peak"] = pa_model.predict(sim_trans[["ratio_TSS", "TSS_cov"]].to_numpy()).tolist()

    # For those CAGE peaks predicted before ("true"), sample their peak length and distance to TSS
    sim_trans = sim_trans.loc[sim_trans["within_CAGE_peak"] == 1].copy()

    # Resample from ECDF TRUE CAGE peaks
    true_unique_pairs = true_peaks_df[["CAGE_peak_length", "distance"]].values
    true_joint_probs = true_peaks_df["ecdf"].values
    true_rand_ind = numpy.random.choice(len(true_unique_pairs), size=len(sim_trans), replace=True, p=true_joint_probs)
    sim_true_len_dist = true_unique_pairs[true_rand_ind]

    sim_trans["peak_length"] = sim_true_len_dist[:, 0]
    sim_trans["dist2TSS"] = sim_true_len_dist[:, 1]
    # Genomic coord of center CAGE peak
    sim_trans["CAGE_center"] = sim_trans["TSS_genomic_coord"] + numpy.where(sim_trans["strand"] == "+", - sim_trans["dist2TSS"], sim_trans["dist2TSS"])

    # Create BED file
    cond1 = (sim_trans["dist2TSS"] > 0) & (sim_trans["strand"] == "+")
    cond2 = (sim_trans["dist2TSS"] < 0) & (sim_trans["strand"] == "-")
    offset_start = numpy.where(cond1 | cond2,
                    numpy.floor(sim_trans["peak_length"] / 2),
                    numpy.ceil(sim_trans["peak_length"] / 2))
    
    offset_end = numpy.where(cond1 | cond2,
                    numpy.ceil(sim_trans["peak_length"] / 2),
                    numpy.floor(sim_trans["peak_length"] / 2))

    CAGE_true_sim = pandas.DataFrame({
        "chrom": sim_trans["chrom"],
        "start": (sim_trans["CAGE_center"] - offset_start).astype("int"),
        "end": (sim_trans["CAGE_center"] + offset_end).astype("int") + 1,
        "feature_id": "",
        "score": 0,
        "strand": sim_trans["strand"]
    })
    CAGE_true_sim["feature_id"] = CAGE_true_sim["chrom"] + "_" + CAGE_true_sim["start"].astype(str) + "_" + CAGE_true_sim["end"].astype(str) + "_" + CAGE_true_sim["strand"]
    print("[SQANTI-SIM] Simulated", len(CAGE_true_sim), "true CAGE peaks", file=sys.stdout)

    # Simulate "false" CAGE peaks (not hitting TSS)
    n_false_peaks = round(len(sim_trans) * prop_false_peaks / (1 - prop_false_peaks))
        
    # Generate random variable by resampling from ECDF
    false_unique_pairs = false_peaks_df[["CAGE_peak_length", "distance"]].values
    false_joint_probs = false_peaks_df["ecdf"].values
    false_rand_ind = numpy.random.choice(len(false_unique_pairs), size=n_false_peaks, replace=True, p=false_joint_probs)
    sim_false_len_dist = false_unique_pairs[false_rand_ind]

    false_peak_length = sim_false_len_dist[:, 0]
    false_peak_TSSdist = sim_false_len_dist[:, 1]

    # Decide where to place those CAGE peaks
    false_peak_chrom, false_peak_strand, false_peak_coord, false_peak_TSSdist, false_peak_length = sample_position(chrom_TSS_pos, false_peak_TSSdist, false_peak_length)
    false_peak_strand = numpy.array(false_peak_strand)
    false_peak_coord = numpy.array(false_peak_coord)
    false_peak_TSSdist = numpy.array(false_peak_TSSdist)
    false_peak_length = numpy.array(false_peak_length)

    # Create BED file (pd.Dataframe)
    false_peak_center = false_peak_coord + numpy.where(false_peak_strand == "+", - false_peak_TSSdist, false_peak_TSSdist)
    CAGE_false_sim = pandas.DataFrame({
        "chrom": false_peak_chrom,
        "start": (false_peak_center - numpy.ceil(false_peak_length / 2)).astype("int"),
        "end": (false_peak_center + numpy.floor(false_peak_length / 2)).astype("int") + 1,
        "feature_id": "",
        "score": 0,
        "strand": false_peak_strand
    })
    CAGE_false_sim["feature_id"] = CAGE_false_sim["chrom"] + "_" + CAGE_false_sim["start"].astype(str) + "_" + CAGE_false_sim["end"].astype(str) + "_" + CAGE_false_sim["strand"]
    print("[SQANTI-SIM] Simulated", len(CAGE_false_sim), "false CAGE peaks", file=sys.stdout)

    CAGE_bed_file = pandas.concat([CAGE_true_sim, CAGE_false_sim])

    # CAGE to transcript correspondence
    cage_to_trans = []
    for cageid, transid in zip(CAGE_true_sim["feature_id"], sim_trans["transcript_id"]):
        cage_to_trans.append((cageid, transid))
    for cageid in CAGE_false_sim["feature_id"]:
        cage_to_trans.append((cageid, "NA"))

    return CAGE_bed_file, cage_to_trans

def main():
    parser = argparse.ArgumentParser(prog="cage_sim.py", description="cage_sim.py parse options", )
    subparsers = parser.add_subparsers(dest="mode", description="\t\tTwo modes: train (provide a reconstructed transcriptome and sample-specific short-read RNA-seq and CAGE bed file to generate a custom train model) and sim (simulate CAGE peak bed file for your isoforms given a pre-trained model)")

    parser_t = subparsers.add_parser("train", help="\t\tRun in train mode")
    parser_t.add_argument("--gtf", help='\t\tReconstructed transcriptome (GTF format)')
    parser_t.add_argument("--genome", help='\t\tReference genome (Fasta format)')
    parser_t.add_argument("--CAGE_peak", help='\t\tSample-specific Cage Peaks (BED format)')
    parser_t.add_argument("--short_reads", help='\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.')
    parser_t.add_argument("--SR_bam" , help='\t\t Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome')
    parser_t.add_argument("-d","--dir", help='\t\tDirectory for output files. Default: Directory where the script was run.')
    parser_t.add_argument("-t", "--cores", default=8, type=int, help='\t\tNumber of threads used during alignment by aligners. (default: 8)')

    parser_s = subparsers.add_parser("sim", help="\t\tRun in sim mode")
    parser_s.add_argument("-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM (*_index.tsv)", )
    parser_s.add_argument("--gtf", help='\t\tReference annotation from where transcripts were simulated  (GTF format)')
    parser_s.add_argument("--genome", help='\t\tReference genome (Fasta format)')
    parser_s.add_argument("--short_reads", help='\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.')
    parser_s.add_argument("--SR_bam" , help='\t\t Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome')
    parser_s.add_argument("--CAGE_model", type=str, default=None, help="\t\tDirectory of the pre-trained CAGE model")
    parser_s.add_argument("--falseCAGE_prop", type=float, default=0.2, help="\t\tProportion (0, 1) of simulated CAGE peaks that are not derived from actual TSS locations (default: 0.2)")
    parser_s.add_argument("-d","--dir", help='\t\tDirectory for output files. Default: Directory where the script was run.')
    parser_s.add_argument("-t", "--cores", default=8, type=int, help='\t\tNumber of threads used during alignment by aligners. (default: 8)')
    parser_s.add_argument("-s", "--seed", type=int, default=None, help="\t\tRandomizer seed", )

    args, unknown = parser.parse_known_args()

    if unknown:
        print("[SQANTI-SIM] cage_sim.py train mode unrecognized arguments: {}\n".format(" ".join(unknown)), file=sys.stderr)

    if args.mode == "sim":
        if not 0 <= args.falseCAGE_prop < 1:
            raise ValueError("--falseCAGE_prop: value should be between 0 and 1.")
    
        if args.seed:
            random.seed(args.seed)
            numpy.random.seed(args.seed)

    
    # 1. Prepare short-read alignments
    if  args.SR_bam is not None:
        print("[SQANTI-SIM] Using provided BAM files for calculating TSS ratio", file=sys.stdout)
        if os.path.isdir(args.SR_bam):
            bams = []
            for files in os.listdir(args.SR_bam):
                if files.endswith(".bam"):
                    bams.append(args.SR_bam + "/" + files)
        else:
            b = open(args.SR_bam , "r")
            bams = []
            for file in b:
                bams.append(file.rstrip())
        chr_order = get_bam_header(bams[0])

    elif args.short_reads is not None:
            print("[SQANTI-SIM] Starting STAR mapping. We need this for calculating TSS ratio. It may take some time...")
            star_out, star_index = star(args.genome, args.short_reads, args.dir, args.cores)
            chr_order = star_index + "/chrNameLength.txt"
            bams=[]
            for filename in os.listdir(star_out):
                if filename.endswith(".bam"):
                    bams.append(star_out + "/" + filename)
    print('[SQANTI-SIM] BAM files identified: '+str(bams))

    # 2. Calculate TSS ratio for each isoform
    print("[SQANTI-SIM] Calculating TSS ratio", file=sys.stdout)
    inside_bed, outside_bed = get_TSS_bed(args.gtf, chr_order)
    ratio_rep_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order)
    
    # 3. Calculate TSS coverage (20 bp downstream)
    print("[SQANTI-SIM] Calculating TSS coverage", file=sys.stdout)
    TSS_cov_dict = get_TSS_cov(args.gtf, bams, chr_order)

    if args.mode == "train":
        print("[SQANTI-SIM] Characterizing sample-specific CAGE data")
        # 4. Characterize sample-specific CAGE data
        cage_peak_obj = CAGEPeak(args.CAGE_peak)

        matched_CAGE = set()
        TSS_by_chrom = defaultdict(lambda: [])
        CAGE_char = []

        # Isoforms with and without CAGE support
        limit_info = dict(gff_type=["transcript"])
        in_handle=open(args.gtf)
        for rec in BCBio_GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
            chrom=rec.id
            trans_id=rec.features[0].qualifiers["transcript_id"][0]
            loc=str(rec.features[0].location)
            loc=re.split("[\(\)\[\]\:]",loc)
            loc=list(filter(None,loc))
            strand=str(loc[2])
            if strand == "+":
                TSS_coord = int(loc[0])
            else:
                TSS_coord = int(loc[1])

            # Check if there's a CAGE peak supporting the transcript TSS
            within_CAGE, dist_CAGE, CAGE_hit, CAGE_matches = cage_peak_obj.find(chrom, strand, TSS_coord) 
            matched_CAGE = matched_CAGE.union(CAGE_matches)
            CAGE_char.append([trans_id, chrom, TSS_coord, ratio_rep_dict[trans_id]["max_ratio_TSS"], TSS_cov_dict[trans_id]["TSS_cov"], dist_CAGE, cage_peak_obj.cage_length[CAGE_hit], within_CAGE, True, CAGE_hit])
            TSS_by_chrom[(chrom, strand)].append(TSS_coord)
        CAGE_char = pandas.DataFrame(CAGE_char)
        CAGE_char.columns=["id", "chrom", "pos", "ratio_TSS", "TSS_cov", "distance", "CAGE_peak_length", "CAGE", "TSS", "CAGE_id"]

        # Add CAGE peaks not hitting any isoform TSS
        CAGE_char2 = []
        for line in open(args.CAGE_peak):
                raw = line.strip().split()
                chrom = raw[0]
                start0 = int(raw[1])
                end1 = int(raw[2])
                cageid = str(raw[3])
                strand = raw[5]
                tss0 = int((start0+end1)/2)
                cage_length = (abs(end1-start0))-1

                if cageid not in matched_CAGE and (chrom, strand) in TSS_by_chrom:
                    dist_CAGE = (min(TSS_by_chrom[(chrom, strand)], key=lambda x:abs(x-tss0)) - tss0) * (-1 if strand=="-" else +1)
                    CAGE_char2.append([cageid, chrom, tss0, "NA", "NA", dist_CAGE, cage_length, True, False, cageid])
        CAGE_char2 = pandas.DataFrame(CAGE_char2)
        CAGE_char2.columns=["id", "chrom", "pos", "ratio_TSS", "TSS_cov", "distance", "CAGE_peak_length", "CAGE", "TSS", "CAGE_id"]

        CAGE_char = pandas.concat([CAGE_char, CAGE_char2])
        CAGE_char.to_csv(os.path.join(args.dir, "tmp_CAGE_train.tsv"), sep="\t", header=True, index=False, na_rep="NA")
        CAGE_char= pandas.read_csv(os.path.join(args.dir, "tmp_CAGE_train.tsv"), sep="\t")

        # 5. Presence/Absence model fitting
        # Transform variables
        iso_peaks = CAGE_char[CAGE_char["TSS"] == True].copy()
        iso_peaks["ratio_TSS"] = numpy.log(iso_peaks["ratio_TSS"].astype('float')) # log ratio TSS
        iso_peaks["TSS_cov"] = numpy.log((iso_peaks["TSS_cov"].astype('float') + 1) / (iso_peaks["TSS_cov"].astype('float') + 1).sum()) # log of smoothed (+1) proportion

        X = iso_peaks[["ratio_TSS", "TSS_cov"]].to_numpy()
        y = iso_peaks.loc[iso_peaks["TSS"], "CAGE"].to_numpy().ravel().astype(int)
        pa_model = LogisticRegression(solver="newton-cg")
        pa_model.fit(X, y)

        # 6. Save sample-specific CAGE model
        with open(os.path.join(args.dir, "CAGE_pa_model.pkl"), "wb") as model_file:
            pickle.dump(pa_model, model_file)
        model_file.close()

        # 7. Compute bivariate ECDF for (PeakLength, TSSdistance) for peaks hitting isoforms TSS (true) and those not hitting (false)
        true_cage_values = CAGE_char.loc[(CAGE_char["TSS"] == True) & (CAGE_char["CAGE"] == True), ["CAGE_peak_length", "distance"]].values
        false_cage_values = CAGE_char.loc[(CAGE_char["TSS"] == False) & (CAGE_char["CAGE"] == True), ["CAGE_peak_length", "distance"]].values

        # Compute unique pairs and their counts for observed data
        true_unique_pairs, true_counts = numpy.unique(true_cage_values, axis=0, return_counts=True)
        false_unique_pairs, false_counts = numpy.unique(false_cage_values, axis=0, return_counts=True)

        # Calculate empirical joint probabilities for observed data
        true_joint_probs = true_counts / numpy.sum(true_counts)
        false_joint_probs = false_counts / numpy.sum(false_counts)

        # 8. Save ECDF to file
        true_cage_values_df = pandas.DataFrame({"CAGE_peak_length": true_unique_pairs[:, 0], "distance": true_unique_pairs[:, 1], "ecdf": true_joint_probs})
        false_cage_values_df = pandas.DataFrame({"CAGE_peak_length": false_unique_pairs[:, 0], "distance": false_unique_pairs[:, 1], "ecdf": false_joint_probs})
        
        true_cage_values_df.to_csv(os.path.join(args.dir, "cage_true_peaks_ecdf.tsv"), index=False, sep="\t")
        false_cage_values_df.to_csv(os.path.join(args.dir, "cage_false_peaks_ecdf.tsv"), index=False, sep="\t")
        
        print("[SQANTI-SIM] cage_sim.py train done!")

    elif args.mode == "sim":
        print("[SQANTI-SIM] Loading pre-trained model")
        # 4. Prepare input
        if not args.CAGE_model:
            src_dir = os.path.dirname(os.path.realpath(__file__))
            models = os.path.join(src_dir, "../pre-trained_models/")
            model_name = "human_WTC11_PacBio_FLAIR_ssCAGE"
            args.CAGE_model = models + model_name + "/"

            # Untar if necessary
            if not os.path.exists(args.CAGE_model):
                print("[SQANTI-SIM] Untar CAGEsim model")
                cwd = os.getcwd()
                os.chdir(models)
                sys.stdout.flush()
                res = subprocess.run(["tar", "-xzf", model_name + ".tar.gz"])
                os.chdir(cwd)
                if res.returncode != 0:
                    print("[SQANTI-SIM] ERROR: Unpacking CAGEsim pre-trained model failed", file=sys.stderr)
                    sys.exit(1)

        # Load pre-trained model
        with open(os.path.join(args.CAGE_model, "CAGE_pa_model.pkl"), "rb") as model_file:
            pa_model = pickle.load(model_file)
        model_file.close()

        true_cage_values_df = pandas.read_csv(os.path.join(args.CAGE_model,"cage_true_peaks_ecdf.tsv"), sep="\t")
        false_cage_values_df = pandas.read_csv(os.path.join(args.CAGE_model,"cage_false_peaks_ecdf.tsv"), sep="\t")

        # Simulated isoforms short-read support
        trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
        sim_trans = trans_index.loc[trans_index["sim_counts"] >= 1, "transcript_id"].tolist()
        
        iso_sr = []
        limit_info = dict(gff_type=["transcript"])
        in_handle=open(args.gtf)
        for rec in BCBio_GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
            chrom=rec.id
            trans_id=rec.features[0].qualifiers["transcript_id"][0]
            if trans_id not in sim_trans:  # Only simulated transcripts
                continue
            loc=str(rec.features[0].location)
            loc=re.split("[\(\)\[\]\:]",loc)
            loc=list(filter(None,loc))
            strand=str(loc[2])
            if strand == "+":
                TSS_coord = int(loc[0])
            else:
                TSS_coord = int(loc[1])

            iso_sr.append([trans_id, chrom, strand, TSS_coord, ratio_rep_dict[trans_id]["max_ratio_TSS"], TSS_cov_dict[trans_id]["TSS_cov"]])
        iso_sr = pandas.DataFrame(iso_sr)
        iso_sr.columns=["transcript_id", "chrom", "strand", "TSS_genomic_coord", "ratio_TSS", "TSS_cov"]
        iso_sr.to_csv(os.path.join(args.dir, "tmp_CAGE_sim.tsv"), sep="\t", header=True, index=False, na_rep="NA")
        #iso_sr= pandas.read_csv(os.path.join(args.dir, "tmp_CAGE_sim.tsv"), sep="\t")

        # 5. Simulate CAGE bed file
        print("[SQANTI-SIM] Simulating CAGE peaks")
        sim_df_bed, cage_to_trans = CAGEsim_bed(iso_sr, pa_model, true_cage_values_df, false_cage_values_df, args.falseCAGE_prop)

        sim_bed_file = pybedtools.BedTool.from_dataframe(sim_df_bed)
        sim_bed_file_sorted = sim_bed_file.sort(g=chr_order)
        sim_bed_file_sorted.saveas(os.path.join(args.dir, "CAGE_simulated.bed"))

        with open(os.path.join(args.dir, "CAGE_to_transcript.tsv"), "w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(cage_to_trans)
        file.close()
        print("[SQANTI-SIM] cage_sim.py sim done!")

main()
