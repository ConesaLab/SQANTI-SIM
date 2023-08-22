#!/usr/bin/env python3
"""
simulate_reads.py

Simulate data: PacBio (PBSIM3, IsoSeqSim), ONT (NanoSim),
Illumina (Polyester) and CAGE-seq

Author: Jorge Mestre Tomas (jorge.mestre.tomas@csic.es)
"""

import numpy
import os
import pandas
import pysam
import random
import re
import subprocess
import sys
from collections import defaultdict

from Bio import SeqIO
from subprocess import call


def pbsim_simulation(args):
    """Simulate PacBio reads using the PBSIM3
    Written by Tianyuan Liu
    """
    def counts_to_index(row):
        return id_counts[row["transcript_id"]]

    def maf_to_fasta_and_mapping(maf_file, fasta_file, mapping_file):
        """Convert MAF file to FASTA and mapping file."""
        with open(maf_file, 'r') as maf, open(fasta_file, 'w') as fasta, open(mapping_file, 'w') as mapping:
            id, number, sequence = None, None, None
            for line in maf:
                if line.startswith('s'):
                    parts = line.split()
                    if id is None:
                        id = parts[1]
                    else:
                        number = parts[1].split('_')[1]
                        sequence = parts[6].replace('-', '')
                        fasta.write(f'>{id}_PBSIM_simulated_read_{number}\n{sequence}\n')
                        mapping.write(f'{id}_PBSIM_simulated_read_{number}\n{id}\n')
                        id, number, sequence = None, None, None

    def create_transcript_file(index_file, fasta_file, output_file, long_count = None):
        # Load the index file into a pandas DataFrame
        index_df = pandas.read_csv(index_file, sep='\t')

        # Load the fasta file using Biopython's SeqIO
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

        # Initialize an empty list to store the transcript data
        transcripts = []

        # Iterate over the rows in the index DataFrame
        for _, row in index_df.iterrows():
            # Extract the necessary information
            transcript_id = row['transcript_id']
            if not long_count:
                sense_count = row['requested_counts']  # Assuming this is the sense count
            else:
                sense_count = int(round((row['requested_tpm']*long_count)/1000000))
 
            antisense_count = 0

            # Fetch the sequence from the fasta file
            sequence = str(fasta_dict[transcript_id].seq)

            # Check if the sequence is empty
            if not sequence:
                print(
                    f"No sequence found for transcript {transcript_id}")
            # Append a tuple with the transcript data to the list
            transcripts.append((transcript_id, sense_count, antisense_count, sequence))

        # Convert the list of transcripts to a DataFrame
        transcript_df = pandas.DataFrame(transcripts,
                                         columns=['transcript_id', 'sense_count', 'antisense_count', 'sequence'])

        # Write the DataFrame to a tab-delimited file
        transcript_df.to_csv(output_file, sep='\t', index=False)

    # convert gtf to gdp

    src_dir = os.path.dirname(os.path.realpath(__file__))
    pbsim = os.path.join(src_dir, "pbsim3")
    if os.path.isdir(args.dir):
        print("WARNING: output directory already exists. Overwriting!", file=sys.stderr)
    else:
        os.makedirs(args.dir)
    # temporary directory
    tmp_dir = os.path.join(args.dir, "temp_pbsim")
    if os.path.isdir(tmp_dir):
        print("WARNING: temporary directory already exists. Overwriting!", file=sys.stderr)
    else:
        os.makedirs(tmp_dir)

    udir = src_dir + "/IsoSeqSim/" + "utilities/"
    cmd_gtf2gdp = udir + "/py_isoseqsim_gtf2gpd.py -i " + args.gtf + " -o " + tmp_dir + "/normal_annotation.gpd"
    call(cmd_gtf2gdp.split())
    sys.stdout.flush()

    cmd_gpd2fa = udir + "/py_isoseqsim_gpd2fa_normal.py -a " + args.genome + " -g " + tmp_dir + "/normal_annotation.gpd" + " -o " + tmp_dir + "/normal_transcriptome.fa"
    call(cmd_gpd2fa.split())

    # Create transcript file for PBSIM3
    if not args.long_count:
        create_transcript_file(index_file = args.trans_index, fasta_file= tmp_dir + "/normal_transcriptome.fa", output_file=os.path.join(tmp_dir, "sample.transcript"))
    else:
        create_transcript_file(index_file = args.trans_index, fasta_file= tmp_dir + "/normal_transcriptome.fa", output_file=os.path.join(tmp_dir, "sample.transcript"), long_count=args.long_count)

    # PBSIM3 simulation
    print("[SQANTI-SIM] Simulating PacBio reads with PBSIM3")

    if args.pass_num == 1:
        cmd = [
            "pbsim",
            "--strategy", "trans",
            "--method", "qshmm",
            "--qshmm", pbsim + "/data/QSHMM-RSII.model",
            "--transcript", tmp_dir + "/sample.transcript",
            "--accuracy-mean", "0.95",
            "--seed", str(args.seed)
        ]

        cmd = " ".join(cmd)
        sys.stdout.flush()
        if subprocess.check_call(cmd, shell=True) != 0:
            print("ERROR running PBSIM3: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)

        # Move output files to tmp directory
        os.rename("sd.fastq", os.path.join(tmp_dir, "sd.fastq"))
        os.rename("sd.maf", os.path.join(tmp_dir, "sd.maf"))
        maf_file = os.path.join(tmp_dir, "sd.maf")
        pbsm_fasta = os.path.join(args.dir, "PBSIM_simulated.fasta")
        read_to_iso = os.path.join(args.dir, "PBSIM_simulated.read_to_isoform.tsv")

        # Remove unaligned reads from the MAF file
        print("[SQANTI-SIM] Creating read_to_isoform file")
        # Define a function to convert the MAF file to a FASTA file
        maf_to_fasta_and_mapping(maf_file, pbsm_fasta, read_to_iso)

    elif args.pass_num > 1:
        cmd = [
            "pbsim",
            "--strategy", "trans",
            "--method", "qshmm",
            "--qshmm", pbsim + "/data/QSHMM-RSII.model",
            "--transcript", tmp_dir + "/sample.transcript",
            "--accuracy-mean", "0.95",
            "--pass-num", str(args.pass_num),
            "--seed", str(args.seed)
        ]

        cmd = " ".join(cmd)
        sys.stdout.flush()
        if subprocess.check_call(cmd, shell=True) != 0:
            print("ERROR running PBSIM3: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)

        # Move output files to tmp directory
        os.rename("sd.sam", os.path.join(tmp_dir, "sd.sam"))
        os.rename("sd.maf", os.path.join(tmp_dir, "sd.maf"))
        # convert sam file to bam file
        cmd = "samtools view -bS " + tmp_dir + "/sd.sam > " + tmp_dir + "/sd.bam"
        if subprocess.check_call(cmd, shell=True) != 0:
            print("ERROR running samtools: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)

        # use the ccs program to generate consensus reads (HiFi reads)
        cmd = "ccs " + tmp_dir + "/sd.bam " + tmp_dir + "/sd.ccs.bam"
        if subprocess.check_call(cmd, shell=True) != 0:
            print("ERROR running ccs: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)

        # use samtools convert bam to fasta
        cmd = "samtools fasta " + tmp_dir + "/sd.ccs.bam > " + tmp_dir + "/sd.ccs.fasta"
        if subprocess.check_call(cmd, shell=True) != 0:
            print("ERROR running samtools: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)

        def rename_fasta_headers(transcript_table_path, fasta_path, output_path):
            # Parse the transcript table
            transcript_table = {}
            total_sense_count = 0
            with open(transcript_table_path, 'r') as f:
                next(f)  # Skip header line
                for line in f:
                    transcript_id, sense_count, antisense_count, _ = line.split('\t')
                    sense_count = int(sense_count)
                    transcript_table[total_sense_count + 1] = (transcript_id, sense_count)
                    total_sense_count += sense_count

            # Rename the headers in the FASTA file
            with open(fasta_path, 'r') as in_f, open(output_path, 'w') as out_f:
                for line in in_f:
                    if line.startswith('>'):
                        read_num = int(line.split('/')[1])
                        for start, (transcript_id, sense_count) in transcript_table.items():
                            if start <= read_num < start + sense_count:
                                read_id = read_num - start + 1
                                out_f.write(f'>{transcript_id}_PBSIM_simulated_read_{read_id}\n')
                                break
                    else:
                        out_f.write(line)

        pbsm_fasta = os.path.join(args.dir, "PBSIM3_simulated.fasta")
        read_to_iso = os.path.join(args.dir, "PBSIM3_simulated.read_to_isoform.tsv")
        rename_fasta_headers(tmp_dir + "/sample.transcript", tmp_dir + "/sd.ccs.fasta", pbsm_fasta)


    else:
        print("ERROR: pass_num must be >= 1", file=sys.stderr)
        sys.exit(1)

    print("[SQANTI-SIM] Counting PBSIM3 reads")

    output_read_info = open(read_to_iso, "w")
    id_counts = defaultdict(lambda: 0)
    with open (pbsm_fasta, "r") as sim_fasta:
        for line in sim_fasta:
            if line.startswith(">"):
                line = line.lstrip(">").rstrip()
                line_split = line.split("_")
                trans_id = "_".join(line_split[:-4])
                output_read_info.write(line + "\t" + trans_id + "\n")
                id_counts[trans_id] += 1

    sim_fasta.close()
    output_read_info.close()

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    trans_index["sim_counts"] = trans_index.apply(counts_to_index, axis=1)
    trans_index["sim_counts"] = trans_index["sim_counts"].fillna(0)
    trans_index.to_csv(
        args.trans_index, sep="\t", header=True, index=False, na_rep="NA"
    )
    print("[SQANTI-SIM] PBSIM3 simulation done")
    return

def isoseqsim_simulation(args):
    """Simulate PacBio reads using the IsoSeqSim pipeline"""

    def counts_to_index(row):
        return id_counts[row["transcript_id"]]
    
    # Generate IsoSeqSim template expression file
    expr_f = os.path.join(
        os.path.dirname(os.path.abspath(args.trans_index)),
        "tmp_expression.tsv",
    )
    index_file_requested_counts = 0
    f_out = open(expr_f, "w")
    f_out.write("target_id\test_counts\ttpm\n")
    with open(args.trans_index, "r") as idx:
        header_names = idx.readline()
        header_names = header_names.split()
        i = header_names.index("requested_counts")
        j = header_names.index("requested_tpm")
        for line in idx:
            line = line.split()
            if int(line[i]) == 0:
                continue
            f_out.write(line[0] + "\t" + line[i] + "\t" + line[j] + "\n")
            index_file_requested_counts += int(line[i])
    idx.close()
    f_out.close()

    if not args.long_count:
        args.long_count = index_file_requested_counts

    if os.path.isdir(args.dir):
        print("WARNING: output direcory already exists. Overwritting!", file=sys.stderr)
    else:
        os.makedirs(args.dir)

    # PacBio Sequel simulation -> error rates from IsoSeqSim GitHub
    print("[SQANTI-SIM] Simulating PacBio reads with IsoSeqSim")
    src_dir = os.path.dirname(os.path.realpath(__file__))
    isoseqsim = os.path.join(src_dir, "IsoSeqSim/bin/isoseqsim")
    util_dir = os.path.join(src_dir, "IsoSeqSim/utilities/")
    cmd = [
        isoseqsim,
        "-g",
        str(args.genome),
        "-a",
        str(args.gtf),
        "--expr",
        str(expr_f),
        "--c5",
        os.path.join(util_dir, "5_end_completeness.PacBio-Sequel.tab"),
        "--c3",
        os.path.join(util_dir, "3_end_completeness.PacBio-Sequel.tab"),
        "-o",
        os.path.join(args.dir, "IsoSeqSim_simulated"),
        "-t",
        os.path.join(args.dir, "IsoSeqSim_simulated.tsv"),
        "--es",
        "0.01731",
        "--ed",
        "0.01090",
        "--ei",
        "0.02204",
        "-n",
        str(args.long_count),
        "-m",
        "normal",
        "--cpu",
        str(args.cores),
        "--tempdir",
        os.path.join(args.dir, "temp_isoseqsim"),
        "--seed",
        str(args.seed),
    ]

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("ERROR running IsoSeqSim: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)
    os.remove(expr_f)

    print("[SQANTI-SIM] Counting PacBio reads")
    read_to_iso = os.path.join(args.dir, "IsoSeqSim_simulated.read_to_isoform.tsv")
    output_read_info = open(read_to_iso, "w")
    id_counts = defaultdict(lambda: 0)
    isoseqsim_fasta = os.path.join(args.dir, "IsoSeqSim_simulated.fasta")
    with open(isoseqsim_fasta, "r") as sim_fasta:
        for line in sim_fasta:
            if line.startswith(">"):
                line = line.lstrip(">")
                line_split = line.split("_")
                trans_id = "_".join(line_split[:-4])
                output_read_info.write(line + "\t" + trans_id + "\n")
                id_counts[trans_id] += 1
    sim_fasta.close()
    output_read_info.close()

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    trans_index["sim_counts"] = trans_index.apply(counts_to_index, axis=1)
    trans_index["sim_counts"] = trans_index["sim_counts"].fillna(0)
    trans_index.to_csv(
        args.trans_index, sep="\t", header=True, index=False, na_rep="NA"
    )

    print("[SQANTI-SIM] IsoSeqSim simulation done")
    return


def ont_simulation(args):
    """Simulate ONT reads using the NanoSim pipeline"""

    def counts_to_index(row):
        return id_counts[row["transcript_id"]]

    # Generate NanoSim template expression file
    expr_f = os.path.join(os.path.dirname(os.path.abspath(args.trans_index)), "tmp_expression.tsv")
    index_file_requested_counts = 0
    f_out = open(expr_f, "w")
    f_out.write("target_id\test_counts\ttpm\n")
    with open(args.trans_index, "r") as idx:
        header_names = idx.readline()
        header_names = header_names.split()
        i = header_names.index("requested_counts")
        j = header_names.index("requested_tpm")
        for line in idx:
            line = line.split()
            if int(line[i]) == 0:
                continue
            f_out.write(line[0] + "\t" + line[i] + "\t" + line[j] + "\n")
            index_file_requested_counts += int(line[i])
    idx.close()
    f_out.close()

    if not args.long_count:
        args.long_count = index_file_requested_counts

    if os.path.isdir(args.dir):
        print("WARNING: output direcory already exists. Overwritting!")
    else:
        os.makedirs(args.dir)

    src_dir = os.path.dirname(os.path.realpath(__file__))
    nanosim = os.path.join(src_dir, "NanoSim/src/simulator.py")

    if not args.nanosim_model:
        models = os.path.join(src_dir, "../pre-trained_models/")
        if args.read_type == "dRNA":
            model_name = "human_WTC11_dRNA_guppy_NanoSim"
        elif args.read_type == "cDNA":
            model_name = "human_WTC11_cDNA_guppy_NanoSim"
        prefix = "training"
        model_dir = models + model_name + "/"
        model_complete = models + model_name + "/" + prefix
    else:
        model_complete = args.nanosim_model
        model_dir = os.path.dirname(model_complete)
        prefix = os.path.basename(model_complete)
        models = os.path.dirname(os.path.normpath(model_dir))
        model_name = os.path.basename(os.path.normpath(model_dir))
    
    if not os.path.exists(model_dir):
        print("[SQANTI-SIM] Untar NanoSim model")
        cwd = os.getcwd()
        os.chdir(models)
        sys.stdout.flush()
        if args.read_type == "cDNA":
            if os.path.isfile(model_name + ".tar.gz.partaa") and os.path.isfile(model_name + ".tar.gz.partab"):
                os.system(" ".join(["cat", model_name + ".tar.gz.part*", ">", model_name + ".tar.gz"]))
            else:
                print("[SQANTI-SIM] ERROR: NanoSim cDNA pre-trained model is missing. Please check if %s exists." % (model_name + ".tar.gz.part*"), file=sys.stderr)
                sys.exit(1)
        res = subprocess.run(["tar", "-xzf", model_name + ".tar.gz"])
        os.chdir(cwd)
        if res.returncode != 0:
            print("[SQANTI-SIM] ERROR: Unpacking NanoSim pre-trained model failed", file=sys.stderr)
            sys.exit(1)

    if args.read_type == "dRNA":
        r_type = "dRNA"
        uracil = True
    elif args.read_type == "cDNA":
        r_type = "cDNA_1D2"
        uracil = False
    else:
        print("[SQANTI-SIM] ERROR not valid read_type value %s" % (args.read_type), file=sys.stderr)
        return

    # Extract fasta transcripts
    print("[SQANTI-SIM] Extracting transcript sequences")
    ref_t = os.path.join(os.path.dirname(args.genome), "sqanti-sim.transcripts.fa")
    if os.path.exists(ref_t):
        print("[SQANTI-SIM] WARNING: %s already exists, it will be overwritten" %(ref_t))

    cmd = ["gffread", "-w", str(ref_t), "-g", str(args.genome), str(args.gtf)]
    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running gffread: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)

    print("[SQANTI-SIM] Simulating ONT reads with NanoSim")
    cmd = [
        nanosim,
        "transcriptome",
        "-rt",
        str(ref_t),
        "-rg",
        str(args.genome),
        "-e",
        str(expr_f),
        "-c",
        str(model_complete),
        "-o",
        os.path.join(args.dir, "ONT_simulated"),
        "-n",
        str(args.long_count),
        "-r",
        r_type,
        "-b",
        "guppy",
        "-t",
        str(args.cores),
        "--seed",
        str(args.seed),
        "--fastq",
        "--no_model_ir",
    ]

    if uracil:
        cmd.append("--uracil")

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("ERROR running NanoSim: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)
    os.remove(expr_f)

    print("[SQANTI-SIM] Counting and renaming ONT reads")
    fastqs = [
        os.path.join(args.dir, "ONT_simulated_aligned_reads.fastq"),
        os.path.join(args.dir, "ONT_simulated_unaligned_reads.fastq"),
    ]

    n_read = 0
    pair_id = []
    id_counts = defaultdict(lambda: 0)
    f_name = os.path.join(args.dir, "ONT_simulated.fastq")
    f_out = open(f_name, "w")

    for f in fastqs:
        f_in = open(f, "r")
        for line in f_in:
            if line.startswith("@"):
                line = line.lstrip("@")
                trans_id = "_".join(line.split("_")[:-7])
                id_counts[trans_id] += 1
                read_id = trans_id + "_ONT_simulated_read_" + str(n_read)
                n_read += 1
                pair_id.append((read_id, trans_id))

                f_out.write("@{}\n".format(read_id))

            else:
                f_out.write(line)
    f_in.close()
    f_out.close()

    # Saving counts and read-to-isoform files
    f_name = os.path.join(args.dir, "ONT_simulated.read_to_isoform.tsv")
    f_out = open(f_name, "w")
    for pair in pair_id:
        f_out.write(str(pair[0]) + "\t" + str(pair[1]) + "\n")
    f_out.close()

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    trans_index["sim_counts"] = trans_index.apply(counts_to_index, axis=1)
    trans_index["sim_counts"] = trans_index["sim_counts"].fillna(0)
    trans_index.to_csv(args.trans_index, sep="\t", header=True, index=False, na_rep="NA")

    print("[SQANTI-SIM] NanoSim simulation done")
    return


def illumina_simulation(args):
    """Simulate Illumina reads using the Polyester pipeline"""

    def counts_to_index(row):
        return id_counts[row["transcript_id"]]

    # Extract fasta transcripts
    print("[SQANTI-SIM] Extracting transcript sequences")
    ref_t = os.path.join(os.path.dirname(args.genome), "sqanti-sim.transcripts.fa")
    if os.path.exists(ref_t):
        print("[SQANTI-SIM] WARNING: %s already exists, it will be overwritten" %(ref_t))

    cmd = ["gffread", "-w", str(ref_t), "-g", str(args.genome), str(args.gtf)]
    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running gffread: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)

    # Generate Polyester template expression file
    count_d = defaultdict(float)
    n = 0
    with open(args.trans_index, "r") as idx:
        header_names = idx.readline()
        header_names = header_names.split()
        i = header_names.index("requested_counts")
        j = header_names.index("requested_tpm")
        for line in idx:
            line = line.split()
            if int(line[i]) == 0:
                continue
            count_d[line[0]] = float(line[j])
            n += int(line[i])
    idx.close()

    if not args.short_count:
        args.short_count = n

    for k in count_d:
        count_d[k] = round((count_d[k] * args.short_count) / 1000000)

    expr_f = os.path.join(args.dir, "tmp_expression.tsv")
    f_out = open(expr_f, "w")
    with open(ref_t, "r") as rt:
        for line in rt:
            if line.startswith(">"):
                line = line[1:].split()
                line = line[0]
                f_out.write(str(count_d[line]) + "\n")
    rt.close()
    f_out.close()

    print("[SQANTI-SIM] Simulating Illumina reads with Polyester")
    src_dir = os.path.dirname(os.path.realpath(__file__))

    cmd = [
        "Rscript",
        os.path.join(src_dir, "polyester_sim.R"),
        ref_t,
        expr_f,
        args.dir,
        str(args.seed),
    ]

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print(
            "ERROR simulatin with Polyester: {0}".format(cmd), file=sys.stderr
        )
        sys.exit(1)
    os.remove(expr_f)

    print("[SQANTI-SIM] Counting Illumina reads")
    os.rename(
        os.path.join(args.dir, "sample_01_1.fasta"),
        os.path.join(args.dir, "Illumina_simulated_1.fasta"),
    )
    os.rename(
        os.path.join(args.dir, "sample_01_2.fasta"),
        os.path.join(args.dir, "Illumina_simulated_2.fasta"),
    )

    id_counts = defaultdict(lambda: 0)
    with open(
        os.path.join(args.dir, "Illumina_simulated_1.fasta"), "r"
    ) as illumina_sim:
        for line in illumina_sim:
            if line.startswith(">"):
                line = line[1:].strip()
                line = re.split("/|;|\s+|\n",line)
                line = line[1]
                id_counts[line] += 1
    illumina_sim.close()

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    trans_index["illumina_counts"] = trans_index.apply(counts_to_index, axis=1)
    trans_index["illumina_counts"] = trans_index["illumina_counts"].fillna(0)
    trans_index.to_csv(args.trans_index, sep="\t", header=True, index=False, na_rep="NA")

    print("[SQANTI-SIM] Polyester simulation done")
    return

def CAGE_simulation(args):
    """Simulate a sample-specific CAGE peak BED file"""

    print("[SQANTI-SIM] Simulating CAGE peaks with cage_sim.py")
    src_dir = os.path.dirname(os.path.realpath(__file__))
    cagesim = os.path.join(src_dir, "cage_sim.py")

    sr1=os.path.join(args.dir, "Illumina_simulated_1.fasta")
    sr2=os.path.join(args.dir, "Illumina_simulated_2.fasta")
    srfofn=os.path.join(args.dir, "Illumina_simulated.fofn")
    with open(srfofn, "w") as fofn:
        fofn.write(sr1 + " " + sr2)
    fofn.close()

    cmd = [
        cagesim,
        "sim",
        "--trans_index",
        str(args.trans_index),
        "--gtf",
        str(args.gtf),
        "--genome",
        str(args.genome),
        "--short_reads",
        str(srfofn),
        "--falseCAGE_prop",
        str(args.falseCAGE_prop),
        "--dir",
        str(args.dir),
        "-t",
        str(args.cores)
    ]

    if args.CAGE_model:
        cmd.extend(["--CAGE_model", str(args.CAGE_model)])
    if args.seed:
        cmd.extend(["--seed", str(args.seed)])


    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("ERROR running cage_sim.py: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)

    return
