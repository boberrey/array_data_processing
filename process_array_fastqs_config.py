# Configuration file 

import os
import sys
import glob

# Input data

FASTQ_DIR = 'input/fastqs'
METADATA_FILE = 'input/metadata.txt'

# Resources
GENOME_BUILD = "HB27"
REFERENCE_FILE = '/share/PI/wjg/lab/genomes/HB27/HB27'
BLACKLIST = None


# EXE_DIR is required for certain software (e.g. fastqc, etc.) How many of these do we want to keep anyway?
EXE_DIR = '/share/PI/wjg/lab/bin'
PICARD_JAR = '/share/PI/wjg/lab/bin/picard.jar'
ARRAY_TOOLS = '/home/users/boberrey/git_clones/array_data_processing'



# metadata file

def make_meta(filename):
    """
    Generate a metadata file with the names of all the samples to be processed.
    Sample names are inferred from fastq files.
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    r1_files = list(map(os.path.abspath, glob.glob(os.path.join(FASTQ_DIR,"*_R1*.fastq*"))))
    if (len(r1_files) < 1):
        sys.exit("No fastqs with _R1 found.")
    r2_files = [os.path.join(os.path.dirname(r1_file), 
        os.path.basename(r1_file).replace('R1', 'R2')) for r1_file in r1_files]
    if  all([os.path.isfile(r2_file) for r2_file in r2_files]) is False:
        sys.exit("Not all matching _R2 files found.")
    sample_labels = [os.path.basename(r1_file).split("_R1")[0] for r1_file in r1_files]
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name","Read1","Read2"]) + "\n")
        for sample_label, r1_file, r2_file in zip(sample_labels, r1_files, r2_files):
            if len(sample_label) > 30:
                sample_label = sample_label[:20] + "..." + sample_label[-10:]
            outfile.write("\t".join([sample_label, r1_file, r2_file]) + "\n")
