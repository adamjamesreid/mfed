#!/usr/bin/env python3

# Convert nf-core/chipseq design file into a DiffBind samplesheet suitable for mfed.nf

import sys
import argparse

parser = argparse.ArgumentParser(description='Convert nf-core/chipseq design file into a DiffBind samplesheet suitable for mfed.nf')
parser.add_argument('-i', '--input', help="Design file", required=True)
parser.add_argument('-d', '--bamdir', help="BAM file directory", required=True)
#parser.add_argument('-O', '--outfmt', help="Output format - bed or gtf [gtf]")
#parser.add_argument('file_path', help="Fasta file of reference genome sequence")
args = parser.parse_args()

# Name of fragment bed file generated by mfed.nf
filtered_fragment_file = "filtered_fragments.bed"

# Print header for DiffBind samplesheet
print("SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller")

# Open design file
with open(args.input) as i:
    for x in i.readlines():
        x = x.rstrip()
        # Skip design file header
        if x.startswith('group,'):
            continue
        # group,replicate,fastq_1,fastq_2,antibody,control
        group,replicate,fastq_1,fastq_2,antibody,control = x.split(',')
        # Print out correct format for a DiffBind samplesheet
        print("{}_{},TEST,TEST,TEST,{},{},{}/{}_R{}.mLb.clN.sorted.bam,,,{},bed".format(group, replicate, group, replicate, args.bamdir, group, replicate, filtered_fragment_file))