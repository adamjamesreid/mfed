#!/usr/bin/env python3
# Filter GATC fragments on size and read count using featureCount results, prior to running DiffBind

import sys
import argparse

parser = argparse.ArgumentParser(description='Filter featureCount results on size and read count')
parser.add_argument('-i', '--input', help="featureCounts output file")
parser.add_argument('-m', '--minsize', help="Minimum fragment size [300]")
parser.add_argument('-r', '--mincount', help="Minimum read count, in at least one sample [10]")
args = parser.parse_args()

if args.input:
    i = args.input
else:
    print ("featureCounts results required (-i)")
    exit()
if args.minsize:
    min_size = int(args.minsize)
else:
    min_size = 300
if args.mincount:
    min_read_count = int(args.mincount)
else:
    min_read_count = 10

with open(i) as file:
    # featureCounts output format
    # Geneid	Chr	Start	End	Strand	Length	DAM_R1.mLb.clN.sorted.bam	DAM_R2.mLb.clN.sorted.bam	DAM_R3.mLb.clN.sorted.bam	DAM_R4.mLb.clN.sorted.bam	Suv39_R1.mLb.clN.sorted.bam	Suv39_R2.mLb.clN.sorted.bam	Suv39_R3.mLb.clN.sorted.bam	Suv39_R4.mLb.clN.sorted.bam
    for x in file.readlines():
        x = x.rstrip()

        if x.startswith('#') or x.startswith('Geneid'):
            continue

        v = x.split('\t')
        # filter small fragments
        if int(v[5]) < min_size:
            continue
        max_count = 0
        for i in range(6,len(v)):
            #print(v[i])
            if int(v[i]) > max_count:
                max_count = int(v[i])
        if (max_count >= min_read_count):
            #print(v[1], v[2], v[3], max_count, v[5], sep='\t')
            print(v[1], v[2], v[3], sep='\t')
