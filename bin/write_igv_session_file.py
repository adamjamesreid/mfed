#!/usr/bin/env python3

# Generate an IGV session file for viewing data from the mfed pipeline

# n.b. IGV will not read gtf files with 'ensGene' or 'refGene' in the name, so these need to be changed somehow (https://github.com/igvteam/igv/issues/588)

import argparse
import sys

# Name of genome/fa
# Bigwig files
# GTF annotation
# foldchange bedgraph
# enriched_fragments.bed

# Gather filenames from command line
parser = argparse.ArgumentParser(description='Generate an IGV XML seesion file of data from the Mfed pipeline')
parser.add_argument('-f', '--fasta', help="genome fasta", required=True)
parser.add_argument('-g', '--gtf', help="GTF of annotation", required=True)
parser.add_argument('-fc', '--foldchange', help="Bedgraph of fold changes", required=True)
parser.add_argument('-ef', '--enrichfrag', help="BED file of enriched fragments", required=True)
parser.add_argument('bw', nargs='+', help="BigWig files of mapped reads")
args = parser.parse_args()

# Print out IGV XML session file
print("""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{}" hasGeneTrack="true" hasSequenceTrack="true">
    <Resources>
        <Resource name="log fold change" path="{}" type="bedgraph"/>""".format(args.fasta, args.foldchange))

for x in args.bw:
    print("        <Resource name=\"Coverage {}\" path=\"{}\" type=\"bw\"/>".format(x, x))

print("""        <Resource name="Enriched fragments" path=\"{}\" type=\"bed\"/>        
        <Resource name="Genes" path="{}"/>
    </Resources>
</Session>
""".format(args.enrichfrag, args.gtf))

"""
        <Resource path="damonly_R2.mLb.clN.sorted.bam" type="bam"/>
        <Resource path="significant_fragments.bed" type="bed"/>
        <Resource path="hp1fusion_R1.mLb.clN.sorted.bam" type="bam"/>
        <Resource path="hp1fusion_R1.mLb.clN.sorted.deeptools.bw" type="bw"/>
        <Resource path="hp1fusion_R1.bigWig" type="bigwig"/>
        <Resource path="hp1fusion_R2.mLb.clN.sorted.bam" type="bam"/>
        <Resource path="enriched_fragments.bed" type="bed"/>
        <Resource path="damonly_R1.mLb.clN.sorted.bam" type="bam"/>
        <Resource path="hp1fusion_R1.mLb.clN.sorted.bw" type="bw"/>
        <Resource path="foldchange.bedgraph" type="bedgraph"/>
        <Resource name="Refseq Genes" path="https://s3.amazonaws.com/igv.org.genomes/dm6/ncbiRefSeq.txt.gz" type="refgene"/>
"""
