#!/usr/bin/env python3
# Generate a BED or GTF file of in silico GATC fragments for the mfed DAMid analysis pipeline

#from Bio import SeqIO
import re
import sys
import argparse

pattern = 'GATC'

parser = argparse.ArgumentParser(description='Generate a BED file of in silico GATC fragments for the mfed DAMid analysis pipeline')
parser.add_argument('-p', '--pattern', help="Patern at which to cut [GATC]")
parser.add_argument('-O', '--outfmt', help="Output format - bed or gtf [gtf]")
parser.add_argument('file_path', help="Fasta file of reference genome sequence")
args = parser.parse_args()

if args.pattern:
    pattern = args.pattern
else:
    pattern = 'GATC'

if args.file_path:
    file_path = args.file_path

if args.outfmt == 'bed':
    outfmt = 'bed'
elif args.outfmt == 'gtf' or not args.outfmt:
    outfmt = 'gtf'
else:
    sys.stderr.write("Output format '{}' not recognised, try 'bed' or 'gtf'".format(args.outfmt))
    exit()

def read_fasta (fa_file):
    '''Generator to parse fasta file and return header, sequence tuples'''
    seq = ''
    seqid = ''
    header = ''
    # Open fasta file
    with open(fa_file, 'r') as f:
        # Read in a line
        for x in f.readlines():
            x = x.rstrip()
            # If it is a header line, check whether we have already stored a sequence (if so yield it)
            if x.startswith('>'):
                if header != '':
                    yield((seqid[0], seq))
                # Store this header
                header = x.replace('>', '')
                seqid = header.split(' ')
                # reset sequence
                seq = ''
            # Store sequence unless it is an empty line
            elif x != '':
                seq = seq + x
    # Yield final sequence
    yield((seqid[0], seq))

def get_dpnI_fragments(pattern, file_path):
    ''' Generator to get in silico fragments of a genome sequence'''
    for (header, seq) in read_fasta(file_path):
        previous = 0
        for match in re.finditer(pattern, str(seq).upper()):
            mid_pos = match.start() + 2 # +1 for the right start position, +2 for cut in middle of GATC

            if previous > 0:
                yield(header, previous, mid_pos)
            previous = mid_pos
        yield(header, previous, mid_pos)

for (c, p, m) in get_dpnI_fragments(pattern, file_path):
    if outfmt == 'bed':
        print ('\t'.join([c, str(p), str(m)]))
    elif outfmt == 'gtf':
        frag_id = 'frag_id="{}:{}_{}"'.format(c, p, m)
        print ('\t'.join([c, pattern+'_frag', 'region', str(p), str(m), '.', '.', '.', frag_id]))
        

#cat gatc_frags.bed | perl -ne 'chomp;@a=split/\t/;print "$a[0]\tGATC_frag\tregion\t$a[1]\t$a[2]\t.\t.\t.\tfrag_id=\"$a[0]\:$a[1]\_$a[2]\"\n"' > gatc_frags.gtf
