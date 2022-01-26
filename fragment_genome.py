from Bio import SeqIO
import re
import sys


pattern = 'GATC'

file_path = sys.argv[1]

def get_dpnI_fragments(pattern, file_path):
    for record in SeqIO.parse(open(file_path, "r"), "fasta"):
        chrom = record.id
        previous = 0
        for match in re.finditer(pattern, str(record.seq).upper()):
            mid_pos = match.start() + 2 # +1 for the right start position, +2 for cut in middle of GATC

            if previous > 0:
                print (chrom, previous, mid_pos, sep='\t')
            previous = mid_pos

get_dpnI_fragments(pattern, file_path)
