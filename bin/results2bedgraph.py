# Take results from the mfed pipeline and make a bedgraph file of log fold changes

import sys

# Results file i.e. DiffBind results written out with write.table (sep="\t", quote=FALSE, row.names=FALSE)
rf = sys.argv[1]

# Print bedgraph header line
print('track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20')

with open(rf, 'r') as r:
    for x in r.readlines():
        x = x.rstrip()
        v = x.split('\t')
        # Skip header
        if v[0] == 'seqnames':
            continue
        # Print out seqname, start, end, Fold change
        print(" ".join([v[0], v[1], v[2], v[7]]))
