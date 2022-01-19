library(DiffBind)
library(parallel)

## ARGS
# 1. samplesheet
# 2. Control condition (as defined in the samplesheet)
# 3. Treatment condition (as defined in the samplesheet)
# 4. Fold change cutoff
# 5. FDR cutoff

args <- commandArgs(TRUE)

ss <- args[1]
control <- args[2]
treatment <- args[3]
fc_cut <- args[4]
fdr_cut <- args[5]

db <- dba(sampleSheet=ss)

# minOverlap is the proportion of samples in which the peak must be found to be included (if 1 or more confusingly it becomes the number of samples)
# fragmentSize is used as the length of the reads and they will be extended, size from BAM file used if set to 0
# summits - TRUE means that peaksets are unaffected by use of summit calculations
db <- dba.count(db, minOverlap=0.5, fragmentSize = 0, summits = TRUE)

info <- dba.show(db)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

pdf("plot.pdf")
plot(db)
dev.off()

db <- dba.normalize(db)

norm <- dba.normalize(db, bRetrieve=TRUE)
norm

normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs

db <- dba.contrast(db, contrast=c("Treatment", treatment, control))
db

db <- dba.analyze(db)
dba.show(db, bContrasts=TRUE)

#dev.off()
pdf("sample_heatmap_post_contrast.pdf")
plot(db, contrast=1)
dev.off()

db.DB <- dba.report(db)
db.DB

hist(db.DB$Fold, breaks=50)

# Filter for significance (fold change and FDR)
db.DB.conf <- subset(db.DB, ((db.DB$Fold >= fc_cut | db.DB$Fold <= fc_cut) & db.DB$FDR <= fdr_cut))
write.table(db.DB.conf, file="results.txt")

# Write out significant peaks 
db.DB.conf.df <- as.data.frame(db.DB.conf)  
out_table <- db.DB.conf.df[, c("seqnames", "start", "end")]
write.table(out_table, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, file="significant_fragments.bed")

# Write out consensus peak set
dba.peakset(db, peak.format="bed", bRetrieve=TRUE, writeFile="filtered_fragments_diffbind.bed")


