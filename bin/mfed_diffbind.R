# Check where we are loading libraries from
.libPaths()

library(DiffBind)
library(parallel)
library(ChIPseeker)
library(GenomicFeatures)

## ARGS
# 1. samplesheet
# 2. Control condition (as defined in the samplesheet)
# 3. Treatment condition (as defined in the samplesheet)
# 4. Fold change cutoff
# 5. FDR cutoff
# 6. GTF file

args <- commandArgs(TRUE)

ss <- args[1]
control <- args[2]
treatment <- args[3]
fc_cut <- args[4]
fdr_cut <- args[5]
gtf_file <- args[6]
annotation_level <- args[7]
tss_region_start <- args[8]
tss_region_end <- args[9]

# Comma-separated file listing priority of annotations, default used if not supplied
#e.g. "Exon,Intron,5UTR,3UTR,Promoter,Downstream,Intergenic"
genomicAnnotationPriorityfile <- args[10] # Comma-separated file listing priority of annotations, default used if not supplied

ss
control
treatment
fc_cut
fdr_cut
gtf_file
annotation_level
tss_region_start
tss_region_end


# Read in samplesheet
db <- dba(sampleSheet=ss)

# minOverlap is the proportion of samples in which the peak must be found to be included (if 1 or more confusingly it becomes the number of samples)
# fragmentSize is used as the length of the reads and they will be extended, size from BAM file used if set to 0
# summits - TRUE means that peaksets are unaffected by use of summit calculations
db <- dba.count(db, minOverlap=0.5, fragmentSize = 0, summits = TRUE)

info <- dba.show(db)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

# Plot a heatmap of samples
pdf("sample_heatmap_plot.pdf")
dba.plotHeatmap(db)
dev.off()

# Plot a PCA plot of samples
pdf("sample_pca_plot.pdf")
dba.plotPCA(db)
dev.off()

# Normalise data
db <- dba.normalize(db)

norm <- dba.normalize(db, bRetrieve=TRUE)
norm

normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs

# Compare fusion and dam-only
db <- dba.contrast(db, contrast=c("Treatment", treatment, control))
db

db <- dba.analyze(db, bBlacklist=FALSE, bGreylist=FALSE)
dba.show(db, bContrasts=TRUE)

# Output an MA plot
pdf("MA_plot.pdf")
dba.plotMA(db, th=as.double(fdr_cut))
dev.off()

#Output a volcano plot
pdf("volcano_plot.pdf")
dba.plotVolcano(db, th=as.double(fdr_cut))
dev.off()

# Generate report
db.DB <- dba.report(db)
db.DB

# Filter for significance (fold change and FDR)
db.DB.conf <- subset(db.DB, ((db.DB$Fold >= as.double(fc_cut) | db.DB$Fold <= as.double(fc_cut)) & db.DB$FDR <= as.double(fdr_cut)))

# Write out all significant results - enriched in fusion or control
write.table(db.DB.conf, file="results_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Write out peaks significantly enriched in fusion
db.DB.conf.enrich.df <- as.data.frame(subset(db.DB.conf, db.DB.conf$Fold > 0))
out_table <- db.DB.conf.enrich.df[, c("seqnames", "start", "end")]
write.table(out_table, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, file="enriched_fragments.bed")

# Write out results significantly enriched in fusion
write.table(subset(db.DB.conf, db.DB.conf$Fold > 0), file="results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Add gene annotation for fragments
# Make a custom TxDB object
Dm_custom_TXDb <- makeTxDbFromGFF(gtf_file, format="gtf")

#Read in annotation priorities
if(file.exists(genomicAnnotationPriorityfile)) {
  apf_df = read.csv(genomicAnnotationPriorityfile, header=FALSE)

  # Annotate fragments
  peakAnno.custom.edb <- annotatePeak(subset(db.DB.conf, db.DB.conf$Fold > 0), tssRegion=c(as.numeric(tss_region_start), as.numeric(tss_region_end)),
                                      TxDb=Dm_custom_TXDb,
                                      genomicAnnotationPriority = as.character(apf_df[1,]),
                                      level=annotation_level)
} else {
  # Annotate fragments
  peakAnno.custom.edb <- annotatePeak(subset(db.DB.conf, db.DB.conf$Fold > 0), tssRegion=c(as.numeric(tss_region_start), as.numeric(tss_region_end)),
                                      TxDb=Dm_custom_TXDb,
                                      level=annotation_level)

}

# Outputs summary of annotation
peakAnno.custom.edb

peakAnno.custom.df <- as.data.frame(peakAnno.custom.edb)
write.table(peakAnno.custom.edb, file="results_annotated.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Save an image so e.g. figures can be reproduced as needed
save.image(file = "mfed.RData")
