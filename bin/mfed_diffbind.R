# Check where we are loading libraries from
.libPaths()

library(DiffBind)
library(parallel)

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

ss
control
treatment
fc_cut
fdr_cut
gtf_file
annotation_level

# Comma-separated file listing priority of annotations, default used if not supplied
#e.g. "Exon,Intron,5UTR,3UTR,Promoter,Downstream,Intergenic"
#genomicAnnotationPriorityfile <- args[8] # Comma-separated file listing priority of annotations, default used if not supplied
genomicAnnotationPriorityfile  <- 'mfed/annotation_priority.csv'

tssRegion = c(-100, 100)

db <- dba(sampleSheet=ss)

# minOverlap is the proportion of samples in which the peak must be found to be included (if 1 or more confusingly it becomes the number of samples)
# fragmentSize is used as the length of the reads and they will be extended, size from BAM file used if set to 0
# summits - TRUE means that peaksets are unaffected by use of summit calculations
db <- dba.count(db, minOverlap=0.5, fragmentSize = 0, summits = TRUE)

info <- dba.show(db)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

pdf("sample_heatmap_plot.pdf")
dba.plotHeatmap(db)
dev.off()

pdf("sample_pca_plot.pdf")
dba.plotPCA(db)
dev.off()

db <- dba.normalize(db)

norm <- dba.normalize(db, bRetrieve=TRUE)
norm

normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs

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

#dev.off()
#pdf("sample_heatmap_post_contrast.pdf")
#plot(db, contrast=1)
#dev.off()

db.DB <- dba.report(db)
db.DB

hist(db.DB$Fold, breaks=50)

# Write out all results
# write.table(db.DB, file="results_all.txt", sep="\t", quote=FALSE, row.names=FALSE)

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

# Testing to add gene annotation for fragments
#library(TxDb.Dmelanogaster.UCSC.dm6.ensGene) # Which mitochondrial sequence does this annotation use?
library(ChIPseeker)

# Try making a custom TxDB object
library(GenomicFeatures)
Dm_custom_TXDb <- makeTxDbFromGFF(gtf_file, format="gtf")

#Read in annotation priorities
if(file.exists(genomicAnnotationPriorityfile)) {
  apf_df = read.csv(genomicAnnotationPriorityfile, header=FALSE)

  # Annotate fragments
  peakAnno.custom.edb <- annotatePeak(subset(db.DB.conf, db.DB.conf$Fold > 0), tssRegion=tssRegion,
                                      TxDb=Dm_custom_TXDb, #,annoDb="org.Hs.eg.db",
                                      genomicAnnotationPriority = as.character(apf_df[1,]),
                                      level=annotation_level)
} else {
  # Annotate fragments
  peakAnno.custom.edb <- annotatePeak(db.DB.conf, tssRegion=tssRegion,
                                      TxDb=Dm_custom_TXDb, #,annoDb="org.Hs.eg.db",
                                      level=annotation_level)

}

# Outputs summary of annotation
peakAnno.custom.edb

peakAnno.custom.df <- as.data.frame(peakAnno.custom.edb)
write.table(peakAnno.custom.edb, file="results_annotated.tsv", sep="\t", quote=FALSE, row.names=FALSE)



