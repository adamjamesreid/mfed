# mfed
Pipeline for analysing DamID data

## Begin by running nf-core/chipseq with --macs_gsize 0 so that it doesn't run MACS2 and fall over because it can't calculate t for damid data
`nextflow run nf-core/chipseq -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome BDGP6 --input design.csv --macs_gsize 0`

## Then run mfed nextflow script (here i use the UCSC version of the ensembl gene set - which has 'chr' prepended to the sequence names, consistent with the BDGP6 genome version used above for nf-core/chipseq - however I think they use different mitochondrial genomes and so the genome should be explicitly defined above for nf-core/chipseq)
`nextflow run ~/code_development/mfed/mfed.nf --ss suv39_samplesheet.csv --treatment Suv39 --control DAM --frags ~/code_development/mfed/gatc_frags.gtf --outdir outdir --anngtf peaks2genes_test/dm6.ensGene.gtf --annpriority annotation_priority.txt -c /mnt/home3/nextflow/gurdon.config -with-singularity mfed_cruk.sif`

## Files
annotation_priority.csv - example file with priority list for fragment annotations using ChIPseeker

## Singularity image
*mfed_cruk.def* is in development as a Singularity definition file capturing the dependencies for mfed.nf

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cruk.def`

##To Do
Genome fragmentation script current outputs BED, when it should do GTF

I have run this to convert bed to GTF:

`cat gatc_frags.bed | perl -ne 'chomp;@a=split/\t/;print "$a[0]\tGATC_frag\tregion\t$a[1]\t$a[2]\t.\t.\t.\tfrag_id=\"$a[0]\:$a[1]\_$a[2]\"\n"' > gatc_frags.gtf`
