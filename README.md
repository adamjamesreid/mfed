# mfed
Pipeline for analysing DamID data

## Begin by running nf-core/chipseq with --macs_gsize 0 so that it doesn't run MACS2 and fall over because it can't calculate t for damid data
`nextflow run nf-core/chipseq -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome BDGP6 --input design.csv --macs_gsize 0`

## Then run mfed nextflow script (Using --genome BDGP6 uses the older mitochondrial sequence (not ISO1 MT - the larger one) and gives chromosome names which are not suitable for how I'm using ChIPSeeker for annotation. Here I explicitly supply the reference from Hansong.)
`nextflow run ~/code_development/mfed/mfed.nf --ss suv39_samplesheet.csv --treatment Suv39 --control DAM --frags ~/code_development/mfed/gatc_frags.gtf --outdir outdir`

## Singularity image
*mfed_cruk.def* is in development as a Singularity definition file capturing the dependencies for mfed.nf

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cru.def`

n.b. when using the singularity image auto mounts need to be enabled as in the gurdon.config file

e.g. `nextflow run ~/code_development/mfed/mfed.nf --ss suv39_samplesheet.csv --treatment Suv39 --control DAM --frags ~/code_development/mfed/gatc_frags.gtf --outdir outdir -with-singularity mfed_cruk.sif -c /mnt/home3/nextflow/gurdon.config`


##To Do
Genome fragmentation script current outputs BED, when it shoud do GTF

I have run this to convert bed to GTF:

`cat gatc_frags.bed | perl -ne 'chomp;@a=split/\t/;print "$a[0]\tGATC_frag\tregion\t$a[1]\t$a[2]\t.\t.\t.\tfrag_id=\"$a[0]\:$a[1]\_$a[2]\"\n"' > gatc_frags.gtf`
