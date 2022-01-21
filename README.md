# mfed
Pipeline for analysing DamID data

## Begin by running nf-core/chipseq with --macs_gsize 0 so that it doesn't run MACS2 and fall over because it can't calculate t for damid data
`nextflow run nf-core/chipseq -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome BDGP6 --input design.csv --macs_gsize 0`

## Then run mfed nextflow script
`nextflow run ~/code_development/mfed/mfed.nf --ss suv39_samplesheet.csv --treatment Suv39 --control DAM --frags ~/code_development/mfed/gatc_frags.gtf --outdir outdir`

## Singularity image
*mfed_cruk.def* is in development as a Singularity definition file capturing the dependencies for mfed.nf

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cru.def`
