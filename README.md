# mfed
Pipeline for analysing DamID data

## Setup
Install Miniconda (for installing various packages)

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh`

Install Nextflow if you do not have it already (for running pipelines)

`conda install -c bioconda nextflow`

Install Git if you do not have it already (For downloading code repository)

`conda install -c anaconda git`

Pull the mfed git repository

`git clone https://github.com/adamjamesreid/mfed.git`

## Run the pipeline 
### Map the reads with nf-core/chipseq 
1. Set up a design file describing the data as described [here](https://nf-co.re/chipseq/1.2.2/usage), also see the example 'design.csv' in this repository

2. Run nf-core/chipseq (--macs_gsize is set to 0 so that it doesn't run MACS2 and fall over because it can't calculate t for single-end damid reads)

`nextflow run nf-core/chipseq -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome BDGP6 --input design.csv --macs_gsize 0`

### Run mfed nextflow script 

`nextflow run ~/code_development/mfed/mfed.nf --ss mfed_samplesheet_test.csv --treatment hp1fusion --control damonly --frags gatc_frags.gtf --outdir outdir --anngtf dm6.ensGene.gtf --annpriority annotation_priority.txt -c /mnt/home3/nextflow/gurdon.config -with-singularity /mnt/home3/nextflow/mfed/mfed_cruk.sif`

n.b. here i use the UCSC version of the ensembl gene set - which has 'chr' prepended to the sequence names, consistent with the BDGP6 genome version used above for nf-core/chipseq - however I think they use different mitochondrial genomes and so the genome should be explicitly defined above for nf-core/chipseq

## Input files
annotation_priority.csv - example file with priority list for fragment annotations using ChIPseeker
dm6.ensGene.gtf - GFT file of genome annotation which works nicely with the BDGP6 reference
mapping_design_test.csv - example input file for nf-core/chipseq describing the experiment
gatc_frags.gtf - GATC fragments file for fly genome
mfed_cruk.sif - Singularity image required for mfed.nf (this is currently available here /mnt/home3/nextflow/mfed/mfed_cruk.sif)
mfed_samplesheet_test.csv - Example mfed samplesheet for running mfed.nf

Test fastq files are currently located on the Gurdon cluster here: /mnt/bioinfo_sharing/sharing/brand/mfed/

## Output files

## Singularity image
*mfed_cruk.def* is in development as a Singularity definition file capturing the dependencies for mfed.nf

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cruk.def`

## To Do
Genome fragmentation script current outputs BED, when it should do GTF

I have run this to convert bed to GTF:

`cat gatc_frags.bed | perl -ne 'chomp;@a=split/\t/;print "$a[0]\tGATC_frag\tregion\t$a[1]\t$a[2]\t.\t.\t.\tfrag_id=\"$a[0]\:$a[1]\_$a[2]\"\n"' > gatc_frags.gtf`
