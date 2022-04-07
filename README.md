# mfed
The mfed pipeline is designed for analysing DamID-seq data. Initially the reads are mapped using the nf-core/chipseq pipeline, then all in silico GATC fragments are filtered based on size and number of mapping reads. The filtered fragments are used as a peakset in DiffBind and DESeq2 is used to look for fragments enriched between DAM-fusion samples and DAM-only samples.

## Setup
1. Install Miniconda if you have not already (for installing various packages)

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`bash Miniconda3-latest-Linux-x86_64.sh`

2. Install Nextflow if you do not have it already (for running pipelines)

`conda install -c bioconda nextflow`

3. Install Git if you do not have it already (For downloading code repository)

`conda install -c anaconda git`

4. Pull the mfed git repository

`git clone https://github.com/adamjamesreid/mfed.git`

## Run the pipeline with test data
### Map the reads with nf-core/chipseq 
1. Set up a design file describing the data as described [here](https://nf-co.re/chipseq/1.2.2/usage), for testing use 'mapping_design_test.csv' from this repository

2. Copy fastq files locally

`cp /mnt/bioinfo_sharing/sharing/brand/mfed/*gz .`

3. Run nf-core/chipseq (--macs_gsize is set to 0 so that it doesn't run MACS2 and fall over because it can't calculate t for single-end damid reads)

`nextflow run nf-core/chipseq -r 1.2.2 -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome BDGP6 --input mfed/mapping_design_test.csv --macs_gsize 0`

### Run mfed nextflow script 

1. Copy bam files to current directory

`cp results/bwa/mergedLibrary/*bam .`

2. Set up a samplesheet for mfed describing the data (this is formatted a bit differently to the one for nf-core/chipseq and conforms to the format of DiffBind samplesheets described in the "Reading in the peakset" section [here](https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf). The test example file is called mfed_samplesheet_test.csv in this repository.)

3. Run mfed Nextflow pipeline

`nextflow run mfed/mfed.nf --ss mfed/mfed_samplesheet_test.csv --treatment hp1fusion --control damonly --frags mfed/gatc_frags.gtf --outdir outdir --anngtf mfed/dm6.ensGene.gtf --annpriority mfed/annotation_priority.csv -c /mnt/home3/nextflow/gurdon.config -with-singularity /mnt/home3/nextflow/mfed/mfed_cruk.sif`

n.b. here i use the UCSC version of the ensembl gene set - which has 'chr' prepended to the sequence names, consistent with the BDGP6 genome version used above for nf-core/chipseq - however I think they use different mitochondrial genomes and so the genome should be explicitly defined above for nf-core/chipseq

## Input files
annotation_priority.csv - example file with priority list for fragment annotations using ChIPseeker

dm6.ensGene.gtf - GFT file of genome annotation which works nicely with the BDGP6 reference

gatc_frags.gtf - GATC fragments file for fly genome

mapping_design_test.csv - example input file for nf-core/chipseq describing the experiment

mfed_cruk.sif - Singularity image required for mfed.nf (this is currently available here /mnt/home3/nextflow/mfed/mfed_cruk.sif)

mfed_samplesheet_test.csv - Example mfed samplesheet for running mfed.nf

Test fastq files are currently located on the Gurdon cluster here: /mnt/bioinfo_sharing/sharing/brand/mfed/

## Output files

1. From nf-core/chipseq 
* MultiQC results - copy to local machine and view in a web browser

results/multiqc/broadPeak/multiqc_report.html

* BAM files of mapped reads

results/bwa/mergedLibrary/*bam

2. From mfed pipeline (in outdir/ directory)
* feature_counts.txt - lists every GATC fragment in the genome, describing the size and number of reads mapping to it in each sample

* filtered_fragments.bed - All the fragments which have passed initial filtering and will go into the DiffBind analysis

* filtered_fragments_diffbind.bed - Concensus fragment set from DiffBind
 
* foldchange.bedgraph - fold changes for significant fragments across the genome (view in IGV)

* MA_plot.pdf - plot of average abundance versus fold change for each fragment

* plot.pdf - heatmap comparing samples, before normalisation

* results_annotated.tsv - Details of significant fragments, with fold changes and FDRs, with nearest gene features

* results.txt - Details of significant fragments, with fold changes and FDRs

* sample_heatmap_post_contrast.pdf - heatmap comparing samples, after normalisation

* significant_fragments.bed - Fragments which pass the thresholds for significant enrichment

* volcano_plot.pdf - volcano plot of fold changes against p values

## More options

You can make your own GATC fragment file using the script mfed/bin/fragment_genome.py

When running mfed.nf:
* *--min_reads* sets the minimum number of reads mapping to a fragment across all samples added together for it to pass the inital filtering (default = 10
* *--min_length* sets the minimum length in of a fragment for it to pass initial filtering (default = 300)
* *--control* is the name of the control condition as specified in the samplesheet e.g. damonly
* *--treatment* is the name of the treatment/fusion condition specified in the samplesheet e.g. hp1fusion
* *--fc_cut* is the fold change cutoff to consider a fragment as significant (default = 2)
* *--fdr_cut* is the False Discovery Rate (FDR) cutoff to consider a fragment as significant (default = 0.01)
* *--outdir* is the directory for the output files (default = outdir)
* *--anngtf* is a GTF format file of genome annotation for determining the nearest gene
* *--annlevel* used to determine whether to use gene or transcript features - can be 'gene' or 'transcript' (default = gene)
* *--annpriority* file used to determine the order of priority of annotations e.g. are you most interested in Promoter or Exon features (default = use the file provided in the repository - annotation_priority.csv


## Singularity image
*mfed_cruk.def* is in development as a Singularity definition file capturing the dependencies for mfed.nf

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cruk.def`

## To Do
Genome fragmentation script current outputs BED, when it should do GTF

I have run this to convert bed to GTF:

`cat gatc_frags.bed | perl -ne 'chomp;@a=split/\t/;print "$a[0]\tGATC_frag\tregion\t$a[1]\t$a[2]\t.\t.\t.\tfrag_id=\"$a[0]\:$a[1]\_$a[2]\"\n"' > gatc_frags.gtf`

Generate fragments on the fly in mfed.nf

Try to make it a requirement for only one samplesheet

