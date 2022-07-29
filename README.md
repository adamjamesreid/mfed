# mfed v1.0
The mfed pipeline is designed for analysing DamID-seq data. Initially the reads are mapped using the nf-core/chipseq pipeline, then all in silico GATC fragments are filtered based on size and number of mapping reads. The filtered fragments are used as a peakset in DiffBind and DESeq2 is used to look for fragments enriched between Dam-fusion samples and Dam-only samples.

## Contents

* [Setup](#Setup)

* [Run the pipeline with test data](#Run-the-pipeline-with-test-data)

  * [Map the reads with nf-core/chipseq](#Map-the-reads-with-nf-core/chipseq)
  * [Run mfed nextflow script](#Run-mfed-nextflow-script)

* [Input files](#Input-files)

* [Output files](#Output-files)

* [More options](#More-options)

* [Load an IGV sesssion](#Load-an-IGV-sesssion)

* [Singularity image](#Singularity-image)

* [FAQ](#FAQ)


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

N.b. this is designed to work nicely on the Gurdon compute cluster, but will need adapting elsewhere

### Map the reads with nf-core/chipseq 
1. Copy test data to the current directory:

`cp /mnt/bioinfo_sharing/sharing/brand/mfed/*gz .`

2. Set up a design file describing the data as described [here](https://nf-co.re/chipseq/1.2.2/usage). For testing use 'mapping_design_test.csv' from this repository. Make sure you use the full paths to the fastq files, or that they are in the current directory.

3. Run nf-core/chipseq (--macs_gsize is set to 0 so that it doesn't run MACS2 and fall over because it can't calculate t for single-end damid reads)

n.b. If you have a Nexflow version <=21, leave out '-dsl1'

`nextflow run nf-core/chipseq -r 1.2.2 -profile singularity -c /mnt/home3/nextflow/gurdon.config --single_end --genome dm6 --input mfed/mapping_design_test.csv --macs_gsize 0 -dsl1`

n.b. this step helpfully removes duplicates and multimapping reads

### Run mfed nextflow pipeline 

1. Generate a samplesheet for mfed describing the data 

This is formatted a bit differently to the one for nf-core/chipseq and conforms to the format of DiffBind samplesheets described in the "Reading in the peakset" section [here](https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf). The test example file is called mfed_samplesheet_test.csv in this repository. Note that *filtered_fragments.bed* is the name of a file generated by the mfed.nf pipeline. You can use this script to make a DiffBind/mfed samplesheet from the nf-core/chipseq one (the -d argument is the full path to the BAM files generated by nf-core/chipseq):

`mfed/bin/design_to_samplesheet.py -i mfed/mapping_design_test.csv -d $PWD/results/bwa/mergedLibrary/ > mfed_samplesheet_test.csv`

2. Run mfed Nextflow pipeline

n.b. If you have a Nexflow version <=21, leave out '-dsl1'

`nextflow run mfed/mfed.nf --ss mfed_samplesheet_test.csv --treatment hp1fusion --control damonly --frags mfed/gatc_frags.gtf --outdir outdir --anngtf mfed/dm6.ensGene.gtf --annpriority mfed/annotation_priority.csv --genome results/genome/genome.fa -c /mnt/home3/nextflow/gurdon.config -with-singularity /mnt/home3/nextflow/mfed/mfed_cruk.sif -dsl1`

n.b. here i use the UCSC version of the ensembl gene set - which has 'chr' prepended to the sequence names, consistent with the dm6 genome version used above for nf-core/chipseq

n.b. When you run the pipeline with your own data, make sure to change '--treatment' and '--control' to match entries in the 'group' column of the samplesheet

## Input files

**annotation_priority.csv** - example file with priority list for fragment annotations using ChIPseeker

**dm6.ensGene.gtf** - GFT file of genome annotation which works nicely with the dm6 reference

**gatc_frags.gtf** - GATC fragments file for fly genome

**mapping_design_test.csv** - example input file for nf-core/chipseq describing the experiment

**mfed_cruk.sif** - Singularity image required for mfed.nf (this is currently available here /mnt/home3/nextflow/mfed/mfed_cruk.sif)

**mfed_samplesheet_test.csv** - Example mfed samplesheet for running mfed.nf

**results/genome/genome.fa** - This is the genome sequence (as used in the nf-core/chipseq pipeline). This is required for making a working IGV tarball

**gurdon.config** - nextflow config file for running jobs on Slurm at the Gurdon Institute

Test fastq files are currently located on the Gurdon cluster here: /mnt/bioinfo_sharing/sharing/brand/mfed/

## Output files

1. From nf-core/chipseq 
* MultiQC results - copy to local machine and view in a web browser

**results/multiqc/broadPeak/multiqc_report.html**

* BAM files of mapped reads

**results/bwa/mergedLibrary/*bam**

2. From mfed pipeline (in outdir/ directory)

* ***.bw** - BigWig files generated from BAM files

* **enriched_fragments.bed** - BED file of fragments enriched in fusion versus dam-only

* **feature_counts.txt** - lists every GATC fragment in the genome, describing the size and number of reads mapping to it in each sample

* **filtered_fragments.bed** - All the fragments which have passed initial filtering and will go into the DiffBind analysis
 
* **foldchange.bedgraph**  - fold changes for significant fragments across the genome (view in IGV)

* **MA_plot.pdf** - plot of average abundance versus fold change for each fragment

* **mfed_diffbind.Rout** - STDOUT from R session. Has some useful results from ChIPSeeker not available elsewhere

* **mfed.RData** - R data object, allowing reanalysis/redrawing of figures etc. (see FAQ for details)

* **mfed_results_for_igv.tar.gz** - tarball of files to load an IGV session

* **results_all.tsv** - Details of significant fragments, enriched in both fusion and dam-only, with fold changes and FDRs

* **results_annotated.tsv** - Details of significant fragments, enriched in fusion versus dam-only, with fold changes and FDRs and nearest gene features

* **results.tsv** - Details of significant fragments, enriched in fusion versus dam-only, with fold changes and FDRs

* **sample_heatmap_plot.pdf** - heatmap comparing samples

* **sample_pca_plot.pdf** - PCA plot comparing samples

* **volcano_plot.pdf** - volcano plot of fold changes against p values

## More options

You can make your own GATC fragment file using the script mfed/bin/fragment_genome.py, e.g.

`mfed/bin/fragment_genome.py results/genome/genome.fa > gatc_frags.gtf`

When running mfed.nf:
* **--min_reads** sets the minimum number of reads mapping to a fragment across all samples added together for it to pass the inital filtering (default = 10
* **--min_length** sets the minimum length in of a fragment for it to pass initial filtering (default = 150)
* **--control** is the name of the control condition as specified in the samplesheet e.g. damonly
* **--treatment** is the name of the treatment/fusion condition specified in the samplesheet e.g. hp1fusion
* **--fc_cut** is the log2 fold change cutoff to consider a fragment as significant (default = 2)
* **--fdr_cut** is the False Discovery Rate (FDR) cutoff to consider a fragment as significant (default = 0.01)
* **--outdir** is the directory for the output files (default = outdir)
* **--anngtf** is a GTF format file of genome annotation for determining the nearest gene
* **--annlevel** used to determine whether to use gene or transcript features - can be 'gene' or 'transcript' (default = gene)
* **--annpriority** file used to determine the order of priority of annotations e.g. are you most interested in Promoter or Exon features (default = use the file provided in the repository - annotation_priority.csv
* **--tss_region_start** TSS region start used in ChIPSeeker analysis to determine location of significant peaks with respect to genes (default = -1000)
* **--tss_region_end** TSS region end (default = 1000)

## Load an IGV sesssion

The mfed pipeline produces a tarball of files (mfed_results_for_igv.tar.gz), including a IGV session file to browse the results in the IGV genome browser

Copy the tarball to where you can use IGV - for those at Gurdon, this might be your laptop. On your laptop run:

```
mkdir igv_files

cd igv_files

scp <user>@cb-milan1.gurdon.private.cam.ac.uk:<path>/mfed_results_for_igv.tar.gz .

tar -xzvf mfed_results_for_igv.tar.gz
```

Install IGV if you need to from [here](https://software.broadinstitute.org/software/igv/download)

Open IGV, 'File' -> 'Open Session', select *mfed_results_for_igv.tar.gz*

## Singularity image
The software dependencies for mfed.nf are captured in a Singularity image. This can be generated using this Singularity definition file: *mfed_cruk.def*.

It makes use of the DiffBind install in the CRUK DiffBind workshop Singularity image described [here](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core/software/diffbind-tool-for-chip-seq-and-atac-seq-analysis)

Build like this: `sudo singularity build mfed_cruk.sif mfed_cruk.def`

## FAQ

* *Why does conda not run?*

It may not be in your path. Try specifying the full path e.g. `~/miniconda3/bin/conda`. Ideally add `~/miniconda3/bin/` in your $PATH environment variable. See this [tutorial:]( https://riptutorial.com/bash/example/19613/add-a-path-to-the-path-environment-variable).

* *Can I have underlying data files for PCA, volcano plot, heatmap, MA plot etc.?*

Mfed now saves an R image so that you can revisit the analysis. The file is called ‘mfed.RData’. In R, do:

```
library(DiffBind)
load(‘mfed.RData’)
```

* *Can I get the R object from DiffBind?*

See above

* *Can I process lots of different fusion and dam-only samples together?*

You can process any number of different fusion and dam-only conditions together in the mapping step. However you have to run the mfed analysis step once for each fusion versus dam-only comparison. You can do this in two ways. Firstly, if you use the whole samplesheet, with all the samples in it, they will get normalised together. If some conditions have lots of variation between replicates this will increase the estimate of variance across all conditions and will affect all comparisons. Alternatively you can generate the mfed samplesheet with `design_to_samplesheet.py` and then split out the samples you need for each fusion into separate samplesheets, running each one completely separately.

* *Can I change the TSS size for ChipSeeker?*

You can! Recently arguments were added to mfed.nf - `--tss_region_start` and `--tss_region_end`.

* *The pipeline fails at the DiffBind stage. How can I tell what is wrong?*

The nextflow output should tell you which folder the DiffBind analysis was run in. Something like ‘work/0b/71e25020b58f2b7192dabce04f4931/’. The error message ought to tell you which directory, or you could look in the output for a line like this:

`[13/b2f5aa] process > diffbind           [100%] 1 of 1 ✔`

In square brackets is the start of the name of the folder where diffbind was run. In that folder will be a file called `mfed_diffbind.Rout`. Towards the bottom of that file will be some error messages from R about what went wrong. This file is now also saved to the mfed output directory.

One thing that sometimes causes an error is not having ‘treatment’ and ‘control’ arguments which match the names in the ‘group’ column of the mfed samplesheet. Worth checking!

* *How can I see the output of DiffBind and ChIPSeeker?*

The output (STDOUT) from the R session is saved in the file ‘mfed_dffbind.R’ in the mfed output directory

* *Why does Nextflow use so much disk space?*

Nextflow saves all the intermediate files it uses so that failed runs can be troubleshooted and then resumed. If you are happy that none of your runs need to be troubleshooted or resumed, you can remove the ‘work/’ directory and make your sys admin happy again.


## To Do

Generate fragments on the fly in mfed.nf

Try to make it a requirement for only one samplesheet

Generate fragment file with cut sites annotated for viewing in IGV

