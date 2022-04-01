#!/usr/bin/env nextflow

// input is a samplesheet defining the experiment
params.ss = 'samplesheet.csv'
params.frags = 'gatc_frags.gtf'
params.min_reads = 10
params.min_length = 300
params.control = 'DAM'
params.treatment = ''
params.fc_cut = 2
params.fdr_cut = 0.01
params.outdir = "$baseDir/outdir"
params.anngtf = ""
params.annlevel = 'gene' // gene or transcript
params.annpriority = ''// File of element priority for fragment annotation

log.info """\
        METHYLATION FRAGMENT ENRICHMENT FOR DAMID (mfed)
        ==================================================
   
        Samplesheet: 		${params.ss}
        Genome annotations GTF: ${params.anngtf}
        Annotation level:	${params.annlevel}
        Annotation priority file:	${params.annpriority}
        GATC fragment GTF:	${params.frags}
        Minimum frag length:	${params.min_length}
        Minimum reads (in biggest samples):	${params.min_reads}
        Control:		${params.control}
        Treatment:		${params.treatment}
        Fold change cutoff:	${params.fc_cut}
        FDR cutoff:		${params.fdr_cut}	
        Output dir:		${params.outdir}
        """
        .stripIndent()

////
// Set up input channels

// File of GATC fragments which will be the things we are looking for enrichment of
if (params.frags) { 
  Channel
    .fromPath(params.frags, checkIfExists: true)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{frag_ch}
} else { exit 1, 'GATC fragments GTF file not specified!' }

// Set up Annotation GTF channel (used by ChIPseeker in the mfed_diffbind.R script)
if (params.anngtf) { 
  Channel
    .fromPath(params.anngtf, checkIfExists: true)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{anngtf_ch}
} else { exit 1, 'Genome annotation GTF file not specified!' }

// set up optional annotation priority file channel (used by ChIPseeker in the mfed_diffbind.R script)
if (params.annpriority) {
  Channel
    .fromPath(params.annpriority, checkIfExists: false)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{annpriority_ch}
}

// Read in samplesheet and get bams (twice! Not sure how to duplicate channel in this syntax)
if (params.ss) {
  Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bamReads)}
    .set{ bams_ch1 }
} else { exit 1, 'Samplesheet file not specified!' }

Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bamReads)}
    .set{ bams_ch2 }

Channel
    .fromPath(params.ss)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{ss_ch}


// use featureCounts to get counts for each GATC fragment
process featureCounts {
    publishDir params.outdir, mode:'copy' 

    input:
    path frags from frag_ch
    path bams from bams_ch1.collect()

    output:
    path "feature_counts.txt" into counts_ch

    script:
    """
    featureCounts -a ${frags} -o feature_counts.txt -t region -g frag_id ${bams.join(' ')}
    """
}

// filter fragments by length and coverage
process filter_frags {
    publishDir params.outdir, mode:'copy'

    input:
    path fc_result from counts_ch

    output:
    path "filtered_fragments.bed" into filt_frag_ch

    script:
    """
    ${baseDir}/bin/filter_fragments.py -i $fc_result -r ${params.min_reads} -m ${params.min_length} > filtered_fragments.bed
    """
}

// Run DiffBind R script to get enriched fragments - report BED file of all peaks, BED file of significant peaks and Excel file of fold change/FDRs
process diffbind {
    publishDir params.outdir, mode:'copy'

    input:
    path ss from ss_ch
    path bams from bams_ch2.collect()
    path anngtf from anngtf_ch
    path annpriority from annpriority_ch
    path "filtered_fragments.bed" from filt_frag_ch

    output:
    path "plot.pdf" into results_ch1
    path "sample_heatmap_post_contrast.pdf" into results_ch5
    path "volcano_plot.pdf" into results_ch7
    path "MA_plot.pdf" into results_ch8
    path "significant_fragments.bed" into results_ch2
    path "filtered_fragments_diffbind.bed" into results_ch3
    path "results.txt" into results_ch4
    path "results_annotated.tsv" into results_ch6

    script:
    """
    R CMD BATCH --no-save --no-restore \"--args ${ss} ${params.control} ${params.treatment} ${params.fc_cut} ${params.fdr_cut} ${anngtf} ${params.annlevel} ${annpriority}\" ${baseDir}/bin/mfed_diffbind.R .mfed_diffbind.Rout
    """
}

// Make a bedGraph file of log fold changes
process fc_bedgraph
{
    publishDir params.outdir, mode:'copy'

    input:
    path res from results_ch6

    output:
    path "foldchange.bedgraph" into bedgraph_ch

    """
    ${baseDir}/bin/results2bedgraph.py ${res} > foldchange.bedgraph 
    """
}
