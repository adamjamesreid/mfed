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

log.info """\
        METHYLATION FRAGMENT ENRICHMENT FOR DAMID (mfed)
        ==================================================
   
        Samplesheet: 		${params.ss}
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

// Set up input channels
Channel
    .fromPath(params.frags)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{frag_ch}

// Read in samplesheet and get bams (twice! Not sure how to duplicate channel in this syntax)
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bamReads)}
    .set{ bams_ch1 }

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
    path "filtered_fragments.bed" from filt_frag_ch

    output:
    path "plot.pdf" into results_ch1
    path "sample_heatmap_post_contrast.pdf" into results_ch5
    path "significant_fragments.bed" into results_ch2
    path "filtered_fragments_diffbind.bed" into results_ch3
    path "results.txt" into results_ch4

    script:
    """
    R CMD BATCH --no-save --no-restore \"--args ${ss} ${params.control} ${params.treatment} ${params.fc_cut} ${params.fdr_cut}\" ${baseDir}/bin/mfed_diffbind.R .mfed_diffbind.Rout
    """
}
