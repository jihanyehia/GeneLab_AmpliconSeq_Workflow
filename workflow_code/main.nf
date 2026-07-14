nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

def prefix = params.output_prefix ?: ""
params.cleaned_prefix = (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix

validateParameters()

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";

/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:

         Amplicon target region : ${params.target_region}
         GLDS or OSD accession : ${params.accession}
         Input csv file : ${params.input_file}
         Output directory: ${params.outdir}
         Database Store Directory: ${params.database_store_path}
         Genelab Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         Trim Primers: ${params.trim_primers}

         Cutadapt Parameters:
         Forward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum Trimmed Reads length: ${params.min_cutadapt_len}
         Primers Are linked: ${params.primers_linked}
         Primers Are Anchored: ${params.anchored_primers}
         Discard Untrimmed Reads: ${params.discard_untrimmed}
 
         Dada2 Parameters:
         Truncate left: ${params.left_trunc}bp
         Truncate right: ${params.right_trunc}bp
         Max error left: ${params.left_maxEE}
         Max error right: ${params.right_maxEE}
         Concatenate Reads: ${params.concatenate_reads_only}
         
         Diversity and Differential abundance Parameters:
         Method: ${params.diff_abund_method}
         Rarefaction Depth: ${params.rarefaction_depth}
         Remove Structural Zeros: ${params.remove_struc_zeros}
         Remove Rare Taxa and Samples: ${params.remove_rare}
         Taxa Prevalence Cut Off: ${params.prevalence_cutoff}
         Sample Library Cut Off: ${params.library_cutoff}
         Groups to Comapre Column: ${params.group}
         Samples Column: ${params.samples_column}

         Debugging Options:
         Limit Samples for Testing: ${params.limit_samples_to}
         Force Processing Single-End: ${params.force_single_end}

         Test Mode:
         Test Mode Active: ${params.test_mode}
         Test Level: ${params.test_level}
         Test Primer n-reads: ${params.test_primer_nreads}
         Test truncLen left: ${params.test_trunc_left}
         Test truncLen right: ${params.test_trunc_right}
         Test maxEE left: ${params.test_maxEE_left}
         Test maxEE right: ${params.test_maxEE_right}

         General Pipeline Settings:
         Nextflow Directory publishing mode: ${params.publish_dir_mode}
         MultiQC configuration file: ${params.multiqc_config}
         Nextflow Error strategy: ${params.errorStrategy}

         Conda Environments:
         dp_tools: ${params.conda_dp_tools}
         fastqc: ${params.conda_fastqc}
         multiqc: ${params.conda_multiqc}
         cutadapt: ${params.conda_cutadapt}
         R: ${params.conda_R}
         Diversity and Differential abundance : ${params.conda_diversity}
         zip: ${params.conda_zip}
         wget: ${params.conda_wget}
         ${c_reset}"""
}

// Stage analysis setup (inputs, and raw reads)
include { STAGE_ANALYSIS } from './subworkflows/stage_analysis.nf'

// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from './modules/quality_assessment.nf'

// Test mode subworkflow
include { TEST_MODE } from './subworkflows/test_mode.nf'

// Production mode subworkflow
include { PRODUCTION_MODE } from './subworkflows/production_mode.nf'

ch_dp_tools_plugin = params.dp_tools_plugin ? channel.value(file(params.dp_tools_plugin)) : channel.value(file("$projectDir/bin/dp_tools__NF_AmpIllumina_${params.target_region}"))

ch_input_file = params.input_file ? channel.fromPath(params.input_file) : null
ch_isa_archive = params.isa_archive ? channel.fromPath(params.isa_archive) : null

// A function to delete white spaces from an input string and covert it to lower case
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}

workflow {
    main:

    //  ---------------------  Sanity Checks ------------------------------------- //
    // Test input requirement
    if (!params.accession &&  !params.input_file){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply either an accession (OSD or Genelab number) or an input CSV file
              by passing either to the --accession or --input_file parameter, respectively.
              ${c_reset}""")
    } 
    
    // Test input csv file
    if(params.input_file){
        // Test primers
        if(!params.F_primer || !params.R_primer){

            error("""${c_back_bright_red}PRIMER ERROR! 
                  When using a csv file as input (--input_file) to this workflow you must provide 
                  forward and reverse primer sequences. Please provide your forward 
                  and reverse primer sequences as arguements to the --F_primer 
                  and --R_primer parameters, respectively.
                  ${c_reset}""")
         }
    }

    // Test ISA archive and accession
    if (params.isa_archive && !params.accession) {
        error """${c_back_bright_red}INPUT ERROR!
            --isa_archive requires --accession to resolve OSD/GLDS accessions
            for the ISA-to-runsheet conversion.${c_reset}"""
    }

    // Test mode sanity check
    if (params.test_mode) {
        def valid_levels = ['primers', 'filter', 'full']
        if (!(params.test_level in valid_levels)) {
            error("""${c_back_bright_red}TEST MODE ERROR!
                  --test_level must be one of: primers, filter, full
                  Got: ${params.test_level}
                  ${c_reset}""")
        }

        if (params.limit_samples_to) {
            log.warn """${c_back_bright_red}TEST MODE WARNING!
                  --limit_samples_to is set but test mode stages the full dataset
                  and selects representative samples automatically via quality ranking.
                  limit_samples_to will be IGNORED in test mode.${c_reset}"""
        }

        if (params.input_file) {
            log.info """${c_bright_green}
            TEST MODE ACTIVE  —  level: ${params.test_level}
            Runs test mode on samples provided in ${params.input_file}, 
            ignoring automatic sample selection and quality ranking.
         ${c_reset}"""
        }
        else {
            log.info """${c_bright_green}
            TEST MODE ACTIVE  —  level: ${params.test_level}
            Stages full dataset, selects ${String.valueOf(params.test_n_samples ?: 3)} representative samples to run test mode on.
         ${c_reset}"""
        }
    }

    software_versions_ch = channel.empty()
        
    // Stage analysis setup (inputs, and raw reads)
    STAGE_ANALYSIS(
        params.accession,
        params.target_region,
        ch_input_file,
        ch_isa_archive,
        params.api_url,
        ch_dp_tools_plugin
    )
    staged_reads_ch = STAGE_ANALYSIS.out.staged_reads
    runsheet_ch = STAGE_ANALYSIS.out.runsheet
    isa_archive_ch = STAGE_ANALYSIS.out.isa_archive
    gl_file_ch = STAGE_ANALYSIS.out.gl_file
    primers_ch = STAGE_ANALYSIS.out.primers
    
    STAGE_ANALYSIS.out.software_versions
        | mix(software_versions_ch)
        | set { software_versions_ch }


    // Read quality check and trimming
    RAW_FASTQC(staged_reads_ch)
    raw_fastqc_files = RAW_FASTQC.out.fastqc.flatten().collect()
    
    RAW_MULTIQC("raw", params.multiqc_config,raw_fastqc_files)

    RAW_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    RAW_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // ─────────────────────────────────────────────────────────────────────────
    // TEST MODE or FULL PIPELINE ROUTING
    // ─────────────────────────────────────────────────────────────────────────
    test_results = params.test_mode ?
        TEST_MODE(staged_reads_ch, runsheet_ch, RAW_MULTIQC.out.data, primers_ch) :
        null

    full_results = !params.test_mode ?
        PRODUCTION_MODE(staged_reads_ch, runsheet_ch, isa_archive_ch, gl_file_ch, 
                         primers_ch, software_versions_ch) :
        null

publish:
    // Metadata
    runsheet = runsheet_ch
    isa_archive = isa_archive_ch
    gl_file = gl_file_ch

    // Test mode report
    test_mode_report = params.test_mode ? test_results.report : channel.empty()

    // Raw reads
    raw_reads = staged_reads_ch

    // Trimmed reads
    trimmed_reads = !params.test_mode ? full_results.trimmed_reads : channel.empty()
    trimmed_count = !params.test_mode ? full_results.trimmed_count : channel.empty()
    cutadapt_logs = !params.test_mode ? full_results.cutadapt_logs : channel.empty()
    
    // Filtered reads
    filtered_reads = !params.test_mode ? full_results.filtered_reads : channel.empty()
    filtered_count = !params.test_mode ? full_results.filtered_count : channel.empty()

    // FastQC
    raw_fastqc = RAW_FASTQC.out.fastqc
    filtered_fastqc = !params.test_mode ? full_results.filtered_fastqc : channel.empty()

    // MultiQC
    zip_multiqc_raw = RAW_MULTIQC.out.zipped_data
    html_multiqc_raw = RAW_MULTIQC.out.html
    zip_multiqc_filtered = !params.test_mode ? full_results.zip_multiqc_filtered : channel.empty()
    html_multiqc_filtered = !params.test_mode ? full_results.html_multiqc_filtered : channel.empty()

    // Dada2 outputs
    asv = !params.test_mode ? full_results.asv : channel.empty()
    counts = !params.test_mode ? full_results.counts : channel.empty()
    taxonomy = !params.test_mode ? full_results.taxonomy : channel.empty()
    taxonomy_counts = !params.test_mode ? full_results.taxonomy_counts : channel.empty()
    biom_zip = !params.test_mode ? full_results.biom_zip : channel.empty()
    read_count_tracking = !params.test_mode ? full_results.read_count_tracking : channel.empty()

    // Alpha and beta diversity outputs
    alpha_diversity = !params.test_mode ? full_results.alpha_diversity : channel.empty()
    zip_alpha_plots = !params.test_mode ? full_results.zip_alpha_plots : channel.empty()
    beta_diversity = !params.test_mode ? full_results.beta_diversity : channel.empty()
    zip_beta_euclidean_plots = !params.test_mode ? full_results.zip_beta_euclidean_plots : channel.empty()
    zip_beta_bray_plots = !params.test_mode ? full_results.zip_beta_bray_plots : channel.empty()

    // Taxonomy plots
    taxonomy_plots = !params.test_mode ? full_results.taxonomy_plots : channel.empty()
    zip_taxonomy_samples = !params.test_mode ? full_results.zip_taxonomy_samples : channel.empty()
    zip_taxonomy_groups = !params.test_mode ? full_results.zip_taxonomy_groups : channel.empty()

    // Differential abundance outputs
    da_contrasts = !params.test_mode ? full_results.da_contrasts : channel.empty()
    da_sampleTable = !params.test_mode ? full_results.da_sampleTable : channel.empty()
    ancombc1 = !params.test_mode ? full_results.ancombc1 : channel.empty()
    zip_ancombc1 = !params.test_mode ? full_results.zip_ancombc1 : channel.empty()
    ancombc2 = !params.test_mode ? full_results.ancombc2 : channel.empty()
    zip_ancombc2 = !params.test_mode ? full_results.zip_ancombc2 : channel.empty()
    deseq2 = !params.test_mode ? full_results.deseq2 : channel.empty()
    zip_deseq2 = !params.test_mode ? full_results.zip_deseq2 : channel.empty()

    // GeneLab
    software_versions = !params.test_mode ? full_results.software_versions : channel.empty()

}

output {
    // Metadata
    runsheet { path "Metadata" }
    isa_archive { path "Metadata" }
    gl_file { path "Metadata" }

    // Test mode report
    test_mode_report { path "Test_Mode_Report" }

    // Raw reads
    raw_reads { path "Raw_Sequence_Data" }

    // Trimmed reads
    trimmed_reads { path "Trimmed_Sequence_Data" }
    trimmed_count { path "Trimmed_Sequence_Data" }
    cutadapt_logs { path "Trimmed_Sequence_Data" }
    
    // Filtered reads
    filtered_reads { path "Filtered_Sequence_Data" }
    filtered_count { path "Filtered_Sequence_Data" }

    // FastQC
    raw_fastqc { path {html, zip -> "Raw_Sequence_Data/FastQC_Outputs" } }
    filtered_fastqc { path {html, zip -> "Filtered_Sequence_Data/FastQC_Outputs" } }

    // MultiQC
    zip_multiqc_raw { path "Raw_Sequence_Data/MultiQC_Reports" }
    html_multiqc_raw { path "Raw_Sequence_Data/MultiQC_Reports" }

    zip_multiqc_filtered { path "Filtered_Sequence_Data/MultiQC_Reports" }
    html_multiqc_filtered { path "Filtered_Sequence_Data/MultiQC_Reports" }

    // Dada2 outputs
    asv { path "Final_Outputs" }
    counts { path "Final_Outputs" }
    taxonomy { path "Final_Outputs" }
    taxonomy_counts { path "Final_Outputs" }
    biom_zip { path "Final_Outputs" }
    read_count_tracking { path "Final_Outputs" }

    // Alpha and beta diversity outputs
    alpha_diversity { path "Final_Outputs" }
    zip_alpha_plots { path "Final_Outputs/alpha_diversity" }

    beta_diversity { path "Final_Outputs" }
    zip_beta_euclidean_plots { path "Final_Outputs/beta_diversity" }
    zip_beta_bray_plots { path "Final_Outputs/beta_diversity" }

    // Taxonomy plots
    taxonomy_plots { path "Final_Outputs" }
    zip_taxonomy_samples { path "Final_Outputs/taxonomy_plots" }
    zip_taxonomy_groups { path "Final_Outputs/taxonomy_plots" }

    // Differential abundance outputs
    da_contrasts { path "Final_Outputs" }
    da_sampleTable { path "Final_Outputs" }

    ancombc1 { path "Final_Outputs" }
    zip_ancombc1 { path "Final_Outputs/differential_abundance/ancombc1" }

    ancombc2 { path "Final_Outputs" }
    zip_ancombc2 { path "Final_Outputs/differential_abundance/ancombc2" }

    deseq2 { path "Final_Outputs" }
    zip_deseq2 { path "Final_Outputs/differential_abundance/deseq2" }

    // GeneLab
    software_versions { path "GeneLab" }
}