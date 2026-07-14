// ═════════════════════════════════════════════════════════════════════════════
// PRODUCTION_MODE: Encapsulates full production logic:
// Routes through trimming, filtering, diversity analysis, and differential abundance
// Emits all outputs needed by top-level workflow publish block
// ═════════════════════════════════════════════════════════════════════════════

// Trim primers if requested
include { CUTADAPT ; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from '../modules/cutadapt.nf'

// Cluster ASVs
include { DOWNLOAD_DATABASE } from '../modules/download_database.nf'
include { RUN_DADA2 } from '../modules/run_dada.nf'

// Filtered quality check
include { FASTQC as FILTERED_FASTQC ; MULTIQC as FILTERED_MULTIQC  } from '../modules/quality_assessment.nf'

// Diversity, differential abundance and visualizations
include { ALPHA_DIVERSITY; BETA_DIVERSITY } from '../modules/diversity.nf'
include { PLOT_TAXONOMY } from '../modules/taxonomy_plots.nf'
include { ZIP as ZIP_BIOM; ZIP as ZIP_ALPHA; ZIP as ZIP_BETA_EUCLIDEAN; ZIP as ZIP_BETA_BRAY; ZIP as ZIP_TAXONOMY_SAMPLES; ZIP as ZIP_TAXONOMY_GROUPS } from '../modules/zip.nf'
include { ANCOMBC as ANCOMBC1 } from '../modules/ancombc.nf'
include { ANCOMBC as ANCOMBC2 } from '../modules/ancombc.nf'
include { DESEQ } from '../modules/deseq.nf'
include { ZIP as ZIP_DA; ZIP as ZIP_ANCOMBC1; ZIP as ZIP_ANCOMBC2; ZIP as ZIP_DESEQ2 } from '../modules/zip.nf'

include { SOFTWARE_VERSIONS } from '../modules/utils.nf'

workflow PRODUCTION_MODE {
    take:
    staged_reads_ch
    runsheet_ch
    isa_archive_ch
    gl_file_ch
    primers_ch
    software_versions_ch

    main:
    
    // Download reference database for taxonomic classification
    def db_config = [
            "16S": ["SILVA_SSU_r138_2_v2.RData", "https://api.figshare.com/v2/file/download/64078939"],
            "ITS": ["UNITE_v2025.RData", "https://api.figshare.com/v2/file/download/64079011"],
            "18S": ["PR2_v4_13_March2021.RData", "https://api.figshare.com/v2/file/download/46241917"]
        ]
    target_region_ch = Channel.value(params.target_region)
        .map { region -> tuple(region, db_config[region][0], db_config[region][1]) }
    
    DOWNLOAD_DATABASE(target_region_ch)

    trimmed_reads_ch = channel.empty()
    trimmed_reads_counts = channel.empty()
    cutadapt_logs = channel.empty()
    if(params.trim_primers){

        //if(!params.accession) primers_ch = channel.value([params.F_primer, params.R_primer]) // to be removed once stage analysis workflow is implemented
        CUTADAPT(staged_reads_ch, primers_ch)
        logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
        counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
        trimmed_reads_ch = CUTADAPT.out.reads.map{ 
                                            sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                                            }.flatten().collect()

        COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(counts, logs, runsheet_ch)
        trimmed_reads_counts = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts
        cutadapt_logs = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.logs

        isPaired_ch = CUTADAPT.out.reads.map{ 
                                            sample_id, reads, isPaired -> isPaired
                                            }.first()

        samples_ch = runsheet_ch.first()
                    .concat(isPaired_ch)
                    .collate(2)
        
        
        // Run dada2
        RUN_DADA2(samples_ch, trimmed_reads_ch, trimmed_reads_counts, DOWNLOAD_DATABASE.out.database)

        CUTADAPT.out.version | mix(software_versions_ch) | set{software_versions_ch}
    }else{
        raw_reads_ch = staged_reads_ch.map{
                        sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                        }.flatten().collect()

        isPaired_ch = staged_reads_ch.map{sample_id, reads, isPaired -> isPaired}.first()
        samples_ch = runsheet_ch.first()
                    .concat(isPaired_ch)
                    .collate(2)
        
        // Run dada2 without primer trimming
        RUN_DADA2(samples_ch, raw_reads_ch, file("NO_FILE"), DOWNLOAD_DATABASE.out.database)
    }

    dada_counts = RUN_DADA2.out.counts
    dada_taxonomy = RUN_DADA2.out.taxonomy
    dada_biom = RUN_DADA2.out.biom
    filtered_count = RUN_DADA2.out.filtered_count

    filtered_reads_ch = RUN_DADA2.out.reads
            .flatten()
            .map { file ->
                    // derive sample_id from filename
                    def sample_id
                    if (file.name.endsWith("${params.assay_suffix}_R1_filtered.fastq.gz")) {
                            sample_id = file.name.replace("${params.assay_suffix}_R1_filtered.fastq.gz", "")
                    } else if (file.name.endsWith("${params.assay_suffix}_R2_filtered.fastq.gz")) {
                            sample_id = file.name.replace("${params.assay_suffix}_R2_filtered.fastq.gz", "")
                    }

                    tuple(sample_id, file)
            }
            .groupTuple(by:0)  // group R1/R2 by sample_id
            .map { sample_id, files ->
                    def pathFiles = files.collect { it instanceof String ? file(it) : it }  // ensure Path objects
                    def isPaired = pathFiles.size() > 1
                    tuple(sample_id, pathFiles, isPaired)
            }

    FILTERED_FASTQC(filtered_reads_ch)
        filtered_fastqc_files = FILTERED_FASTQC.out.fastqc.flatten().collect()

    FILTERED_MULTIQC("filtered", params.multiqc_config, filtered_fastqc_files)

    RUN_DADA2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Zip biom file
    dada_biom
        .map { biom -> tuple("taxonomy-and-counts", biom) }
        | ZIP_BIOM

    ZIP_BIOM.out.version | mix(software_versions_ch) | set{software_versions_ch}


    
    // Diversity, differential abundance testing and their corresponding visualizations
    if(params.accession){

        values = ["samples": "Sample Name",
                "group" : "groups",
                "depth" : params.rarefaction_depth,
                "assay_suffix" : params.assay_suffix,
                "output_prefix" : params.cleaned_prefix,
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" : params.remove_rare ? "--remove-rare" : "",
                "struc_zero": params.remove_struc_zeros ? "--remove-structural-zeros" : ""
                ]
    }else{

        values = ["samples": params.samples_column,
                "group" : params.group,
                "depth" : params.rarefaction_depth,
                "assay_suffix" : params.assay_suffix,
                "output_prefix" : params.cleaned_prefix,
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" :  params.remove_rare ? "--remove-rare" : "",
                "struc_zero": params.remove_struc_zeros ? "--remove-structural-zeros" : ""
                ]
    }
    meta  = channel.of(values)
    metadata = runsheet_ch
    
    // Diversity analysis
    ALPHA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)
    BETA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping diversity plots
    // Alpha diversity (if rarefaction succeeded)
    ALPHA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.endsWith('.png') }
            pngs ? tuple(
                "alpha_diversity_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_ALPHA

    // Beta diversity - euclidean distance
    BETA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('euclidean') && it.name.endsWith('.png') }
            pngs ? tuple(
                "euclidean_distance_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_BETA_EUCLIDEAN

    // Beta diversity - bray curtis (if rarefaction succeeded)
    BETA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('bray') && it.name.endsWith('.png') }
            pngs ? tuple(
                "bray_curtis_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_BETA_BRAY

    // Taxonomy plotting
    PLOT_TAXONOMY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping taxonomy plots
// Sample plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('samples') && it.name.endsWith('.png') }
            pngs ? tuple(
                "sample_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_SAMPLES

// Group taxonomy plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('groups') && it.name.endsWith('.png') }
            pngs ? tuple(
                "group_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_GROUPS
    
    ALPHA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    BETA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    PLOT_TAXONOMY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
    // Differential abundance testing
    ancombc1_ch = channel.empty()
    zip_ancombc1_ch = channel.empty()
    ancombc2_ch = channel.empty()
    zip_ancombc2_ch = channel.empty()
    deseq2_ch = channel.empty()
    zip_deseq2_ch = channel.empty()
    method = channel.of(params.diff_abund_method)
    if (params.diff_abund_method == "deseq2"){
    
        DESEQ(meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        deseq2_ch = DESEQ.out.output_dir
        da_contrasts_ch = DESEQ.out.contrasts_file
        da_sampleTable_ch = DESEQ.out.sample_table_file
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping DESeq2 plots
        DESEQ.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip
    
    }else if (params.diff_abund_method == "ancombc1"){
    
        ANCOMBC1(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC1 plots
        ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip

    }else if (params.diff_abund_method == "ancombc2"){

        ANCOMBC2(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc2_ch = ANCOMBC2.out.output_dir
        da_contrasts_ch = ANCOMBC2.out.contrasts_file
        da_sampleTable_ch = ANCOMBC2.out.sample_table_file
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC2 plots
        ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip

    }else{

        ANCOMBC1("ancombc1", meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}

        ANCOMBC2("ancombc2", meta, dada_counts, dada_taxonomy, metadata, ANCOMBC1.out.output_dir)
        ancombc2_ch = ANCOMBC2.out.output_dir
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}

        DESEQ(meta, dada_counts, dada_taxonomy, metadata, ANCOMBC2.out.output_dir)
        deseq2_ch = DESEQ.out.output_dir
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}

        // Zipping DA plots
        //ANCOMBC1
        ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip
        //ANCOMBC2
        ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip
        // DESeq2
        DESEQ.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip
    }
    

    // Software Version Capturing - combining all captured software versions
    nf_version = "Nextflow Version ".concat("${nextflow.version}")
    nextflow_version_ch = channel.value(nf_version)

    // Collect software versions and emit for publishing
    collected_versions_ch = software_versions_ch | map { it.text.strip() }
                            | unique
                            | mix(nextflow_version_ch)
                            | collectFile({it -> it}, newLine: true, cache: false)
                            | SOFTWARE_VERSIONS

    emit:
    trimmed_reads = trimmed_reads_ch
    trimmed_count = trimmed_reads_counts
    cutadapt_logs = cutadapt_logs
    filtered_reads = filtered_reads_ch
    filtered_count = filtered_count
    filtered_fastqc = FILTERED_FASTQC.out.fastqc
    zip_multiqc_filtered = FILTERED_MULTIQC.out.zipped_data
    html_multiqc_filtered = FILTERED_MULTIQC.out.html
    asv = RUN_DADA2.out.fasta
    counts = RUN_DADA2.out.counts
    taxonomy = RUN_DADA2.out.taxonomy
    taxonomy_counts = RUN_DADA2.out.taxonomy_count
    biom_zip = ZIP_BIOM.out.zip
    read_count_tracking = RUN_DADA2.out.read_count
    alpha_diversity = ALPHA_DIVERSITY.out.output_dir
    zip_alpha_plots = ZIP_ALPHA.out.zip
    beta_diversity = BETA_DIVERSITY.out.output_dir
    zip_beta_euclidean_plots = ZIP_BETA_EUCLIDEAN.out.zip
    zip_beta_bray_plots = ZIP_BETA_BRAY.out.zip
    taxonomy_plots = PLOT_TAXONOMY.out.output_dir
    zip_taxonomy_samples = ZIP_TAXONOMY_SAMPLES.out.zip
    zip_taxonomy_groups = ZIP_TAXONOMY_GROUPS.out.zip
    da_contrasts = da_contrasts_ch
    da_sampleTable = da_sampleTable_ch
    ancombc1 = ancombc1_ch
    zip_ancombc1 = zip_ancombc1_ch
    ancombc2 = ancombc2_ch
    zip_ancombc2 = zip_ancombc2_ch
    deseq2 = deseq2_ch
    zip_deseq2 = zip_deseq2_ch
    software_versions = collected_versions_ch
}