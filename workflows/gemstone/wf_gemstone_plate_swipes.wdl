version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf
import "../utilities/wf_metaspades_assembly.wdl" as metaspades_assembly_wf
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken_task
import "../../tasks/taxon_id/task_metawrap.wdl" as metawrap_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/gene_typing/task_mob_suite.wdl" as mob_suite_task
import "../../tasks/gene_typing/task_bakta.wdl" as bakta_task
import "../../tasks/gene_typing/task_amrfinderplus.wdl" as amrfinderplus
import "../standalone_modules/wf_strainge_pe.wdl" as strainge_wf

workflow gemstone_plate_swipes {
  meta {
    author: "Marco Teixeira"
    email: "mcarvalh@broadinstitute.org"
    description: "Reference-based consensus calling or de novo assembly for metagenomic sequencing data with the TheiaMeta Illumina PE workflow from Theiagen PHB. Strain-level detection, genome binning, taxonomy assignment, and AMR genotyping."
  }
  input {
    File read1
    File read2
    String samplename
    File? reference
    String lab_determined_genus
    Boolean output_additional_files = false
    # StrainGE options
    Boolean call_strainge = false
    File strainge_db_config
    Int strainge_max_strains = 5
    Int strainge_disk_size = 100
    Int strainge_cpus = 4
    Int strainge_memory = 128
    Int strainge_db_kmer_size = 23
    # Kraken options
    Boolean call_kraken2 = false
    File kraken2_db
    Int bracken_read_len = 150
    String bracken_classification_level = "G"
    Int kraken2_mem = 32
    Int kraken2_cpu = 4
    Int kraken2_disk_size = 256
    # metaWRAP options
    Boolean call_metawrap = false
    File? metawrap_checkm_db
    Int metawrap_completion = 80
    Int metawrap_contamination = 10
    Int metawrap_min_contig_length = 1000
    String binning_flags = "--metabat2 --maxbin2 --concoct"
  }
  call read_qc_wf.read_QC_trim_pe as read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1,
      read2_raw = read2,
      workflow_series = "theiameta"
  }
  if (call_kraken2) {
    call kraken_task.kraken2_standalone as kraken2_clean {
      input:
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        kraken2_db = kraken2_db,
        kraken2_args = "",
        classified_out = "classified#.fastq",
        unclassified_out = "unclassified#.fastq",
        mem = kraken2_mem,
        cpu = kraken2_cpu,
        disk_size = kraken2_disk_size,
        bracken_read_len = bracken_read_len,
        bracken_classification_level = bracken_classification_level
    }
  }
  call metaspades_assembly_wf.metaspades_assembly_pe as metaspades {
  input:
    read1 = read_QC_trim.read1_clean,
    read2 = read_QC_trim.read2_clean,
    samplename = samplename
  }
  # if reference is provided, perform mapping of assembled contigs to 
  # reference with minimap2, and extract those as final assembly
  if (defined(reference)){
    call minimap2_task.minimap2 as minimap2_assembly {
      input:
        query1 = metaspades.assembly_fasta,
        reference = select_first([reference]),
        samplename = samplename,
        mode = "asm20",
        output_sam = false
    }
    call parse_mapping_task.retrieve_aligned_contig_paf {
      input:
        paf = minimap2_assembly.minimap2_out,
        assembly = metaspades.assembly_fasta,
        samplename = samplename
    }
    call parse_mapping_task.calculate_coverage_paf {
      input:
        paf = minimap2_assembly.minimap2_out
    }
  }
  call quast_task.quast {
    input:
      assembly = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta]),
      samplename = samplename,
      min_contig_len = 1
  }
  if (output_additional_files){
    call minimap2_task.minimap2 as minimap2_reads {
      input:
        query1 = read_QC_trim.read1_clean,
        query2 = read_QC_trim.read2_clean, 
        reference = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta]),
        samplename = samplename,
        mode = "sr",
        output_sam = true
    }
    call parse_mapping_task.sam_to_sorted_bam {
      input:
        sam = minimap2_reads.minimap2_out,
        samplename = samplename
    }
    call parse_mapping_task.calculate_coverage {
      input:
        bam = sam_to_sorted_bam.bam,
        bai = sam_to_sorted_bam.bai
    }
    call parse_mapping_task.retrieve_pe_reads_bam as retrieve_unaligned_pe_reads_sam {
      input:
        bam = sam_to_sorted_bam.bam,
        samplename = samplename,
        prefix = "unassembled",
        sam_flag = 4
    }
    call parse_mapping_task.retrieve_pe_reads_bam as retrieve_aligned_pe_reads_sam {
      input:
        bam = sam_to_sorted_bam.bam,
        samplename = samplename,
        sam_flag = 2,
        prefix = "assembled"
    }
    call parse_mapping_task.assembled_reads_percent {
      input:
        bam = sam_to_sorted_bam.bam,
    } 
  }  
  call bakta_task.bakta {
    input:
      assembly = metaspades.assembly_fasta,
      samplename = samplename
  }
  call mob_suite_task.mob_recon {
    input:
      assembly = metaspades.assembly_fasta,
      samplename = samplename
  }
  call amrfinderplus.amrfinderplus_nuc as amrfinderplus_task {
    input:
      assembly = metaspades.assembly_fasta,
      samplename = samplename
  }
  if (call_strainge) {
    call strainge_wf.strainge_pe_wf as strainge {
      input:
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        strainge_db_kmer_size = strainge_db_kmer_size,
        strainge_disk_size = strainge_disk_size,
        strainge_cpus = strainge_cpus,
        strainge_memory = strainge_memory,
        strainge_max_strains = strainge_max_strains,
        strainge_db_config = strainge_db_config,
        predicted_taxonomy = lab_determined_genus
    }
  }
  if (call_metawrap) {
    call metawrap_task.metawrap_binning as metawrap {
      input:
        samplename = samplename,
        assembly_fasta = metaspades.assembly_fasta,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        metawrap_completion = metawrap_completion,
        metawrap_contamination = metawrap_contamination,
        metawrap_min_contig_length = metawrap_min_contig_length,
        checkm_database = select_first([metawrap_checkm_db]),
        binning_flags = binning_flags
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version capture
    String gemstone_wf_version = version_capture.wf_version
    String analysis_date = version_capture.date
    # Kraken2 outputs
    String? kraken2_version = kraken2_clean.kraken2_version
    String? kraken2_docker = kraken2_clean.kraken2_docker
    File? kraken2_report = kraken2_clean.kraken2_report
    Float? kraken2_percent_human = kraken2_clean.kraken2_percent_human
    File? bracken_report = kraken2_clean.bracken_report
    String? bracken_version = kraken2_clean.bracken_version
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # Read QC - fastq_scan outputs
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    String fastq_scan_docker = read_QC_trim.fastq_scan_docker
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - bbduk outputs
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - Read stats
    Float? average_read_length = read_QC_trim.average_read_length
    # Assembly - metaspades 
    File assembly_fasta = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta])
    String metaspades_version = metaspades.metaspades_version
    String metaspades_docker = metaspades.metaspades_docker
    # Assembly - minimap2
    String minimap2_version = metaspades.minimap2_version
    String minimap2_docker = metaspades.minimap2_docker
    # Assembly - samtools
    String samtools_version = metaspades.samtools_version
    String samtools_docker = metaspades.samtools_docker
    # Assembly - pilon
    String pilon_version = metaspades.pilon_version
    String pilon_docker = metaspades.pilon_docker
    # Assembly QC - quast
    Int assembly_length = quast.genome_length
    Int contig_number = quast.number_contigs
    Int largest_contig = quast.largest_contig
    String quast_version = quast.version
    String quast_docker = quast.quast_docker
    # Assembly QC - minimap2
    Float? percent_coverage = calculate_coverage_paf.percent_coverage
    # Assembly QC - bedtools
    Float? assembly_mean_coverage = calculate_coverage.mean_depth_coverage
    String? bedtools_version = calculate_coverage.bedtools_version
    String? bedtools_docker = calculate_coverage.bedtools_docker
    # Read retrieval
    File? read1_unmapped = retrieve_unaligned_pe_reads_sam.read1
    File? read2_unmapped = retrieve_unaligned_pe_reads_sam.read2
    File? read1_mapped = retrieve_aligned_pe_reads_sam.read1 
    File? read2_mapped = retrieve_aligned_pe_reads_sam.read2
    # Assembly stats
    Float? percentage_mapped_reads = assembled_reads_percent.percentage_mapped
    # StrainGE outputs
    Array[File]? straingst_kmerized_reads = strainge.straingst_kmerized_reads
    Array[File]? straingst_selected_db = strainge.straingst_selected_db
    Boolean? straingst_found_db = strainge.straingst_found_db
    Array[File]? straingst_strains = strainge.straingst_strains
    Array[File]? straingst_statistics = strainge.straingst_statistics
    Array[File]? straingr_concat_fasta = strainge.straingr_concat_fasta
    Array[File]? straingr_read_alignment = strainge.straingr_read_alignment
    Array[File]? straingr_variants = strainge.straingr_variants
    Array[File]? straingr_report = strainge.straingr_report
    Array[String]? strainge_docker = strainge.strainge_docker
    Array[String]? strainge_version = strainge.strainge_version
    # Bakta outputs
    File bakta_gbff = bakta.bakta_gbff
    File bakta_gff3 = bakta.bakta_gff3
    File bakta_tsv = bakta.bakta_tsv
    File bakta_summary = bakta.bakta_txt
    String bakta_version = bakta.bakta_version
    # MOB-recon Results
    File mob_recon_results = mob_recon.mob_recon_results
    File mob_typer_results = mob_recon.mob_typer_results
    File mob_recon_chromosome_fasta = mob_recon.chromosome_fasta
    Array[File] mob_recon_plasmid_fastas = mob_recon.plasmid_fastas
    String mob_recon_docker = mob_recon.mob_recon_docker
    String mob_recon_version = mob_recon.mob_recon_version
    # NCBI-AMRFinderPlus Outputs
    File amrfinderplus_all_report = amrfinderplus_task.amrfinderplus_all_report
    File amrfinderplus_amr_report = amrfinderplus_task.amrfinderplus_amr_report
    File amrfinderplus_stress_report = amrfinderplus_task.amrfinderplus_stress_report
    File amrfinderplus_virulence_report = amrfinderplus_task.amrfinderplus_virulence_report
    String amrfinderplus_amr_core_genes = amrfinderplus_task.amrfinderplus_amr_core_genes
    String amrfinderplus_amr_plus_genes = amrfinderplus_task.amrfinderplus_amr_plus_genes
    String amrfinderplus_stress_genes = amrfinderplus_task.amrfinderplus_stress_genes
    String amrfinderplus_virulence_genes = amrfinderplus_task.amrfinderplus_virulence_genes
    String amrfinderplus_amr_classes = amrfinderplus_task.amrfinderplus_amr_classes
    String amrfinderplus_amr_subclasses = amrfinderplus_task.amrfinderplus_amr_subclasses
    String amrfinderplus_version = amrfinderplus_task.amrfinderplus_version
    String amrfinderplus_db_version = amrfinderplus_task.amrfinderplus_db_version
    # MetaWRAP outputs
    String? metawrap_docker = metawrap.metawrap_docker
    String? metawrap_version = metawrap.metawrap_version
    File? metawrap_stats = metawrap.metawrap_stats
    Int? metawrap_n_bins = metawrap.metawrap_n_bins
    String? metawrap_binning_flags = metawrap.metawrap_binning_flags
    Array[File]? metawrap_fasta = metawrap.metawrap_fasta
    File? metawrap_contigs = metawrap.metawrap_contigs
  }
}