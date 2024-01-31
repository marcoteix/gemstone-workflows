version 1.0

import "../../tasks/taxon_id/task_metawrap.wdl" as metawrap_task

workflow metawrap_pe_wf {
  meta {
    description: "Assembly binning with metaWRAP."
  }
  input {
    File assembly_fasta
    File read1
    File read2
    String samplename
    Int metawrap_completion = 80
    Int metawrap_contamination = 10
    Int metawrap_min_contig_length = 1000
    File checkm_database
    String binning_flags = "--metabat2 --maxbin2 --concoct"
  }
  call metawrap_task.metawrap_binning as metawrap {
    input:
      assembly_fasta = assembly_fasta,
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      metawrap_completion = metawrap_completion,
      metawrap_contamination = metawrap_contamination,
      metawrap_min_contig_length = metawrap_min_contig_length,
      checkm_database = checkm_database,
      binning_flags = binning_flags
  }
  output {
    String metawrap_docker = metawrap.metawrap_docker
    String metawrap_version = metawrap.metawrap_version
    String metawrap_analysis_date = metawrap.analysis_date
    File metawrap_stats = metawrap.metawrap_stats
    Int metawrap_n_bins = metawrap.metawrap_n_bins
    String metawrap_binning_flags = metawrap.metawrap_binning_flags
    Array[File] metawrap_fasta = metawrap.metawrap_fasta
    File metawrap_contigs = metawrap.metawrap_contigs
  }
}