version 1.0

import "../../tasks/taxon_id/task_metawrap.wdl" as metawrap_task

workflow strainge_pe_wf {
  meta {
    description: "Assembly binning and taxon assignment with metaWRAP."
  }
  input {
    File assembly_fasta
    File read1
    File read2
    String samplename
    Int metawrap_completion = 80
    Int metawrap_contamination = 10
    Int metawrap_min_contig_length = 1000
    File ncbi_nt_database
    File ncbi_taxonomy_database
    File checkm_database
  }
  call metawrap_task.metawrap {
    input:
      assembly_fasta = assembly_fasta,
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      metawrap_completion = metawrap_completion,
      metawrap_contamination = metawrap_contamination,
      metawrap_min_contig_length = metawrap_min_contig_length,
      ncbi_nt_database = ncbi_nt_database,
      ncbi_taxonomy_database = ncbi_taxonomy_database,
      checkm_database = checkm_database
  }
  output {
    String metawrap_docker = metawrap.metawrap_docker
    String metawrap_version = metawrap.metawrap_version
    String metawrap_analysis_date = metawrap.analysis_date
    File metawrap_insert_sizes = metawrap.metawrap_insert_sizes
    File metawrap_stats = metawrap.metawrap_stats
    Int metawrap_n_bins = metawrap.metawrap_n_bins
    File metawrap_taxonomy = metawrap.metawrap_taxonomy 
  }
}