version 1.0

import "../../../tasks/species_typing/task_strainge.wdl" as strainge_task
import "../../../tasks/species_typing/task_select_strainge_db.wdl" as select_strainge_db_task

workflow test_strainge {
  input {
    String gambit_taxon = 'Escherichia coli'
    File read1
    File read2
    String samplename
    Array[File] strainge_dbs
    Int strainge_kmer_size = 23
  } 
  call select_strainge_db_task.select_reference_db {
    input:
      samplename = samplename,
      gambit_taxonomy = gambit_taxon,
      strainge_dbs = strainge_dbs
  }
  call strainge_task.strainge {
    input:
      samplename = samplename,
      reads_1 = read1,
      reads_2 = read2,
      kmer_size = strainge_kmer_size,
      strainge_db = strainge_dbs[select_reference_db.selected_db],
      disk_size = 16,
      memory = 8,
      cpus = 2
  } 
  output {
    String? strainge_genus_used = select_reference_db.selected_db
    File? straingst_kmerized_reads = strainge.straingst_kmerized_reads
    String? straingst_reference_db = strainge.straingst_reference_db_used
    File? straingst_strains = strainge.straingst_strains
    File? straingst_statistics = strainge.straingst_statistics
    File? straingr_concat_fasta = strainge.straingr_concat_fasta
    File? straingr_read_alignment = strainge.straingr_read_alignment
    File? straingr_variants = strainge.straingr_variants
    File? straingr_report = strainge.straingr_report
    String? strainge_docker = strainge.strainge_docker
    String? strainge_version = strainge.strainge_version
  }
}