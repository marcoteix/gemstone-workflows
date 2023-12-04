version 1.0

import "../../tasks/species_typing/task_strainge.wdl" as strainge
import "../../tasks/species_typing/task_select_strainge_db.wdl" as select_db

workflow strainge_isolates_pe_wf {
  meta {
    description: "Isolate strain-level detection with StrainGE."
  }
  input {
    String samplename
    File read1
    File read2
    Array[File] strainge_dbs
    String predicted_taxonomy
    Int db_kmer_size = 23
    Int strainge_disk_size = 100
    Int strainge_cpus = 4
    Int strainge_memory = 128
  }
  call select_db.select_reference_db {
    input:
        samplename = samplename,
        gambit_taxonomy = predicted_taxonomy,
        strainge_dbs = strainge_dbs
  }
  call strainge.strainge as strainge_isolate{
    input:
        samplename = samplename,
        reads_1 = read1,
        reads_2 = read2,
        kmer_size = db_kmer_size,
        strainge_db = strainge_dbs[select_reference_db.selected_db],
        disk_size = strainge_disk_size,
        cpus = strainge_cpus,
        memory = strainge_memory
  }
  output {
    File straingst_kmerized_reads = strainge_isolate.straingst_kmerized_reads
    String straingst_reference_db_used = strainge_isolate.straingst_reference_db_used
    File straingst_strains = strainge_isolate.straingst_strains
    File straingst_statistics = strainge_isolate.straingst_statistics
    File straingr_concat_fasta = strainge_isolate.straingr_concat_fasta
    File straingr_read_alignment = strainge_isolate.straingr_read_alignment
    File straingr_variants = strainge_isolate.straingr_variants
    File straingr_report = strainge_isolate.straingr_report
    String strainge_docker = strainge_isolate.strainge_docker
  }
}
