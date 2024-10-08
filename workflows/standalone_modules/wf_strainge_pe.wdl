version 1.0

import "../../tasks/taxon_id/task_strainge.wdl" as strainge
import "../../tasks/taxon_id/task_select_strainge_db.wdl" as select_db
import "../../tasks/task_versioning.wdl" as versioning

workflow strainge_pe_wf {
  meta {
    author: "Marco Teixeira"
    email: "mcarvalh@broadinstitute.org"
    description: "Strain-level detection with StrainGE."
  }
  input {
    String samplename
    File read1
    File read2
    File strainge_db_config
    String predicted_taxonomy
    Int strainge_db_kmer_size = 23
    Int strainge_disk_size = 100
    Int strainge_cpus = 4
    Int strainge_memory = 128
    Int strainge_max_strains = 5
    Boolean strainge_prepare_straingr = false
    Float min_coverage = 0.8
  }
  call select_db.select_reference_db_lite as select_reference_db {
    input:
        samplename = samplename,
        predicted_taxonomy = predicted_taxonomy,
        strainge_db_config = strainge_db_config
  }
  if (select_reference_db.found_db) {
    scatter (strainge_db in select_reference_db.selected_db) {
      call strainge.strainge as strainge_isolate{
        input:
            samplename = samplename,
            reads_1 = read1,
            reads_2 = read2,
            kmer_size = strainge_db_kmer_size,
            strainge_db = strainge_db,
            disk_size = strainge_disk_size,
            cpus = strainge_cpus,
            memory = strainge_memory,
            max_strains = strainge_max_strains,
            prepare_straingr = strainge_prepare_straingr
      }
    }
    call strainge.top_strain {
      input:
        straingst_strains = strainge_isolate.straingst_strains,
        min_coverage = min_coverage
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    String strainge_pe_wf_version = version_capture.wf_version
    String strainge_pe_wf_analysis_date = version_capture.date
    Array[File]? straingst_kmerized_reads = strainge_isolate.straingst_kmerized_reads
    Array[String] straingst_selected_db = select_reference_db.selected_db
    Boolean straingst_found_db = select_reference_db.found_db
    Array[File]? straingst_strains = strainge_isolate.straingst_strains
    Array[File]? straingst_statistics = strainge_isolate.straingst_statistics
    Array[File?]? straingr_concat_fasta = strainge_isolate.straingr_concat_fasta
    Array[File?]? straingr_read_alignment = strainge_isolate.straingr_read_alignment
    Array[File?]? straingr_variants = strainge_isolate.straingr_variants
    Array[File?]? straingr_report = strainge_isolate.straingr_report
    Array[String]? strainge_docker = strainge_isolate.strainge_docker
    Array[String]? strainge_version = strainge_isolate.strainge_version
    String? straingst_top_strain = top_strain.straingst_top_strain
  }
}
