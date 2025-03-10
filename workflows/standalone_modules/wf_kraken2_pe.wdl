version 1.0

import "../../tasks/taxon_id/task_kraken2.wdl" as kraken2
import "../../tasks/task_versioning.wdl" as versioning

workflow kraken2_pe_wf {
  meta {
    author: "Marco Teixeira"
    email: "mcarvalh@broadinstitute.org"
    description: "Classify paired-end reads using Kraken2. Estimate abundances with Bracken."
  }
  input {
    String samplename
    File read1
    File read2
    File kraken2_db
    Int bracken_read_len = 150
    String bracken_classification_level = "G"
    Int bracken_min_reads = 10
    Int kraken2_mem = 32
    Int kraken2_cpu = 4
    Int kraken2_disk_size = 256
  }
  call kraken2.kraken2_standalone as kraken2_pe {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = kraken2_db,
      mem = kraken2_mem,
      cpu = kraken2_cpu,
      disk_size = kraken2_disk_size,
      bracken_read_len = bracken_read_len,
      bracken_classification_level = bracken_classification_level,
      bracken_min_reads = bracken_min_reads
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Captures
    String kraken2_pe_wf_version = version_capture.wf_version
    String kraken2_pe_wf_analysis_date = version_capture.date
    # Kraken2
    String kraken2_version = kraken2_pe.kraken2_version
    String kraken2_docker = kraken2_pe.kraken2_docker
    File kraken2_report = kraken2_pe.kraken2_report
    File kraken2_classified_report = kraken2_pe.kraken2_classified_report
    File kraken2_unclassified_read1 = kraken2_pe.kraken2_unclassified_read1
    File? kraken2_unclassified_read2 = kraken2_pe.kraken2_unclassified_read2
    File kraken2_classified_read1 = kraken2_pe.kraken2_classified_read1
    File? kraken2_classified_read2 = kraken2_pe.kraken2_classified_read2
    Float kraken2_percent_human = kraken2_pe.kraken2_percent_human
    File bracken_report = kraken2_pe.bracken_report
    String bracken_version = kraken2_pe.bracken_version
  }
}
