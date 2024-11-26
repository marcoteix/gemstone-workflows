version 1.0

import "../epidemiology/task_straingr_matrix.wdl" as straingr_matrix

task straingr_compare {
  input {
    String samplename_1
    String samplename_2
    File variants_1
    File variants_2
    String docker = "marcoteix/strainge:1.0.0"
    Int disk_size = 8
    Int memory = 8
  }
  command <<<
    base="/opt/conda/envs/strainge/bin"

    echo "Comparing variants in samples ~{samplename_1} and ~{samplename_2}..."
    $base/straingr compare ~{variants_1} ~{variants_2} \
        -o ~{samplename_1}.vs.~{samplename_2}.summary.tsv \
        -d ~{samplename_1}.vs.~{samplename_2}.details.tsv
  >>>
  output {
    File straingr_summary = "~{samplename_1}.vs.~{samplename_2}.summary.tsv"
    File straingr_details = "~{samplename_1}.vs.~{samplename_2}.details.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

workflow straingr_compare_query {
  meta {
    description: "Uses StrainGR to compare one sample against a set of samples."
  }
  input {
    String reference_name
    Array[String] query_names
    File reference_variants
    Array[File] query_variants
    Int disk_size = 8
    Int memory = 16
  }
  scatter (i in range(length(query_names))) {
    call straingr_compare {
      input:
        samplename_1 = reference_name,
        samplename_2 = query_names[i],
        variants_1 = reference_variants,
        variants_2 = query_variants[i],
        disk_size = disk_size,
        memory = memory
    }
  }
  call straingr_matrix.straingr_gather_query {
    input:
      reference_name = reference_name,
      query_names = query_names,
      straingr_summaries = straingr_compare.straingr_summary
  }
  output {
    File straingr_summary = straingr_gather_query.straingr_summary
  }
}