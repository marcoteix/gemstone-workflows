version 1.0

task qc_flags_row {
  input {
    String? raw_read_screen
    String? clean_read_screen
    Float? est_coverage_clean
    String? species
    String? gambit_predicted_taxon
    Float? checkm2_completeness
    Float? checkm2_contamination
    Float min_coverage = 40.0
    Float max_contamination = 2.0
    Float min_completeness = 80.0
    String docker = "marcoteix/gemstone-qc:1.0.0"
    Int disk_size = 10
  }
  command <<<

    cd /scripts/qc/bin
    python row.py \
      -rs ~{raw_read_screen} \
      -cs ~{clean_read_screen} \
      -x ~{est_coverage_clean} \
      -ls ~{species} \
      -g ~{gambit_predicted_taxon} \
      -c ~{checkm2_contamination} \
      -C ~{checkm2_completeness} \
      -mx ~{min_coverage} \
      -Mc ~{max_contamination} \
      -mc ~{min_completeness} \
      -o "~/qc_results"

  >>>
  output {
    String qc_check = read_string("~/qc_results/qc_check")
    String qc_note = read_string("~/qc_results/qc_note")
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

task qc_flags_table {
  input {
    String table
    String workspace
    Float min_coverage = 40.0
    Float max_contamination = 2.0
    Float min_completeness = 80.0
    String docker = "marcoteix/gemstone-qc:1.0.0"
    Int disk_size = 10
  }
  command <<<

    cd /scripts/qc/bin
    python table.py \
      -t ~{table} \
      -w ~{workspace} \
      -mx ~{min_coverage} \
      -Mc ~{max_contamination} \
      -mc ~{min_completeness}

  >>>
  output {
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}