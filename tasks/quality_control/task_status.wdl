version 1.0

task qc_flags_table {
  input {
    String table
    String workspace
    String docker = "marcoteix/gemstone-qc:1.0.0"
    Int disk_size = 10
  }
  command <<<

    cd /scripts/qc/bin
    python status.py \
      -t ~{table} \
      -w ~{workspace}

    echo "Completed! Wrote the status flags to the table ~{table} in the workspace ~{workspace}."
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