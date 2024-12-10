version 1.0

task string_to_array {
  input {
    String string
  }
  command <<<
    echo ~{sep='\n' string} > ARRAY.txt
  >>>  
  output {
    Array[String] array = read_lines("ARRAY.txt")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    disk: "10 GB"
    maxRetries: 0
    preemptible: 0
  }
}