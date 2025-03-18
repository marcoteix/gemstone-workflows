version 1.0

task filter_variants {
  input {
    File vcf
    String samplename
    String docker = "marcoteix/cleansweep:0.2"
  }
  command <<<

    mkdir -p ~{samplename}

    bcftools view  \
        ~{vcf} \
        -o ~{samplename}/~{samplename}.variants.vcf \
        -O v \
        -e 'INFO/AC = 0' \
        -f ".,PASS"

  >>>
  output {
    File filtered_vcf = "~{samplename}/~{samplename}.variants.vcf"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
    disks:  "local-disk 16" + " SSD"
    disk: "16 GB"
    preemptible: 0
  }
}