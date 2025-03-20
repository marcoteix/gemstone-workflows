version 1.0

task filter_variants {
  input {
    Array[File] vcfs
    String samplename
  }
  command <<<

    # Write a list of input VCF paths 
    echo ~{sep='\n' vcfs} > filelist.txt

    bcftools 
    ./vcf2phylip.py 

  >>>
  output {
    File filtered_vcf = "~{samplename}/~{samplename}.variants.vcf"
  }
  runtime {
    docker: "marcoteix/gemstone-utils:1.0.0"
    memory: "8 GB"
    cpu: 1
    disks:  "local-disk 16" + " SSD"
    disk: "16 GB"
    preemptible: 0
  }
}