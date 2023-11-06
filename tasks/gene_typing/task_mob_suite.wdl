task mob_recon {
  input {
    File assembly
    String samplename
    String docker = "kbessonov/mob_suite:3.0.3"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # version capture
    mob_recon --version | cut -d ' ' -f2 > VERSION.txt

    mkdir mob_recon

    # unzip assembly FASTA
    recompress=false
    assembly=~{assembly}
    if [[ ~{assembly} == *.gz ]]; then 
      gzip -d ~{assembly}
      recompress=true
      assembly=${assembly: 0:-3}
    fi

    # run mob-recon
    mob_recon --infile $assembly --outdir mob_recon/~{samplename}

    # If the assembly FASTA was originally compressed, compress it again
    if $recompress; then
      gzip $assembly

    # Compress output FASTAs
    for fasta in mob_recon/~{samplename}/*.fasta; do 
      gzip $fasta
    done
  >>>
  output {
    File mob_recon_results = "mob_recon/~{samplename}/contig_report.txt"
    File chromosome_fasta = "mob_recon/~{samplename}/chromosome.fasta.gz"
    String mob_recon_version = read_string("VERSION.txt")
    String mob_recon_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}