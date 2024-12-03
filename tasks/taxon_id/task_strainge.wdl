version 1.0

task strainge {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    File strainge_db
    String docker = "marcoteix/strainge:1.0.1"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 16
    Int max_strains = 5
    Boolean prepare_straingr = false
  }
  command <<<
    # Decompress reference database
    mkdir database
    tar --zstd -C ./database/ -xf ~{strainge_db}
    fastas_dir="database/FASTA"
    similarities="database/similarities.tsv"
    hdf5=$(dir database/*.hdf5)
    echo $hdf5 > ~{samplename}_reference_used.txt
    strainge --version > VERSION.txt
    straingst kmerize -k ~{kmer_size} -o ~{samplename}_kmerized_reads.hdf5 ~{reads_1} ~{reads_2}
    echo run -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    straingst run -i ~{max_strains} -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    if [ "~{prepare_straingr}" = true ]; then
    straingr prepare-ref -s ~{samplename}_straingst_results.strains.tsv -p "$fastas_dir/{ref}.fa.gz" \
        -S "$similarities" -o ~{samplename}_refs_concat.fasta
    bwa index ~{samplename}_refs_concat.fasta
    bwa mem -I 300 -t 2 ~{samplename}_refs_concat.fasta ~{reads_1} ~{reads_2} | samtools sort -@ 2 -O BAM -o ~{samplename}_straingr_alignment.bam
    samtools index ~{samplename}_straingr_alignment.bam
    straingr call ~{samplename}_refs_concat.fasta ~{samplename}_straingr_alignment.bam --hdf5-out \
        ~{samplename}_straingr_variants.hdf5 --summary ~{samplename}_straingr.tsv --tracks all
    fi
  >>>
  output {
    File straingst_kmerized_reads = "~{samplename}_kmerized_reads.hdf5"
    String straingst_reference_db_used = "~{strainge_db}"
    File straingst_strains = "~{samplename}_straingst_results.strains.tsv"
    File straingst_statistics = "~{samplename}_straingst_results.stats.tsv"
    File? straingr_concat_fasta = "~{samplename}_refs_concat.fasta"
    File? straingr_read_alignment = "~{samplename}_straingr_alignment.bam"
    File? straingr_variants = "~{samplename}_straingr_variants.hdf5"
    File? straingr_report = "~{samplename}_straingr.tsv"
    String strainge_docker = "~{docker}"    
    String strainge_version = read_string("VERSION.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpus
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task top_strain {
  input {
    Array[File] straingst_strains
    Float min_coverage = 0.8
    String docker = "marcoteix/gemstone-qc:1.0.0"
  }
  command <<<

    python3 << EOF
    import pandas as pd 
    paths = '~{sep=" " straingst_strains}'.split(" ")
    strains = pd.concat([pd.read_table(x, index_col=0) for x in paths])
    top_strain = strains[strains["cov"].ge(float("~{min_coverage}"))] \
      .sort_values("rapct", ascending=False).strain.iloc[0]
    with open("TOP_STRAIN.txt", "w") as file: file.write(top_strain)
    EOF

  >>>
  output {
    String straingst_top_strain = read_string("TOP_STRAIN.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    disk: "10 GB"
    maxRetries: 0
    preemptible: 0
  }
}