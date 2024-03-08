version 1.0

task strainge {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    File strainge_db
    String docker = "marcoteix/strainge:1.0.0"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 16
    Int max_strains = 5
  }
  command <<<
    # Decompress reference database
    mkdir database
    tar --zstd -C ./database/ -xf ~{strainge_db}
    fastas_dir="database/FASTA"
    similarities="database/similarities.tsv"
    hdf5=$(dir database/*.hdf5)
    echo $hdf5 > ~{samplename}_reference_used.txt
    /opt/conda/envs/strainge/bin/strainge --version > VERSION.txt
    /opt/conda/envs/strainge/bin/straingst kmerize -k ~{kmer_size} -o ~{samplename}_kmerized_reads.hdf5 ~{reads_1} ~{reads_2}
    echo run -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingst run -i ~{max_strains} -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingr prepare-ref -s ~{samplename}_straingst_results.strains.tsv -p "$fastas_dir/{ref}.fa.gz" \
        -S "$similarities" -o ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa index ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa mem -I 300 -t 2 ~{samplename}_refs_concat.fasta ~{reads_1} ~{reads_2} | /opt/conda/envs/strainge/bin/samtools sort -@ 2 -O BAM -o ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/samtools index ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/straingr call ~{samplename}_refs_concat.fasta ~{samplename}_straingr_alignment.bam --hdf5-out \
        ~{samplename}_straingr_variants.hdf5 --summary ~{samplename}_straingr.tsv --tracks all
  >>>
  output {
    File straingst_kmerized_reads = "~{samplename}_kmerized_reads.hdf5"
    File straingst_reference_db_used = "~{strainge_db}"
    File straingst_strains = "~{samplename}_straingst_results.strains.tsv"
    File straingst_statistics = "~{samplename}_straingst_results.stats.tsv"
    File straingr_concat_fasta = "~{samplename}_refs_concat.fasta"
    File straingr_read_alignment = "~{samplename}_straingr_alignment.bam"
    File straingr_variants = "~{samplename}_straingr_variants.hdf5"
    File straingr_report = "~{samplename}_straingr.tsv"
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

