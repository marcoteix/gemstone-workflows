version 1.0

struct Strainge_db_{
  String genus
  File hdf5
  Array[File] fastas # Path to all FASTA files. Forces Cromwell to mount these files when running the task
  String fastas_dir # Directory containing all reference FASTA files
  File similarities
}

task StrainGE_PE {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    String genus
    String escherichia_db_genus
    File escherichia_db_hdf5
    String escherichia_db_fastas_dir
    File escherichia_db_similarities
    Array[File] escherichia_db_fastas
    String pseudomonas_db_genus
    File pseudomonas_db_hdf5
    String pseudomonas_db_fastas_dir
    File pseudomonas_db_similarities
    Array[File] pseudomonas_db_fastas
    String proteus_db_genus
    File proteus_db_hdf5
    String proteus_db_fastas_dir
    File proteus_db_similarities
    Array[File] proteus_db_fastas
    String klebsiella_db_genus
    File klebsiella_db_hdf5
    String klebsiella_db_fastas_dir
    File klebsiella_db_similarities
    Array[File] klebsiella_db_fastas
    String staphylococcus_db_genus
    File staphylococcus_db_hdf5
    String staphylococcus_db_fastas_dir
    File staphylococcus_db_similarities
    Array[File] staphylococcus_db_fastas
    String acinetobacter_db_genus
    File acinetobacter_db_hdf5
    String acinetobacter_db_fastas_dir
    File acinetobacter_db_similarities
    Array[File] acinetobacter_db_fastas
    String enterococcus_db_genus
    File enterococcus_db_hdf5
    String enterococcus_db_fastas_dir
    File enterococcus_db_similarities
    Array[File] enterococcus_db_fastas
    String enterobacter_db_genus
    File enterobacter_db_hdf5
    String enterobacter_db_fastas_dir
    File enterobacter_db_similarities
    Array[File] enterobacter_db_fastas
    String docker = "marcoteix/strainge:0.0.2"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 16
  }
  command <<<
    match=false
    if [ ~{genus} == "Escherichia" ] || [ ~{genus} == "Shigella" ]; then 
      hdf5=~{escherichia_db_hdf5}
      similarities=~{escherichia_db_similarities}
      fastas=(~{sep=" " escherichia_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{pseudomonas_db_genus} ]; then 
      hdf5=~{pseudomonas_db_hdf5}
      fastas_dir=~{pseudomonas_db_fastas_dir}
      similarities=~{pseudomonas_db_similarities}
      fastas=(~{sep=" " pseudomonas_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{proteus_db_genus} ]; then 
      hdf5=~{proteus_db_hdf5}
      fastas_dir=~{proteus_db_fastas_dir}
      similarities=~{proteus_db_similarities}
      fastas=(~{sep=" " proteus_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{klebsiella_db_genus} ]; then 
      hdf5=~{klebsiella_db_hdf5}
      fastas_dir=~{klebsiella_db_fastas_dir}
      similarities=~{klebsiella_db_similarities}
      fastas=(~{sep=" " klebsiella_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{staphylococcus_db_genus} ]; then 
      hdf5=~{staphylococcus_db_hdf5}
      fastas_dir=~{staphylococcus_db_fastas_dir}
      similarities=~{staphylococcus_db_similarities}
      fastas=(~{sep=" " staphylococcus_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{acinetobacter_db_genus} ]; then 
      hdf5=~{acinetobacter_db_hdf5}
      fastas_dir=~{acinetobacter_db_fastas_dir}
      similarities=~{acinetobacter_db_similarities}
      fastas=(~{sep=" " acinetobacter_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{enterococcus_db_genus} ]; then 
      hdf5=~{enterococcus_db_hdf5}
      fastas_dir=~{enterococcus_db_fastas_dir}
      similarities=~{enterococcus_db_similarities}
      fastas=(~{sep=" " enterococcus_db_fastas})
      match=true
    fi
    if [ ~{genus} == ~{enterobacter_db_genus} ]; then 
      hdf5=~{enterobacter_db_hdf5}
      fastas_dir=~{enterobacter_db_fastas_dir}
      similarities=~{enterobacter_db_similarities}
      fastas=(~{sep=" " enterobacter_db_fastas})
      match=true
    fi
    if [ match = false ]; then
      exit 0
    fi
    fastas_dir=$(echo $(dirname ${fastas[0]}))
    echo $hdf5 > ~{samplename}_reference_used.txt
    /opt/conda/envs/strainge/bin/strainge --version > VERSION.txt
    /opt/conda/envs/strainge/bin/straingst kmerize -k ~{kmer_size} -o ~{samplename}_kmerized_reads.hdf5 ~{reads_1} ~{reads_2}
    /opt/conda/envs/strainge/bin/straingst run -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingr prepare-ref -s ~{samplename}_straingst_results.strains.tsv -p "$fastas_dir/{ref}" \
        -S "$similarities" -o ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa index ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa mem -I 300 -t 2 ~{samplename}_refs_concat.fasta ~{reads_1} ~{reads_2} | /opt/conda/envs/strainge/bin/samtools sort -@ 2 -O BAM -o ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/samtools index ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/straingr call ~{samplename}_refs_concat.fasta ~{samplename}_straingr_alignment.bam --hdf5-out \
        ~{samplename}_straingr_variants.hdf5 --summary ~{samplename}_straingr.tsv --tracks all
  >>>
  output {
    File straingst_kmerized_reads = "~{samplename}_kmerized_reads.hdf5"
    String straingst_reference_db_used = read_string("~{samplename}_reference_used.txt")
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

task strainge {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    File strainge_db
    String docker = "marcoteix/strainge:0.0.2"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 16
  }
  command <<<
    # Decompress reference database
    mkdir database
    tar -C ./database/ -xzvf ~{strainge_db}
    fastas_dir="database/fasta"
    similarities="database/similarities.tsv"
    hdf5=$(dir database/*.hdf5)
    echo $hdf5 > ~{samplename}_reference_used.txt
    /opt/conda/envs/strainge/bin/strainge --version > VERSION.txt
    /opt/conda/envs/strainge/bin/straingst kmerize -k ~{kmer_size} -o ~{samplename}_kmerized_reads.hdf5 ~{reads_1} ~{reads_2}
    echo run -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingst run -O -o ~{samplename}_straingst_results $hdf5 ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingr prepare-ref -s ~{samplename}_straingst_results.strains.tsv -p "$fastas_dir/{ref}" \
        -S "$similarities" -o ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa index ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa mem -I 300 -t 2 ~{samplename}_refs_concat.fasta ~{reads_1} ~{reads_2} | /opt/conda/envs/strainge/bin/samtools sort -@ 2 -O BAM -o ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/samtools index ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/straingr call ~{samplename}_refs_concat.fasta ~{samplename}_straingr_alignment.bam --hdf5-out \
        ~{samplename}_straingr_variants.hdf5 --summary ~{samplename}_straingr.tsv --tracks all
  >>>
  output {
    File straingst_kmerized_reads = "~{samplename}_kmerized_reads.hdf5"
    String straingst_reference_db_used = read_string("~{samplename}_reference_used.txt")
    File straingst_strains = "~{samplename}_straingst_results.strains.tsv"
    File straingst_statistics = "~{samplename}_straingst_results.stats.tsv"
    File straingr_concat_fasta = "~{samplename}_refs_concat.fasta"
    File straingr_read_alignment = "~{samplename}_straingr_alignment.bam"
    File straingr_variants = "~{samplename}_straingr_variants.hdf5"
    File straingr_report = "~{samplename}_straingr.tsv"
    String strainge_docker = "~{docker}"    String strainge_version = read_string("VERSION.txt")
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