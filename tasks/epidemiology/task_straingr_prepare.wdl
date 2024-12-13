version 1.0

task straingr_prepare {
  input {
    String samplename
    File strainge_db
    File straingst_strains
    String docker = "marcoteix/strainge:1.0.1"
    Int disk_size = 32
    Int cpus = 1
    Int memory = 32
  }
  command <<<
    # Decompress reference database
    mkdir database
    tar --zstd -C ./database/ -xf ~{strainge_db}
    fastas_dir="database/FASTA"
    similarities="database/similarities.tsv"

    strainge --version > VERSION.txt

    # Prepare the reference for alignment
    echo "Preparing the reference sequence for alignment..."
    straingr prepare-ref \
        -s ~{straingst_strains} \
        -p "$fastas_dir/{ref}.fa.gz" \
        -S "$similarities" \
        -o ~{samplename}_refs_concat.fasta
  >>>
  output {
    File straingr_concat_fasta = "~{samplename}_refs_concat.fasta"
    File straingr_repetitiveness = "~{samplename}_refs_concat.meta.json"
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

task straingr_concatenate_references {
  input {
    Array[File] reference_fastas
    Int disk_size = 32
    Int cpus = 1
    Int memory = 32
  }
  command <<<
    # Concatenate input FASTAs
    echo "Concatenating input FASTAs into \"concatenated_ref.fasta\"..."
    cat ~{sep=' ' reference_fastas} > concatenated_ref.fasta
  >>>  
  output {
    File straingr_concat_fasta = "concatenated_ref.fasta"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "~{memory} GB"
    cpu: cpus
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task straingr_call {
  input {
    String samplename
    File reads_1
    File reads_2
    File reference_fasta
    Int insert_size = 300
    String docker = "marcoteix/strainge:1.0.1"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 64
  }
  command <<<

    strainge --version > VERSION.txt

    # Index the reference FASTA
    echo "Indexing the reference FASTA..."
    bwa index ~{reference_fasta}
    # Align
    echo "Aligning reads with bwa mem..."
    bwa mem \
        -I ~{insert_size} \
        -t ~{cpus} \
        ~{reference_fasta} \
        ~{reads_1} ~{reads_2} | \
        samtools sort -@ 2 -O BAM \
        -o ~{samplename}_straingr_alignment.bam

    # Index BAM file
    echo "Indexing BAM..."
    samtools index ~{samplename}_straingr_alignment.bam

    # Call variants
    echo "Calling variants..."
    straingr call ~{reference_fasta} \
        ~{samplename}_straingr_alignment.bam \
        --hdf5-out ~{samplename}_straingr_variants.hdf5 \
        --summary ~{samplename}_straingr.tsv
    
  >>>
  output {
    File straingr_variants_hdf5 = "~{samplename}_straingr_variants.hdf5"
    File straingr_summary = "~{samplename}_straingr.tsv"
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