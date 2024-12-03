version 1.0

task straingr_prepare {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    File strainge_db
    File straingst_strains
    Int insert_size = 300
    String docker = "marcoteix/strainge:1.0.1"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 64
  }
  command <<<
    # Decompress reference database
    mkdir database
    tar --zstd -C ./database/ -xf ~{strainge_db}
    fastas_dir="database/FASTA"
    similarities="database/similarities.tsv"

    base="/opt/conda/envs/strainge/bin"

    $base/strainge --version > VERSION.txt

    # Prepare the reference for alignment
    echo "Preparing the reference sequence for alignment..."
    $base/straingr prepare-ref \
        -s ~{straingst_strains} \
        -p "$fastas_dir/{ref}.fa.gz" \
        -S "$similarities" \
        -o ~{samplename}_refs_concat.fasta
    # Index the reference FASTA
    echo "Indexing the reference FASTA..."
    $base/bwa index ~{samplename}_refs_concat.fasta
    # Align
    echo "Aligning reads with bwa mem..."
    $base/bwa mem \
        -I ~{insert_size} \
        -t ~{cpus} \
        ~{samplename}_refs_concat.fasta \
        ~{reads_1} ~{reads_2} | \
        $base/samtools sort -@ 2 -O BAM \
        -o ~{samplename}_straingr_alignment.bam
    # Index BAM file
    echo "Indexing BAM..."
    $base/samtools index ~{samplename}_straingr_alignment.bam
    # Call variants
    echo "Calling variants..."
    $base/straingr call ~{samplename}_refs_concat.fasta \
        ~{samplename}_straingr_alignment.bam \
        --hdf5-out ~{samplename}_straingr_variants.hdf5 \
        --summary ~{samplename}_straingr.tsv
    
  >>>
  output {
    File straingr_concat_fasta = "~{samplename}_refs_concat.fasta"
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