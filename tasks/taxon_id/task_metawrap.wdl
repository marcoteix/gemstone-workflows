version 1.0

task metawrap {
  input {
    File assembly_fasta
    File read1
    File? read2
    String samplename
    String docker = "quay.io/biocontainers/metawrap:1.2--hdfd78af_2"
    Int metawrap_completion = 80
    Int metawrap_contamination = 10
    Int metawrap_min_contig_length = 1000
    File ncbi_nt_database
    File ncbi_taxonomy_database
    File checkm_database
    Int mem = 64
    Int cpu = 8
    Int disk_size = 100
  }
  command <<<
    echo $(metawrap --version) | sed 's/^.*metaWRAP v=//;s/ .*$//' | tee METAWRAP_VERSION
    date | tee DATE

    entrypoint=$(pwd)

    # Configure databases
    mkdir databases && cd databases
    mkdir checkm
    echo "$(date) - Decompressing the CheckM database..."
    tar -C ./checkm/ -xzvf ~{checkm_database}
    echo "$(date) - Setting CheckM data root to $(pwd)..."
    checkm data setRoot $(pwd)

    mkdir ~/NCBI_TAX_DB
    echo "$(date) - Decompressing the NCBI taxonomy database..."
    tar -C ~/NCBI_TAX_DB/ -xzvf ~{ncbi_taxonomy_database}
    
    mkdir ~/NCBI_NT_DB
    echo "$(date) - Decompressing the NCBI nt BLAST database..."
    tar -C ~/NCBI_NT_DB/ -xzvf ~{ncbi_nt_database}
    for a in nt.*.tar.gz; do tar -xzf $a; done
    echo "$(date) - Finished setting up databases. Now to the interesting part..."

    cd $entrypoint

    mkdir ~{samplename} && cd ~{samplename}
    binning_out=~{samplename}/metawrap_binning
    refinement_out=~{samplename}/metawrap_refinement
    classification_out=~{samplename}/metawrap_classification

    echo "$(date) - Creating initial bins with metabat2, maxbin2, and concoct. This may take a while."
    metawrap binning -o $binning_out -t ~{cpu} -m ~{mem} -a ~{assembly_fasta} --metabat2 --maxbin2 --concoct \
      -l ~{metawrap_min_contig_length} ~{read1} ~{read2}
    echo "$(date) - Refining bins..."
    metawrap bin_refinement -o $refinement_out -t ~{cpu} -m ~{mem} -A $binning_out/metabat2_bins/ \
      -B $binning_out/maxbin2_bins/ -C $binning_out/concoct_bins/ -c ~{metawrap_completion} \
      -x ~{metawrap_contamination}

    # Get number of bins meeting the contamination and completeness thresholds
    echo $(cat $refinement_out/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins.stats | \
      awk '$2>~{metawrap_completion} && $3<~{metawrap_contamination}' | wc -l) > N_BINS

    echo "$(date) - Classifying bins..."
    metawrap classify_bins -b $refinement_out/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins \
      -o $classification_out -t ~{cpu}

    echo "$(date) - Done!"

  >>>
  output {
    String metawrap_docker = '~{docker}'
    String metawrap_version = '~{docker}'
    String analysis_date = read_string("DATE")
    File metawrap_insert_sizes = "~{samplename}/metawrap_binning/insert_sizes.txt"
    File metawrap_stats = "~{samplename}/metawrap_refinement/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins.stats"
    Int metawrap_n_bins = read_string("~{samplename}/N_BINS")
    File metawrap_taxonomy = "~{samplename}/metawrap_classification/bin_taxonomy.tab"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{mem} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      preemptible: 0
  }
}

task metawrap_binning {
  input {
    File assembly_fasta
    File read1
    File? read2
    String samplename
    String docker = "quay.io/biocontainers/metawrap:1.2--hdfd78af_2"
    Int metawrap_completion = 80
    Int metawrap_contamination = 10
    Int metawrap_min_contig_length = 1000
    File checkm_database
    Int mem = 32
    Int cpu = 16
    Int disk_size = 100
    # Tools used in metawrap
    String binning_flags = "--metabat2 --maxbin2 --concoct"
  }
  command <<<
    echo $(metawrap --version) | sed 's/^.*metaWRAP v=//;s/ .*$//' | tee METAWRAP_VERSION
    date | tee DATE

    entrypoint=$(pwd)

    # Configure databases
    mkdir databases && cd databases
    mkdir checkm
    echo "$(date) - Decompressing the CheckM database..."
    tar -C ./checkm/ -xzvf ~{checkm_database}
    echo "$(date) - Setting CheckM data root to $(pwd)..."
    # Run the command twice because it fails if you do it just once, for some reason
    echo checkm | checkm data setRoot
    echo checkm | checkm data setRoot
    echo "$(date) - Finished setting up databases. Unzipping reads..."
    # Decompress reads. Rename it to end with _R1.fastq and _R2.fastq
    cd $entrypoint && mkdir ~{samplename} && mkdir ~{samplename}/reads
    r1=~{samplename}/reads/~{samplename}_1.fastq
    r2=~{samplename}/reads/~{samplename}_2.fastq
    gzip -dc ~{read1} > $r1
    gzip -dc ~{read2} > $r2

    cd $entrypoint

    echo "CWD: $(pwd)"
    binning_out=$entrypoint/~{samplename}/metawrap_binning
    refinement_out=$entrypoint/~{samplename}/metawrap_refinement

    echo "$(date) - Creating initial bins with metabat2, maxbin2, and concoct. This may take a while."
    metawrap binning -o $binning_out -t ~{cpu} -m ~{mem} -a ~{assembly_fasta} ~{binning_flags} \
      -l ~{metawrap_min_contig_length} $r1 $r2
    echo "$(date) - Refining bins..."
    metawrap bin_refinement -o $refinement_out -t ~{cpu} -m ~{mem} -A $binning_out/metabat2_bins/ \
      -B $binning_out/maxbin2_bins/ -C $binning_out/concoct_bins/ -c ~{metawrap_completion} \
      -x ~{metawrap_contamination}

    # Get number of bins meeting the contamination and completeness thresholds
    echo $(cat $refinement_out/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins.stats | \
      awk '$2>~{metawrap_completion} && $3<~{metawrap_contamination}' | wc -l) > ~{samplename}/N_BINS

    # Write output bin FASTA locations
    for fasta in $refinement_out/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins/*; do
      echo $fasta >> ~{samplename}/BINS
    done

    echo $binning_out : $(ls $binning_out)
    echo $refinement_out : $(ls $refinement_out)
    echo metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins : $(ls $refinement_out/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins)

    echo "$(date) - Done!"

  >>>
  output {
    String metawrap_docker = '~{docker}'
    String metawrap_version = read_string("METAWRAP_VERSION")
    String analysis_date = read_string("DATE")
    String metawrap_binning_flags = "~{binning_flags}"
    File metawrap_stats = "~{samplename}/metawrap_refinement/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins.stats"
    Int metawrap_n_bins = read_string("~{samplename}/N_BINS")
    Array[File] metawrap_fasta = glob("~{samplename}/metawrap_refinement/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins/*")
    File metawrap_contigs = "~{samplename}/metawrap_refinement/metawrap_~{metawrap_completion}_~{metawrap_contamination}_bins.contigs"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{mem} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      preemptible: 0
  }
}