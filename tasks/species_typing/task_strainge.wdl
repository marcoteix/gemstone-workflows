version 1.0

struct Strainge_db = {
  String genus,
  File hdf5,
  Array[File] fastas, # Path to all FASTA files. Forces Cromwell to mount these files when running the task
  File fastas_dir, # Directory containing all reference FASTA files
  File similarities
}

task StrainGE_PE {
  input {
    String samplename
    File reads_1
    File reads_2
    Int kmer_size
    Strainge_db straingst_reference_db
    String docker = "marcoteix/strainge:0.0.2"
    Int disk_size = 100
    Int cpus = 4
    Int memory = 16
  }
  parameter_meta {
    samplename: "Sample ID."
    reads_1: "Input file containing clean reads."
    reads_2: "Input file containing clean reads."
    kmer_size: "K-mer sizes used to k-merize the input reads. Should match the value used in the construction of the reference database."
    straingst_reference_db: "Strainge_db struct with the target genus name, the path to the HDF5 file, an array of reference FASTA files, the directory containing the reference FASTAS, and the similarities TSV."
    straingst_kmerized_reads: "HDF5 file containing the k-merized input reads."
    straingst_reference_db_used: "HDF5 file containing the StrainGST reference database used."
    straingst_strains: "TSV file with the strains detected by StrainGST."
    straingst_statistics: "TSV file with StrainGST sample statistics."
    straingr_concat_fasta: "Concatenated FASTA file of all representative sequences in the StrainGST reference database."
    straingr_read_alignment: "BAM file with reads aligned to the closest reference."
    straingr_variants: "HDF5 file with variants detected by StrainGR."
    straingr_report: "Human readable TSV file with a summary of StrainGR results."
    strainge_docker: "StrainGE docker image."
    strainge_version: "StrainGE version."
    straingst_reference_similarities: "TSV with similarities between the sequences in the StrainGST reference database."
  }
  command <<<
    /opt/conda/envs/strainge/bin/strainge --version > VERSION.txt
    /opt/conda/envs/strainge/bin/straingst kmerize -k ~{kmer_size} -o ~{samplename}_kmerized_reads.hdf5 ~{reads_1} ~{reads_2}
    /opt/conda/envs/strainge/bin/straingst run -O -o ~{samplename}_straingst_results ~{straingst_reference_db.hdf5} ~{samplename}_kmerized_reads.hdf5
    /opt/conda/envs/strainge/bin/straingr prepare-ref -s ~{samplename}_straingst_results.strains.tsv -p "~{straingst_reference_db.fastas_dir}/{ref}" \
        -S ~{straingst_reference_db.similarities} -o ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa index ~{samplename}_refs_concat.fasta
    /opt/conda/envs/strainge/bin/bwa mem -I 300 -t 2 ~{samplename}_refs_concat.fasta ~{reads_1} ~{reads_2} | /opt/conda/envs/strainge/bin/samtools sort -@ 2 -O BAM -o ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/samtools index ~{samplename}_straingr_alignment.bam
    /opt/conda/envs/strainge/bin/straingr call ~{samplename}_refs_concat.fasta ~{samplename}_straingr_alignment.bam --hdf5-out \
        ~{samplename}_straingr_variants.hdf5 --summary ~{samplename}_straingr.tsv --tracks all
  >>>
  output {
    File straingst_kmerized_reads = "~{samplename}_kmerized_reads.hdf5"
    File straingst_reference_db_used = "~{straingst_reference_db.hdf5}"
    String straingst_reference_genus_used = "~{straingst_reference_db.genus}" 
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
    maxRetries: 1
    preemptible: 0
  }
}

# Return the Strainge_db struct that should be used with StrainGE
task select_reference_db {
  input {
    String samplename
    String gambit_taxonomy
    Strainge_db escherichia_db
    Strainge_db pseudomonas_db
    Strainge_db proteus_db
    Strainge_db klebsiella_db
    Strainge_db staphylococcus_db
    Strainge_db acinetobacter_db
    Strainge_db enterococcus_db
    Strainge_db enterobacter_db
  }
  command <<<

    mkdir ~{samplename}

    python3 <<CODE

    gambit_taxonomy = ~{gambit_taxonomy}
    gambit_taxonomy = gambit_taxonomy.replace('_', ' ').split(' ')[0].capitalize()       
    # Write to TXT file
    with open(f"~{samplename}/GAMBIT_taxonomy.txt", "w") as f: f.write(gambit_taxonomy)

    CODE

    gambit_taxonomy=$(cat ~{samplename}/GAMBIT_taxonomy.txt)

    if [ gambit_taxonomy=='Escherichia' ]; then
      json_file=write_json(~{escherichia_db})
    fi
    if [ gambit_taxonomy=='Shigella' ]; then
      json_file=write_json(~{escherichia_db})
    fi
    if [ gambit_taxonomy=='Pseudomonas' ]; then
      json_file=write_json(~{pseudomonas_db})
    fi
    if [ gambit_taxonomy=='Proteus' ]; then
      json_file=write_json(~{proteus_db})
    fi
    if [ gambit_taxonomy=='Klebsiella' ]; then
      json_file=write_json(~{klebsiella_db})
    fi
    if [ gambit_taxonomy=='Staphylococcus' ]; then
      json_file=write_json(~{staphylococcus_db})
    fi
    if [ gambit_taxonomy=='Klebsiella' ]; then
      json_file=write_json(~{klebsiella_db})
    fi
    if [ gambit_taxonomy=='Acinetobacter' ]; then
      json_file=write_json(~{acinetobacter_db})
    fi
    if [ gambit_taxonomy=='Enterobacter' ]; then
      json_file=write_json(~{enterobacter_db})
    fi
    if [ gambit_taxonomy=='Enterococcus' ]; then
      json_file=write_json(~{enterococcus_db})
    fi
    echo json_file > "~{samplename}/JSON_PATH.txt"
  >>>
  output {
    String selected_strainge_db_json_path = read_string("~{samplename}/JSON_PATH.txt")
    Strainge_db selected_strainge_db = read_json(read_string("~{samplename}/JSON_PATH.txt"))
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}