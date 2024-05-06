version 1.0

task select_reference_db {
  input {
    String samplename
    String gambit_taxonomy
    Array [File] strainge_dbs  # Must be filenames of the type [genus].tar.gz
    Int disk_size = 128
  }
  command <<<

    mkdir ~{samplename} && cd ~{samplename}
    mkdir databases

    counter=0

    for db in ~{sep=" " strainge_dbs}; do
      # Decompress the database
      echo $(date) - Decompressing the StrainGE database at $db...
      tar -C ./databases/ -xzvf $db

      # Compare reference db name against the assigned genus
      echo $(date) - Comparing database genus with the sample assigned genus...
      python3 << EOF

  import os 

  gambit_taxonomy = "~{gambit_taxonomy}"
  gambit_taxonomy = gambit_taxonomy.replace('_', ' ').split(' ')[0].capitalize()       
  # Write to TXT file
  with open("GAMBIT_taxonomy.txt", "w") as f: f.write(gambit_taxonomy)

  with open("databases/genus.txt") as f:
    db_genus = f.read().replace('_', ' ').replace('\n', '').split('/')
  db_genus = [x.capitalize() for x in db_genus]
  with open("SEARCH_RESULT", "w") as f:
    f.write(f'{gambit_taxonomy in db_genus}')
  EOF
      
      # Check comparison results
      result=$(cat SEARCH_RESULT)
      if [ $result = True ]; then
        echo $counter > SELECTED_DB.txt
        echo $(date) - Found matching database.
        break
      fi

      echo $(date) - Removing intermediate files...
      rm databases/fasta/*
      rmdir databases/fasta
      rm databases/*
      echo $(date) - Analyzing the next database...
      counter=$(($counter + 1))
    done

    echo $(date) - Done!
    cd ..
  >>>
  output {
    Int selected_db = read_int("~{samplename}/SELECTED_DB.txt")
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

task select_reference_db_lite {
  input {
    String samplename
    String predicted_taxonomy
    File strainge_db_config  # TSV containing database names and paths
    Int disk_size = 4
  }
  command <<<

    mkdir ~{samplename} && cd ~{samplename}

    # Parse the predicted genus
    echo $(date) - Parsing the predicted genus...

    python3 << EOF

  import os 

  taxonomy = "~{predicted_taxonomy}"
  # If there are multiple taxa separated by /, include all
  taxonomy = [x.replace('_', ' ').split(' ')[0].capitalize() for x in taxonomy.split("/")]
  # Write to TXT file
  with open("TAXON", "w") as f: f.write("\n".join(taxonomy))

  # Iterate over the TSV with the strainge_dbs and search for a match
  matched = False
  with open("~{strainge_db_config}") as f:
    with open("DATABASE", "w") as out_f: out_f.write("")
    for line in f:
        db_taxon, db = line.replace("\n", "").split("\t")
        if db_taxon in taxonomy:
          with open("DATABASE", "a") as out_f: out_f.write(db+"\n")
          matched = True
  
  # Check if there was a match
  with open("MATCH", "w") as f: f.write(str(matched).lower())

  EOF
  >>>
  output {
    Array[String] selected_db = read_lines("~{samplename}/DATABASE")
    Boolean found_db = read_boolean("~{samplename}/MATCH")
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