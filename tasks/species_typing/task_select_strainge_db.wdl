version 1.0

import "task_strainge.wdl"

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
    Int disk_size = 20
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