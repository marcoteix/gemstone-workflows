version 1.0

task select_straingr_inputs {
  input {
    String samplename
    String genus
    File strainge_db_config  # TSV containing database names and paths
    String straingst_selected_db
    String straingst_strains_tsvs
    Int disk_size = 4
  }
  command <<<

    mkdir ~{samplename} && cd ~{samplename}

    # Parse the predicted genus
    echo $(date) - Searching for the ~{genus} database...

    python3 << EOF

  import os 

  # Iterate over the TSV with the strainge_dbs and search for a match
  matched = False
  with open("~{strainge_db_config}") as f:
    with open("DATABASE", "w") as out_f: out_f.write("")
    for line in f:
        db_taxon, db = line.replace("\n", "").split("\t")
        if db_taxon.lower() == "~{genus}".lower():
          with open("DATABASE", "a") as out_f: out_f.write(db+"\n")
          matched = True
          break
  
  if not matched: print(f"Failed to find a matching database for ~{genus}.")
  # Check if there was a match
  with open("MATCH", "w") as f: f.write(str(matched).lower())

  # Select the StrainGST report file matching the desired genus
  if matched:
    found_report = False
    straingst_dbs = [x.removepreffix("[").removesuffix("]") 
      for x in "~{straingst_selected_db}".split(" ")]
    strains_tsvs = [x.removepreffix("[").removesuffix("]") 
      for x in "~{straingst_strains}".split(" ")]
    for a, b in zip(straingst_dbs, strains_tsvs):
      print(a, b)
      if a == db: 
        with open("STRAINGST", "a") as out_f: out_f.write(b+"\n")
        found_report = True
        break

    if not found_report: print(f"Failed to find a matching StrainGST report.")
  EOF
  >>>
  output {
    String selected_db = read_string("~{samplename}/DATABASE")
    Boolean found_db = read_boolean("~{samplename}/MATCH")
    String straingst_strains = read_string("~{samplename}/STRAINGST")
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