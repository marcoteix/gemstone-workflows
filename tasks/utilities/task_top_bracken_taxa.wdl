version 1.0

task top_bracken_taxa {
  input {
    File bracken_report
    Float min_abundance = 0.01
    Int memory = 8
    Int disk_size = 16
  }
  command <<<
        
        python3 <<CODE

        import pandas as pd

        bracken_file = "~{bracken_report}"
        min_abundance = float("~{min_abundance}")

        # Read Bracken report
        bracken = pd.read_csv(
            bracken_file,
            sep = "\t"
        )

        # Get species with an abundance > min_abundance
        species = bracken[ 
            bracken.taxonomy_lvl.eq("S") & \
            bracken.fraction_total_reads.ge(min_abundance)
        ].name

        species = "/".join(species.values)

        # Get the most abundant species
        top_species = bracken[ 
            bracken.taxonomy_lvl.eq("S")
        ].sort_values(
            "fraction_total_reads",
            ascending = False
        ).name.iloc[0]

        # Now do the same for genera
        genera_abundances = bracken.assign(
            genus = bracken.name.apply(
                lambda x: x.split(" ")[0]
            )
        ).groupby(
            "genus",
            as_index = False
        ).fraction_total_reads \
        .sum()

        genera = genera_abundances.loc[ 
            genera_abundances.fraction_total_reads.ge(min_abundance),
            "genus"
        ]

        genera = "/".join(genera.values)

        top_genus = genera_abundances.sort_values(
            "fraction_total_reads", 
            ascending = False 
        ).genus.iloc[0]

        for file, content in zip(
            ["species.txt", "genera.txt", "top_species.txt", "top_genus.txt"],
            [species, genera, top_species, top_genus]
        ):
            with open(file, "w") as f:
                f.write(content)
        
        CODE

  >>>
  output {
    String bracken_most_abundant_species = read_string("top_species.txt")
    String bracken_most_abundant_genus = read_string("top_genera.txt")
    String bracken_genera = read_string("genera.txt")
    String bracken_species = read_string("species.txt")
  }
  runtime {
    docker: "marcoteix/gemstone-utils:1.0.0"
    memory: memory + " GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}