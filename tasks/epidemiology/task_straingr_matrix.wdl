version 1.0

task straingr_gather_query {
  meta {
    description: "Concatenates a set of StrainGR summary files."
  }
  input {
    String reference_name
    Array[String] query_names
    Array[File] straingr_summaries
    Int disk_size = 8
    Int memory = 8
  }
  command <<<

    mkdir ~{reference_name} && cd ~{reference_name}

    python3 << EOF

  import pandas as pd

  queries = "~{sep=' ' query_names}".split(" ")
  summaries = "~{sep=' ' straingr_summaries}".split(" ")

  # Concatenate summary files, excluding the TOTAL line
  X = pd.concat(
    {
      k: pd.read_csv(v, index_col=0, sep="\t").assign(query_id=k).drop("TOTAL").reset_index()
      for k,v in zip(queries, summaries)
    }
  ).assign(reference_id="~{reference_name}")

  print("Writing the concatenated StrainGR summary to ~{reference_name}/concat_summary.tsv")
  X.to_csv("concat_summary.tsv", sep="\t")

  EOF
  >>>
  output {
    File straingr_summary = "~{reference_name}/concat_summary.tsv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "~{memory} GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task straingr_gather_pairwise {
  meta {
    description: "Concatenates a set of concatenated StrainGR summary files."
  }
  input {
    Array[File] straingr_summaries
    Int disk_size = 8
    Int memory = 8
  }
  command <<<

    python3 << EOF

  import pandas as pd

  files = "~{sep=' ' straingr_summaries}".split(" ")

  # Concatenate summary files
  X = pd.concat(
    [pd.read_csv(file, sep="\t") for file in files]
  )

  print("Writing the concatenated StrainGR summary to \"concat_summary.tsv\"...")
  X.to_csv("concat_summary.tsv", sep="\t")

  EOF
  >>>
  output {
    File straingr_summary = "concat_summary.tsv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "~{memory} GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}