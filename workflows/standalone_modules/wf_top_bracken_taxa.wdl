version 1.0

import "../../tasks/utilities/task_top_bracken_taxa.wdl" as top_bracken_taxa_task

workflow top_bracken_taxa_wf {
    input {
        File bracken_report
        Float min_abundance = 0.01
        Int memory = 8
        Int disk_size = 16
    }
    call top_bracken_taxa_task.top_bracken_taxa {
        input:
            bracken_report = bracken_report,
            min_abundance = min_abundance,
            memory = memory,
            disk_size = disk_size
    }
    output {
        String bracken_most_abundant_species = top_bracken_taxa.bracken_most_abundant_species
        String bracken_most_abundant_genus = top_bracken_taxa.bracken_most_abundant_genus
        String bracken_genera = top_bracken_taxa.bracken_genera
        String bracken_species = top_bracken_taxa.bracken_species
    }
}