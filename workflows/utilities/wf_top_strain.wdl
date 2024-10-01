version 1.0

import "../../../tasks/taxon_id/task_strainge.wdl" as strainge_task

workflow top_strain {
    meta {
        description: "Takes an array of StrainGST output TSV files and finds the strain with highest relative abundance."
    }
    input {
        Array[File] straingst_strains
        Float min_coverage = .8
    }
    call strainge_task.top_strain {
        input:
          straingst_strains = straingst_strains,
          min_coverage = min_coverage
    }
    output {
        String straingst_top_strain = top_strain.straingst_top_strain
    }
}