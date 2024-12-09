version 1.0

import "../../tasks/epidemiology/task_straingr_compare.wdl" as straingr_compare_task
import "../../tasks/epidemiology/task_straingr_matrix.wdl" as straingr_matrix

workflow straingr_compare {
    meta {
        description: "Uses StrainGR to compare one sample against a set of samples given HDF5 files from StrainGR prepare."
    }
    input {
        Array[String] samplenames
        Array[File] straingr_variants_hdf5
        Int disk_size = 8
        Int memory = 16
    }
    scatter (i in range(length(samplenames))) {
        call straingr_compare_task.straingr_compare_query {
            input:
                reference_name = samplenames[i],
                query_names = samplenames,
                reference_variants = straingr_variants_hdf5[i],
                query_variants = straingr_variants_hdf5,
                disk_size = disk_size,
                memory = memory
        }
    }
    call straingr_matrix.straingr_gather_pairwise {
        input:
            straingr_summaries = straingr_compare_query.straingr_summary
    }
    output {
        File straingr_summary = straingr_gather_pairwise.straingr_summary
    }
}