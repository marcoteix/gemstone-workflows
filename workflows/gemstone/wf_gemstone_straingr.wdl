version 1.0

import "../standalone_modules/wf_straingr_prepare.wdl" as straingr_prepare_wf
import "../../tasks/epidemiology/task_straingr_compare.wdl" as straingr_compare
import "../../tasks/epidemiology/task_straingr_matrix.wdl" as straingr_matrix


workflow straingr {
    input {
        Array[String] samplenames
        String genus
        Array[File] reads_1
        Array[File] reads_2
        File strainge_db_config
        Array[String] straingst_selected_dbs
        Array[String] straingst_strains
        Int kmer_size
        Int insert_size = 300
        Int straingr_prepare_memory = 64
        Int straingr_prepare_cpus = 4
        Int straingr_prepare_disk_size = 100
        Int straingr_compare_memory = 8
        Int straingr_compare_disk_size = 16
    }
    scatter (i in range(length(samplenames))) {
        call straingr_prepare_wf.straingr_prepare {
            input:
                samplename = samplenames[i],
                genus = genus,
                reads_1 = reads_1[i],
                reads_2 = reads_2[i],
                strainge_db_config = strainge_db_config,
                straingst_selected_db = straingst_selected_dbs[i],
                straingst_strains = straingst_strains[i],
                kmer_size = kmer_size,
                insert_size = insert_size,
                memory = straingr_prepare_memory,
                cpus = straingr_prepare_cpus,
                disk_size = straingr_prepare_disk_size
        }
    }
    scatter (i in range(length(samplenames))) {
        call straingr_compare.straingr_compare_query {
            input:
                reference_name = samplenames[i],
                query_names = samplenames,
                reference_variants = straingr_prepare.straingr_variants_hdf5[i],
                query_variants = straingr_prepare.straingr_variants_hdf5,
                disk_size = straingr_compare_disk_size,
                memory = straingr_compare_memory
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