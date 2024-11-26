version 1.0

import "../../tasks/epidemiology/task_straingr_inputs.wdl" as straingr_inputs_task
import "../../tasks/epidemiology/task_straingr_prepare.wdl" as straingr_prepare_task


workflow straingr_prepare {
    input {
        String samplename
        String genus
        File reads_1
        File reads_2
        File strainge_db_config  # TSV containing database names and paths
        Array[String] straingst_selected_db
        Array[String] straingst_strains
        Int kmer_size
        Int insert_size = 300
        Int memory = 64
        Int cpus = 4
        Int disk_size = 100
    }
    call straingr_inputs_task.select_straingr_inputs {
        input:
            samplename = samplename,
            genus = genus,
            strainge_db_config = strainge_db_config,
            straingst_selected_db = straingst_selected_db,
            straingst_strains = straingst_strains
    }
    call straingr_prepare_task.straingr_prepare {
        input:
            samplename = samplename,
            reads_1 = reads_1,
            reads_2 = reads_2,
            kmer_size = kmer_size,
            strainge_db = select_straingr_inputs.selected_db,
            straingst_strains = select_straingr_inputs.straingst_strains,
            insert_size = insert_size,
            disk_size = disk_size,
            cpus = cpus,
            memory = memory
    }
    output {
        File straingr_variants_hdf5 = straingr_prepare.straingr_variants_hdf5
        File straingr_variants_vcf = straingr_prepare.straingr_variants_vcf
        String strainge_version = straingr_prepare.strainge_version
    }
}