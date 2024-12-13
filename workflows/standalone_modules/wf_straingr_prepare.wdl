version 1.0

import "../../tasks/epidemiology/task_straingr_prepare.wdl" as straingr_prepare_task


workflow straingr_prepare {
    input {
        String samplename
        File reads_1
        File reads_2
        Array[File] straingst_strains
        Array[File] straingst_selected_dbs
        Int insert_size = 300
        Int memory = 64
        Int cpus = 4
        Int disk_size = 100
    }
    scatter (i in range(length(straingst_strains))) {
        call straingr_prepare_task.straingr_prepare as prepare_reference {
            input:
                samplename = samplename,
                strainge_db = straingst_selected_dbs[i],
                straingst_strains = straingst_strains[i]
        }
    }
    call straingr_prepare_task.straingr_concatenate_references {
        input:
            reference_fastas = prepare_reference.straingr_concat_fasta
            
    }
    call straingr_prepare_task.straingr_call {
        input:
            samplename = samplename,
            reads_1 = reads_1,
            reads_2 = reads_2,
            reference_fasta = straingr_concatenate_references.straingr_concat_fasta,
            insert_size = insert_size,
            disk_size = disk_size,
            cpus = cpus,
            memory = memory
    }
    output {
        Array[File] straingr_repetitiveness = prepare_reference.straingr_repetitiveness
        File straingr_variants_hdf5 = straingr_call.straingr_variants_hdf5
        File straingr_summary = straingr_call.straingr_summary
        String strainge_version = straingr_call.strainge_version 
        String strainge_docker = straingr_call.strainge_docker
    }
}