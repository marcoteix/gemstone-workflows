version 1.0

import "../../../tasks/quality_control/task_taxonomy_qc.wdl" as taxonomy_qc_task
import "../../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../../tasks/quality_control/task_checkm2.wdl" as checkm2_task

workflow test_qc {
    input {
        String samplename = "test_sample"
        String lab_determined_genus = "Escherichia"
        String gambit_taxonomy = "Escherichia coli"
        String read_qc_flag = "PASS"
        File assembly
        File? qc_check_table
        Int mem = 8
        Int cpu = 2
        Int disk_size = 16
        Float contamination_threshold = 0.01
        Float checkm2_contamination = 0.001
    }
    if (read_qc_flag == "PASS") {
        call taxonomy_qc_task.taxonomy_qc_check {
            input:
                samplename = samplename,
                gambit_taxonomy = gambit_taxonomy,
                checkm2_contamination = checkm2_contamination,
                lab_determined_genus = lab_determined_genus,
                contamination_threshold = contamination_threshold
        }
        if (defined(qc_check_table)) {
            call qc_check.qc_check_phb as qc_check_task {
                input:
                    disk_size = disk_size
            }
        }
    } 
    call qc_check.qc_flags_check {
        input:
            qc_taxonomy_flag = taxonomy_qc_check.qc_check,
            qc_check_table_flag = qc_check_task.qc_check,
            qc_clean_reads_flag = read_qc_flag
    }
    output {
        # Taxon QC Results
        String? qc_taxonomy_check = taxonomy_qc_check.qc_check
        File? qc_taxonomy_report = taxonomy_qc_check.qc_report
        # QC_Check Results
        String? read_and_table_qc_check = qc_check_task.qc_check
        File? read_and_table_qc_standard = qc_check_task.qc_standard
        String? qc_check = qc_flags_check.qc_check
        String? qc_note = qc_flags_check.qc_note
    }
}