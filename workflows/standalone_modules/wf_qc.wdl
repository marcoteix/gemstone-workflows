version 1.0

import "../../../tasks/quality_control/task_qc_flags.wdl" as task_qc_flags

workflow test_row_qc {
    meta {
        description: "Generates a qc_check flag and qc_note description."
    }
    input {
        String? raw_read_screen
        String? clean_read_screen
        Float? est_coverage_clean
        String? species
        String? gambit_predicted_taxon
        Float? checkm2_completeness
        Float? checkm2_contamination
        Float min_coverage = 40.0
        Float max_contamination = 5.0
        Float min_completeness = 80.0
    }
    call task_qc_flags.qc_flags_row {
        input:
            raw_read_screen = raw_read_screen,
            clean_read_screen = clean_read_screen,
            est_coverage_clean = est_coverage_clean,
            species = species,
            gambit_predicted_taxon = gambit_predicted_taxon,
            checkm2_completeness = checkm2_completeness,
            checkm2_contamination = checkm2_contamination,
            min_coverage = min_coverage,
            max_contamination = max_contamination,
            min_completeness = min_completeness
    }
    output {
        String qc_check = qc_flags_row.qc_check
        String qc_note = qc_flags_row.qc_note
    }
}