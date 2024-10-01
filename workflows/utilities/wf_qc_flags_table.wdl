version 1.0

import "../../tasks/quality_control/task_qc_flags.wdl" as task_qc_flags

workflow qc_flags_table {
    meta {
        description: "Re-generates QC flags for all rows in an entity data table."
    }
    input {
        String table
        String workspace
        Float min_coverage = 40.0
        Float max_contamination = 5.0
        Float min_completeness = 80.0
    }
    call task_qc_flags.qc_flags_table {
        input:
            table = table,
            workspace = workspace,
            min_coverage = min_coverage,
            max_contamination = max_contamination,
            min_completeness = min_completeness
    }
    output {
    }
}