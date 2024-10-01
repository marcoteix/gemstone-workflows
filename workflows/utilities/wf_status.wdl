version 1.0

import "../../tasks/quality_control/task_status.wdl" as task_status

workflow status {
    meta {
        description: "Generates a status (preferred/duplicate/uci_check/resequence/fail) flag for each isolate, taking a data table as input."
    }
    input {
        String table
        String workspace
    }
    call task_status.status {
        input:
            table = table,
            workspace = workspace
    }
    output {
    }
}