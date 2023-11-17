version 1.0

import "../../../tasks/taxon_id/task_kraken2.wdl" as kraken2

workflow test_kraken {
    input {
        String kraken2_db
        File read1_raw
        File read2_raw
        String samplename = "test_sample"
    }
    call kraken2.kraken2_standalone as kraken2 {
        input:
            samplename = samplename,
            read1 = read1_raw,
            read2 = read2_raw,
            kraken2_db = kraken2_db,
            disk_size = 10,
            cpu = 2
    }
    output {
        String? kraken2_version = kraken2.kraken2_version
        String? kraken2_docker = kraken2.kraken2_docker
        String? kraken2_analysis_date = kraken2.analysis_date
        File? kraken2_report = kraken2.kraken2_report
        File? kraken2_classified_report = kraken2.kraken2_classified_report
        File? kraken2_unclassified_read1 = kraken2.kraken2_unclassified_read1
        File? kraken2_unclassified_read2 = kraken2.kraken2_unclassified_read2
        File? kraken2_classified_read1 = kraken2.kraken2_classified_read1
        File? kraken2_classified_read2 = kraken2.kraken2_classified_read2
        Float? kraken2_percent_human = kraken2.kraken2_percent_human
    }
}