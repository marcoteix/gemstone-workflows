version 1.0

task checkm2 {
    input {
        File assembly
        String samplename
        String docker = "marcoteix/checkm2:0.0.1"
        Int cpu = 8
        Int memory = 32
        Int disk_size = 100
    }
    command <<<
        checkm2 predict --threads 4 --tmpdir /tmp/ --input ~{assembly} --output-directory ~{samplename} --lowmem

        # Parse report
        while read -r name completeness contamination others
        do
            echo $completeness > ~{samplename}/completeness.txt
            echo $contamination > ~{samplename}/contamination.txt
        done < ~{samplename}/quality_report.tsv

        # versioning
        checkm2 --version > ~{samplename}/VERSION.txt
    >>>
    output {
        File report = "~{samplename}/quality_report.tsv"
        Float completeness = read_float("~{samplename}/completeness.txt")
        Float contamination = read_float("~{samplename}/contamination.txt")
        String checkm2_docker = docker
        String checkm2_version = read_string("~{samplename}/VERSION.txt")
    }
    runtime {
        docker: "~{docker}"
        memory: "~{memory} GB"
        cpu: cpu
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB" # TES
        preemptible: 0
        maxRetries: 0
    }
}