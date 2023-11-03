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
        checkm2 predict --threads 30 --input ~{assembly} --output-directory ~{samplename}

        # Parse report
        while read -r name completedness contamination model translation_table coding_density \
        n50 avg_gene_len genome_size gc total_cds notes
        do
            echo $completedness > ~{samplename}/completedness.txt
            echo $contamination > ~{samplename}/contamination.txt
            echo $coding_density > ~{samplename}/coding_density.txt
            echo $n50 > ~{samplename}/contig_n50.txt
            echo $avg_gene_len > ~{samplename}/avg_gene_len.txt
            echo $genome_size > ~{samplename}/genome_size.txt
            echo $gc > ~{samplename}/gc_content.txt
            echo $total_cds > ~{samplename}/total_cds.txt
        done

        # versioning
        checkm2 --version > ~{samplename}/VERSION.txt
    >>>
    output {
        File report = "~{samplename}/quality_report.tsv"
        Float completedness = read_float("~{samplename}/completedness.txt")
        Float contamination = read_float("~{samplename}/contamination.txt")
        Float coding_density = read_float("~{samplename}/coding_density.txt")
        Int contig_n50 = read_int("~{samplename}/contig_n50.txt")
        Float avg_gene_len = read_float("~{samplename}/avg_gene_len.txt")
        Int genome_size = read_int("~{samplename}/genome_size.txt")
        Float gc_content = read_float("~{samplename}/gc_content.txt")
        Int total_cds = read_int("~{samplename}/total_cds.txt")
        String checkm2_docker = docker
        String version = read_string("~{samplename}/VERSION.txt")
    }
    runtime {
        docker: "~{docker}"
        memory: "~{memory} GB"
        cpu: cpu
        disks:  "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB" # TES
        preemptible: 0
        maxRetries: 3
    }
}