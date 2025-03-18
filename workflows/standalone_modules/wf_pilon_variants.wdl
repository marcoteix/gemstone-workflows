version 1.0

import "../../tasks/alignment/task_bwa.wdl" as bwa
import "../../tasks/quality_control/task_pilon.wdl" as pilon
import "../../tasks/utilities/task_filter_variants.wdl" as filter_variants

workflow pilon_variants {
    meta {
        author: "Marco Teixeira"
        email: "mcarvalh@broadinstitute.org"
        description: "Pipeline for variant calling with Pilon."
    }
    input {
        String samplename
        File reads_1
        File reads_2
        File reference_fasta
        Int alignment_cpu = 6
        Int alignment_disk_size = 100
        Int pilon_memory = 32
    }
    call bwa.bwa {
        input:
            read1 = reads_1,
            read2 = reads_2,
            samplename = samplename,
            reference_genome = reference_fasta,
            cpu = alignment_cpu,
            disk_size = alignment_disk_size
    }
    call pilon.pilon {
        input:
            assembly = reference_fasta,
            bam = bwa.sorted_bam,
            bai = bwa.sorted_bai,
            samplename = samplename,
            fix = "bases",
            memory = pilon_memory
    }
    call filter_variants.filter_variants {
        input:
            vcf = pilon.vcf,
            samplename = samplename
    }
    output {
        File variants_vcf = pilon.vcf
        String pilon_version = pilon.pilon_version
        String pilon_docker = pilon.pilon_docker
        String bwa_version = bwa.bwa_version
    }
}