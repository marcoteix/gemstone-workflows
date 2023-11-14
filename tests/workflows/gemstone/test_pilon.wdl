version 1.0

import "../../../tasks/quality_control/task_pilon.wdl" as pilon
import "../../../tasks/alignment/task_bwa.wdl" as bwa

workflow test_pilon {
    input {
        String samplename = "test_sample"
        File read1_raw
        File read2_raw
        File reference_genome
        Int mem = 8
        Int cpu = 2
        Int disk_size = 16
    }
    call bwa.bwa as bwa_pre_pilon{
        input:
            read1 = read1_raw,
            read2 = read2_raw,
            reference_genome = reference_genome,
            samplename = samplename,
            cpu = cpu,
            disk_size = disk_size
    }
    call pilon.pilon {
        input:
          assembly = reference_genome,
          bam = bwa_pre_pilon.sorted_bam,
          bai = bwa_pre_pilon.sorted_bai,
          samplename = samplename,
          cpu = cpu,
          memory = mem,
          disk_size = disk_size
    }
    call bwa.bwa as bwa_post_pilon{
        input:
            read1 = read1_raw,
            read2 = read2_raw,
            reference_genome = pilon.assembly_fasta,
            samplename = samplename,
            disk_size = disk_size,
            cpu = cpu
    }
    output {
        # Assembly - pilon outputs
        File? assembly_fasta = pilon.assembly_fasta
        File? pilon_changes = pilon.changes
        File? pilon_vcf = pilon.vcf
        String? pilon_version = pilon.pilon_version
        # Alignment - bwa pre-pilon outputs
        String? bwa_version = bwa_pre_pilon.bwa_version
        String? sam_version = bwa_pre_pilon.sam_version
        File? raw_assembly_alignment_bam = bwa_pre_pilon.sorted_bam
        File? raw_assembly_alignment_index = bwa_pre_pilon.sorted_bai
        File? raw_assembly_alignment_read1 = bwa_pre_pilon.read1_aligned
        File? raw_assembly_alignment_read2 = bwa_pre_pilon.read2_aligned
        # Alignment - bwa post-pilon outputs
        File? assembly_alignment_bam = bwa_post_pilon.sorted_bam
        File? assembly_alignment_index = bwa_post_pilon.sorted_bai
        File? assembly_alignment_read1 = bwa_post_pilon.read1_aligned
        File? assembly_alignment_read2 = bwa_post_pilon.read2_aligned
    }
}