version 1.0

import "../../../tasks/gene_typing/task_mob_suite.wdl" as mob_suite_task

workflow test_pilon {
    input {
        String samplename = "test_sample"
        File read1_raw
        File read2_raw
        File assembly
        Int mem = 8
        Int cpu = 2
        Int disk_size = 16
    }
    call mob_suite_task.mob_recon {
        input:
          assembly = assembly,
          samplename = samplename,
          cpu = cpu,
          memory = mem,
          disk_size = disk_size
    }
    output {
        # MOB-recon Results
        File? mob_recon_results = mob_recon.mob_recon_results
        File? mob_typer_results = mob_recon.mob_typer_results
        File? mob_recon_chromosome_fasta = mob_recon.chromosome_fasta
        Array[File]? mob_recon_plasmid_fastas = mob_recon.plasmid_fastas
        String? mob_recon_docker = mob_recon.mob_recon_docker
        String? mob_recon_version = mob_recon.mob_recon_version
    }
}