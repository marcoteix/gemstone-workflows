version 1.0

import "../../tasks/cleansweep/task_cleansweep_filter.wdl" as cleansweep_filter
import "../../tasks/cleansweep/task_cleansweep_prepare.wdl" as cleansweep_prepare
import "../../tasks/utilities/task_cleansweep_find_straingst_references.wdl" as find_references
import "../../tasks/alignment/task_bwa.wdl" as bwa
import "../../tasks/quality_control/task_pilon.wdl" as pilon

workflow cleansweep {
    meta {
        author: "Marco Teixeira"
        email: "mcarvalh@broadinstitute.org"
        description: "Pipeline for strain-specific variant calling from plate swipe data with CleanSweep."
    }
    input {
        String samplename
        File reads_1
        File reads_2
        # Find input FASTAs based on StrainGST outputs
        String query_name
        Array[File] straingst_strains
        String fasta_location 
        String? fasta_extension = ".fa"
        # CleanSweep prepare options
        Float cleansweep_max_identity = 0.95
        Int cleansweep_min_length = 150 
        String cleansweep_docker = "marcoteix/cleansweep:main"
        # Alignment options
        Int alignment_cpu = 6
        Int alignment_disk_size = 32
        # Cleansweep filter options
        Int cleansweep_min_depth = 10
        Int cleansweep_min_alt_bc = 10
        Int cleansweep_min_ref_bc = 0
        Int cleansweep_num_variants_estimator = 200
        Int cleansweep_num_variants_coverage = 100000
        Float cleansweep_max_overdispersion = 0.55
        Int cleansweep_overdispersion_bias = 500
        Int cleansweep_random_state = 23
        Int cleansweep_num_chains = 5
        Int cleansweep_num_draws = 100000
        Int cleansweep_num_burnin = 1000
        Int cleansweep_cpu = 5
    }

    call find_references.find_straingst_references {
        input:
            query_strain = query_name,
            straingst_strains = straingst_strains,
            fasta_location = fasta_location,
            fasta_extension = fasta_extension
    }
    call cleansweep_prepare.cleansweep_prepare as prepare_straingst {
        input:
            samplename = samplename,
            query_reference = find_straingst_references.query_fasta,
            background_references = find_straingst_references.background_fasta,
            max_identity = cleansweep_max_identity,
            min_length = cleansweep_min_length,
            docker = cleansweep_docker
    }
    call bwa.bwa {
        input:
            read1 = reads_1,
            read2 = reads_2,
            samplename = samplename,
            reference_genome = prepare_straingst.cleansweep_reference_fasta,
            cpu = alignment_cpu,
            disk_size = alignment_disk_size,
            open = 36,
            extend = 36,
            clip = 40,
            unpaired = 54,
            mismatch = 24
    }
    call pilon.pilon {
        input:
            assembly = prepare_straingst.cleansweep_reference_fasta,
            bam = bwa.sorted_bam,
            bai = bwa.sorted_bai,
            samplename = samplename,
            fix = "bases",
            extra_options = "--nostrays --duplicates"
    }
    call cleansweep_filter.cleansweep_filter {
        input:
            samplename = samplename,
            variants_vcf = pilon.vcf,
            cleansweep_prepare_swp = prepare_straingst.cleansweep_prepare_swp,
            min_depth = cleansweep_min_depth,
            min_alt_bc = cleansweep_min_alt_bc,
            min_ref_bc = cleansweep_min_ref_bc,
            num_variants_estimator = cleansweep_num_variants_estimator,
            num_variants_coverage = cleansweep_num_variants_coverage,
            max_overdispersion = cleansweep_max_overdispersion,
            overdispersion_bias = overdispersion_bias,
            random_state = cleansweep_random_state,
            num_chains = cleansweep_num_chains,
            num_draws = cleansweep_num_draws,
            num_burnin = cleansweep_num_burnin,
            docker = cleansweep_docker,
            cpu = cleansweep_cpu
    }
    output {
        File cleansweep_reference_fasta = prepare_straingst.cleansweep_reference_fasta
        File cleansweep_prepare_swp = prepare_straingst.cleansweep_prepare_swp
        String cleansweep_version = prepare_straingst.cleansweep_version
        File cleansweep_alignment = bwa.sorted_bam
        File cleansweep_pilon_vcf = pilon.vcf
        File cleansweep_variants = cleansweep_filter.cleansweep_variants
        File cleansweep_filter_swp = cleansweep_filter.cleansweep_filter
        File cleansweep_report = cleansweep_filter.cleansweep_report
        File cleansweep_allele_depths_plot = cleansweep_filter.cleansweep_allele_depths_plot
        File cleansweep_query_depths_plot = cleansweep_filter.cleansweep_query_depths_plot
    }
}