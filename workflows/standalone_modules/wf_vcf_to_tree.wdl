version 1.0

import "../../tasks/utilities/task_vcf_to_msa.wdl" as vcf_to_msa_task
import "../../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins_task
import "../../tasks/phylogenetic_inference/task_iqtree2.wdl" as iqtree2_task
import "../../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists_task

workflow vcf_to_tree {
    input {
        String collection_name
        Array[String] samplenames
        Array[File] variants_vcfs
        String vcf_filters = "PASS,."
        Int min_samples = 1
        Boolean use_gubbins = true
    }
    # Convert a set of VCFs to a multisequence alignment FASTA
    call vcf_to_msa_task.vcf_to_msa {
        input:
            samplenames = samplenames,
            vcfs = variants_vcfs,
            collection_name = collection_name,
            filters = vcf_filters,
            min_samples = min_samples
    }
    # Mask recombinant regions with Gubbins
    if (use_gubbins) {
        call gubbins_task.gubbins {
            input:
                alignment = vcf_to_msa.msa,
                cluster_name = collection_name,
                tree_builder = "iqtree"
        }
    }
    # Create a ML tree with IQTree2
    call iqtree2_task.iqtree2 {
        input:
            alignment = select_first([gubbins.gubbins_polymorphic_fasta, vcf_to_msa.msa]),
            cluster_name = collection_name
    }
    # Create a pairwise SNP matrix
    call snp_dists_task.snp_dists {
        input:
            alignment = select_first([gubbins.gubbins_polymorphic_fasta, vcf_to_msa.msa]),
            cluster_name = collection_name
    }
    output {

        File msa_fasta = select_first([gubbins.gubbins_polymorphic_fasta, vcf_to_msa.msa])
        File vcf_to_tree_final_tree = iqtree2.ml_tree

        # Gubbins outputs
        String? gubbins_docker = gubbins.gubbins_docker
        String? gubbins_version = gubbins.gubbins_version
        File? gubbins_recombination_gff = gubbins.gubbins_recombination_gff
        File? gubbins_branch_stats = gubbins.gubbins_branch_stats

        # IQTree2 outputs
        String iqtree2_version = iqtree2.iqtree2_version
        String iqtree2_docker = iqtree2.iqtree2_docker
        String iqtree2_model_used = iqtree2.iqtree2_model_used
        String vcf_to_tree_date = iqtree2.date

        # SNP dists outputs
        String snp_dists_version = snp_dists.snp_dists_version
        String snp_dists_docker = snp_dists.snp_dists_docker
        File vcf_to_tree_snp_matrix = snp_dists.snp_matrix
    }
}