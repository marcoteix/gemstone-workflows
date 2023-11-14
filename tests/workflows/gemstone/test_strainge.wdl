version 1.0

import "../../../tasks/species_typing/task_strainge.wdl" as strainge_task
import "../../../tasks/species_typing/task_select_strainge_db.wdl" as select_strainge_db_task

workflow test_strainge {
    input {
        String gambit_taxon = 'Escherichia coli'
        File read1
        File read2
        String samplename
        Strainge_db strainge_escherichia_db
        Strainge_db strainge_pseudomonas_db
        Strainge_db strainge_proteus_db
        Strainge_db strainge_klebsiella_db
        Strainge_db strainge_staphylococcus_db
        Strainge_db strainge_acinetobacter_db
        Strainge_db strainge_enterococcus_db
        Strainge_db strainge_enterobacter_db
        Int strainge_kmer_size = 23
    } 
    call select_strainge_db_task.select_reference_db {
        input:
          samplename = samplename,
          gambit_taxonomy = gambit_taxon,
          escherichia_db = strainge_escherichia_db,
          pseudomonas_db = strainge_pseudomonas_db,
          proteus_db = strainge_proteus_db,
          klebsiella_db = strainge_klebsiella_db,
          staphylococcus_db = strainge_staphylococcus_db,
          acinetobacter_db = strainge_acinetobacter_db,
          enterococcus_db = strainge_enterococcus_db,
          enterobacter_db = strainge_enterobacter_db
    }
    call strainge_task.StrainGE_PE as strainge {
    input:
        samplename = samplename,
        reads_1 = read1,
        reads_2 = read2,
        kmer_size = strainge_kmer_size,
        straingst_reference_db = select_reference_db.selected_strainge_db
    } 
    output {
        # StrainGE Results
        File? straingst_kmerized_reads = strainge.straingst_kmerized_reads
        File? straingst_reference_db = strainge.straingst_reference_db_used
        File? straingst_strains = strainge.straingst_strains
        String? straingst_reference_genus = strainge.straingst_reference_genus_used
        File? straingst_statistics = strainge.straingst_statistics
        File? straingr_concat_fasta = strainge.straingr_concat_fasta
        File? straingr_read_alignment = strainge.straingr_read_alignment
        File? straingr_variants = strainge.straingr_variants
        File? straingr_report = strainge.straingr_report
        String? strainge_docker = strainge.strainge_docker
        String? strainge_version = strainge.strainge_version
    }
}