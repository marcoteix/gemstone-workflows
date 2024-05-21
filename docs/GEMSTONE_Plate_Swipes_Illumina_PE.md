# GEMSTONE_Plate_Swipes_Illumina_PE

<span style="font-size:200%; color:grey"> `v1.0.3` </span>

This workflow processes paired-end Illumina reads from plate swipes/plate bacterial metagenomes. It could also be used for broader bacterial metagenomic analysis. It removes human reads, performs QA and QC, and metagenomic assembly. Furthermore, it estimates taxa abundances with Kraken2 and Bracken, does strain-level identification with StrainGE, bins MAGs, performs AMR genotyping, and identifies plasmid contigs.

This workflow was based on the [PHB v1.0.0 TheiaMeta workflow](https://theiagen.notion.site/TheiaMeta-7e73c010559148b380144d7c3b3ec4e6) from Theiagen.

![GEMSTONE_Plate_Swipes_Illumina_PE workflow diagram](https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE_Plate_Swipes_Illumina_PE.svg "GEMSTONE_Plate_Swipes_Illumina_PE workflow diagram")

## Inputs
<details open>
<summary>Click to open or hide</summary>


### Workflow options

<details open>
<summary>Click to open or hide</summary>

#### qc_only
`Boolean` `Optionl` `Default = true` If `true`, performs QC and estimates taxa abundances with Kraken2/Braken.
If `false`, runs all other tasks in the pipeline, depending on the `call_strainge`, `call_kraken`, and `call_metawrap` options. Set to `false` to assemble genomes.

#### output_aditional_files
`Boolean` `Optional` `Default = false` If `true`, aligns reads to the assembly, computes coverage, and retrieves aligned and unaligned reads.

#### reference
`File` `Optional` Reference file for consensus calling, in FASTA format. If provided, performs consensus assembly and calling.

#### call_strainge
`Boolean` `Optional` `Default = false` If `true`, performs strain-level detection with StrainGE.

#### call_kraken
`Boolean` `Optional` `Default = false` If `true`, identifies taxa and estimates their abundance with Kraken2/Bracken.

#### call_metawrap
`Boolean` `Optional` `Default = false` If `true`, performs genome binning and refinement of bins with metaWRAP.

</details>

### Sample data

<details open>
<summary>Click to open or hide</summary>
    
#### read1 
`File` `Required` FASTQ file with forward raw reads. Must be Illumina paired-end.

#### read2
`File` `Required` FASTQ file with reverse raw reads. Must be Illumina paired-end.

#### samplename 
`String` `Required` Name or ID of the sample.

#### lab_determined_genus 
`String` `Required` Genus or species name, as determined in the lab or by GAMBIT. Must be written in full, with whitespaces (e.g., *Escherichia coli* and not *E. coli* nor *Escherichia_coli*). Used to select a StrainGE reference database by matching the genus passed in this parameter with those in each database. Multiple genera are supported through a "/" separator, such as "Escherichia/Klebsiella" - this will call StrainGE twice, once using an Escherichia database, and once with a Klebsiella database.

</details>

### Strain-level detection

<details open>
<summary>Click to open or hide</summary>

#### strainge_db_config
`File` `Required` TSV configuration file for StrainGE databases. It should be a table with two columns: one with the database genus name (e.g., *Escherichia* or *Proteus*), and another with the path to the tar archive with the StrainGE database for that genus. An example of this table is:

<table>
    <tr> 
        <th> Escherichia </th>
        <td> gs://fc-secure-uuid/databases/strainge/escherichia_shigella.tar.gz </td> 
    </tr>
    <tr> 
        <th> Shigella </th>
        <td> gs://fc-secure-uuid/databases/strainge/escherichia_shigella.tar.gz </td> 
    </tr>
    <tr> 
        <th> Pseudomonas </th>
        <td> gs://fc-secure-uuid/databases/strainge/pseudomonas.tar.gz </td> 
    </tr>	
    <tr> 
        <th> Staphylococcus </th>
        <td> gs://fc-secure-uuid/databases/strainge/staphylococcus.tar.gz </td> 
    </tr>
    <tr> 
        <th> Proteus </th>
        <td> gs://fc-secure-uuid/databases/strainge/proteus.tar.gz </td> 
    </tr>	
    <tr> 
        <th> Klebsiella </th>
        <td> gs://fc-secure-uuid/databases/strainge/klebsiella.tar.gz </td> 
    </tr>	
    <tr> 
        <th> Acinetobacter </th>
        <td> gs://fc-secure-uuid/databases/strainge/acinetobacter.tar.gz </td> 
    </tr>	
    <tr> 
        <th> Enterobacter </th>
        <td> gs://fc-secure-uuid/databases/strainge/enterobacter.tar.gz </td> 
    </tr>	
    <tr> 
        <th> Enterococcus </th>
        <td> gs://fc-secure-uuid/databases/strainge/enterococcus.tar.gz </td> 
    </tr>	
</table>

#### strainge_max_strains
`Int` `Optional` `Default = 5` Maximum number of strains searched by StrainGST.

#### strainge_disk_size
`Int` `Optional` `Default = 100` Disk size in Gb for the StrainGE task.

#### strainge_cpus
`Int` `Optional` `Default = 4` Number of CPUs used in the StrainGE task.

#### strainge_memory
`Int` `Optional` `Default = 128` RAM in Gb for the StrainGE task.

#### strainge_db_kmer_size
`Int` `Optional` `Default = 23` K-mer sized used when creating the StrainGST databases.

</details>

### Taxonomy assignment

<details open>
<summary>Click to open or hide</summary>

#### kraken2_db
`File` `Optional` Compressed Kraken2/Bracken database as a tar archive. Required if `call_kraken` is `true`. Please make sure that the archive contains all the files needed to run Bracken (refer to the [Bracken docs](https://github.com/jenniferlu717/Bracken)).

#### bracken_read_len
`Int` `Optional` `Default = 150` Input read length.

#### bracken_classification_level
`String` `Optional` `Default = "G"` Taxonomic level for Bracken abundance estimation. Defaults to genus level (G).  Other possible options are K (kingdom level), P (phylum), C (class), O (order), F (family), and S (species)

#### kraken2_disk_size
`Int` `Optional` `Default = 256` Disk size in Gb for the Kraken2/Bracken task.

#### kraken2_cpu
`Int` `Optional` `Default = 4` Number of CPUs used in the Kraken2/Bracken task task.

#### kraken2_mem
`Int` `Optional` `Default = 32` RAM in Gb for the Kraken2/Bracken task task.

</details>

### MAG binning

<details open>
<summary>Click to open or hide</summary>

#### metawrap_checkm_db
`File` `Optional` Compressed checkM database as a tar archive to be used in metaWRAP. Required if `call_metawrap` is `true`.

#### metawrap_completion
`Int` `Optional` `Default = 80` Minimum completion of MAG bins, as determined by checkM, in percentage (i.e., a value of `80` means only bins with comlpetion greater or equal to 80% will be returned).

#### metawrap_contamination
`Int` `Optional` `Default = 10` Maximum contamination of MAG bins, as determined by checkM, in percentage (i.e., a value of `10` means only bins with contamination less than or equal to 10% will be returned).

#### metawrap_min_contig_length
`Int` `Optional` `Default = 1000` Minimum length in bp of contigs included in the MAG bins returned by metaWRAP.

#### binning_flags
`String` `Optional` `Default = "--metabat2 --maxbin2 --concoct"` Contig binning tools used in metaWRAP, as flags. These flags are used in the metaWRAP *binning* command.

</details>
</details>

## Outputs

<details open>
<summary>Click to open or hide</summary>

#### theiameta_illumina_pe_version
`String` Version of the TheiaMeta workflow used.

#### analysis_date
`String` Analysis date.

### Taxon identification

<details open>
<summary>Click to open or hide</summary>

#### kraken2_version
`String` `Optional`
 
Kraken2 version used.

#### kraken2_docker
`String` `Optional` Name of the Kraken2 and Bracken Docker image used.

#### kraken2_report
`File` `Optional` Kraken2 report with taxa classifications.

#### kraken2_percent_human
`File` `Optional` Percentage of human-classified reads in the sample as determined by Kracken2.

#### bracken_report
`File` `Optional` Bracken report with estimated abundances per taxon.

#### bracken_version
`String` `Optional` Bracken version used.

</details>

### Read QC

<details open>
<summary>Click to open or hide</summary>

#### read1_dehosted
`File` `Optional`
 
FASTQ file of forward reads with human reads removed.

#### read2_dehosted
`File` `Optional`
 
FASTQ file of reverse reads with human reads removed.

#### ncbi_scrub_docker
`String` `Optional` Name of the NCBI Human Scrubber Docker image used.

#### num_reads_raw1
`Int` Number of reads in `read1`.

#### num_reads_raw2
`Int` Number of reads in `read2`.

#### num_reads_raw_pairs
`Int` Number of pairs of reads in `read1` and `read2` (raw reads).

#### fastq_scan_version
`String` Version of fastq_scan used.

#### fastq_scan_docker
`String` `Optional` Name of the fastq_scan Docker image used.

#### num_reads_clean1
`Int` Number of reads after QC in `read1_clean`.

#### num_reads_clean2
`Int` Number of reads after QC in `read2_clean`.

#### num_reads_clean_pairs
`Int` Number of read pairs after QC in `read1_clean` and `read2_clean` (clean reads).

#### trimmomatic_version
`String` `Optional` Version of trimmomatic used.

#### trimmomatic_docker
`String` `Optional` Name of the trimmomatic Docker image used.

#### read1_clean
`File` FASTQ file with forward cleaned reads, after QC and de-hosting.

#### read2_clean
`File` FASTQ file with reverse cleaned reads, after QC and de-hosting.

#### bbduk_docker
`String` Name of the BBDuk Docker image used.

#### average_read_length
`Float` Average length in bp of clean reads (`read1_clean` and `read2_clean`).

</details>

### MAG and assembly QC

<details open>
<summary>Click to open or hide</summary>

#### assembly_fasta
`File` FASTA file of the final assembly (MAG). If no `reference` is used, it is the output from metaSPAdes + Pilon.

#### metaspades_version
`String` Version of metaSPAdes used.

#### metaspades_docker
`String` Name of the metaSPAdes Docker image used.

#### minimap2_version
`String` Version of minimap2 used.

#### minimap2_docker
`String` Name of the minimap2 Docker image used.

#### samtools_version
`String` Version of samtools used.

#### samtools_docker
`String` Name of the samtools Docker image used.

#### pilon_version
`String` Version of Pilon used.

#### pilon_docker
`String` Name of the Pilon Docker image used.

#### assembly_length
`Int` Assembly (MAG) length in bp.

#### contig_number
`Int` Number of MAG contigs.

#### largest_contig
`Int` Length of the largest MAG contig in bp.

#### quast_version
`String` Version of Quast used for assembly QC.

#### quast_docker
`String` Name of the Quast Docker image used for assembly QC.

#### percent_coverage
`Float` `Optional` Percentage coverage of the reference genome (if provided).

#### assembly_mean_coverage
`Float` `Optional` Mean depth of coverage of the final assembly. Returned only if `output_additional_files` is `true`.

#### bedtools_version
`String` `Optional` Version of bedtools used for assembly QC. Returned only if `output_additional_files` is `true`.

#### bedtools_docker
`String` `Optional` Name of the bedtools Docker image used for assembly QC. Returned only if `output_additional_files` is `true`.

#### read1_unmapped
`File` `Optional` Unmapped forwards reads to the assembly. Returned only if `output_additional_files` is `true`.

#### read2_unmapped
`File` `Optional` Unmapped reverse reads to the assembly. Returned only if `output_additional_files` is `true`.

#### read1_unmapped
`File` `Optional` Mapped forwards reads to the assembly. Returned only if `output_additional_files` is `true`.

#### read2_unmapped
`File` `Optional` Mapped reverse reads to the assembly. Returned only if `output_additional_files` is `true`.

#### percentage_mapped_reads
`Float` `Optional` Percentage of mapped reads to the assembly. Returned only if `output_additional_files` is `true`.

</details>

### Strain-level detection

<details open>
<summary>Click to open or hide</summary>

#### straingst_kmerized_reads
`Array[File]` `Optional` Files with k-merized input reads. The size of the array depends on how many genera are assigned in `lab_determined_genus`, but the contents of each file should be the same. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingst_selected_db
`Array[String]` `Optional` StrainGST databases used in each call to StrainGE. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingst_found_db
`Boolean` `Optional` Whether a StrainGST database matching the genera in `lab_determined_genus` was found. Returned only if `call_strainge` is `true`.

#### straingst_strains
`Array[File]` `Optional` Text files of strains found by StrainGST when using each database. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingst_statistics
`Array[File]` `Optional` Reports with StrainGST statistics, including strain relative abundances, with each databased used. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingr_concat_fasta
`Array[File]` `Optional` Concatenated references as FASTA files needed for downstream StrainGR analysis. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingr_read_alignment
`Array[File]` `Optional` Indexed BAM file of clean reads mapped to the concatenated references (in `straingr_concat_fasta`). The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingr_variants
`Array[File]` `Optional` HDF5 files with variant calling results from StrainGR. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### straingr_report
`Array[File]` `Optional` StrainGR reports with variant calling statistics. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` and `straingst_found_db` are `true`.

#### strainge_docker
`Array[String]` `Optional` Name of the StrainGE Docker image used for strain-level detection. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` is `true`.

#### strainge_version
`Array[String]` `Optional` Version of StrainGE used for strain-level detection. The size of the array depends on how many genera are assigned in `lab_determined_genus`. Returned only if `call_strainge` is `true`.

</details>

### Gene annotation

<details open>
<summary>Click to open or hide</summary>

#### bakta_gbff
`File` MAG gene annotations from Bakta in GenBank format.

#### bakta_gff3
`File` MAG gene annotations from Bakta in GFF3 format.

#### bakta_tsv
`File` MAG gene annotations from Bakta in TSV format.

#### bakta_summary
`File` Summary report of MAG gene annotation from Bakta.

#### bakta_version
`String` Version of Bakta used for MAG gene annotation.

</details>

### Plasmid identification

<details open>
<summary>Click to open or hide</summary>

#### mob_recon_results
`File` TSV file with plasmid/chromosome classification of contigs from MOB-recon.

#### mob_typer_results
`File` TSV file with plasmid typing results from MOB-typer.

#### mob_recon_chromosome_fasta
`File` FASTA file of chromosomal contigs in the MAG.

#### mob_recon_plasmid_fastas
`File` FASTA file of plasmid contigs in the MAG.

#### mob_recon_docker
`String` Name of the MOB-recon/MOB-suite Docker image used for plasmid identification.

#### mob_recon_version
`String` Version of MOB-recon/MOB-suite used for plasmid identification.

</details>

### AMR genotyping

<details open>
<summary>Click to open or hide</summary>

#### amrfinderplus_all_report
`File` Report of all genes (virulence, stress, and AMR) found by AMRFinderPlus, as a TSV file.

#### amrfinderplus_amr_report
`File` Report of AMR genes found by AMRFinderPlus, as a TSV file.

#### amrfinderplus_stress_report
`File` Report of stress genes found by AMRFinderPlus, as a TSV file.

#### amrfinderplus_virulence_report
`File` Report of virulence genes found by AMRFinderPlus, as a TSV file.

#### amrfinderplus_amr_core_genes
`String` Comma separated list of core AMR genes found by AMRFinderPlus.

#### amrfinderplus_amr_plus_genes
`String` Comma separated list of plus AMR genes found by AMRFinderPlus.

#### amrfinderplus_stress_genes
`String` Comma separated list of stress genes found by AMRFinderPlus.

#### amrfinderplus_virulence_genes
`String` Comma separated list of virulence genes found by AMRFinderPlus.

#### amrfinderplus_amr_classes
`String` Comma separated list of classes of antimicrobial drugs for which AMR genes were found by AMRFinderPlus.

#### amrfinderplus_amr_subclasses
`String` Comma separated list of subclasses of antimicrobial drugs for which AMR genes were found by AMRFinderPlus.

#### amrfinderplus_version
`String` Version of AMRFinderPlus used for AMR genotyping.

#### amrfinderplus_db_version
`String` Version of AMRFinderPlus database used for AMR genotyping.

</details>

### MAG binning

<details open>
<summary>Click to open or hide</summary>

#### metawrap_docker
`String` `Optional` Name of the metaWRAP Docker image used for MAG binning and bin refinement. Returned only if `call_metawrap` is `true`.

#### metawrap_version
`String` `Optional` Version of metaWRAP used for MAG binning and bin refinement. Returned only if `call_metawrap` is `true`.

#### metawrap_stats
`File` `Optional` TSV with MAG bins statistics (including contamination and completeness). Returned only if `call_metawrap` is `true`.

#### metawrap_n_bins
`Int` `Optional` Number of MAG bins found by metaWRAP. Returned only if `call_metawrap` is `true`.

#### metawrap_binning_flags
`String` `Optional` Contig binning tools used in metaWRAP, as flags. These flags were used in the metaWRAP *binning* command and provided as an input parameter. They are returned as outputs for record keeping. Returned only if `call_metawrap` is `true`.

#### metawrap_fasta
`Array[File]` `Optional` MAG bins found by metaWRAP as FASTA files. Returned only if `call_metawrap` is `true`.

#### metawrap_contigs
`File` `Optional` Bin assignments to each contig in the MAG, as a TSV file. Returned only if `call_metawrap` is `true`.

</details>
</details>

## Authors

__Marco Teixeira__

<span style="font-size:80%; margin-top:0px"> mcarvalh@broadinstitute.org </span>

__Colin Worby__

<span style="font-size:80%; margin-top:0px"> cworby@broadinstitute.org </span>

<img src="https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE%20logo.png" alt="GEMSTONE logo" width=30%>