# GEMSTONE_Isolates_Illumina_PE

<span style="font-size:200%; color:grey"> `v1.0.3` </span>

This workflow processes paired-end Illumina reads from bacterial isolates. It performs QA and QC, flagging low quality and contaminated samples. Samples for which reads fail QC are not subject to further analysis. Performs *de novo* assembly with SPAdes and refines it with Pilon. Finds AMR, virulence, and stress genes with AMRFinderPlus, identifies and types plasmid contigs with MOB-recon, types sequences with TS-MLST, annotates them with Bakta, and infers taxonomy with GAMBIT. Also performs taxa-specific analysis - please refer to [TheiaProk documentation](https://theiagen.notion.site/TheiaProk-Workflow-Series-1296f3df6a7c41c7937053a03f664a5a) for those tasks. Optionally, it estimates taxa abundances with Kraken2 and Bracken and does strain-level identification with StrainGE (based on the GAMBIT predicted taxon).

This workflow was based on the [PHB v1.0.0 TheiaProk workflow](https://theiagen.notion.site/TheiaProk-Workflow-Series-1296f3df6a7c41c7937053a03f664a5a) from Theiagen.

![GEMSTONE_Isolates_Illumina_PE workflow diagram](https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE_Isolates_Illumina_PE.svg "GEMSTONE_Isolates_Illumina_PE workflow diagram")

## Inputs
<details open>
<summary>Click to open or hide</summary>

### Workflow options

<details open>
<summary> Click to open or hide </summary>

#### call_kmerfinder
`Boolean` `Optional` `Default = false` If `true`, predicts species based on k-mer similarity (with k = 16) between the assebly and genomes in a [Theiagen bacterial database](https://theiagen.notion.site/TheiaProk-Workflow-Series-1296f3df6a7c41c7937053a03f664a5a?p=27c5fc0f5dc94df9a01523d5da1dc22b&pm=s).

#### call_kraken
`Boolean` `Optional` `Default = false` If `true`, identifies taxa and estimates their abundance with Kraken2/Bracken.

#### call_strainge
`Boolean` `Optional` `Default = false` If `true`, performs strain-level detection with StrainGE.

#### expected_taxon
`String` `Optional` If set, overrides the GAMBIT predicted species when setting a species in AMRFinderPlus and when comparing QC metrics agains the thresholds defined in `qc_check_table`. Useful for when GAMBIT predictions are incorrect. If not set, uses GAMBIT the prediction.

</details>

### Sample data

<details open>
<summary> Click to open or hide </summary>

#### read1_raw
`File` `Required` FASTQ file with forward raw reads. Must be Illumina paired-end.

#### read2_raw
`File` `Required` FASTQ file with reverse raw reads. Must be Illumina paired-end.

#### samplename 
`String` `Required` Name or ID of the sample.

#### lab_determined_genus 
`String` `Optional` Genus or species name, as determined in the lab. Must be written in full, with whitespaces (e.g., *Escherichia coli* and not *E. coli* nor *Escherichia_coli*). It is compared against the GAMBIT predicted taxonomy to derive a *taxonomy QC flag*, corresponding to the `qc_taxonomy_check` output. If the lab predicted genus mathces the GAMBIT prediction, then `qc_taxonomy_check` is set to `PASS`; otherwise, it is set to `ALERT`.

</details>

### Global QA/QC

<details open>
<summary> Click to open or hide </summary>

#### contamination_threshold 
`Float` `Optional` `Default = 2` Maximum contamination as a percentage (as determined by checkM2) allowed for a sample to pass taxonomy QC. The default of 2 means that the contamination threshold is 2%.

#### qc_check_table
`File` `Optional` User-defined, taxa-specific, thresholds for QC metrics as a TSV file. If all QC metrics meet the threshold, the `qc_check` output variable will read `QC_PASS`. Otherwise, the output will read `QC_NA` if the task could not proceed or `QC_ALERT` followed by a string indicating what metric failed. Each row in the table should be a species or genus, written in full, with underscores instead of whitespaces (matching the format from GAMBIT). Column names should be *taxon*, followed by the QC metric name. The sample taxa is taken from the `gambit_predicted_taxon` value inferred by GAMBIT or from a user-defined `expected_taxon`. Example of a `qc_check_table`:

<table><tbody><tr><th scope="col" style="font-weight: 500; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">taxon</div></div></th><th scope="col" style="font-weight: 500; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">est_coverage_raw</div></div></th><th scope="col" style="font-weight: 500; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">est_coverage_clean</div></div></th><th scope="col" style="font-weight: 500; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">assembly_length_min</div></div></th><th scope="col" style="font-weight: 500; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">assembly_length_max</div></div></th></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Listeria_monocytogenes</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">20</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">2800000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">3200000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Escherichia_coli</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">40</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4900000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">6000000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Shigella</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">40</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4200000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4900000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Salmonella</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">30</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4400000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">5700000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Campylobacter</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">20</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">1400000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">2200000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Vibrio_cholerae</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">40</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">3800000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4300000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Vibrio_parahaemolyticus</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">40</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4900000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">5500000</div></div></td></tr><tr><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">Vibrio_vulnificus</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">40</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 120px; max-width: 240px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px; min-height: 1em; color: rgb(55, 53, 47); -webkit-text-fill-color: rgba(55, 53, 47, 0.5);"></div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 170px; max-width: 170px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">4700000</div></div></td><td style="color: inherit; fill: inherit; border: 1px solid rgb(233, 233, 231); position: relative; vertical-align: top; text-align: start; min-width: 166px; max-width: 166px; min-height: 32px;"><div><div style="max-width: 100%; width: 100%; ; padding: 7px 9px; font-size: 14px; line-height: 20px;">5300000</div></div></td></tr></tbody></table>

</details>

### Read QA/QC

<details open>
<summary> Click to open or hide </summary>

#### genome_size
`Int` `Optional` Expected genome size in bp. Used during read and assembly QC. If not provided, the workflow uses the genome length estimated by QUAST for assembly QC.

#### skip_screen
`Int` `Optional` `Default = false` If `true`, skips read QA/QC.

#### min_reads
`Int` `Optional` `Default = 7472` Minimum number of reads (raw and clean) needed to pass QC.

#### min_basepairs
`Int` `Optional` `Default = 2241820` Minimum number of total bp in reads (raw and clean) needed to pass QC.

#### min_genome_size
`Int` `Optional` `Default = 7472` Minimum estimated genome size in bp needed to pass QC (both for raw and clean reads).

#### max_genome_size
`Int` `Optional` `Default = 18040666` Maximum estimated genome size in bp needed to pass QC (both for raw and clean reads).

#### min_coverage
`Int` `Optional` `Default = 10` Minimum estimated genome coverage needed to pass QC (both for raw and clean reads).

#### min_proportion
`Int` `Optional` `Default = 40` The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than these proportion of basepairs are in either the forward or reverse reads files.

</details>

### Read trimming

<details open>
<summary> Click to open or hide </summary>

#### trim_minlen
`Int` `Optional` `Default = 75` Minimum read length in bp required after trimming for it to be included in downstream analyses.

#### trim_quality_trim_score
`Int` `Optional` `Default = 20` Average quality of bases in a sliding window needed for those bases to be kept.

#### trim_window_size
`Int` `Optional` `Default = 10` Length in bp of the window used for trimming.

</details>

### Taxon identification with Kracken2/Bracken

<details open>
<summary> Click to open or hide </summary>

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

### Strain-level detection with StrainGE

<details open>
<summary> Click to open or hide </summary>

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

#### strainge_db_kmer_size
`Int` `Optional` `Default = 23` K-mer sized used when creating the StrainGST databases.

</details>
</details>

## Outputs

<details open>
<summary>Click to open or hide</summary>

#### theiameta_illumina_pe_version
`String` Version of the TheiaProk workflow used.

#### analysis_date
`String` Analysis date.

### Global QA/QC

<details open>
<summary>Click to open or hide</summary>

#### raw_read_screen
`String` "PASS" or "FAIL" result from raw read screening. If the result is "FAIL", the flag is accompanied by the reason for failure.

#### clean_read_screen
`String` `Optional` "PASS" or "FAIL" result from clean read screening. If the result is "FAIL", the flag is accompanied by the reason for failure. If the raw reads did not pass QC, `clean_read_screen` will not be returned.

#### qc_taxonomy_check
`String` `Optional` "QC_PASS" or "QC_ALERT" flag for taxon identification. If the `lab_predicted_genus` and the taxon prediction from GAMBIT match **and** the estimated genome contamination is less than `contamination_threshold`, this flag is set to "QC_PASS"; otherwise, it is set to "QC_ALERT". If the `lab_predicted_genus` is unknown, this flag only depends on the estimated contamination.

#### qc_taxonomy_report
`File` `Optional` Text file containing detailed flags for taxon identification.
- "Lab genus QC" holds the results from the comparison between the `lab_predicted_genus` and the taxon prediction from GAMBIT match: if these match or `lab_predicted_genus` is unknown, this flag is set to "QC_PASS"; otherwise, it is set to "QC_ALERT". 
- "Contamination QC" is set to "QC_PASS" if the estimated contamination by checkM is less than or equal to `contamination_threshold`, and to "QC_ALERT" otherwise.
- "Global QC" holds the global QC flag  for taxon identification: "QC_PASS" if both "Lab genus QC" and "Contamination QC" are "QC_PASS", or "QC_ALERT" otherwise. It has the same value as `qc_taxonomy_check`. 

#### read_and_table_qc_check
`String` `Optional` "QC_PASS"/"QC_ALERT" flag resulting from the comparison between QC metrics and the user-defined thresholds in the `qc_check_table`. Returned only if `qc_check_table` is provided.

#### read_and_table_qc_standard
`File` `Optional` The `qc_check_table` with the user-defined thresholds. If no `qc_check_table` is provided, `read_and_table_qc_check` is not returned.

#### qc_check
`String` `Optional` Global "QC_PASS"/"QC_ALERT" flag. It is set to "QC_ALERT" if and only if at least one of `read_and_table_qc_check`, or `qc_taxonomy_check` are "QC_ALERT". Otherwise, it is set to "QC_PASS".

#### qc_note
`String` `Optional` Reason for a `qc_check` "QC_ALERT" flag.

</details>

### Read QA/QC

<details open>
<summary>Click to open or hide</summary>

#### num_reads_raw1
`Int` `Optional` Number of reads in `read1`.

#### num_reads_raw2
`Int` `Optional` Number of reads in `read2`.

#### num_reads_raw_pairs
`Int` `Optional` Number of pairs of reads in `read1` and `read2` (raw reads).

#### fastq_scan_version
`String` `Optional` Version of fastq_scan used.

#### num_reads_clean1
`Int` `Optional` Number of reads after QC in `read1_clean`.

#### num_reads_clean2
`Int` `Optional` Number of reads after QC in `read2_clean`.

#### num_reads_clean_pairs
`Int` `Optional` Number of read pairs after QC in `read1_clean` and `read2_clean` (clean reads).

#### trimmomatic_version
`String` `Optional` Version of trimmomatic used.

#### read1_clean
`File` `Optional` FASTQ file with forward cleaned reads, after QC and de-hosting.

#### read2_clean
`File` `Optional` FASTQ file with reverse cleaned reads, after QC and de-hosting.

#### bbduk_docker
`String` `Optional` Name of the BBDuk Docker image used.

#### r1_mean_q_raw
`Float` `Optional` Mean quality score of forward raw reads.

#### r2_mean_q_raw
`Float` `Optional` Mean quality score of reverse raw reads.

#### combined_mean_q_raw
`Float` `Optional` Mean quality score of forward and reverse raw reads.

#### combined_mean_q_clean
`Float` `Optional` Mean quality score of forward and reverse clean reads.

#### r1_mean_readlength_raw
`Float` `Optional` Mean read length in bp of forward raw reads.

#### r2_mean_readlength_raw
`Float` `Optional` Mean read length in bp of reverse raw reads.

#### combined_mean_readlength_raw
`Float` `Optional` Mean read length in bp of forward and reverse raw reads.

#### combined_mean_readlength_clean
`Float` `Optional` Mean read length in bp of forward and reverse clean reads.

#### est_coverage_raw
`Float` `Optional` Estimated coverage of raw reads, given the estimated genome length.

#### cg_pipeline_report_clean
`File` `Optional` TSV file with read metrics from clean reads, including average read length, number of reads, and estimated genome coverage.

#### est_coverage_clean
`Float` `Optional` Estimated coverage of clean reads, given the estimated genome length.

#### cg_pipeline_report_raw
`File` `Optional` TSV file with read metrics from raw reads, including average read length, number of reads, and estimated genome coverage.

#### cg_pipeline_docker
`String` `Optional` Name of the Docker image used in the CG pipeline (used to get raw and clean reads QC metrics).

</details>

### Assembly

<details open>
<summary>Click to open or hide</summary>

#### assembly_fasta
`File` `Optional` Assembly FASTA file. This file is the refined assembly output from Pilon. See the [shovill pipeline documentation](https://github.com/tseemann/shovill) for more information on how the assembly is generated.

#### raw_assembly_gfa
`File` `Optional` Assembly graph as a GFA file. This is the assembly graph from SPAdes, **prior to the Pilon refinement**. See the [shovill pipeline documentation](https://github.com/tseemann/shovill) for more information on how the assembly is generated.

#### shovill_pe_version
`String` `Optional` Version of the shovill pipeline used for assembly.

#### quast_report
`File` `Optional` Assembly QC report from QUAST as a text file.

#### quast_version
`String` `Optional` QUAST version used for assembly QC.

#### assembly_length
`Int` `Optional` Total contig length in bp.

#### number_contigs
`Int` `Optional` Total number of contigs in the assembly.

#### n50_value
`Int` `Optional` Assembly N50 value (minimum contig length in bp of the largest contigs containing 50% of the total assembly length) as computed by QUAST.

#### quast_gc_percent
`Float` `Optional` Assembly GC percentage as computed by QUAST.

#### checkm2_report
`File` `Optional` TSV report from checkM2, including the estimated assembly completeness and contamination.

#### checkm2_completeness
`Float` `Optional` Estimated genome completeness, as a percentage (i.e., a value of 100 means 100% completeness), estimated by checkM2.

#### checkm2_contamination
`Float` `Optional` Estimated genome contamination, as a percentage (i.e., a value of 1 means 1% contamination), estimated by checkM2.

#### checkm2_version
`String` `Optional` Version of checkM2 used to estimate genome contamination and completeness.

</details>

### Taxon identification

<details open>
<summary>Click to open or hide</summary>

**GAMBIT**

<details open>
<summary> Click to open or hide </summary>

#### gambit_report
`File` `Optional` Report from GAMBIT as a text file, including the predicted species.

#### gambit_closest_genomes
`File` `Optional` CSV file listing genomes in the GAMBIT database that are most similar to the assembly.

#### gambit_predicted_taxon
`String` `Optional` Predicted taxon of the assembly as estimated by GAMBIT. 

#### gambit_predicted_taxon_rank
`String` `Optional` Rank of the predicted taxon of the assembly as estimated by GAMBIT (e.g. species, genus...).

#### gambit_docker
`String` `Optional` Name of the Docker image used in the GAMBIT task (for taxon identification).

#### gambit_version
`String` `Optional` Version of GAMBIT used for taxon identification.

#### gambit_db_version
`String` `Optional` Version of the GAMBIT database used for taxon identification.

</details>

**k-mer similarity-based taxonomy identification**

<details open>
<summary> Click to open or hide </summary>

#### kmerfinder_results_tsv
`File` `Optional` Results of the k-mer similarity-based taxonomy identification, as a TSV file. Returned only if `call_kmerfinder` is `true`.

#### kmerfinder_top_hit
`String` `Optional` Top hit species of the k-mer similarity-based taxonomy identification. Returned only if `call_kmerfinder` is `true`.

#### kmerfinder_query_coverage
`String` `Optional` Query coverage of the top hit result of the k-mer similarity-based taxonomy identification. Returned only if `call_kmerfinder` is `true`.

#### kmerfinder_template_coverage
`String` `Optional` Template coverage of the top hit result of the k-mer similarity-based taxonomy identification. Returned only if `call_kmerfinder` is `true`.

#### kmerfinder_docker
`String` `Optional` Name of the Docker image used for k-mer similarity-based taxonomy identification. Returned only if `call_kmerfinder` is `true`.

#### kmerfinder_database
`String` `Optional` Reference database used for k-mer similarity-based taxonomy identification. Returned only if `call_kmerfinder` is `true`.

</details>

**MLST**
<details open>
<summary> Click to open or hide </summary>

#### ts_mlst_results
`File` `Optional` TSV report with detailed MLST profile, including [missing data symbols](https://github.com/tseemann/mlst#missing-data).

#### ts_mlst_predicted_st
`String` `Optional` Predicted sequence type.

#### ts_mlst_pubmlst_scheme
`String` `Optional` PubMLST scheme used to infer sequence type.

#### ts_mlst_allelic_profile
`File` `Optional` Allelic profile detected when infering sequence type.

#### ts_mlst_novel_alleles
`String` `Optional` FASTA file containing nucleotide sequences of any alleles that are not in the MLST database.

#### ts_mlst_version
`String` `Optional` Version of [MLST](https://github.com/tseemann/mlst) used to infer sequence type.

#### ts_mlst_docker
`String` `Optional` Name of the Docker image used in the MLST task (for sequence typing).

</details>
</details>

### AMR genotyping

<details open>
<summary> Click to open or hide </summary>

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

### Gene annotation

<details open>
<summary> Click to open or hide </summary>

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
<summary> Click to open or hide </summary>

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

### Strain-level detection

<details open>
<summary> Click to open or hide </summary>

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
</details>

## Authors

__Marco Teixeira__

<span style="font-size:80%; margin-top:0px"> mcarvalh@broadinstitute.org </span>

__Colin Worby__

<span style="font-size:80%; margin-top:0px"> cworby@broadinstitute.org </span>

<img src="https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE%20logo.png" alt="GEMSTONE logo" width=30%>