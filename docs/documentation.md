# GEMSTONE Workflows

## Isolate workflow

`GEMSTONE_Isolates_Illumina_PE`

- [Source code](https://github.com/marcoteix/gemstone-workflows/blob/main/workflows/gemstone/wf_gemstone_isolates.wdl)
- [Detailed documentation](GEMSTONE_Isolates_Illumina_PE.md)
- [Dockstore](https://dockstore.org/workflows/github.com/marcoteix/gemstone-workflows/GEMSTONE_Isolates_Illumina_PE:main?tab=info)
- [Launch in Terra](https://app.terra.bio/#import-tool/dockstore/github.com/marcoteix/gemstone-workflows/GEMSTONE_Isolates_Illumina_PE:main)

This workflow processes paired-end Illumina reads from bacterial isolates. It performs QA and QC, flagging low quality and contaminated samples. Performs *de novo* assembly with the *shovill* pipeline. Finds AMR, virulence, and stress genes with AMRFinderPlus, identifies and types plasmid contigs with MOB-recon, types sequences with TS-MLST, annotates them with Bakta, and infers taxonomy with GAMBIT. Optionally, it estimates taxa abundances with Kraken2 and Bracken and does strain-level identification with StrainGE.

![GEMSTONE_Isolates_Illumina_PE workflow diagram](https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE_Isolates_Illumina_PE.svg "GEMSTONE_Isolates_Illumina_PE workflow diagram")


## Plate swipe workflow

`GEMSTONE_Plate_Swipes_Illumina_PE`

- [Source code](https://github.com/marcoteix/gemstone-workflows/blob/main/workflows/gemstone/wf_gemstone_plate_swipes.wdl)
- [Detailed documentation](GEMSTONE_Plate_Swipes_Illumina_PE.md)
- [Dockstore](https://dockstore.org/workflows/github.com/marcoteix/gemstone-workflows/GEMSTONE_Plate_Swipes_Illumina_PE:main?tab=info)
- [Launch in Terra](https://app.terra.bio/#import-tool/dockstore/github.com/marcoteix/gemstone-workflows/GEMSTONE_Plate_Swipes_Illumina_PE:main)

This workflow processes paired-end Illumina reads from plate swipes/plate bacterial metagenomes. It could also be used for broader bacterial metagenomic analysis. It removes human reads, performs QA and QC, and metagenomic assembly. Furthermore, it estimates taxa abundances with Kraken2 and Bracken, does strain-level identification with StrainGE, bins MAGs, performs AMR genotyping, and identifies plasmid contigs.

![GEMSTONE_Plate_Swipes_Illumina_PE workflow diagram](https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE_Plate_Swipes_Illumina_PE.svg "GEMSTONE_Plate_Swipes_Illumina_PE workflow diagram")

## Authors

__Marco Teixeira__

<span style="font-size:80%; margin-top:0px"> mcarvalh@broadinstitute.org </span>

__Colin Worby__

<span style="font-size:80%; margin-top:0px"> cworby@broadinstitute.org </span>

<img src="https://github.com/marcoteix/gemstone-workflows/blob/main/docs/figures/GEMSTONE%20logo.png" alt="GEMSTONE logo" width=30%>
