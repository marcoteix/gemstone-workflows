version 1.0

task vcf_to_msa {
  input {
    Array[String] samplenames
    Array[File] vcfs
    String collection_name = "variants"
    String filters = "PASS,."
    Int min_samples = 1
    Int memory = 8
    Int disk_size = 16
  }
  command <<<

    mkdir ./vcfs

    # Filter variants based on the CleanSweep filter. Iterate over the
    # sample names and VCF files simultaneously
    names=( ~{sep=' ' samplenames} )
    vcfs=( ~{sep=' ' vcfs} )

    # Keep track of sample names so we can change the sample names in 
    # the merged VCF
    touch samplenames.txt
    # Make a file with cleansweep VCF files
    touch filelist.txt

    for i in "${!names[@]}"; do

        name="${names[i]}"
        vcf="${vcfs[i]}"

        echo "Filtering and indexing $vcf..."

        bcftools view $vcf \
            -f ~{filters} \
            -o ./vcfs/$name.pass.vcf.gz \
            -O b

        # Index
        bcftools index ./vcfs/$name.pass.vcf.gz

        echo $name >> samplenames.txt
        echo $(pwd)/vcfs/$name.pass.vcf.gz >> filelist.txt

    done

    echo "First lines of filelist.txt:"
    echo $(head filelist.txt)

    echo "Merging vcfs..."

    bcftools merge \
        -l filelist.txt \
        -o ~{collection_name}.merged.vcf.gz \
        -O b \
        -f ~{filters} \
        --force-samples \
        --missing-to-ref

    # Change sample names in the merged VCF
    echo "Replacing sample names in the merged VCF..."

    bcftools view ~{collection_name}.merged.vcf.gz | \
        bcftools reheader -s samplenames.txt | \
        bcftools view -o ~{collection_name}.merged.vcf.gz -O b

    echo "Generating MSA..."

    python /tmp/scripts/vcf2phylip.py \
        -i ~{collection_name}.merged.vcf.gz \
        --output-folder "msa" \
        --output-prefix ~{collection_name} \
        -f -p -r -m ~{min_samples}

    # Get versions
    bcftools --version | head -1 > bcftools_version.txt
    python /tmp/scripts/vcf2phylip.py --version > vcf2phylip_version.txt

  >>>
  output {
    File msa = "msa/~{collection_name}.min~{min_samples}.fasta"
    String bcftools_version = read_string("bcftools_version.txt")
    String vcf2phylip_version = read_string("vcf2phylip_version.txt")
  }
  runtime {
    docker: "marcoteix/gemstone-utils:1.0.0"
    memory: memory + " GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}