version 1.0

task vcf_to_msa {
  input {
    Array[File] vcfs
    String collection_name = "variants"
    String filters = "PASS,."
    Int min_samples = 1
    Int memory = 8
    Int disk_size = 16
  }
  command <<<

    mkdir vcfs

    # Filter variants based on the CleanSweep filter
    for vcf in ~{sep=' ' vcfs}; do

        echo "Filtering and indexing $vcf..."

        bcftools view $vcf \
            -f ~{filters} \
            -o vcfs/${vcf%.vcf}.pass.vcf.gz \
            -O b

        # Index
        bcftools index vcfs/${vcf%.vcf}.pass.vcf.gz

    done

    # Make a file with cleansweep VCF files
    dir vcfs/*.pass.vcf.gz > filelist.txt

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

    echo "Generating MSA..."

    python /tmp/scripts/vcf2phylip.py -i ~{collection_name}.merged.vcf.gz \
        --output-folder "msa" \
        --output-prefix ~{collection_name} \
        -f -p -m ~{min_samples}

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