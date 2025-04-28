import argparse
import numpy as np
import pandas as pd 

parser = argparse.ArgumentParser(
    "Filter SNPs in MSA"
)

parser.add_argument(
    "input",
    type = str, help = "Input multiple sequence alignment in FASTA format."
)

parser.add_argument(
    "--min-n-samples", "-m",
    type = int, default = 2,
    help = "Minimum number of samples with the same base in a locus for it to be kept. Defaults to %(default)s."
)

parser.add_argument(
    "--output", "-o",
    type = str, required = False,
    help = "Output FASTA file. If not set, prints to stdout."
)

args = parser.parse_args()

# Read input FASTA
with open(args.input) as file:
    msa = file.read()

# Check how many samples there are
n_samples = msa.count("\n>")

# Correct the minimum number of samples if needed 
min_n_samples = np.minimum(
    args.min_n_samples,
    n_samples
)

# Iterate over samples, extract sequences
sequences = {
    x.split("\n")[0]: "\n".join(x.split("\n")[1:])
    for x in msa.split("\n>")
}

# Holds filtered sequences
sequences_out = {
    k: ""
    for k in sequences.keys()
} 

# Get the length of the MSA
n_bases = len(
    list(sequences.values())[0]
)

for i in range(n_bases):

    # Get the nucleotide at this position for all samples
    nucleotides = pd.Series(
        {
            k: v[i]
            for k, v in sequences.items()
        },
        name = "nucleotides"
    )

    # Count the number of occurences for each nucleotide
    counts = nucleotides.value_counts()

    # Replace nucleotides appearing < min-n-samples times with Ns
    for sample in sequences.keys():
        sequences_out[sample] += (
            "N"
            if counts.loc[sequences[sample][i]] < min_n_samples
            else sequences[sample][i]
        )

# Write to output
output_msa = "\n".join(
    [ 
        ">" + k + "\n" + v 
        for k,v in sequences_out.items()
    ]
)

if args.output:
    
    with open(args.output, "w") as file:
        file.write(output_msa)

else:

    print(output_msa)