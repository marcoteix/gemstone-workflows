# %%
import cli
import pandas as pd
from terra import get_entity, upload_entity

args = cli.CLI().qc_table()

# Load table
X = get_entity(args.table, args.workspace)

# Ignore samples without a raw_read_screen (no QC analysis)
X = X.dropna(subset=["raw_read_screen"])

qc_colnames = ["qc_check", "qc_note", "qc_taxonomy_check"]
qc = pd.DataFrame(
    [["PASS", "", ""]]*len(X), 
    columns=qc_colnames, 
    index=X.index
)
# If the isolate fails the raw or clean read QC, it should fail global QC
qc.loc[(X.raw_read_screen.ne("PASS") | X.clean_read_screen.ne("PASS")) 
    & qc.qc_check.eq("PASS"), qc_colnames] = ["FAIL", "Low yield/quality", ""]
# If the isolate does not meet the minimum coverage, it should fail global QC
qc.loc[X.est_coverage_clean.lt(args.min_coverage) 
    & qc.qc_check.eq("PASS"), qc_colnames] = ["FAIL", "Low coverage", ""]

# Extract the genus from the species names
genera = X.species.apply(lambda x: x.split(" ")[0]) \
    .replace({"Klebsiealla": "Klebsiella"})

# Set contaminated samples to "FAIL", but ignore C. auris
qc.loc[~genera.eq("Candida") & X.checkm2_contamination.astype(float).gt(args.max_contamination)
    & qc.qc_check.eq("PASS"), qc_colnames] = ["FAIL", "Contamination", ""]

# Set incomplete samples to "FAIL", but ignore C. auris
qc.loc[~genera.eq("Candida") & X.checkm2_completeness.astype(float).lt(args.min_completeness)
    & qc.qc_check.eq("PASS"), qc_colnames] = ["FAIL", "Low completeness", ""]

# Compare the lab species to the GAMBIT predicted taxon
# Extract the genus from the GAMBIT predictions
gambit = X.gambit_predicted_taxon.apply(lambda x: x.split(" ")[0] 
                                    if not pd.isna(x) else pd.NA)
# Ignore missing lab predictions and Candida; in case of lab and GAMBIT predictions
# not matching, set the qc_flag to ALERT and set the qc_taxonomy_flag to QC_ALERT
qc.loc[~genera.isin(["unknown", "Candida"]) & gambit.ne(genera, fill_value="nan") 
    & qc.qc_check.eq("PASS"), qc_colnames] = ["ALERT", "Taxonomic mismatch", "QC_ALERT"]

# Samples that pass QC should not have a QC note
qc.loc[qc.qc_check.eq("PASS"), "qc_note"] = pd.NA

# [DANGER ZONE] Upload QC rows back to the data table
upload_entity(qc, args.table, args.workspace)
# %%
