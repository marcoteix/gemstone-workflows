#%%
from terra import get_entity, upload_entity
import pandas as pd
import cli 

args = cli.CLI().status()

# Load table
X = get_entity(args.table, args.workspace)

#%% Find the preferred sample (passed all QC and has the highest N50) for each stock
# and strain
preferred = X[
    X.raw_read_screen.eq("PASS", fill_value="-") &
    X.clean_read_screen.eq("PASS", fill_value="-") &
    X.qc_check.eq("PASS", fill_value="-")
].groupby(["stock_id", "straingst_top_strain"]).n50_value.idxmax() \
    .rename("preferred_sample_id")

X = X.drop(columns="preferred_sample_id", errors="ignore") \
    .join(preferred, on=["stock_id", "straingst_top_strain"])

# If a sample fails either raw or clean read QC, set status to fail or resequence
fail_yield = ~(X.raw_read_screen.eq("PASS") & X.clean_read_screen.eq("PASS"))
X.loc[fail_yield & X.preferred_sample_id.isna(), "status"] = "resequence"
X.loc[fail_yield & ~X.preferred_sample_id.isna(), "status"] = "fail"

# Otherwise, if it fails the taxonomy/completeness/contamination check, set status
# to fail or resequence
fail_contam = (~fail_yield & ~X.qc_check.eq("PASS"))
X.loc[fail_contam & X.preferred_sample_id.isna(), "status"] = "uci_check"
X.loc[fail_contam & ~X.preferred_sample_id.isna(), "status"] = "fail"

# If it passes all QC flags and it is not the preferred sample, set status to 
# duplicate
X.loc[~fail_yield & ~fail_contam & ~X.preferred_sample_id.eq(
    X.index.to_series(), fill_value="-"), "status"] = "duplicate"
X.loc[preferred.values, "status"] = "preferred"

# Ignore non-processed samples
X = X[~X.analysis_date.isna()]

# Processed samples 
#%% [DANGER ZONE] Upload QC rows back to the data table
upload_entity(X[["status", "preferred_sample_id"]], args.table, args.workspace)


# %%
