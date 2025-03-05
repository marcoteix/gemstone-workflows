#%%
from pathlib import Path
import cli

args = cli.CLI().qc_row()

qc_check, qc_note = "PASS", ""
if args.raw_read_screen != "PASS" or args.clean_read_screen != "PASS":
    # If the isolate fails the raw or clean read QC, it should fail global QC
    qc_check = "FAIL"
    qc_note = "Low yield/quality"
elif args.coverage < args.min_coverage:
    # If the isolate does not meet the minimum coverage, it should fail global QC
    qc_check = "FAIL"
    qc_note = "Low coverage"
else:
    # Extract the genus from the lab predicted taxon
    lab_genus = args.lab_species.split(" ")[0] if args.lab_species is not None else "unknown"
    if lab_genus == "Klebsiealla": lab_genus="Klebsiella"
    
    # Ignore Candida
    if lab_genus != "Candida": 
        # Set contaminated samples to "FAIL"
        if args.contamination > args.max_contamination:
            qc_check = "FAIL"
            qc_note = "Contamination"
        elif args.completeness < args.min_completeness:
            # Set incomplete samples to "FAIL"
            qc_check = "FAIL"
            qc_note = "Low completeness"
        elif not lab_genus in ["unknown", ""]:
            # Compare the lab species to the GAMBIT predicted taxon
            # Extract the genus from the GAMBIT predictions
            gambit = args.gambit_taxon.split(" ")[0]
            # If the lab and GAMBIT predictions do not match, set the qc_flag to ALERT
            if gambit != lab_genus:
                qc_check = "ALERT"
                qc_note = "Taxonomic mismatch"

if qc_check == "PASS":
    qc_note = ""

# Write outputs
outdir = Path(args.output)
outdir.mkdir(parents=True, exist_ok=True)
outdir.joinpath("qc_check").write_text(qc_check)
outdir.joinpath("qc_note").write_text(qc_note)

# %%
