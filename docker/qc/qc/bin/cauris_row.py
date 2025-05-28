from pathlib import Path
import cli

args = cli.CLI().qc_row_cauris()

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
    if args.contamination > args.max_contamination:
        qc_check = "FAIL"
        qc_note = "Contamination"
    elif args.completeness < args.min_completeness:
        # Set incomplete samples to "FAIL"
        qc_check = "FAIL"
        qc_note = "Low completeness"    
    elif (
        args.gambit_taxon != "Candidozyma auris" or \
        args.kraken2_taxon != "Candidozyma auris"
    ): 
        # If it is not Candidozyma auris, fail QC
        qc_check = "ALERT"
        qc_note = "Taxonomic mismatch"

# Write outputs
outdir = Path(args.output)
outdir.mkdir(parents=True, exist_ok=True)
outdir.joinpath("qc_check").write_text(qc_check)
outdir.joinpath("qc_note").write_text(qc_note)
