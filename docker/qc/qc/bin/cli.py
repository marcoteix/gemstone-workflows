import argparse
from typing import Literal

class CLI:
    """
    Class for parsing arguments from the command line.
    """

    def __init__(self) -> None: pass

    def qc_row(self) -> argparse.Namespace:
        self.__parser(
            name="QC flags - input as a row",
            description="Takes QC metrics as input, generates QC flags, and writes two files to \
the output directory: qc_check (containing a final QC flag) and qc_note (reason for failure)."
        )
        self.__qc_row_args()
        self.__qc_args()
        return self.__cast_qc_row(self.parser.parse_args())

    def qc_table(self) -> argparse.Namespace:

        self.__parser(
            name="QC flags - input as a table",
            description="Takes an entity data table as input, generates QC flags for all samples, and writes the \
output back into the input data table in the columns qc_check, qc_note, and qc_taxonomy_check."
        )
        self.__table_args()
        self.__qc_args()

        return self.parser.parse_args()
    
    def status(self) -> argparse.Namespace:
        
        self.__parser(
            name="Status",
            description="Takes an entity data table as input, generates status flags, and writes the \
output back into the input data table in the columns status and preferred_sample_id."
        )
        self.__table_args()
        return self.parser.parse_args()

    def __parser(self, name: str, description: str, **kwargs):
        self.parser = argparse.ArgumentParser(name, description=description, **kwargs)

    def __table_args(self):

        self.io_group = self.parser.add_argument_group("Input/Output")
        self.io_group.add_argument("--table", "-t", required=True, type=str, 
            help="Name of the input entity data table to use as input. The results will also be \
written to this data table, overwritting existing values.")
        self.io_group.add_argument("--workspace", "-w", required=True, type=str, 
            choices=["DSSC", "test"], help="Workspace containing the input/output entity data table.")

    def __qc_row_args(self):

        self.io_group = self.parser.add_argument_group("Input/Output")
        self.io_group.add_argument("--raw-read-screen", "-rs", type=str, required=True,
help="Column \"raw_read_screen\".")
        self.io_group.add_argument("--clean-read-screen", "-cs", type=str, required=True,
help="Column \"clean_read_screen\".")
        self.io_group.add_argument("--coverage", "-x", type=str, required=True,
help="Column \"est_coverage_clean\".")
        self.io_group.add_argument("--lab-species", "-ls", type=str, required=True,
help="Column \"species\".")
        self.io_group.add_argument("--gambit-taxon", "-g", type=str, required=True,
help="Column \"gambit_predicted_taxon\".")
        self.io_group.add_argument("--contamination", "-c", type=str, required=True,
help="Column \"checkm2_contamination\".")
        self.io_group.add_argument("--completeness", "-C", type=str, required=True,
help="Column \"checkm2_completeness\".")
        self.io_group.add_argument("--output", "-o", type=str, required=True,
help="Output directory.")
    
    def __qc_args(self):
        self.qc_group = self.parser.add_argument_group("QC options")
        self.qc_group.add_argument("--min-coverage", "-mx", type=float, default=40.0,
            help="Minimum estimated coverage of clean reads for a sample to pass QC. The default is 40.0 (40x).")
        self.qc_group.add_argument("--max-contamination", "-Mc", type=float, default=5.0,
            help="Maximum contamination as estimated by checkM2, as a percentage, for a sample to pass QC. The \
default is 5.0 (meaning 5 pct).")
        self.qc_group.add_argument("--min-completeness", "-mc", type=float, default=80.0,
            help="Minimum genome completeness as estimated by checkM2, as a percentage, for a sample to pass QC. \
The default is 80.0 (meaning 80 pct).")
        
    def __cast(self, val, type, if_fail):

        try: return type(val)
        except Exception: return if_fail

    def __cast_qc_row(self, args: argparse.Namespace) -> argparse.Namespace:

        args.coverage = self.__cast(args.coverage, float, 0.0)
        args.contamination = self.__cast(args.contamination, float, 100.0)
        args.completeness = self.__cast(args.completeness, float, 0.0)
        return args