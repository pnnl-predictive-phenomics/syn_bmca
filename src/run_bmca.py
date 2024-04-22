"""Run script for SynBMCA."""

from pathlib import Path

from pymc_model import SynBMCA

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent.resolve()
MODEL = ROOT.joinpath("models/iJB785_w_sucrose_transport.xml")
DATA = ROOT.joinpath("data/processed_data")
METAB = DATA.joinpath("cleaned_metabolomics.csv")
PROT = DATA.joinpath("cleaned_enzyme.csv")
VSTAR = DATA.joinpath("v_star.csv")


def main():
    """Runs main function."""
    ref_state = "Se_axen_d6_1"

    SynBMCA(
        MODEL,
        VSTAR,
        METAB,
        PROT,
        ref_state,
    )


if __name__ == "__main__":
    main()
