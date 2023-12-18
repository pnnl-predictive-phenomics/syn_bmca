"""File to process Omics data collected for Synechococcus elongatus."""

import logging
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

HERE = Path(__file__).parent.resolve()
RAW_DATA = HERE.joinpath("raw_data")
ZIP_DATA = RAW_DATA.joinpath("Selon_omics_Pavlo.zip")
UNZIPPED = RAW_DATA.joinpath("Selon_omics_Pavlo")
OUTPUT = HERE.joinpath("processed_data")
METAB = UNZIPPED.joinpath("Metabolomics")
TRANS = UNZIPPED.joinpath("Transcript")

logging.basicConfig(level=logging.INFO)


def _unzip_data(zip_path: Path, out_dir: Path) -> list:
    """Unzip a zip file to the specified output directory.

    Parameters
    ----------
    zip_path: Path
        Path to the zip file.
    out_dir: Path
        Path to the directory where extracted files should be saved.

    Returns
    -------
        Nothing
    """
    with ZipFile(zip_path, "r") as zip_obj:
        zip_obj.extractall(out_dir)


def _load_metabolomics() -> pd.DataFrame:
    """Loads metabolomic data."""
    logging.info("Parsing Metabolomics Data...")
    df = pd.concat(
        [
            pd.read_csv(f, index_col=0)
            for f in METAB.glob("*.csv")
            if ("EMSL" in f.name or "JHU" in f.name) and "metadata" not in f.name
        ]
    )
    df.index = [label.split(")", 1)[-1] for label in df.index]
    logging.info("Loaded Metabolomics Data")
    return df


def _load_transcriptomics() -> pd.DataFrame:
    """Loads transcriptomic data."""
    logging.info("Parsing Transcriptomics Data...")
    df = pd.read_excel(next(TRANS.glob("*.xlsx")), index_col=0)
    logging.info("Loaded Transcriptomics Data")
    return df


def main():
    """Returns a placeholder."""
    if UNZIPPED.exists():
        logging.info("Directory '%s' found, proceeding...", UNZIPPED.name)
    else:
        logging.info("Directory '%s' not found, unzipping...", UNZIPPED.name)
        _unzip_data(ZIP_DATA, RAW_DATA)
        logging.info("Unzipping complete, proceeding...")

    metabolomics = _load_metabolomics()
    transcriptomics = _load_transcriptomics()

    metabolomics.to_csv(OUTPUT.joinpath("metabolomics.csv"))
    transcriptomics.to_csv(OUTPUT.joinpath("transcriptomics.csv"))


if __name__ == "__main__":
    main()
