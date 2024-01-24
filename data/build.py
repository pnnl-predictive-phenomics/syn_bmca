"""File to process Omics data collected for Synechococcus elongatus."""

import logging
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import numpy as np
import json

HERE = Path(__file__).parent.resolve()
RAW_DATA = HERE.joinpath("raw_data")
KEGG_TO_BIGG_MAP = RAW_DATA.joinpath("KEGG_to_BIGG_mapper_for_Syn_Rt.json")
ZIP_DATA = RAW_DATA.joinpath("Selon_omics_Pavlo.zip")
UNZIPPED = RAW_DATA.joinpath("Selon_omics_Pavlo")
OUTPUT = HERE.joinpath("processed_data")
METAB = UNZIPPED.joinpath("Metabolomics")
TRANS = UNZIPPED.joinpath("Transcript")
EXP_ID_MAP = METAB.joinpath("Metaboites_Se_Rt_IDs_9day.xlsx")


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
    # df.index = [label.split(")", 1)[-1] for label in df.index]
    logging.info("Loaded Metabolomics Data")
    return df


def _load_transcriptomics() -> pd.DataFrame:
    """Loads transcriptomic data."""
    logging.info("Parsing Transcriptomics Data...")
    df = pd.read_excel(next(TRANS.glob("*.xlsx")), index_col=0)
    logging.info("Loaded Transcriptomics Data")
    return df

# Function to extract time points from Transcriptomics data
def _extract_timepoints(this_df: pd.DataFrame) -> list:
    """Extracts the time points from the dataframe columsn
    Note: expects X_T_N format, 
    where, 
        X is the prefix to the experimental condition and measurement
        T is the time point
        N is the measurement replicate at time point T
    """
    col_timepts = np.sort(np.unique(np.array([int(cn.split("_")[-2]) for cn in list(this_df.columns.values)])))
    return list(col_timepts)


def _clean_metabolomics(metabolomics_df: pd.DataFrame, transcriptomics_timepts: list) -> pd.DataFrame:
    """Reduces the metabolomics data

    We remove metabolomics data based on the criteria:
    1. Unnamed metabolites from data
    2. Time points of data that don't match transcriptomics data

    We also rename rows and columns using the following criteria:
    1. Rename rows using BIGG-IDs
    2. Rename transctiptomics columns to match metabolomics column names that are kept [done in `_clean_transciptomics()` function]
    """
    
    ## Remove columns that don't match transcriptomics data time points 
    # Get time points in hours from metabolomics data column names
    metab_col_timepts = [24.0*int(cn.split("_")[-2][1]) for cn in list(metabolomics_df.columns.values)]

    # Create mapping between transcriptomics time points and the closest metabolomics time points
    trans_to_metab_time_map = dict(zip(transcriptomics_timepts, [min(metab_col_timepts, key=lambda x:abs(x-y)) for y in transcriptomics_timepts]))

    # Keep metabolomics data columns closest to transcriptomics time points
    col_keep_idx = [(x in trans_to_metab_time_map.values()) for x in metab_col_timepts]
    metabolomics_df = metabolomics_df.iloc[:, col_keep_idx]

   
    ## Rename rows using BIGG-IDs
    # Get KEGG_to_BIGG mapping
    with open(KEGG_TO_BIGG_MAP) as f:
        kegg_to_bigg_map = json.load(f)
    
    # Create mapping from short names in experiment to BIGG IDs
    exp_map = pd.read_excel(EXP_ID_MAP, skiprows=1, header=0, index_col=0)

    exp_to_bigg_map = {}
    count_missing_IDs = 0
    for idx in exp_map.index:
        # Only include known metabolites with higher certainty
        if not('*' in idx) and not('unk-' in idx):
            kid = exp_map.loc[idx,'Kegg']
            bid = exp_map.loc[idx,'BiGG']
            # shortname = exp_map.loc[idx,'Abbr name']
            if not(kid=='' or pd.isna(kid)):
                if kid in kegg_to_bigg_map:
                    # exp_to_bigg_map[shortname] = kegg_to_bigg_map[kid]
                    exp_to_bigg_map[idx] = kegg_to_bigg_map[kid]
                elif not(bid=='' or pd.isna(bid)):
                    # exp_to_bigg_map[shortname] = bid
                    exp_to_bigg_map[idx] = bid
                else:
                    print(f"{kid} not found in any mappings")
                    count_missing_IDs+=1
    
    print(f"Missing IDs for {count_missing_IDs} metabolites. These rows will be removed for now")
    
    # Cleanup metabolics df by renaming index is KEGG and BIGG IDs are available, drop row otherwise
    
    print(metabolomics_df)
    
    print('---')
    for k,v in exp_to_bigg_map.items():
        print(k,v)

    print('---')
    # new_index = []
    # for idx in metabolomics_df.index:
    #     if idx in exp_to_bigg_map:
    #         new_index.append(exp_to_bigg_map[idx])
    #     else:
    #         metabolomics_df.drop(idx)
    
    keep_row_map = {k:v for k,v in exp_to_bigg_map.items() if k in metabolomics_df.index}
    print(keep_row_map)
    print('---')
    for idx in metabolomics_df.index:
        if idx in keep_row_map:
            print(idx, keep_row_map[idx])
    print('---')
    
    print(metabolomics_df.index)
    print('---')

    keep_rows = [k for k in exp_to_bigg_map.keys() if k in metabolomics_df.index]
    print("Keep rows")
    print(keep_rows)

    metabolomics_df = metabolomics_df.loc[[k for k in exp_to_bigg_map.keys() if k in metabolomics_df.index]]
    metabolomics_df = metabolomics_df.rename(index=exp_to_bigg_map)



    return metabolomics_df



def _clean_transcriptomics(transcriptomics_df: pd.DataFrame, exp_prefix) -> pd.DataFrame:
    """Reduces the transcriptomics data

    We remove transcriptomics data based on the criteria:
    1. Remove rows with unknown transcript
    """

    ## Remove rows with unknown transcript
    

    return reduced_trans_df

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

    print("Metabolomics data:")
    print(metabolomics)

    print("Transcriptomics data:")
    print(transcriptomics)
    
    
    metabolomics.to_csv(OUTPUT.joinpath("metabolomics.csv"))
    transcriptomics.to_csv(OUTPUT.joinpath("transcriptomics.csv"))

    # Get time points from Transcriptomics data
    transcript_timepts = _extract_timepoints(transcriptomics)

    print("Timepoints from Transcriptomics data:")
    print(transcript_timepts)

    # Clean Metabolomics data
    reduced_metabolomics = _clean_metabolomics(metabolomics, transcript_timepts)
    
    print("Reduced Metabolomics data:")
    print(reduced_metabolomics)

if __name__ == "__main__":
    main()
