# Script to generate metabolic abundance rates from time series data

# Use the original metabolomics data to better approximate rates
#   Note1:  Use the version processed in build.py that combines data 
#           from multiple sources into one file
#   Note2:  In order to calculate rate values as local as possible, 
#           we take average differences between metabolite abundances 
#           when possible, and take just a difference of consecutive 
#           measurements at the first/last time steps accordingly. 
#           This is calculated for each replicate in the data separately

from pathlib import Path
import pandas as pd
import numpy as np

HERE = Path(__file__).parent.resolve()
OUTPUT = HERE.joinpath("processed_data")
METAB = OUTPUT.joinpath("metabolomics.csv")
RATES = OUTPUT.joinpath("calculated_metabolomic_abundance_rates.csv")

def _load_metabolomics_data() -> pd.DataFrame:
    return pd.read_csv(METAB, index_col='Sample')


def _calculate_rates(metab_df) -> pd.DataFrame:
    rates_df = pd.DataFrame(index=metab_df.index, columns=metab_df.columns)

    # Calculate rates for day 1: f(d2) - f(d1)
    day1_cols = [c for c in metab_df.columns if 'd1' in c]
    day2_cols = [c for c in metab_df.columns if 'd2' in c]
    for rep in range(1,4): # We know there are 3 replicates per day
        day1_c = [c for c in day1_cols if c[-1]==str(rep)]
        day2_c = [c for c in day2_cols if c[-1]==str(rep)]
        for i in range(len(day1_c)):
            rates_df.loc[:,day1_c[i]] = metab_df.loc[:,day2_c[i]] - metab_df.loc[:,day1_c[i]]
    
    # Calculate rates for day 2-8: Mean( f(dB) - f(dA), f(dC) - f(dB) )
    for day in range(2,9):
        dayA_cols = [c for c in metab_df.columns if f'd{day-1}' in c]
        dayB_cols = [c for c in metab_df.columns if f'd{day}' in c]
        dayC_cols = [c for c in metab_df.columns if f'd{day+1}' in c]
        for rep in range(1,4): # We know there are 3 replicates per day
            dayA_c = [c for c in dayA_cols if c[-1]==str(rep)]
            dayB_c = [c for c in dayB_cols if c[-1]==str(rep)]
            dayC_c = [c for c in dayC_cols if c[-1]==str(rep)]
            for i in range(len(dayB_c)):
                this_rate_df = pd.DataFrame({'Pre': metab_df.loc[:,dayB_c[i]] - metab_df.loc[:,dayA_c[i]],
                                             'Post': metab_df.loc[:,dayC_c[i]] - metab_df.loc[:,dayB_c[i]]})
                rates_df.loc[:,dayB_c[i]] = this_rate_df.mean(axis=1)
    
    # Calculate rates for day 9: f(d9) - f(d8)
    day8_cols = [c for c in metab_df.columns if 'd8' in c]
    day9_cols = [c for c in metab_df.columns if 'd9' in c]
    for rep in range(1,4): # We know there are 3 replicates per day
        day8_c = [c for c in day8_cols if c[-1]==str(rep)]
        day9_c = [c for c in day9_cols if c[-1]==str(rep)]
        for i in range(len(day9_c)):
            rates_df.loc[:,day9_c[i]] = metab_df.loc[:,day9_c[i]] - metab_df.loc[:,day8_c[i]]
    
    return rates_df

def main():
    # Load metabolomics data
    metabolomics_df = _load_metabolomics_data()
    # calcualte rates
    metab_rates_df = _calculate_rates(metabolomics_df)

    print(metab_rates_df)

    # save to file
    metab_rates_df.to_csv(RATES)


if __name__ == "__main__":
    main()