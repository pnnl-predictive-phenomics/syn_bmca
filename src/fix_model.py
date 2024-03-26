"""Model to modify cobra model for use with YMCA."""

from pathlib import Path

import cobra
import pandas as pd

"""
This code must accomplish the following objectives:
___________________________________________________

1. Load in a cobra model
2. Use that cobra model to run model.optimize()
3. Use the calculated fluxes to drop 0 fluxes from model
4. Use the calculated fluxes to swap the direction of reactions.
5. Rerun the optimize to ensure that all fluxes are positive non-zero
"""

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent.resolve()
MODEL = ROOT.joinpath("models/syn_elong.xml")


def calculate_fluxes(model: cobra.core.model.Model) -> pd.DataFrame:
    """Loads in model and calculates fluxes."""
    solution = model.optimize()
    fluxes = solution.fluxes.to_frame()

    return fluxes


def drop_zero_fluxes(model: cobra.core.model.Model, fluxes: pd.DataFrame) -> cobra.core.model.Model:
    """Drops zero fluxes from the model."""
    zero_fluxes = fluxes[fluxes.any(axis=1) == 0].index.to_list()
    model.remove_reactions(zero_fluxes)

    return model


def correct_negative_fluxes():
    """Sets direction of reactions in model to flip negative."""
    pass


def main():
    """Runs script."""
    model = cobra.io.read_sbml_model(MODEL)
    fluxes = calculate_fluxes(model)
    model = drop_zero_fluxes(model, fluxes)
    new_fluxes = calculate_fluxes(model)

    return new_fluxes


if __name__ == "__main__":
    new_fluxes = main()
