"""Model to modify cobra model for use with YMCA."""

from pathlib import Path

import cobra
import numpy as np
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
MODEL = ROOT.joinpath("models/iJB785_w_sucrose_transport.xml")


def calculate_fluxes(model: cobra.core.model.Model) -> pd.DataFrame:
    """Loads in model and calculates fluxes."""
    solution = model.optimize()
    fluxes = solution.fluxes.to_frame()

    return fluxes


def drop_zero_fluxes(model, fluxes: pd.DataFrame) -> cobra.core.model.Model:
    """Finds zero fluxes from the model."""
    zero_fluxes = fluxes[fluxes.any(axis=1) == 0].index.to_list()
    model.remove_reactions(zero_fluxes)

    return model


def flip_reaction(model: cobra.core.model.Model, reaction_id: str) -> None:
    """Flips the reaction direction for a given model and reaction id."""
    reaction = model.reactions.get_by_id(reaction_id)
    if "<=>" in reaction.reaction:
        left, right = reaction.reaction.split(" <=> ")
        reaction.reaction = right + " <=> " + left
    elif "-->" in reaction.reaction:
        left, right = reaction.reaction.split(" --> ")
        reaction.reaction = right + " <-- " + left
    elif "<--" in reaction.reaction:
        left, right = reaction.reaction.split(" <-- ")
        reaction.reaction = right + " --> " + left

    model.repair()
    reaction.upper_bound = max(np.abs(reaction.lower_bound), reaction.upper_bound)
    reaction.lower_bound = 0
    reaction.update_variable_bounds()


# def create_slack_variables(model: cobra.core.model.Model, reactions: list[str]):
#     """Creates slack variables and adds them to objective."""
#     for reaction in reactions:
#         rxn = model.reactions.get_by_id(reaction)
#         slack_var = model.problem.Variable(f"{rxn.id}_slack", lb=0)
#         model.add_cons_vars(slack_var)

#         constraint = model.problem.Constraint(
#             rxn.flux_expression - slack_var,
#             ub=0,
#         )
#         model.add_cons_vars(constraint)
#     return model


def main() -> None:
    """Runs script."""
    model = cobra.io.read_sbml_model(MODEL)
    fluxes = calculate_fluxes(model)
    model = drop_zero_fluxes(model, fluxes)
    # model = create_slack_variables(model=model, reactions=fluxes)
    new_fluxes = calculate_fluxes(model)
    cobra.io.write_sbml_model(model, ROOT / "models/iJB785_no_zero_flux.xml")
    new_fluxes.to_csv("v_star_iJB785.csv")


if __name__ == "__main__":
    main()
