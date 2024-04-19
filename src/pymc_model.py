"""Script to generate PyMC results for Synechococcus."""

import gzip
from pathlib import Path

import cloudpickle
import cobra
import emll
import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt

HERE = Path(__file__).parent.resolve()
ROOT = HERE.parent.resolve()
MODEL = ROOT.joinpath("models/iJB785_no_zero_flux_flipped.xml")
DATA = ROOT.joinpath("data/processed_data")
METAB = DATA.joinpath("cleaned_metabolomics.csv")
PROT = DATA.joinpath("reduced_enzyme.csv")
VSTAR = DATA.joinpath("v_star_iJB785.csv")


class SynBMCA:
    """Class to run BMCA for the Synechococcus elongatus model."""

    def __init__(
        self,
        model_path,
        v_star_path,
        metabolite_concentrations_path,
        enzyme_measurements_path,
        reference_state,
        run_inference=True,
        # boundary_fluxes_path,
    ):
        """Placeholder Docstring."""
        self.model = cobra.io.read_sbml_model(model_path)
        self.v_star = pd.read_csv(v_star_path, header=None, index_col=0)[1]
        self.x = pd.read_csv(metabolite_concentrations_path, index_col=0)
        # self.v = pd.read_csv(boundary_fluxes_path, index_col=0)
        self.e = pd.read_csv(enzyme_measurements_path, index_col=0)

        self.ref_state = reference_state
        self.preprocess_data()
        self.build_pymc_model()

        # If only building PyMC model, set run_inference to False
        if self.run_inference:
            self.approx, self.hist = self.run_emll()
            self.save_results(self.approx, self.hist)

    def preprocess_data(self):
        """Read in cobra model as components."""
        # Establish compartments for reactions and metabolites
        self.r_compartments = [
            r.compartments if "e" not in r.compartments else "t" for r in self.model.reactions
        ]
        # TODO: Find why these were included in Hackett
        # self.r_compartments[self.model.reactions.index("SUCCt2r")] = "c"
        # self.r_compartments[self.model.reactions.index("ACt2r")] = "c"
        for rxn in self.model.exchanges:
            self.r_compartments[self.model.reactions.index(rxn)] = "t"
            self.m_compartments = [m.compartment for m in self.model.metabolites]

        # Reindex arrays to have the same column ordering
        to_consider = self.x.columns
        # self.v = self.v.loc[:, to_consider]
        self.x = self.x.loc[:, to_consider]
        self.e = self.e.loc[:, to_consider]

        self.n_exp = len(to_consider) - 1
        # Drop reference state
        # self.vn = (
        #     self.v.drop(self.ref_state)
        #     .divide(self.v.drop(self.ref_state)["flux"], axis=0)
        #     .drop("flux", axis=1)
        #     .T
        # )
        self.xn = (self.x.subtract(self.x[self.ref_state], 0) * np.log(2)).T
        self.en = self.e.T  # (2 ** np.abs(self.e.subtract(self.e[self.ref_state], 0))).T

        self.xn = self.xn.drop(self.ref_state)
        self.en = self.en.drop(self.ref_state)

        # Get indexes for measured values
        self.x_inds = np.array([self.model.metabolites.index(met) for met in self.xn.columns])
        self.e_inds = np.array([self.model.reactions.index(rxn) for rxn in self.en.columns])
        # self.v_inds = np.array([self.model.reactions.index(rxn) for rxn in self.vn.columns])

        self.e_laplace_inds = []
        self.e_zero_inds = []

        for i, rxn in enumerate(self.model.reactions):
            if rxn.id not in self.en.columns:
                if ("e" not in rxn.compartments) and (len(rxn.compartments) == 1):
                    self.e_laplace_inds += [i]
                else:
                    self.e_zero_inds += [i]

        self.e_laplace_inds = np.array(self.e_laplace_inds)
        self.e_zero_inds = np.array(self.e_zero_inds)
        self.e_indexer = np.hstack([self.e_inds, self.e_laplace_inds, self.e_zero_inds]).argsort()

        self.N = cobra.util.create_stoichiometric_matrix(self.model)
        self.Ex = emll.util.create_elasticity_matrix(self.model)
        self.Ey = emll.util.create_Ey_matrix(self.model)

        self.Ex *= 0.1 + 0.8 * np.random.rand(*self.Ex.shape)
        self.v_star = abs(self.v_star)
        self.ll = emll.LinLogLeastNorm(self.N, self.Ex, self.Ey, self.v_star.values, driver="gelsy")

    def build_pymc_model(self):
        """Build the PyMC probabilistic model."""
        with pm.Model() as pymc_model:
            # Priors on elasticity values
            self.Ex_t = pm.Deterministic(
                "Ex",
                emll.util.initialize_elasticity(
                    self.ll.N,
                    b=0.01,
                    sigma=1,
                    alpha=None,
                    m_compartments=self.m_compartments,
                    r_compartments=self.r_compartments,
                ),
            )

            self.Ey_t = pt.as_tensor_variable(self.Ey)

            e_measured = pm.Normal(
                "log_e_measured",
                mu=np.log(self.en),
                sigma=0.2,
                shape=(self.n_exp, len(self.e_inds)),
            )
            e_unmeasured = pm.Laplace(
                "log_e_unmeasured", mu=0, b=0.1, shape=(self.n_exp, len(self.e_laplace_inds))
            )
            log_en_t = pt.concatenate(
                [e_measured, e_unmeasured, pt.zeros((self.n_exp, len(self.e_zero_inds)))], axis=1
            )[:, self.e_indexer]

            pm.Deterministic("log_en_t", log_en_t)

            # Priors on external concentrations
            yn_t = pm.Normal(
                "yn_t",
                mu=0,
                sigma=10,
                shape=(self.n_exp, self.ll.ny),
                initval=0.1 * np.random.randn(self.n_exp, self.ll.ny),
            )

            chi_ss, _ = self.ll.steady_state_pytensor(self.Ex_t, self.Ey_t, pt.exp(log_en_t), yn_t)
            pm.Deterministic("chi_ss", chi_ss)
            # pm.Deterministic("vn_ss", vn_ss)

            # log_vn_ss = pt.log(pt.clip(vn_ss[:, self.v_inds], 1e-8, 1e8))
            # log_vn_ss = pt.clip(log_vn_ss, -1.5, 1.5)

            chi_clip = pt.clip(chi_ss[:, self.x_inds], -1.5, 1.5)

            chi_obs = pm.Normal(
                "chi_obs",
                mu=chi_clip,
                sigma=0.2,
                observed=self.xn.clip(
                    lower=-1.5,
                    upper=1.5,
                ),
            )

            # log_vn_obs = pm.Normal(
            #     "vn_obs",
            #     mu=log_vn_ss,
            #     sigma=0.1,
            #     observed=np.log(self.vn).clip(lower=-1.5, upper=1.5),
            # )

        self.pymc_model = pymc_model

    def run_emll(self):
        """Build linlog model and run inference."""
        with self.pymc_model:
            approx = pm.ADVI()
            hist = approx.fit(
                n=80000,
                obj_optimizer=pm.adagrad_window(learning_rate=0.005),
                total_grad_norm_constraint=100,
            )

            # trace = hist.sample(500)
            # ppc = pm.sample_ppc(trace)

        return approx, hist

    def save_model_data(self):
        """Save model input data and pymc model."""
        with gzip.open("Model_Data.pgz", "wb") as f:
            cloudpickle.dump(
                {
                    "model": self.pymc_model,
                    "cobra": self.model,
                    "xn": self.xn,
                    "en": self.en,
                },
                f,
            )

    def save_results(self, approx, hist):
        """Save ADVI results in cloudpickle."""
        with gzip.open("Syn_80k.pgz", "wb") as f:
            cloudpickle.dump(
                {
                    "model": self.pymc_model,
                    "approx": approx,
                    "hist": hist,
                },
                f,
            )


def main():
    """Run SynBMCA for default case."""
    ref_state = "Se_axen_d8_3"

    SynBMCA(
        MODEL,
        VSTAR,
        METAB,
        PROT,
        ref_state,
        run_inference=False,
    )


if __name__ == "__main__":
    main()