"""Test of fba_utils."""

import numpy as np
import pandas as pd
import pytest
from syn_bmca.fba_utils import prepare_data_for_bmca


@pytest.fixture(scope="module")
def test_prepare_data_for_bmca():
    """Test the data preparation step."""
    # Initialize test data
    all_conditions = ["condition_1", "condition_2", "condition_3"]
    measured_data = pd.DataFrame({
        "gene_1": [1, np.Inf, 3],
        "gene_2": [4, 5, 6],
        "gene_3": [7, 8, 9],
    })
    unmeasured_variables = ["gene_4", "gene_5"]
    unmapped_variables = ["gene_6", "gene_7"]

    # Run function
    result = prepare_data_for_bmca(
        all_conditions, measured_data, unmeasured_variables, unmapped_variables
    )

    # Initialize expected result
    expected_result = pd.DataFrame(
        {
            "gene_1": [1, np.Inf, 3],
            "gene_2": [4, 5, 6],
            "gene_3": [7, 8, 9],
            "gene_4": [np.Inf, np.Inf, np.Inf],
            "gene_5": [np.Inf, np.Inf, np.Inf],
            "gene_6": [np.nan, np.nan, np.nan],
            "gene_7": [np.nan, np.nan, np.nan],
        },
        index=["condition_1", "condition_2", "condition_3"],
    )

    # Compare
    assert result.equals(expected_result)
