"""Check the catalogue-correction sign against real Eros astrometry and Horizons."""

import numpy as np
import pytest

from .validation.eros_catalogue_debiasing import CASES, FIXTURES, run_case, statistics


@pytest.mark.parametrize("case_name", CASES)
def test_eros_catalogue_debiasing_against_horizons(case_name):
    """The reader's correction improves two frozen Eros series; reversing it worsens both."""
    result = run_case(FIXTURES / case_name)
    summary = statistics(result)

    # Check the applied sign and sky-plane units against independently sampled JPL bias maps.
    golden_bias = result[["bias_ra_cosdec_arcsec", "bias_dec_arcsec"]].to_numpy()
    applied = result[["applied_ra_arcsec", "applied_dec_arcsec"]].to_numpy()
    np.testing.assert_allclose(applied, -golden_bias, atol=1e-10, rtol=0)

    # Every measurement is retained; the correction must be large enough to discriminate signs.
    assert len(result) == {"eros_t05_2016": 91, "eros_711_1996": 20}[case_name]
    assert np.isfinite(result.filter(regex="_arcsec$").to_numpy()).all()
    assert summary["correction_to_raw_rms"] > 0.4

    # Independently queried Horizons station angles bound the complete Tudat model error to 1 mas.
    assert summary["horizons_observer_difference_maximum_arcsec"] < 0.001

    # This is an empirical check against Horizons, with no fit or residual-based outlier rejection.
    assert summary["subtract"]["rms_arcsec"] < 0.8 * summary["raw"]["rms_arcsec"]
    assert summary["add"]["rms_arcsec"] > 1.3 * summary["raw"]["rms_arcsec"]

    # Correction application changes observations only: the three predictions must be identical.
    raw = result[["raw_ra_arcsec", "raw_dec_arcsec"]].to_numpy()
    corrected = result[["subtract_ra_arcsec", "subtract_dec_arcsec"]].to_numpy()
    reversed_sign = result[["add_ra_arcsec", "add_dec_arcsec"]].to_numpy()
    np.testing.assert_allclose(corrected - raw, applied, atol=2e-9, rtol=0)
    np.testing.assert_allclose(reversed_sign - raw, -applied, atol=2e-9, rtol=0)
