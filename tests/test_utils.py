"""
Test jwst_gtvt utils.py module
"""

import pytest
from jwst_gtvt.utils import check_jwst_instrument_name


@pytest.mark.parametrize(
    "instrument_name",
    [
        "nircam",
        "nirspec",
        "miri",
        "fgs",
        "v3pa",
        pytest.param("acs", marks=pytest.mark.xfail(reason="Hubble Instrument")),
        pytest.param("cos", marks=pytest.mark.xfail(reason="Hubble Instrument")),
        pytest.param("wfc3", marks=pytest.mark.xfail(reason="Hubble Instrument")),
        pytest.param("stis", marks=pytest.mark.xfail(reason="Hubble Instrument")),
        pytest.param("dsfkjhsdkjfg,/ew3765903", marks=pytest.mark.xfail(reason="Random characters")),
    ],
)
def test_check_jwst_instrument_name(instrument_name):
    check_jwst_instrument_name(instrument_name)
