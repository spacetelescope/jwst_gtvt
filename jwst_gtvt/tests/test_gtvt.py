""" JWST GENERAL AND MOVING TARGET VISIBILITY TOOL TESTS"""
import pytest

from jwst_gtvt.jwst_tvt import Ephemeris
from jwst_gtvt.display_results import display_results
from jwst_gtvt.plotting import plot_visibility

def test_moving_target():
    """Test end-to-end for moving target portion of tool"""
    eph = Ephemeris()
    eph.get_moving_target_positions('Ceres',False)


def test_ra_dec_hour_min_sec():
    """Test end-to-end for jwst_gtvt fixed target given ra and dec in hours:minutes:seconds"""
    eph = Ephemeris()
    ra, dec = '16:52:58.9', '02:24:03'
    eph.get_fixed_target_positions(ra, dec)


def test_ra_dec_sexigesimal():
    """Test end-to-end for jwst_gtvt fixed target given ra and dec in sexigesimal coords"""
    eph = Ephemeris()
    ra, dec = '253.2458', '2.4008'
    eph.get_fixed_target_positions(ra, dec)


@pytest.mark.parametrize("instrument", ['MIRI', 'NIRCAM', 'NIRSPEC', 'NIRISS', 'FGS', 'V3PA'])
@pytest.mark.parametrize("mode", ['fixed', 'moving'])
def test_gtvt_single_instrument(instrument, mode):
    eph = Ephemeris()
    if mode == 'fixed':
        ra, dec = '253.2458', '2.4008'
        eph.get_fixed_target_positions(ra, dec)
    else:
        eph.get_moving_target_positions('Ceres',False)

    plot_visibility(eph, instrument)