""" JWST GENERAL AND MOVING TARGET VISIBILITY TOOL TESTS"""

from jwst_gtvt.find_tgt_info import main as jwst_gtvt_main
from jwst_gtvt.find_tgt_info import get_target_ephemeris

class mimic_parser:
    """ This tool uses argparse. This class lets us mimic the calls to the argparser (i.e. args.variable).
    in the source code. Parsing the arguments makes it difficult to write pytests, this is the easiest near 
    term solution.
    """
    def __init__(self, ra=None, dec=None, desg=None, fixed=True):
        """Class initializer
        
        Parameters
        ----------
        ra : str
            right acsention of target in hours:minutes:seconds or sexigesimal coords
        dec : str
            declination of target in hours:minutes:seconds or sexigesimal coords
        desg : str
            desigantion (target name) in horizons database (for moving targets)
        fixed : bool
            boolean argument for moving target support. True will use the GTVT and False will use MTVT.
        """

        if fixed:
            self.ra = ra
            self.dec = dec

        if not fixed:
            self.desg = desg
            self.smallbody = True


        self.save_table = None
        self.no_verbose = True
        self.start_date = '2021-12-19'
        self.end_date = '2023-12-19'
        self.v3pa = None
        self.instrument = None
        self.name = None
        self.save_plot = 'test'

def test_ra_dec_hour_min_sec():
    """Test end-to-end for jwst_gtvt fixed target given ra and dec in hours:minutes:seconds"""

    args = mimic_parser(ra='16:52:58.9', dec='02:24:03')
    jwst_gtvt_main(args)


def test_ra_dec_sexigesimal():
    """Test end-to-end for jwst_gtvt fixed target given ra and dec in sexigesimal coords"""

    args = mimic_parser(ra='253.2458', dec='2.4008')
    jwst_gtvt_main(args)


def test_moving_target():
    """Test end-to-end for moving target portion of tool"""
    
    args = mimic_parser(fixed=False, desg=['Ceres'])

    name, args.ra, args.dec = get_target_ephemeris(
        ' '.join(args.desg), args.start_date, args.end_date, smallbody=args.smallbody)
    if args.name is None:
        args.name = name

    jwst_gtvt_main(args, fixed=False)
