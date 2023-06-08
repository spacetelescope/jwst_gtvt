usage = """

Main driver for moving target support. This script will display and
generate the figure for 

Usage:
  jwst_mtvt <desg> [--start_date=<obs_start>] [--end_date=<obs_end>] [--instrument=<inst>] [--write_ephemeris=<write_path>] [--write_plot=<plot_path>] [--silent]

Arguments:
  <desg>    Name of moving target. 

Options:
  [--start_date]         Start date for plot
  [--end_date]           End date for plot
  [--instrument]         JWST instrument to plot individually
  [--write_ephemeris]    File name to write ephemeris to
  [--write_plot]         File name to write plot out to
  --silent               Boolean to print results to screen [default: False]
  -h --help              Show this screen.
  --version              Show version.
"""

from astropy.time import Time

from jwst_gtvt.display_results import display_results
from docopt import docopt
from jwst_gtvt.jwst_tvt import Ephemeris
from jwst_gtvt.plotting import plot_visibility


def main(args):
    if args['--start_date'] and args['--end_date']:
        start = Time(args['--start_date'])
        end = Time(args['--end_date'])
        if start > end:
            raise ValueError("Date mismatch, end date is before start date.")
        else:
            eph = Ephemeris(start_date=start, end_date=end)
    else:
        eph = Ephemeris()

    eph.get_moving_target_positions(args['<desg>'])

    if not eph.dataframe['in_FOR'].any():
        in_FOR_msg = ("No position angles in field of regard! "
                      "Check constraints for your target and if it is observable with JWST. \n"
                      "Vist: https://jwst-docs.stsci.edu/jwst-observatory-characteristics/jwst-observatory-coordinate-system-and-field-of-regard for more information")
        raise IndexError(in_FOR_msg)

    if not args['--silent']:
        display_results(eph)

    if args['--write_ephemeris']:
        eph.write_ephemeris(eph.dataframe, args['--write_ephemeris'])

    plot_visibility(eph, args['--instrument'], write_plot=args['--write_plot'])

def driver():
    args = docopt(usage, version='0.6.2') 
    main(args)
