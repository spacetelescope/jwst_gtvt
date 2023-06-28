"""
Utility functions for package
"""

import astropy
import os
import shutil


def delete_cache():
    """Delete astroquery Horizons cache"""

    cache = astropy.config.get_cache_dir()
    horizons_cache = os.path.join(cache, 'astroquery', 'Horizons')

    if os.path.exists(horizons_cache):
        print('Horizons cache located: {}'.format(horizons_cache))
        choice = input("Do you wish to delete this folder? [y/n]: ")

        try:
            if choice == 'y':
                shutil.rmtree(horizons_cache)
                print('Cache deleted..')
            elif choice == 'n':
                print('You chose to not delete the cache, exiting.')
            else:
                raise ValueError
        except ValueError:
            print('Input was not y or n, try again and provide valid input!')
    else:
        print('There is no Horizons cache available, exiting')
