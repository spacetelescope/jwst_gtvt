"""Setup definition for the JWST General Target Visibility Tool
(Based on https://github.com/pypa/sampleproject)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='jwst-gtvt',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.3.1',

    description='JWST General Target Visibility Tool',
    
    long_description=long_description,
    
    long_description_content_type='text/markdown',

    # The project's main homepage.
    url='https://github.com/spacetelescope/jwst_gtvt',

    # Author details
    author='Space Telescope Science Institute',

    # Choose your license
    license='BSD',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

    # What does your project relate to?
    keywords='jwst astro astronomy observation',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=['jwst_gtvt'],
    scripts=['bin/jwst_gtvt', 'bin/jwst_mtvt', 'bin/delete_cache'],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
     install_requires=[
         'numpy',
         'astropy',
         'matplotlib',
         'astroquery',
         'pysiaf',
     ],
     
     python_requires='>=3.8',

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #     'dev': [
    #         'pytest>=2.9.1',
    #         'lxml>=3.5.0',
    #         'check-manifest>=0.31'
    #     ],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'jwst_gtvt': [
            'horizons_EM_jwst_wrt_sun_2021-2024.txt'
        ],
    },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # entry_points={
    #     'gui_scripts': [
    #         'jwst-coronagraph-visibility-gui=jwst_coronagraph_visibility.gui:run',
    #     ],
    # },
    zip_safe=False
)
