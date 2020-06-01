import os
from setuptools import setup
from setuptools.dist import Distribution
import sys

NAME = 'tudatpy'
VERSION = ''
DESCRIPTION = 'TU Delft Astrodynamics Toolbox in Python'
LONG_DESCRIPTION = 'A platform to accelerate Tudat implementation, interface with machine learning frameworks and provide a tool for contemporary education in astrodynamics.'
URL = 'https://github.com/ggarrett13/tudatpy'
AUTHOR = 'The tudatpy development team'
AUTHOR_EMAIL = 'g.h.garrett13@gmail.com'
LICENSE = 'BSD 3-Clause License'
CLASSIFIERS = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
'Development Status :: 3 - Alpha',

    'Operating System :: OS Independent',

    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'License :: BSD 3-Clause License',
    'Programming Language :: Python :: 3'
]
KEYWORDS = 'science math physics optimization ai evolutionary-computing parallel-computing metaheuristics'
INSTALL_REQUIRES = ['numpy', 'pandas', 'jupyterlab', 'tabulate']
PLATFORMS = ['Unix', 'Windows', 'OSX']


class BinaryDistribution(Distribution):

    def has_ext_modules(foo):
        return True


# Setup the list of external dlls.
if os.name == 'nt':
    mingw_wheel_libs = 'mingw_wheel_libs_python{}{}.txt'.format(
        sys.version_info[0], sys.version_info[1])
    l = open(mingw_wheel_libs, 'r').readlines()
    DLL_LIST = [os.path.basename(_[:-1]) for _ in l]

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url=URL,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      classifiers=CLASSIFIERS,
      keywords=KEYWORDS,
      platforms=PLATFORMS,
      install_requires=INSTALL_REQUIRES,
      packages=['tudatpy'],
      # Include pre-compiled extension
      package_data={'tudatpy': ['core.pyd'] + \
                             DLL_LIST if os.name == 'nt' else ['core.so']},
      distclass=BinaryDistribution)
