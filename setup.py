# from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension

from setuptools import find_packages

# try:
#     from Cython.Distutils import build_ext
# except ImportError:
#     use_cython = False
# else:
#     use_cython = True

cmdclass = {}
ext_modules = []

# if use_cython:
#     ext_modules += [
#         Extension("astrotk.simulator.integrators.RK4", ["src/astrotk/simulator/integrators/RK4.pyx"]),
#         Extension("astrotk.simulator.integrators.Euler", ["src/astrotk/simulator/integrators/Euler.pyx"]),
#         Extension("astrotk.simulator.integrators.AB4", ["src/astrotk/simulator/integrators/AB4.pyx"]),
#         Extension("astrotk.simulator.eom.test", ["src/astrotk/simulator/eom/test.pyx"]),
#     ]
#
#     cmdclass.update({'build_ext': build_ext})
#
# else:
#     ext_modules += [
#         Extension("astrotk.simulator.integrators.RK4", ["src/astrotk/simulator/integrators/RK4.c"]),
#         Extension("astrotk.simulator.integrators.Euler", ["src/astrotk/simulator/integrators/Euler.c"]),
#         Extension("astrotk.simulator.integrators.AB4", ["src/astrotk/simulator/integrators/AB4.c"]),
#         Extension("astrotk.simulator.eom.test", ["src/astrotk/simulator/eom/test.c"]),
#     ]

setup(
    name='tudatpy',
    version='0.0.2',
    packages=['tudatpy',
              'tudatpy.tests', 'tudatpy.simulation.environment',
              'tudatpy.constants'] + find_packages(),
    package_dir={'': 'tudatpy'},
    url='',
    license='MIT',
    author='Geoffrey Garrett',
    author_email='g.h.garrett13@gmail.com',
    description='',
    install_requires=[
        'numpy',
        # 'sympy',
        # 'pandas',
        # 'poliastro',
        # 'pytest',
        # 'prettytable',
        # 'Cython'
    ],
    # This is for the AE4878 constants. Will be deprecated.
    include_package_data=True,
    cmdclass=cmdclass,
    ext_modules=ext_modules

)