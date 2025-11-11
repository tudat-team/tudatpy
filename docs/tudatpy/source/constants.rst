.. _constants:

``constants``
=============
This module contains constants that are used in the tudatpy library or commonly used in astrodynamic calculations.

.. References
.. ----------
.. .. [1] NASA SSD. Astrodynamic Parameters, http://ssd.jpl.nasa.gov/?constants#ref, 6th September, 2011, last accessed: 21st February, 2012.
.. .. [2] Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards", in Highlights of Astronomy (I. Appenzeller, ed.), Table 1, Kluwer Academic Publishers, Dordrecht
.. .. [3] IAU 2012, Resolution B2. https://www.iau.org/static/resolutions/IAU2012_English.pdf
.. .. [4] Anderson, J.D. Jr. Hypersonic and High-Temperature Gas Dynamics, Second Edition, p469, 2006
.. .. [5] NIST reference on constants, units and uncertainty. http://physics.nist.gov/cuu/Constants/index.html, last accessed: 11th January, 2013
.. .. [6] Wolfram Research, http://scienceworld.wolfram.com/physics/Stefan-BoltzmannLaw.html, last accessed: 11th January 2013.


Physical Constants
------------------
.. currentmodule:: tudatpy.constants

Fundamental physical constants used in physics and astronomy.

.. autodata:: tudatpy.constants.SEA_LEVEL_GRAVITATIONAL_ACCELERATION
    
    Standard gravitational acceleration at Earth's sea level in :math:`\mathrm{m}/\mathrm{s}^2`.
    
.. autodata:: tudatpy.constants.SPEED_OF_LIGHT

    Speed of light in vacuum in :math:`\mathrm{m}/\mathrm{s}` :cite:p:`standish1995`.

.. autodata:: tudatpy.constants.SPEED_OF_LIGHT_LONG

    Speed of light in vacuum, high-precision (long double) variant in :math:`\mathrm{m}/\mathrm{s}` :cite:p:`standish1995`.

.. autodata:: tudatpy.constants.GRAVITATIONAL_CONSTANT

    Newtonian gravitational constant in :math:`\mathrm{m}^3/(\mathrm{kg} \, \mathrm{s}^2)` :cite:p:`standish1995`.

.. autodata:: tudatpy.constants.ASTRONOMICAL_UNIT

    Astronomical Unit - mean Earth-Sun distance in :math:`\mathrm{m}` :cite:p:`iau2012`.

.. autodata:: tudatpy.constants.SPECIFIC_GAS_CONSTANT_AIR

    Specific gas constant for air in :math:`\mathrm{J}/(\mathrm{kg} \, \mathrm{K})` :cite:p:`anderson2006`.

.. autodata:: tudatpy.constants.MOLAR_GAS_CONSTANT

    Universal molar gas constant in :math:`\mathrm{J}/(\mathrm{mol} \, \mathrm{K})` :cite:p:`nistconstants`. Also known as universal gas constant.

.. autodata:: tudatpy.constants.PLANCK_CONSTANT

    Planck constant in :math:`\mathrm{J} \, \mathrm{s}` :cite:p:`nistconstants`.

.. autodata:: tudatpy.constants.BOLTZMANN_CONSTANT

    Boltzmann constant - gas constant per particle in :math:`\mathrm{m}^{2} \, \mathrm{kg} / ( \mathrm{s}^{2} \, \mathrm{K} )` :cite:p:`nistconstants`.

.. autodata:: tudatpy.constants.STEFAN_BOLTZMANN_CONSTANT

    Stefan-Boltzmann constant for black body radiation in :math:`\mathrm{W}/(\mathrm{m}^{2} \, \mathrm{K}^{4})` :cite:p:`wolfram2013`.

.. autodata:: tudatpy.constants.INVERSE_SQUARE_SPEED_OF_LIGHT

    Precomputed :math:`1/c^2` in :math:`\mathrm{s}^2/\mathrm{m}^2`.

.. autodata:: tudatpy.constants.INVERSE_CUBIC_SPEED_OF_LIGHT

    Precomputed :math:`1/c^3` in :math:`\mathrm{s}^3/\mathrm{m}^3`.

.. autodata:: tudatpy.constants.INVERSE_QUARTIC_SPEED_OF_LIGHT

    Precomputed :math:`1/c^4` in :math:`\mathrm{s}^4/\mathrm{m}^4`.

.. autodata:: tudatpy.constants.INVERSE_QUINTIC_SPEED_OF_LIGHT

    Precomputed :math:`1/c^5` in :math:`\mathrm{s}^5/\mathrm{m}^5`.

.. autodata:: tudatpy.constants.VACUUM_PERMEABILITY

    Magnetic permeability of free space :math:`\mu_0` in :math:`\mathrm{H}/\mathrm{m}`. Defined as :math:`\mu_0 = 4\pi \times 10^{-7}` :math:`\mathrm{H}/\mathrm{m}`.

.. autodata:: tudatpy.constants.VACUUM_PERMITTIVITY

    Electric permittivity of free space :math:`\varepsilon_0` in :math:`\mathrm{F}/\mathrm{m}`. Computed as :math:`\varepsilon_0 = 1/(\mu_0 c^2)`.


Time Constants
--------------

Time conversion factors and relativistic time scale parameters.

.. autodata:: tudatpy.constants.JULIAN_DAY
    
    Number of seconds in a Julian day in :math:`\mathrm{s}` :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_DAY_LONG

    Number of seconds in a Julian day, high-precision (long double) variant in :math:`\mathrm{s}` :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_YEAR_IN_DAYS

    Number of days in a Julian year :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_YEAR_IN_DAYS_LONG

    Number of days in a Julian year, high-precision (long double) variant :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_YEAR

    Number of seconds in a Julian year in :math:`\mathrm{s}`. Result of ``JULIAN_YEAR_IN_DAYS * JULIAN_DAY``.

.. autodata:: tudatpy.constants.SIDEREAL_DAY

    Number of seconds in a sidereal day in :math:`\mathrm{s}` :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.SIDEREAL_YEAR_IN_DAYS

    Number of days in a sidereal year in quasar reference frame :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.SIDEREAL_YEAR

    Number of seconds in a sidereal year in quasar reference frame in :math:`\mathrm{s}`. Result of ``SIDEREAL_YEAR_IN_DAYS * JULIAN_DAY``.

.. autodata:: tudatpy.constants.LG_TIME_RATE_TERM

    Relative time rate difference between Geocentric Coordinate Time (TCG) and Terrestrial Time (TT). Defines the linear drift rate: TCG progresses faster than TT by approximately 0.7 parts per billion. Over one year, TCG gains about 22 milliseconds relative to TT.

.. autodata:: tudatpy.constants.LG_TIME_RATE_TERM_LONG

    Relative time rate difference between TCG and TT, high-precision (long double) variant.


Epoch Constants
---------------

.. autodata:: tudatpy.constants.JULIAN_DAY_ON_J2000

    Julian Day Number at the J2000.0 epoch (2000-01-01 12:00:00 TT).

.. autodata:: tudatpy.constants.JULIAN_DAY_AT_0_MJD

    Julian Day Number at Modified Julian Date zero (1858-11-17 00:00:00).


Mathematical Constants
----------------------

Fundamental mathematical constants.

.. autodata:: tudatpy.constants.E

    Euler's number, base of natural logarithm. Also known as Napier's constant.

.. autodata:: tudatpy.constants.GOLDEN_RATIO

    The golden ratio, also known as the divine proportion, golden mean, or golden section. A number often encountered when taking the ratios of distances in simple geometric figures such as the pentagon, pentagram, decagon and dodecahedron.

.. autodata:: tudatpy.constants.COMPLEX_I

    Imaginary unit :math:`i = \sqrt{-1}`. Independent root of -1.

.. autodata:: tudatpy.constants.PI

    The constant :math:`\pi`, ratio of circle's circumference :math:`C` to its diameter :math:`d = 2r`.

.. autodata:: tudatpy.constants.TUDAT_NAN

    Tudat Not-a-Number (NaN) value indicating undefined/invalid results.


High-Precision Constants
------------------------

The following constants have high-precision ``_LONG`` variants using ``long double`` 
precision instead of standard ``double`` precision. These provide extended precision 
for calculations requiring high accuracy over long time spans.

**When to use high-precision variants:**

- Long-duration orbital propagations (decades to centuries)
- High-accuracy ephemeris calculations
- Relativistic corrections in precise orbit determination
- Time scale transformations requiring sub-nanosecond accuracy

**Available high-precision variants:**

- :data:`JULIAN_DAY_LONG` - For time conversions
- :data:`JULIAN_YEAR_IN_DAYS_LONG` - For long-term propagations
- :data:`SPEED_OF_LIGHT_LONG` - For relativistic corrections
- :data:`LG_TIME_RATE_TERM_LONG` - For TCG-TT transformations
