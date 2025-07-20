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


Constants
---------
.. currentmodule:: tudatpy.constants

.. autodata:: tudatpy.constants.SEA_LEVEL_GRAVITATIONAL_ACCELERATION
    
    Standard gravitational acceleration at sea-level.
    
.. autodata:: tudatpy.constants.JULIAN_DAY
    
    Julian day in seconds :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_YEAR_IN_DAYS

    Julian year in Julian days :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.JULIAN_YEAR

    Julian year in seconds. Result of JULIAN_YEAR_IN_DAYS * JULIAN_DAY.

.. autodata:: tudatpy.constants.SIDEREAL_DAY

    Sidereal day in seconds :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.SIDEREAL_YEAR_IN_DAYS

    Sidereal year in Julian days in quasar reference frame :cite:p:`NASASSD2011`.

.. autodata:: tudatpy.constants.SIDEREAL_YEAR

    Sidereal year in quasar reference frame. Result of SIDEREAL_YEAR_IN_DAYS * JULIAN_DAY.

.. autodata:: tudatpy.constants.SPEED_OF_LIGHT

    Speed of light in meters per second :cite:p:`standish1995`.

.. autodata:: tudatpy.constants.GRAVITATIONAL_CONSTANT

    Gravitational constant in :math:`\mathrm{m}^3/\mathrm{s}^2` :cite:p:`standish1995`.

.. autodata:: tudatpy.constants.ASTRONOMICAL_UNIT

    Astronomical Unit in meters :cite:p:`iau2012`.

.. autodata:: tudatpy.constants.SPECIFIC_GAS_CONSTANT_AIR

    The specific gas constant of air in :math:`\mathrm{J}/(\mathrm{kg} \, \mathrm{K})` :cite:p:`anderson2006`.

.. autodata:: tudatpy.constants.MOLAR_GAS_CONSTANT

    The molar gas constant in :math:`\mathrm{J}/(\mathrm{mol} \, \mathrm{K})` :cite:p:`nistconstants`.
    Also known as universal gas constant.

.. autodata:: tudatpy.constants.PLANCK_CONSTANT

    Planck's constant in :math:`\mathrm{m}^{2} \, \mathrm{kg}/\mathrm{s}` :cite:p:`nistconstants`.

.. autodata:: tudatpy.constants.BOLTZMANN_CONSTANT

    The Boltzmann constant (gas constant per particle) in  :math:`\mathrm{m}^{2} \, \mathrm{kg} / ( \mathrm{s}^{2} \, \mathrm{K} )`, :cite:p:`nistconstants`.

.. autodata:: tudatpy.constants.STEFAN_BOLTZMANN_CONSTANT

    Stefan-Boltzmann constant, used for calculating black body radiation intensity in :math:`\mathrm{J} / (\mathrm{s} \, \mathrm{m}^{2} \, \mathrm{K}^{4} )` :cite:p:`wolfram2013`.

.. autodata:: tudatpy.constants.INVERSE_SQUARE_SPEED_OF_LIGHT

    Precomputed inverse-square of speed of light.

.. autodata:: tudatpy.constants.INVERSE_CUBIC_SPEED_OF_LIGHT

    Precomputed inverse 3rd power of speed of light.

.. autodata:: tudatpy.constants.INVERSE_QUARTIC_SPEED_OF_LIGHT

    Precomputed inverse 4th power of speed of light.

.. autodata:: tudatpy.constants.INVERSE_QUINTIC_SPEED_OF_LIGHT

    Precomputed inverse 5th power of speed of light.

.. autodata:: tudatpy.constants.VACUUM_PERMEABILITY

    Permeability of vacuum.

.. autodata:: tudatpy.constants.VACUUM_PERMITTIVITY

    Permittivity of vacuum.

.. autodata:: tudatpy.constants.LG_TIME_RATE_TERM

    Relative time rate difference between TCG and TT time scales.

.. autodata:: tudatpy.constants.JULIAN_DAY_ON_J2000

    Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT).

.. autodata:: tudatpy.constants.JULIAN_DAY_AT_0_MJD

    Julian day at Modified Julian Date 0, i.e. Nov 17, 1858, 00:00.

.. autodata:: tudatpy.constants.E

    The constant E is base of the natural logarithm, and is also known as Napier's constant.

.. autodata:: tudatpy.constants.GOLDEN_RATIO

    The golden ratio, also known as the divine proportion, golden mean, or golden section, is a number often encountered when taking the ratios of distances in simple geometric figures such as the pentagon, pentagram, decagon and dodecahedron.

.. autodata:: tudatpy.constants.COMPLEX_I

    Independent root of -1, typically denoted :math:`i`.

.. autodata:: tudatpy.constants.PI

    The constant PI, denoted :math:`\pi`, is a real number defined as the ratio of a circle's circumference :math:`C` to its diameter, :math:`d = 2r`.

.. autodata:: tudatpy.constants.TUDAT_NAN

    Tudat Not-A-Number (NaN) value.

