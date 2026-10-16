      SUBROUTINE FCULZD_HPA (LATITUDE, ELLIP_HT, PRESSURE, WVP,
     .                  LAMBDA_UM, FCUL_ZTD, FCUL_ZHD, FCUL_ZWD)
*+
*  - - - - - - - - - - - -
*   F C U L Z D _ H P A 
*  - - - - - - - - - - - - 
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine determines the total zenith delay following (Mendes and Pavlis, 2004). 
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1 model	
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     LATITUDE       d      Geodetic Latitude given in degrees (North Latitude)
*     ELLIP_HT       d      Height above ellipsoid given in meters 
*     PRESSURE       d      Surface pressure given in hPa (mbars) (Note 1)
*     WVP            d      Water vapor pressure in hPa (mbars) (Note 1)
*     LAMBDA_UM      d      Laser wavelength (micrometers)
*
*  Returned:
*     FCUL_ZTD       d      Zenith total delay in meters
*     FCUL_ZHD       d      Zenith hydrostatic delay in meters
*     FCUL_ZWD       d      Zenith non-hydrostatic delay in meters
*
*  Notes:
*  
*  1) The surface pressure provided was converted from inches Hg. 
*     The water vapor pressure was calculated from the surface 
*     temperature (Celsius) and Relative Humidity (% R.H.) at the station.
*
*  Test case:
*     given input: LATITUDE  = 30.67166667D0 degrees (McDonald Observatory)
*                  ELLIP_HT  = 2010.344D0 meters 
*                  PRESSURE  = 798.4188D0 hPa (August 14, 2009)
*                  WVP       = 14.322D0 hPa (August 14, 2009)
*                  LAMBDA_UM = 0.532D0 micrometers (See Mendes et al.)
*     expected output: FCUL_ZTD = 1.935225924846803114D0 m
*                      FCUL_ZHD = 1.932992176591644462D0 m
*                      FCUL_ZWD = 0.2233748255158703871D-02 m
*                      
*  References:
*     Mendes, V.B. and E.C. Pavlis, 2004, 
*     "High-accuracy zenith delay prediction at optical wavelengths,"
*     Geophysical Res. Lett., 31, L14602, doi:10.1029/2004GL020308, 2004
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2000 March  31 V.B. Mendes   Original code
*  2009 August 14 B.E. Stetzler Added header and copyright
*  2009 August 14 B.E. Stetzler Use of DOUBLE PRECISION
*  2009 August 14 B.E. Stetzler Provided test case
*  2009 August 14 B.E. Stetzler Capitalized all variables for FORTRAN 77
*                               compatibility and provided more comments
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION C, PI, ELLIP_HT, LATITUDE, F, PRESSURE, WVP
      DOUBLE PRECISION K0, K1, K2, K3
      DOUBLE PRECISION XC, CORR
      DOUBLE PRECISION SIGMA, W0, W1, W2, W3, FH, FNH
      DOUBLE PRECISION LAMBDA_UM
      DOUBLE PRECISION FCUL_ZHD, FCUL_ZWD, FCUL_ZTD

      PARAMETER ( C = 2.99792458D8)  ! speed of light in vacuum (m/s)
      PARAMETER ( PI = 3.1415926535897932384626433D0 )

*  CO2 content in ppm
      XC = 375.0D0
*         constant values to be used in Equation (20)
*         k1 and k3 are k1* and k3* 
      K0 = 238.0185D0
      K1 = 19990.975D0   
      K2 = 57.362D0
      K3 = 579.55174D0

*         constant values to be used in Equation (32)
      W0 = 295.235D0
      W1 = 2.6422D0
      W2 = -0.032380D0
      W3 = 0.004028D0
         
*  Wave number
      SIGMA = 1/LAMBDA_UM   

*     correction factor - Equation (24)
      F = 1 - 0.00266*COS(2*PI/180*LATITUDE) - 0.00028D-3*ELLIP_HT

*     correction for CO2 content
      CORR = 1.0D0 + 0.534D-6*(XC-450)

*     dispersion equation for the hydrostatic component - Equation (20)
      FH = 0.01D0*CORR*((K1*(K0+SIGMA**2))/((K0-SIGMA**2)**2) +
     .     K3*(K2+SIGMA**2)/((K2-SIGMA**2)**2))

*     computation of the hydrostatic component - Equation (26)
*     caution: pressure in hectoPascal units
      FCUL_ZHD = 2.416579D-3*FH*PRESSURE/F

*     dispersion equation for the non-hydrostatic component - Equation (32)
      FNH = 0.003101D0*(W0+3.0*W1*SIGMA**2 +
     .      5.0*W2*SIGMA**4+7.0*W3*SIGMA**6)

*     computation of the non-hydrostatic component - Equation (38)
*     caution: pressure in hectoPascal units
      FCUL_ZWD = 1.D-4*(5.316D0*FNH-3.759*FH)*WVP/F

*      compute the zenith total delay
      FCUL_ZTD = FCUL_ZHD + FCUL_ZWD

      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      end

