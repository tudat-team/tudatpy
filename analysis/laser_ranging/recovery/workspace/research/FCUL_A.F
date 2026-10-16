      DOUBLE PRECISION FUNCTION FCUL_A (LATITUDE,HEIGHT_M,T_K,ELEV_DEG)
*+
*  - - - - - - - - -
*   F C U L _ A
*  - - - - - - - - -
*
*  This function is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function computes the global total FCULa mapping function (Mendes et al. 2002).
*  It is dependent on latitude, height, and surface temperature.
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
*     LATITUDE       d      Latitude given in degrees (North Latitude)
*     HEIGHT_M       d      Height given in meters (mean sea level)
*     T_K            d      Surface temperature given in Kelvin
*     ELEV_DEG       d      Elevation angle given in degrees (See references)

*  Returned:
*     FCUL_A         d      Mapping function to scale total delay (Note 1)
*
*  Notes:
*
*  1) These coefficients are based on a LS adjustment of 87766 (cleaned)
*     set of traces, based on Ciddor routines to compute refractivity,
*     according to IUGG recommendations (1999).
*
*  Test case:
*     given input: LATITUDE = 30.67166667D0 degrees (McDonald Observatory)
*                  HEIGHT_M = 2075D0 meters (mean sea level)
*                  T_K      = 300.15D0 Kelvin (August 12, 2009)
*                  ELEV_DEG = 15D0 degrees (See Mendes et al.)
*     expected output: FCUL_A = 3.800243667312344087D0
*                      
*  References:
*     Mendes, V.B., G. Prates, E.C. Pavlis, D.E. Pavlis, 
*     and R.B. Langley (2002). "Improved mapping functions for
*     atmospheric refraction correction in SLR", Geophysical 
*     Res. Lett., 29(10), 1414, doi:10.1029/2001GL014394, 2002
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2000 March  31 V.B. Mendes   Original code
*  2009 August 12 B.E. Stetzler Added header and copyright
*  2009 August 13 B.E. Stetzler Use of DOUBLE PRECISION
*  2009 August 13 B.E. Stetzler Provided test case
*  2009 August 13 B.E. Stetzler Capitalized all variables for FORTRAN 77
*                               compatibility and provided more comments
*----------------------------------------------------------------------- 

      IMPLICIT NONE
      DOUBLE PRECISION ELEV_DEG
      DOUBLE PRECISION EPSILON
      DOUBLE PRECISION PI, SINE
      DOUBLE PRECISION A1, A2, A3
      DOUBLE PRECISION MAP_ZEN
      DOUBLE PRECISION HEIGHT_M, LATITUDE, T_K, T_C
      DOUBLE PRECISION COSPHI
      DOUBLE PRECISION A10, A11, A12, A13
      DOUBLE PRECISION A20, A21, A22, A23
      DOUBLE PRECISION A30, A31, A32, A33
 
      PARAMETER ( PI = 3.1415926535897932384626433D0 )

*          Preliminary conversions           *

* Convert elevation angle to radians
      EPSILON = ELEV_DEG * (PI/180D0)
      SINE    = DSIN(EPSILON)
* Convert temperature to Celsius
      T_C     = T_K - 273.15D0
      COSPHI  = DCOS (LATITUDE*(PI/180D0))

* Define coefficients used in the model        

      A10 =  0.121008D-02
      A11 =  0.17295D-05
      A12 =  0.3191D-04
      A13 = -0.18478D-07
*
      A20 =  0.304965D-02
      A21 =  0.2346D-05
      A22 = -0.1035D-03
      A23 = -0.1856D-07

*
      A30 =  0.68777D-01
      A31 =  0.1972D-04
      A32 = -0.3458D-02
      A33 =  0.1060D-06

*     a, b, and c in Marini continued fraction (Eq. 5)
      A1 = A10+A11*T_C+A12*COSPHI+A13*HEIGHT_M 
      A2 = A20+A21*T_C+A22*COSPHI+A23*HEIGHT_M 
      A3 = A30+A31*T_C+A32*COSPHI+A33*HEIGHT_M 

*     numerator in continued fraction
      MAP_ZEN   = (1.0D0 + A1/(1.0D0 + A2/(1.0D0+A3)))
*

* Results
      FCUL_A = MAP_ZEN/(SINE+A1/(SINE+A2/(SINE+A3)))
 
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
      END
 
