      SUBROUTINE ETUTC (YEAR, DELTA)
*+
*  - - - - - - - - -
*   E T U T C
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to compute the difference, delta,
*  between Epheremis Time (ET) and Coordinated Universal Time (UTC).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model	
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
*     year           d      decimal year (Note 1)
*
*  Returned:
*     delta          d      ET - UTC (Note 2)
*
*     :------------------------------------------:
*     :                                          :
*     :                 IMPORTANT                :
*     :                                          :
*     :  A new version of this routine must be   :
*     :  produced whenever a new leap second is  :
*     :  announced.  There are three items to    :
*     :  change on each such occasion:           :
*     :                                          :
*     :  1) Update the nstep variable            :
*     :  2) Update the arrays st and si          :                              
*     :  3) Change date of latest leap second    :
*     :     in 'Latest leap second' and          :
*     :     'Notes: 1)'                          :
*     :                                          :
*     :  Latest leap second:  2016 December 31   :
*     :                                          :
*     :__________________________________________:
*
*  Notes:
*
*  1) This subroutine is valid only from 1700.-until next leap second.
*     Currently, this is up to 2017.0.
* 
*  2) The expression used in given in seconds.
* 
*  3) Leap second table in GAMIT UTC (and UT) is the time most 
*     often used (e.g. in time signals)
*
*  Test case (to within machine precision):
*     given input: year = 2007.0 
*
*     expected output: delta = 65.183999999999997 seconds
*
*     given input: year = 2013.0 
*
*     expected output: delta = 67.183999999999997 seconds
*
*     given input: year = 2016.0 
*
*     expected output: delta = 68.183999999999997 seconds
*
*  References:
*
*     Broucke, R. A., Explanatory Supplement American Ephemeris &
*     Nautical Almanac (cf Comm CACM, 11, 657 (1968) and 15, 918 (1972)),
*     p. 90 (Tables)
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 April 22  B.E. Stetzler   Added header and copyright and
*                                 provided test case
*  2009 August 19 B.E. Stetzler   Capitalized all variables for FORTRAN
*                                 77 compatibility
*  2012 March 13  B.E. Stetzler   Updated for the 30 June 2012 leap
*                                 second
*  2015 April 29  M.A. Davis      Updated for the 30 June 2015 leap
*                                 second
*  2015 May 19    M.A. Davis      Changed DELTA and YEAR from REAL to 
*                                 DOUBLE PRECISION to correspond to
*                                 TDFRPH.F
*  2016 Dec 19    M.A. Davis      Updated for the 31 Dec 2016 leap
*                                 second
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER I,N,NSTEP
      PARAMETER (NSTEP=27)
      REAL D,FRAC,ST,SI,TX,TY
      DOUBLE PRECISION DELTA,YEAR
      DIMENSION D(142),TX(39),TY(39),ST(NSTEP),SI(NSTEP)

* si gives amount of step, at the times given in st

      DATA SI/27*1./
      DATA ST/1972.5,1973.0,1974.0,1975.0,1976.0,1977.0,1978.0,
     .        1979.0,1980.0,1981.5,1982.5,1983.5,1985.5,1988.0,
     .        1990.0,1991.0,1992.5,1993.5,1994.5,1996.0,1997.5,
     .        1999.0,2006.0,2009.0,2012.5,2015.5,2017.0/

      DATA D/ 5.15, 4.64, 5.36, 3.49, 3.27, 2.45, 4.03, 1.76, 3.30,
     .  1.00, 2.42, 0.94, 2.31, 2.27,-0.22, 0.03,-0.05,-0.06,-0.57,
     .  0.03,-0.47, 0.98,-0.86, 2.45, 0.22, 0.37, 2.79, 1.20, 3.52,
     .  1.17, 2.67, 3.06, 2.66, 2.97, 3.28, 3.31, 3.33, 3.23, 3.60,
     .  3.52, 4.27, 2.68, 2.75, 2.67, 1.94, 1.39, 1.66, 0.88, 0.33,
     . -0.17,-1.88,-3.43,-4.05,-5.77,-7.06,-7.36,-7.67,-7.64,-7.93,
     . -7.82,-8.35,-7.91,-8.03,-9.14,-8.18,-7.88,-7.62,-7.17,-8.14,
     . -7.59,-7.17,-7.94,-8.23,-7.88,-7.68,-6.94,-6.89,-7.11,-5.87,
     . -5.04,-3.90,-2.87,-0.58, 0.71, 1.80, 3.08, 4.63, 5.86, 7.21,
     .  8.58,10.50,12.10,12.49,14.41,15.59,15.81,17.52,19.01,18.39,
     . 19.55,20.36,21.01,21.81,21.76,22.35,22.68,22.94,22.93,22.69,
     . 22.94,23.20,23.31,23.63,23.47,23.68,23.62,23.53,23.59,23.99,
     . 23.80,24.20,24.99,24.97,25.72,26.21,26.37,26.89,27.68,28.13,
     . 28.94,29.42,29.66,30.29,30.96,31.09,31.59,32.06,31.82,32.69,
     . 33.05,33.16,33.59/
      DATA TX/61.5,
     .62.     ,62.5     ,63.      ,63.5     ,64.      ,64.5     ,65.   ,
     .65.5    ,66.      ,66.5     ,67.      ,67.5     ,68.      ,68.25 ,
     .68.5    ,68.75    ,69.      ,69.25    ,69.5     ,69.75    ,70.   ,
     .70.25   ,70.5     ,70.75    ,71.      ,71.085   ,71.162   ,71.247,
     .71.329  ,71.414   ,71.496   ,71.581   ,71.666   ,71.748   ,71.833,
     .71.915  ,71.999   ,72.0/
      DATA TY/33.59,
     .34.032  ,34.235   ,34.441   ,34.644   ,34.95    ,35.286   ,35.725,
     .36.16   ,36.498   ,36.968   ,37.444   ,37.913   ,38.39    ,38.526,
     .38.76   ,39.      ,39.238   ,39.472   ,39.707   ,39.946   ,40.185,
     .40.42   ,40.654   ,40.892   ,41.131   ,41.211   ,41.284   ,41.364,
     .41.442  ,41.522   ,41.600   ,41.680   ,41.761   ,41.838   ,41.919,
     .41.996  ,42.184   ,42.184/

*  For the oldest epochs, use approximations

      IF(YEAR.LT.1700.0D0) THEN
        DELTA = 0.0D0
        RETURN
      ENDIF
      IF(YEAR.LT.1785.0D0) THEN
        DELTA = (YEAR-1750.0D0)/5.0D0
        RETURN
      ENDIF
      IF(YEAR.LT.1820.5D0) THEN
        DELTA = 6.0D0
        RETURN
      ENDIF

*  For 1820.5 to 1961.5, data is spaced at yearly intervals

      IF(YEAR.LT.1961.5D0) THEN
         N = YEAR - 1819.5
         FRAC = YEAR - (1819.5 + N)
         DELTA = (D(N+1) - D(N))*FRAC + D(N)
         RETURN
      ENDIF

*  For 1961.5 to 1972.0, interpolate between unequispaced data

      IF(YEAR.LT.1972.0D0) THEN
        DO 150 I = 1,38
           IF(YEAR-1900.0D0.EQ.TX(I)) THEN
              DELTA = TY(I)
              RETURN
           ENDIF
           IF(YEAR-1900.0D0.LT.TX(I)) THEN
              DELTA=TY(I-1) + (TY(I)-TY(I-1))*
     .                ((YEAR-1900.0D0-TX(I-1))/(TX(I)-TX(I-1)))
              RETURN
           ENDIF
150     CONTINUE
      ENDIF

*--------------------------------------------------------------------------*
*   after 1972 et-utc has only step offsets. st is the array of step times,*
*   and si is the step sizes (an added second is +1.)                      *
*--------------------------------------------------------------------------*
      DELTA = 42.184D0
      DO 250 I = 1,NSTEP
         IF(YEAR.GE.ST(I)) DELTA = DELTA + SI(I)
         IF(YEAR.LT.ST(I)) RETURN
250   CONTINUE
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
