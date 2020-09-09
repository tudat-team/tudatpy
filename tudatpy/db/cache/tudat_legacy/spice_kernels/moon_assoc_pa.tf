KPL/FK


   SPICE Lunar PA Reference Frame/Body Association  Kernel
   =====================================================================

   Original file name:                   moon_assoc_pa.tf
   Creation date:                        2007 February 13 17:39
   Created by:                           Nat Bachman  (NAIF/JPL)
   Last updated:                         2008 March 18 22:17
   Purpose of update:
  
      Documentation now refers to DE421 kernels. Deprecated SPICE
      routines and their replacements are noted.

   Overview
   =====================================================================

   In the SPICE system, the default body-fixed reference frame
   associated with the Moon is named

       IAU_MOON

   The IAU_MOON reference frame is implemented via approximate formulas
   provided by the IAU report [1] and is not suitable for high-accuracy
   work.

   This kernel directs the SPICE system to associate the lunar
   "principal axis" reference frame

      MOON_PA

   with the Moon. 

   When this kernel is loaded via FURNSH, the SPICE frame system
   routines CNMFRM and CIDFRM, which identify the reference frame
   associated with a specified body, will indicate that the MOON_PA
   frame is associated with the Moon.  In addition, higher-level SPICE
   geometry routines that rely on the CNMFRM or CIDFRM routines will
   use the MOON_PA frame where applicable. As of the release date of
   this kernel, these SPICE routines are:

      ET2LST
      LSPCN
   
   Any code that calls these routines to obtain results
   involving lunar body-fixed frames are affected.  Within SPICE, the
   only higher-level system that is affected is the dynamic frame
   system.

   The deprecated (as of the N0062 SPICE Toolkit release) routines

      ILLUM
      SRFXPT
      SUBPT
      SUBSOL

   also make use of this kernel; however NAIF recommends that
   users instead call the following routines which, respectively,
   supersede those listed above:

      ILUMIN
      SINCPT
      SUBPNT
      SUBSLR

   The newer routines don't make use of frame association kernels;
   these routines accept the name of the target body-fixed
   frame as an input argument.

   Note:  to direct SPICE to associate the lunar mean Earth/polar 
   axis frame 

      MOON_ME

   with the Moon, load the kernel

      moon_assoc_me.tf

   rather than this one.



   Using this kernel
   =====================================================================

   This kernel must be loaded together with a lunar frame specification
   kernel and a binary lunar PCK.  Below an example meta-kernel that
   loads these files and a small program illustrating use of the
   meta-kernel are shown. The names of the kernels used here are
   current as of the release date of this kernel, but should not be
   assumed to be current at later dates.


      Example meta-kernel
      -------------------

      To use the meta-kernel shown below, the '@' characters must be
      replaced with backslash '\' characters.  Backslashes cannot be
      used in this comment block because they would confuse the SPICE
      text kernel parser.


         KPL/FK

         @begintext

             Kernels to load are:

                Lunar kernels
                -------------
                Binary lunar PCK:          moon_pa_de421_1900-2050.bpc
                Lunar FK:                  moon_080317.tf
                Frame association kernel:  moon_assoc_pa.tf

                Additional kernels to support sub-point computation
                ---------------------------------------------------
                Text PCK for lunar radii:  pck00008.tpc

                Leapseconds kernel (for
                time conversion):          naif0008.tls

                Planetary ephemeris (for
                sub-Earth computation):    de421.bsp

         @begindata

         KERNELS_TO_LOAD = ( 'moon_pa_de421_1900-2050.bpc'
                             'moon_080317.tf'
                             'moon_assoc_pa.tf'    
                             'pck00008.tpc'
                             'naif0008.tls'         
                             'de421.bsp'                   )
         @begintext

         End of kernel


      Example code
      ------------

      Find the geometric (without light time and stellar aberration
      corrections) sub-Earth point on the Moon at a given UTC time,
      using the MOON_PA reference frame. Display the name of the
      body-fixed lunar frame used for the computation.


             PROGRAM EX
             IMPLICIT NONE

             DOUBLE PRECISION      DPR

             INTEGER               FILEN 
             PARAMETER           ( FILEN  = 255 )

             INTEGER               FRNMLN
             PARAMETER           ( FRNMLN = 32 )

             INTEGER               TIMLEN
             PARAMETER           ( TIMLEN = 50 )

             CHARACTER*(FRNMLN)    FRNAME
             CHARACTER*(FILEN)     META
             CHARACTER*(TIMLEN)    TIMSTR

             DOUBLE PRECISION      ALT
             DOUBLE PRECISION      ET
             DOUBLE PRECISION      LAT
             DOUBLE PRECISION      LON
             DOUBLE PRECISION      RADIUS
             DOUBLE PRECISION      SPOINT ( 3 )

             INTEGER               FRCODE

             LOGICAL               FOUND

       C
       C     Obtain name of meta-kernel; load kernel.
       C
             CALL PROMPT ( 'Enter meta-kernel name > ', META   )
             CALL FURNSH ( META )

       C
       C     Obtain input time and convert to seconds past J2000 TDB.
       C
             CALL PROMPT ( 'Enter observation time > ', TIMSTR )
             CALL STR2ET ( TIMSTR, ET )

       C
       C     Find the closest point on the Moon to the center
       C     of the Earth at ET.
       C
             CALL SUBPT  ( 'Near point',  'MOON',  ET,  'NONE', 
            .              'EARTH',       SPOINT,  ALT          )
            .               
       C
       C     Express the sub-observer point in latitudinal
       C     coordinates.
       C
             CALL RECLAT ( SPOINT, RADIUS, LON, LAT )

       C
       C     Look up the name of the lunar body-fixed frame.
       C     
             CALL CNMFRM ( 'MOON', FRCODE, FRNAME, FOUND )

       C
       C     Always check the "found" flag.  Signal an error if we
       C     don't find a frame name.
       C
             IF ( .NOT. FOUND ) THEN
                CALL SETMSG ( 'No body-fixed frame found for the Moon.' )
                CALL SIGERR ( 'SPICE(NOFRAME)'                          )
             END IF

             WRITE(*,*) 'Lunar body-fixed frame is ', FRNAME
             WRITE(*,*) 'Sub-Earth planetocentric longitude (deg):', 
            .            LON*DPR()
             WRITE(*,*) 'Sub-Earth planetocentric latitude  (deg):',
            .            LAT*DPR()
             END


      Example program output
      ----------------------

      Numeric results and output formatting shown below should be
      expected to differ somewhat across different computing platforms.

      When the above example program is run using the example meta-kernel,
      and the (arbitrary) date 2008 Mar 18 00:00:00 UTC is used
      as the observation time, the output will be:

         Lunar body-fixed frame is MOON_PA
         Sub-Earth planetocentric longitude (deg):  5.03631326
         Sub-Earth planetocentric latitude  (deg): -1.63759776


   References
   =====================================================================
   [1]   Seidelmann, P.K., Abalakin, V.K., Bursa, M., Davies, M.E.,
         Bergh, C. de, Lieske, J.H., Oberst, J., Simon, J.L.,
         Standish, E.M., Stooke, P., and Thomas, P.C. (2002).
         "Report of the IAU/IAG Working Group on Cartographic
         Coordinates and Rotational Elements of the Planets and
         Satellites: 2000," Celestial Mechanics and Dynamical
         Astronomy, v.82, Issue 1, pp. 83-111.



   Data
   =====================================================================

   The assignment below directs the SPICE system to associate the MOON_PA 
   reference frame with the Moon.

   For further information, see the Frames Required Reading section
   titled "Connecting an Object to its Body-fixed Frame."

   \begindata

      OBJECT_MOON_FRAME =  'MOON_PA'

   \begintext


   End of kernel
   =====================================================================


