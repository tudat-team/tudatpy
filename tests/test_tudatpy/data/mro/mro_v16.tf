KPL/FK

\beginlabel
PDS_VERSION_ID               = PDS3
RECORD_TYPE                  = STREAM
RECORD_BYTES                 = "N/A"
^SPICE_KERNEL                = "mro_v16.tf"
MISSION_NAME                 = "MARS RECONNAISSANCE ORBITER"
SPACECRAFT_NAME              = "MARS RECONNAISSANCE ORBITER"
DATA_SET_ID                  = "MRO-M-SPICE-6-V1.0"
KERNEL_TYPE_ID               = FK
PRODUCT_ID                   = "mro_v16.tf"
PRODUCT_CREATION_TIME        = 2020-11-24T09:39:52
PRODUCER_ID                  = "NAIF/JPL"
MISSION_PHASE_NAME           = "N/A"
PRODUCT_VERSION_TYPE         = ACTUAL
PLATFORM_OR_MOUNTING_NAME    = "N/A"
START_TIME                   = "N/A"
STOP_TIME                    = "N/A"
SPACECRAFT_CLOCK_START_COUNT = "N/A"
SPACECRAFT_CLOCK_STOP_COUNT  = "N/A"
TARGET_NAME                  = MARS
INSTRUMENT_NAME              = "N/A"
NAIF_INSTRUMENT_ID           = "N/A"
SOURCE_PRODUCT_ID            = "N/A"
NOTE                         = "See comments in the file for details"
OBJECT                       = SPICE_KERNEL
  INTERCHANGE_FORMAT         = ASCII
  KERNEL_TYPE                = FRAMES
  DESCRIPTION                = "MRO SPICE FK file providing the complete
set of frame definitions for the MRO spacecraft, its structures and science
instruments, created by NAIF, JPL. "
END_OBJECT                   = SPICE_KERNEL
\endlabel


Mars Reconnaissance Orbiter Frames Kernel
===============================================================================

   This frame kernel contains complete set of frame definitions for the
   Mars Reconnaissance Orbiter (MRO) spacecraft, its structures and
   science instruments. This frame kernel also contains name - to -
   NAIF ID mappings for MRO science instruments and s/c structures (see
   the last section of the file.)


Version and Date
-------------------------------------------------------------------------------

   Version 1.6 -- November 11, 2020 -- Boris Semenov, NAIF

      Updated MRO_HGA_BASEPLATE and MRO_HGA frame alignments based on
      [22].

   Version 1.5 -- October 22, 2012 -- Boris Semenov, NAIF

      Updated orientations of the MRO_MCS_BASE and
      MRO_MCS_EL_GIMBAL_REF frames to incorporate alignment data
      derived from off-track observations in 2012.

   Version 1.4 -- February 18, 2009 -- Laszlo Keszthelyi, USGS;
                                       Boris Semenov, NAIF

      Adjusted the orientation of the MRO_HIRISE_OPTICAL_AXIS frame to
      compensate for the correction of the optical distortion model
      made in 2008.

      Adjusted the orientation of the MRO_HIRISE_LOOK_DIRECTION frame
      to compensate for the change in the orientation of the
      MRO_HIRISE_OPTICAL_AXIS frame (to keep MRO_HIRISE_LOOK_DIRECTION
      oriented with respect to the spacecraft the same way as it was in
      the versions 0.7-1.3 of the FK).

   Version 1.3 -- January 23, 2009 -- Boris Semenov, NAIF

      Updated orientations of the MRO_MCS_BASE and
      MRO_MCS_EL_GIMBAL_REF frames to incorporate alignment data
      derived from off-track observations in 2008.

   Version 1.2 -- April 24, 2008 -- Boris Semenov, NAIF

      Redefined the MRO_MCS_BASE frame to be with respect to the
      MRO_SPACECRAFT frame with zero offset rotation. "Moved" MCS
      azimuth offset (-0.46 deg about Z) from the MRO_MCS_EL_GIMBAL_REF
      definition to the MRO_MCS_AZ_GIMBAL_REF definition. Added
      MRO_MCS_SOLAR_TARGET frame (ID -74506).

   Version 1.1 -- September 08, 2007 -- Boris Semenov, NAIF

      Re-defined MARCI frames to follow convention used by MSSS.
      Incorporated MARCI alignment provided by MSSS. Changed names for
      MARCI bands in the naif-ID definitions from VIS_BAND1, ...
      UV_BAND2 to VIS_BLUE, ... UV_LONG_UV.

   Version 1.0 -- May 31, 2007 -- Boris Semenov, NAIF

      Changed MCS frame layout based on [17] and incorporated MCS
      misalignment data used by the MCS team, specifically:

         -  renamed MRO_MCS_AZ_GIMBAL_0 to MRO_MCS_AZ_GIMBAL_REF

         -  renamed MRO_MCS_EL_GIMBAL_0 to MRO_MCS_EL_GIMBAL_REF

         -  redefined MRO_MCS frame to be nominally co-aligned with the
            s/c frame in the forward-looking position (az=180,el=90)

         -  incorporated MCS misalignment derived by the MCS team from
            early post MOI observations and used in processing during
            first year of PSP (in the definition of the
            MRO_MCS_EL_GIMBAL_REF frame)
            
      Incorporated final ONC alignment calibrated in flight provided 
      in [15]

      Incorporated CRISM frame layout and misalignment data used by the
      CRISM team ([18]), specifically:

         -  renamed MRO_CRISM_OSU to MRO_CRISM_ART

         -  changed the frame chain diagram and table to show that
            MRO_CRISM_IR is fixed relative to MRO_CRISM_VNIR

         -  replaced previous CRISM frame definition section of the FK
            with sections ``Version and Date'', ``References'',
            ``Contact Information'' and ``CRISM Frame Definitions''
            from [18]

         -  changed MRO_CRISM_VNIR instrument ID from -74012 to -74017
            in the name/ID mapping keywords

         -  changed MRO_CRISM_IR instrument ID from -74013 to -74018
            in the name/ID mapping keywords

         -  deleted MRO_CRISM_OSU/-74011 name/ID mapping keywords

   Version 0.9 -- February 27, 2007 -- Boris Semenov, NAIF

      Fixed comments in the ```MRO Frames'' section.

   Version 0.8 -- April 17, 2006 -- Boris Semenov, NAIF

      Incorporated ONC alignment calibrated in flight (provided by Nick
      Mastrodemos on November 7, 2005.)

      Added a note stating that MRO_MME_OF_DATE frame is the same as
      MME-D frame defined in [16] to the MRO_MME_OF_DATE description
      block.

      Corrected typo in the MRO_MARCI_VIS frame definition (the keyword
      TKFRAME_-74410_AXES was TKFRAME_-74s041_AXES.)

   Version 0.7 -- September 22, 2005 -- Boris Semenov, NAIF

      The following changes were made to make frames defined for HIRISE
      consistent with the terminology used and calibration approach
      proposed by the HIRISE team:

         -  MRO_HIRISE_IF frame was renamed to MRO_HIRISE_OPTICAL_AXIS

         -  MRO_HIRISE frame was renamed to MRO_HIRISE_LOOK_DIRECTION

         -  preliminary HIRISE in-flight calibrated alignment with
            respect to the spacecraft frame was incorporated into the
            MRO_HIRISE_OPTICAL_AXIS frame definition

         -  MRO_HIRISE_LOOK_DIRECTION was redefined to be with respect
            to the MRO_HIRISE_OPTICAL_AXIS frame; the rotation
            incorporated into this definition was computed by combining
            preliminary HIRISE in-flight calibrated alignment for
            MRO_HIRISE_OPTICAL_AXIS frame with the pre-launch alignment
            for the HIRISE boresight.

   Version 0.6 -- August 25, 2005 -- Boris Semenov, NAIF

      Incorporated ground alignment data for HiRISE, CRISM, CTX and
      MSC.

   Version 0.5 -- August 8, 2005 -- Boris Semenov, NAIF

      Added MRO_MME_2000 frame.

   Version 0.4 -- June 2, 2005 -- Boris Semenov, NAIF

      Replaced MRO_MME_OF_DATE placeholder definition with a
      dynamically defined MME of date frame.

   Version 0.3 -- May 16, 2005 -- Boris Semenov, NAIF

      Corrected MRO_HGA frame to align its +X axis with the HGA pattern
      clock angle reference line.

   Version 0.2 -- February 16, 2005 -- Boris Semenov, NAIF

      Changed body ID of the MRO_HIRISE and frame ID of the MRO_HIRISE
      frame from -74600 to -74699. Added MRO_HIRISE_CCD0/-74600 name/ID
      mapping pair. Removed MRO_HIRISE_CCD14/-74614 name/ID mapping pair.

   Version 0.1 -- August 23, 2004 -- Boris Semenov, NAIF

      Added MRO_MME_OF_DATE frame (currently mapped to MARSIAU; when
      SPICE parameterized frames capability is released, this frame
      will be redefined as a dynamic frame.)

   Version 0.0 -- March 15, 2004 -- Boris Semenov, NAIF

      Initial Release.


References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. ``C-Kernel Required Reading''

   4. `Mars Reconnaissance Orbiter. GN&C Hardware Coordinate 
      Frame Definitions and Transformations'', Rev. 3, 11/30/99
  
   5. ``CRISM MICD'', Final Update, Oct 7, 2003

   6. ``CTX ICD'', Final Update, July 8, 2003

   7. ``HIRISE ICD'', Final Update, Oct 17, 2003

   8. ``MARCI ICD'', Final Update, Oct 13, 2003

   9. ``MCS ICD'', Final Update, Oct 13, 2003

  10. ``ONC ICD'', Post-PDR Update, Sep 21, 2002

  11. ``SHARAD ICD'', Final Update, Oct 24, 2003

  12. Misc. PDR/CDR presentations, 2002/2003

  13. E-mail from R. Tung, MRO Telecom, regarding HGA clock reference
      line. May 16, 2005.

  14. Ground Alignment Spreadsheet, ``mro-final-alignment_REV-G.xls''.

  15. ONC-ACS alignment, e-mail from Nick Mastrodemos, ONC Team, 
      May 10, 2007.

  16. "MRO GN&C Hardware Coordinate Frame Definitions and 
      Transformations (LIB-8)", 09/22/04

  17. Suggested changes to the MCS frame layout, e-mail from Steven 
      Gaiser, MCS Team, Aug 4, 2006 

  18. CRISM FK file "MRO_CRISM_FK_0000_000_N_1.TF", created by the CRISM
      Team, 09/14/06.

  19. E-mail from Joe Fahle, MSSS, regarding the MARCI frame definition,
      September 5, 2007.

  20. Analysis of MCS alignment based on 2008 off-track observations,
      by ? (MCS Team, JPL), Fall 2008.

  21. Updated MCS alignment based on 2012 off-track observations,
      e-mails from Dr. John T. Schofield, 07/24/12 & 10/22/12.

  22. E-mail from Paul Salame, LMCO, regarding the HGA frame definitions,
      October 28, 2020.


Contact Information
-------------------------------------------------------------------------------

   Boris V. Semenov, NAIF/JPL, (818)-354-8136, Boris.Semenov@jpl.nasa.gov


Implementation Notes
-------------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must ``load'' the kernel, normally during
   program initialization. The SPICELIB routine FURNSH/furnsh_c loads a
   kernel file into the pool as shown below.

      CALL FURNSH ( 'frame_kernel_name; )
      furnsh_c ( "frame_kernel_name" );

   This file was created and may be updated with a text editor or word
   processor.


MRO Frames
-------------------------------------------------------------------------------

   The following MRO frames are defined in this kernel file:

           Name                  Relative to           Type       NAIF ID
      ======================  =====================  ============   =======

   Non Built-in Mars Frames:
   -------------------------
      MRO_MME_OF_DATE         rel.to J2000           DYNAMIC        -74900
      MRO_MME_2000            rel.to J2000           FIXED          -74901

   Spacecraft frame:
   -----------------
      MRO_SPACECRAFT          rel.to MME_OF_DATE     CK             -74000

   Science Instrument frames:
   --------------------------
      MRO_CRISM_BASE          rel.to SPACECRAFT      FIXED          -74011
      MRO_CRISM_ART           rel.to CRISM_BASE      CK             -74012
      MRO_CRISM_VNIR          rel.to MRO_CRISM_ART   FIXED          -74017
      MRO_CRISM_IR            rel.to MRO_CRISM_VNIR  FIXED          -74018

      MRO_CTX_BASE            rel.to SPACECRAFT      FIXED          -74020
      MRO_CTX                 rel.to CTX_BASE        FIXED          -74021

      MRO_HIRISE_LOOK_DIRECTION  rel.to SPACECRAFT   FIXED          -74699
      MRO_HIRISE_OPTICAL_AXIS    rel.to SPACECRAFT   FIXED          -74690

      MRO_MARCI_BASE          rel.to SPACECRAFT      FIXED          -74400
      MRO_MARCI_VIS           rel.to MARCI_BASE      FIXED          -74410
      MRO_MARCI_UV            rel.to MARCI_BASE      FIXED          -74420

      MRO_MCS_BASE            rel.to SPACECRAFT      FIXED          -74501
      MRO_MCS_AZ_GIMBAL_REF   rel.to MCS_BASE        FIXED          -74502
      MRO_MCS_AZ_GIMBAL       rel.to MCS_AZ_GIMBAL_REF  CK          -74503
      MRO_MCS_EL_GIMBAL_REF   rel.to MCS_AZ_GIMBAL   FIXED          -74504
      MRO_MCS_EL_GIMBAL       rel.to MCS_EL_GIMBAL_REF  CK          -74505
      MRO_MCS                 rel.to MCS_EL_GIMBAL   FIXED          -74500
      MRO_MCS_SOLAR_TARGET    rel.to MCS_AZ_GIMBAL   FIXED          -74506

      MRO_ONC                 rel.to SPACECRAFT      FIXED          -74030

      MRO_SHARAD              rel.to SPACECRAFT      FIXED          -74070

   Antenna frames:
   ---------------
      MRO_HGA_BASEPLATE       rel.to SPACECRAFT      FIXED          -74211
      MRO_HGA_INNER_GIMBAL    rel.to HGA_BASEPLATE   CK             -74212
      MRO_HGA_OUTER_GIMBAL    rel.to HGA_INNER_GIM   CK             -74213
      MRO_HGA                 rel.to HGA_OUTER_GIM   FIXED          -74214
      MRO_LGA1                rel.to HGA             FIXED          -74220
      MRO_LGA2                rel.to HGA             FIXED          -74230
      MRO_UHF                 rel.to SPACECRAFT      FIXED          -74240

   Solar Array frames:
   -------------------
      MRO_SAPX_BASEPLATE      rel.to SPACECRAFT      FIXED          -74311
      MRO_SAPX_INNER_GIMBAL   rel.to SAPX_BASEPLATE  CK             -74312
      MRO_SAPX_OUTER_GIMBAL   rel.to SAPX_INNER_GIM  CK             -74313
      MRO_SAPX                rel.to SAPX_OUTER_GIM  FIXED          -74314
      MRO_SAMX_BASEPLATE      rel.to SPACECRAFT      FIXED          -74321
      MRO_SAMX_INNER_GIMBAL   rel.to SAMX_BASEPLATE  CK             -74322
      MRO_SAMX_OUTER_GIMBAL   rel.to SAMX_INNER_GIM  CK             -74323
      MRO_SAMX                rel.to SAMX_OUTER_GIM  FIXED          -74324


MRO Frames Hierarchy
-------------------------------------------------------------------------------

   The diagram below shows MRO frames hierarchy:


                               "J2000" INERTIAL
        +------------------------------------------------------------+
        |               |              |                             |
        | <--pck        |<--fixed      |<-dynamic                    | <--pck
        |               |              |                             |
        V               |              |                             V
    "IAU_MARS"          V              V                       "IAU_EARTH"
    MARS BFR(*)  "MRO_MME_2000"  "MRO_MME_OF_DATE"             EARTH BFR(*)
    -----------  --------------  -----------------             ------------
                                       |
                                       |
                                       |   "MRO_LGA1"   "MRO_LGA2"
                                       |   ----------   ----------
                                       |     ^                 ^
                                       |     |                 |
                                       |     | <--fixed        | <--fixed
                                       |     |                 |
                                       |     |                 |
                   "MRO_SA*X"          |     |    "MRO_HGA"    |
                   ----------          |     +-----------------+
                        ^              |              ^
                        |              |              |
              fixed-->  |              |              | <--fixed
                        |              |              |
            "MRO_SA*X_OUTER_GIMBAL"    |    "MRO_HGA_OUTER_GIMBAL"
            -----------------------    |    ----------------------
                        ^              |              ^
                        |              |              |
                 ck-->  |              |              | <--ck
                        |              |              |
            "MRO_SA*X_INNER_GIMBAL"    |    "MRO_HGA_INNER_GIMBAL"
            -----------------------    |    ----------------------
                        ^              |              ^
                        |              |              |
                 ck-->  |              |              | <--ck
                        |              |              |
               "MRO_SA*X_BASEPLATE"    |     "MRO_HGA_BASEPLATE"  "MRO_UHF"
               --------------------    |     -------------------  ---------
                        ^              |              ^              ^
                        |              |              |              |
              fixed-->  |              |<--ck         | <--fixed     | <--fdx
                        |              |              |              |
                        |      "MRO_SPACECRAFT"       |              |
         +-----------------------------------------------------------+
         |         |         |         |         |         |         |
         |         |         |         |         |         |         | <--fxd
         |         |         |         |         |         |         |
         |         |         |         |         |         |         V
         |         |         |         |         |         |   "MRO_SHARAD" 
         |         |         |         |         |         |   ------------
         |         |         |         |         |         |
         |         |         |         |         |         | <--fixed
         |         |         |         |         |         |
         |         |         |         |         |         V
         |         |         |         |         |     "MRO_ONC" 
         |         |         |         |         |     ---------
         |         |         |         |         |
         |         |         |         |         | <--fixed
         |         |         |         |         |
         |         |         |         |         V
         |         |         |         |     "MRO_MCS_BASE"
         |         |         |         |     --------------
         |         |         |         |         |
         |         |         |         |         | <--fixed
         |         |         |         |         |
         |         |         |         |         V
         |         |         |         |     "MRO_MCS_AZ_GIMBAL_REF"
         |         |         |         |     -----------------------
         |         |         |         |         |
         |         |         |         |         | <--ck
         |         |         |         |         |
         |         |         |         |         V
         |         |         |         |     "MRO_MCS_AZ_GIMBAL"
         |         |         |         |     ---------------------------+
         |         |         |         |         |                      |
         |         |         |         |         | <--fixed             |
         |         |         |         |         |                      |
         |         |         |         |         V                      |
         |         |         |         |     "MRO_MCS_EL_GIMBAL_REF"    |
         |         |         |         |     -----------------------    |
         |         |         |         |         |                      |
         |         |         |         |         | <--ck                |
         |         |         |         |         |                      |
         |         |         |         |         V                      |
         |         |         |         |     "MRO_MCS_EL_GIMBAL"        |
         |         |         |         |     -------------------        |
         |         |         |         |         |                      |
         |         |         |         |         | <--fixed    fixed--> |
         |         |         |         |         |                      |
         |         |         |         |         V                      V
         |         |         |         |     "MRO_MCS"  "MRO_MCS_SOLAR_TARGET"
         |         |         |         |     ---------  ----------------------
         |         |         |         |
         |         |         |         |
         |         |         |         |
         |         |         |         | <--fixed
         |         |         |         |
         |         |         |         V
         |         |         |    "MRO_MARCI_BASE"
         |         |         |    -------------------------+
         |         |         |         |                   |
         |         |         |         | <--fixed          | <--fixed
         |         |         |         |                   |
         |         |         |         V                   V
         |         |         |    "MRO_MARCI_VIS"    "MRO_MARCI_UV"
         |         |         |    ---------------    --------------
         |         |         | 
         |         |         | 
         |         |         | 
         |         |         +---------------------------------+ 
         |         |         |                                 |
         |         |         | <--fixed                        |<--fixed
         |         |         |                                 |
         |         |         V                                 V
         |         |   "MRO_HIRISE_LOOK_DIRECTION"  "MRO_HIRISE_OPTICAL_AXIS"
         |         |   ---------------------------  -------------------------
         |         | 
         |         | 
         |         | 
         |         | <--fixed
         |         | 
         |         V 
         |     "MRO_CTX_BASE"
         |     --------------
         |         |
         |         | <--fixed
         |         |
         |         V
         |     "MRO_CTX"
         |     ---------
         |
         |
         |
         | <--fixed
         |
         V
    "MRO_CRISM_BASE"
    ----------------
         |
         | <--ck
         |
         V
    "MRO_CRISM_ART"
    ---------------
         |         
         | <--fixed
         |         
         V         
    "MRO_CRISM_VNIR"
    ----------------
         |         
         | <--fixed
         |         
         V         
    "MRO_CRISM_IR"
    --------------


   (*) BFR -- body-fixed rotating frame


Non Built-in Mars Frames:
-------------------------------------------------------------------------------


MME ``Of Date'' Frame
---------------------

   The MRO_MME_OF_DATE frame is based on Mean Mars Equator and IAU
   vector of date computed using IAU 2000 Mars rotation constants. This
   frame is called MME-D in [16]; it is the reference frame of the s/c
   orientation quaternions computed on-board and in the AtArPS program
   and stored in the the MRO s/c CK files.

   In this version of the FK MRO_MME_OF_DATE frame is implemented as as
   Euler frame mathematically identical to the PCK frame IAU_MARS based
   on IAU 2000 Mars rotation constants but without prime meridian
   rotation terms.

   The PCK data defining the IAU_MARS frame are:

      BODY499_POLE_RA          = (  317.68143   -0.1061      0.  )
      BODY499_POLE_DEC         = (   52.88650   -0.0609      0.  )
      BODY499_PM               = (  176.630    350.89198226  0.  )

   These values are from:

      Seidelmann, P.K., Abalakin, V.K., Bursa, M., Davies, M.E., Bergh, C
      de, Lieske, J.H., Oberst, J., Simon, J.L., Standish, E.M., Stooke,
      and Thomas, P.C. (2002). "Report of the IAU/IAG Working Group on
      Cartographic Coordinates and Rotational Elements of the Planets and
      Satellites: 2000," Celestial Mechanics and Dynamical Astronomy, v.8
      Issue 1, pp. 83-111.

   Here pole RA/Dec terms in the PCK are in degrees and degrees/century;
   the rates here have been converted to degrees/sec. Prime meridian
   terms from the PCK are disregarded.

   The 3x3 transformation matrix M defined by the angles is

      M = [    0.0]   [angle_2]   [angle_3]
                   3           1           3

   Vectors are mapped from the J2000 base frame to the MRO_MME_OF_DATE
   frame via left multiplication by M.

   The relationship of these Euler angles to RA/Dec for the
   J2000-to-IAU Mars Mean Equator and IAU vector of date transformation
   is as follows:

      angle_1 is        0.0
      angle_2 is pi/2 - Dec * (radians/degree)
      angle_3 is pi/2 + RA  * (radians/degree), mapped into the
                                                range 0 < angle_3 < 2*pi
                                                        -

   Since when we define the MRO_MME_OF_DATE frame we're defining the
   *inverse* of the above transformation, the angles for our Euler frame
   definition are reversed and the signs negated:

      angle_1 is -pi/2 - RA  * (radians/degree), mapped into the
                                                 range 0 < angle_3 < 2*pi
                                                         -
      angle_2 is -pi/2 + Dec * (radians/degree)
      angle_3 is         0.0

   Then our frame definition is:

   \begindata

      FRAME_MRO_MME_OF_DATE        =  -74900
      FRAME_-74900_NAME            = 'MRO_MME_OF_DATE'
      FRAME_-74900_CLASS           =  5
      FRAME_-74900_CLASS_ID        =  -74900
      FRAME_-74900_CENTER          =  499
      FRAME_-74900_RELATIVE        = 'J2000'
      FRAME_-74900_DEF_STYLE       = 'PARAMETERIZED'
      FRAME_-74900_FAMILY          = 'EULER'
      FRAME_-74900_EPOCH           =  @2000-JAN-1/12:00:00
      FRAME_-74900_AXES            =  ( 3  1  3 )
      FRAME_-74900_UNITS           = 'DEGREES'
      FRAME_-74900_ANGLE_1_COEFFS  = (  -47.68143
                                          0.33621061170684714E-10 )
      FRAME_-74900_ANGLE_2_COEFFS  = (  -37.1135
                                         -0.19298045478743630E-10 )
      FRAME_-74900_ANGLE_3_COEFFS  = (    0.0                     )
      FRAME_-74900_ROTATION_STATE  = 'INERTIAL'

   \begintext

   NOTE 1: The frame definition above will work ONLY with the SPICE
   Toolkits version N0058 or later. Should a need to use this FK with
   an older version of the toolkit (N0057 or earlier) arise, the
   definition above could be replaced with the following keywords:

      FRAME_MRO_MME_OF_DATE        = -74900
      FRAME_-74900_NAME            = 'MRO_MME_OF_DATE'
      FRAME_-74900_CLASS           = 4
      FRAME_-74900_CLASS_ID        = -74900
      FRAME_-74900_CENTER          = -74
      TKFRAME_-74900_SPEC          = 'ANGLES'
      TKFRAME_-74900_RELATIVE      = 'MRO_MME_2000'
      TKFRAME_-74900_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74900_AXES          = ( 1,   2,   3   )
      TKFRAME_-74900_UNITS         = 'DEGREES'

   These keywords simply map the MRO_MME_OF_DATE frame to the
   MRO_MME_2000 frame defined later in this FK. The error introduced by
   such replacement will be about 0.2 milliradian.

   NOTE 2: In order to "freeze" the MRO_MME_OF_DATE frame at an
   arbitrary epoch, in the definition above replace the keyword

      FRAME_-74900_ROTATION_STATE  = 'INERTIAL'

   with the keyword

      FRAME_-74900_FREEZE_EPOCH    = @YYYY-MM-DD/HR:MN:SC.###

   where YYYY-MM-DD/HR:MN:SC.### is the freeze epoch given as ET. For
   example, to freeze this frame at 2006-02-06, which the nominal
   freeze epoch for cruise, provide the keyword with this value:

      FRAME_-74900_FREEZE_EPOCH    = @2006-02-06/00:00:00.000


MME ``2000'' Frame
------------------

   The MRO_MME_2000 frame is the MRO_MME_OF_DATE frame frozen at J2000.
   For computing efficiency reasons this frame is defined as a fixed
   offset frame relative to the J2000 frame. The rotation matrix
   provided in the definition was computed using the following PXFORM
   call:

         CALL PXFORM( 'MRO_MME_OF_DATE', 'J2000', 0.D0, MATRIX )

   \begindata

      FRAME_MRO_MME_2000           = -74901
      FRAME_-74901_NAME            = 'MRO_MME_2000'
      FRAME_-74901_CLASS           = 4
      FRAME_-74901_CLASS_ID        = -74901
      FRAME_-74901_CENTER          = 499
      TKFRAME_-74901_SPEC          = 'MATRIX'
      TKFRAME_-74901_RELATIVE      = 'J2000'
      TKFRAME_-74901_MATRIX        = (
 
         0.6732521982472339       0.7394129276360180       0.0000000000000000
        -0.5896387605430040       0.5368794307891331       0.6033958972853946
         0.4461587269353556      -0.4062376142607541       0.7974417791532832

                                     )

   \begintext


Spacecraft Bus Frame
-------------------------------------------------------------------------------
 

   The spacecraft frame (or AACS control frame) is defined by the s/c design 
   as follows [from 4]:

      -  Z axis is parallel to the nominal HIRISE boresight;
 
      -  Y axis is anti-parallel to the MOI thrust vector;

      -  X axis completes the right hand frame;

      -  the origin of the frame is centered on the launch vehicle
         separation plane.

   (In [4] this frame is designated as "M" frame.)

   These diagrams illustrates the s/c frame:


      -Y view:
      --------                     .o.      HGA
                                 .' | `.
                               .'   |   `.
                           -------------------
                            `.             .'
                              `-._______.-'
                                    o
                               ____/_\____
                              /           \
                             /             \
       Direction            /               \
       of flight            \ +Xsc   +Ysc (into the page)
       <-------              \ <----x      / 
                    ..........o     |     o..........
                 .-'         /|     |     |\         `-.
      SAPX    .-'           / |     V +Zsc| \           `-.     SAMX
           .-\\            /  |-----------|  \            //-.     
        .-'   \\          /   |   |   |   |   \          //   `-.
      -'       \\       ./    .___|   |___.    \.       //       `-  
     \          \\   .-'          .___.          `-.   //          /
      \          \\-'                  HiRISE       `-//          /
       \       .-'                                     `-.       /
        \   .-'                     |                     `-.   /
          -'                        |                        `-
                                    V Nadir

      +Z view:
      --------
                                 . ---- .
                              .'         `. HGA
                            .'             `.
                           /                 \
                          .     .-------.     .
                          |     |   o   |     |
                          .     \       /     .
                           \     \     /     /
                            `.    \   /    .'
                              `.   \ /   .'
      SAPX                      ` --o-- '                      SAMX
      ========================o_____H_____o========================
                              |   / _ \   |
                              |  | '_' | HiRISE
                              |---\___/---|
                              |           |
       Direction              |           |
       of flight              |      +Zsc (out of the page)
       <-------                <----o ____.
                           +Xsc   \_|_/ 
                                   /|\
                                    V 
                                     +Ysc
                                     
                                     

   Since the S/C bus attitude is provided by a C kernel (see [3] for
   more information), this frame is defined as a CK-based frame.

   \begindata

      FRAME_MRO_SPACECRAFT         = -74000
      FRAME_-74000_NAME            = 'MRO_SPACECRAFT'
      FRAME_-74000_CLASS           = 3
      FRAME_-74000_CLASS_ID        = -74000
      FRAME_-74000_CENTER          = -74
      CK_-74000_SCLK               = -74
      CK_-74000_SPK                = -74

   \begintext


MRO Science Instrument Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for MRO science instruments --
   CRISM, CTX, HIRISE, MARCI, MCS, ONC, and SHARAD.


CRISM Frames
------------

   The following frames are defined for CRISM:

      -  CRISM base frame (MRO_CRISM_BASE) -- fixed w.r.t. to the s/c
         frame and nominally has +X axis co-aligned with the s/c +Z
         axis, +Y axis co-aligned with the s/c +X axis, and +Z axis
         co-aligned with the s/c +Y axis.

      -  CRISM articulation frame (MRO_CRISM_ART) -- rotates about +Z
         axis w.r.t. CRISM_BASE frame (and, therefore, defined as a
         CK-based frame) and co-aligned with the CRISM_BASE at "0"
         (nadir) scanner position;

      -  CRISM Visual and Near InfraRed apparent FOV frame
         (MRO_CRISM_VNIR) -- fixed w.r.t. MRO_CRISM_ART and has the +Z
         axis along boresight (the instrument slit center), the -X axis
         along gimbal rotation axis, and the +Y axis completing a
         right-handed frame;

      -  CRISM InfraRed apparent FOV frame (MRO_CRISM_IR) -- fixed
         w.r.t. and defined identically to the MRO_CRISM_VNIR;

   This diagram illustrates CRISM frames for CRISM scanner in "0" (nadir) 
   position:

                                 . ---- .
                              .'         `. HGA
                            .'             `.
                           /                 \
                          .     .-------.     .
                          |     |   o   |     |
                          .     \       /     .
                           \     \     /     /
                            `.    \   /    .'
                              `.   \ /   .'
      SAPX                      ` --o-- '                      SAMX
      ========================o_____H_____o========================
                              |   / _ \   |
                              |  | '      |
                              |---    ^   |
                      +Ycrism_base   .|.+Xcrism_vnir/ir
       Direction    +Ycrism_vnir/ir  |||  |
       of flight              |  <----o   |
       <-------                <----o |___.               +Zsc, +Xcrism_base,
                           +Xsc   \_|_|                   and +Zcrism_vnir/ir 
                                   /|\V +Zcrism_base      are out of the page
                                    V  
                                     +Ysc

   The rest of the comments and frame definitions in this section were
   copied ``as is'' from ``MRO_CRISM_FK_0000_000_N_1.TF'' ([18]).

   ``MRO_CRISM_FK_0000_000_N_1.TF'' Section ``Version and Date''
   -------------------------------------------------------------

      Version 0.1 -- September 14, 2006 -- Lillian Nguyen, JHU/APL
  
         Added alignment information and text.

      Version 0.0 -- April 25, 2006 -- Wen-Jong Shyong, JHU/APL

         Initial Release.


   ``MRO_CRISM_FK_0000_000_N_1.TF'' Section ``References''
   -------------------------------------------------------

      1. CRISM pointing sign conventions, "CALRPT_26_1_V2_PointSign.ppt",
         received from David Humm (JHU/APL).

      2. MRO alignment report, "MRO-final-alignment_REV-G.xls", received
         from David Humm.

      3. "CRISM Alignment Test Report", JHU/APL drawing number 7398-9600.

      4. Discussion between Scott Turner and David Humm regarding CRISM
         alignment.


   ``MRO_CRISM_FK_0000_000_N_1.TF'' Section ``Contact Information''
   ----------------------------------------------------------------

      Lillian Nguyen, JHU/APL, (443)-778-5477, Lillian.Nguyen@jhuapl.edu


   ``MRO_CRISM_FK_0000_000_N_1.TF'' Section ``CRISM Frame Definitions''
   --------------------------------------------------------------------

   The nominal CRISM base frame is defined such that Z is the gimbal axis of
   rotation, X is the projection of the instrument boresight (slit center)
   onto the plane normal to the gimbal axis at 0 degrees, and Y completes
   the right-handed frame. This nominal frame differs from the spacecraft
   frame by the following axis relabelling:

      X                    = Z
       nominal CRISM base     sc

      Y                    = X
       nominal CRISM base     sc

      Z                    = Y
       nominal CRISM base     sc,

   written as a rotation matrix:

     [    ]   [ 0 1 0 ]
     [ R1 ] = [ 0 0 1 ]
     [    ]   [ 1 0 0 ]

   The axes of the nominal CRISM base frame are illustrated in the
   diagram below.


                 ^  instrument slit center
                 |
                _|_                                         X (Z  )
               |   |                                      ^     sc
               |   | gimbal axis into the page            |
           ____|___|____                                  |
          |             |                                 |
          |     .O.     |------> S/C velocity             o------> Y (X  )
          |____/   \____|                            Z (Y  )           sc
          ____/_____\____                                sc
          ///////////////
            spacecraft

   In [2] we are given three alignment matrices from which we can determine
   the rotation matrix taking vectors from the CRISM base frame to vectors in
   the spacecraft frame. The first of those matrices takes vectors in the
   CRISM optical cube frame to vectors in the HiRISE frame:

      [   ]HiRISE          [  0.999917  0.001050  0.012822 ]
      [ A ]             =  [ -0.001056  0.999999  0.000460 ]
      [   ]CRISM           [ -0.012821 -0.000473  0.999918 ]

   where the CRISM frame is defined in [2] as

      Y = Axis of rotation of the instrument.
      Z = Axis perpendicular to Y and lying in the plane formed by the Gimbal
          axis and the Optical axis.
      X axis completes a right-hand-rectangular coordinate frame.

   Note that it is believed that Z was determined using the CRISM optical axis
   (the normal projected from the mirror on the rear of the secondary mirror),
   while the CRISM base frame definition uses the slit center. We will adjust
   for the angular difference between the optical axis and slit center vectors
   later, with matrix [ R2 ].

   Due to circumstances in the alignment tests, the instrument team changed
   the theta Y value (the measure of the gimbal axis rotation) from 0.735
   degree to 0.755 degree [4], resulting in a corrected [ A ] matrix taking
   vectors from the CRISM frame to the HiRISE frame. [2] describes these
   circumstances as follows:

   "Theta Y for CRISM is a measure of the gimbal axis rotation. At the time of
   this measurement, the amount of this rotation was not controlled (CRISM was
   not powered). However, the measured value of 0.735 degree is very close to
   the Pre-environmental measurement of 0.755 degree, which was taken in a
   powered state at zero degrees gimbal rotation."

   The calculations used to determine the corrected [ A ] matrix are
   explained below.

   If we describe the CRISM to HiRISE rotation as

      [   ]HiRISE        [ a d g ]
      [ A ]            = [ b e h ]
      [   ]CRISM         [ c f i ],

   then theta Y is defined in [2] as

      theta Y = atan ( g/i ) (degrees),

   and is equal to 0.755 degrees [4]. To determine the corrected matrix, we
   solve a set of equations defined by the following constraints:

      1) atan(g/i) = 0.755 deg
      2) norm( (g, h, i) ) = 1            (axes are unit length)
      3) dot( (d, e, f), (g, h, i) ) = 0  (orthogonality of axes)

   Note that constraint 3) uses the gimbal axis vector, (d, e, f),  which
   we assume remains fixed.

   Solving for g and i (taking the positive solution to the quadratic equation
   in constraint 2), then readjusting vector (a, b, c) by using the cross
   product to form an orthogonal frame, we get the corrected matrix:

      [   ]HiRISE     [  0.999912630217  0.001049999963  0.013176853036 ]
      [ A ]        =  [ -0.001056141713  0.999998986721  0.000460000010 ]
      [   ]CRISM      [ -0.013176361292 -0.000472999993  0.999913076097 ]

   The second matrix given in [2] takes vectors from the Star Tracker 1
   alignment cube frame to vectors in the HiRISE optical axis frame. The
   alignment report gives the following matrix:

      [   ]HiRISE        [ -0.966071  0.000673 -0.258275 ]
      [ B ]            = [ -0.087542  0.939949  0.329897 ]
      [   ]Star Tr. 1    [  0.242988  0.341314 -0.907999 ]

   The third matrix takes vectors from the Star Tracker 1 alignment cube frame
   to vectors in the spacecraft frame:

      [   ]spacecraft    [ -0.966218  0.000257 -0.257726 ]
      [ C ]            = [ -0.087724  0.939960  0.329818 ]
      [   ]Star Tr. 1    [  0.242337  0.341285 -0.908184 ]

   Finally, we describe the rotation [ R2 ] that takes vectors from the CRISM'
   frame to the CRISM frame defined above, where the CRISM' frame differs
   only in that the Z axis is the projection of the slit center (and not of
   the optical axis as in the CRISM frame) onto the plane perpendicular to
   the gimbal axis. We assume that the gimbal axis as measured in [2] and [3]
   is the same, and use the following measurements given in [3] to determine
   the angle between the optical axis and the slit center:

      gimbal axis:                   [0.0008242301 0.9999940951 0.0033362399]
      slit center at home (0 deg.):  [0.013894654 -0.003154637  0.999898488 ]
      optical axis at home (0 deg.): [0.014603728 -0.003256776  0.999888056 ]

   This rotation is determined by first projecting both the slit center and
   optical axes onto the plane normal to the gimbal axis and calculating the
   angle, alpha, between the two projected vectors, then creating the rotation
   matrix about the gimbal axis (Y in the CRISM frame). The angle between the
   two projected vectors is calculated to be:

      alpha = 0.040635904 deg.

   and the rotation matrix is:

     [    ]CRISM    [ 0.999999748496 0.000000000000 -0.000709230259 ]
     [ R2 ]       = [ 0.000000000000 1.000000000000  0.000000000000 ]
     [    ]CRISM'   [ 0.000709230259 0.000000000000  0.999999748496 ]

   The measured alignment of the CRISM base frame relative to the spacecraft
   frame then is obtained by multiplying the three alignment matrices with the
   two rotation matrices (R1 for axis relabelling and R2 to take into account
   that the measurement in [2] used the CRISM optical axis) as follows (note
   that we are using the corrected [ A ] matrix from above):

      [   ]spacecraft    [   ] [   ]t [   ] [    ] [    ]
      [ R ]            = [ C ] [ B ]  [ A ] [ R2 ] [ R1 ]
      [   ]CRISM base    [   ] [   ]  [   ] [    ] [    ]

   where 't' denotes the matrix transpose. This gives us:

      [   ]spacecraft    [ 0.0117851901010  0.9999301878218  0.0008536843642 ]
      [ R ]            = [ 0.0004938596793 -0.0008595641829  0.9999995086259 ]
      [   ]CRISM base    [ 0.9999304302785 -0.0117847627097 -0.0005039553291 ]

   To review, the sequence of transformations taking vectors from the CRISM
   base frame to the spacecraft frame is as follows:

           CRISM base --[R1]--> CRISM' (slit center as boresight)

               CRISM' --[R2]--> CRISM  (optical axis as boresight)

               CRISM  ---[A]--> HiRISE

               HiRISE ---[B]--> Star Tracker 1

       Star Tracker 1 ---[C]--> spacecraft.

   CRISM Base Frame (MRO_CRISM_BASE):

   \begindata

      FRAME_MRO_CRISM_BASE         = -74011
      FRAME_-74011_NAME            = 'MRO_CRISM_BASE'
      FRAME_-74011_CLASS           = 4
      FRAME_-74011_CLASS_ID        = -74011
      FRAME_-74011_CENTER          = -74
      TKFRAME_-74011_SPEC          = 'MATRIX'
      TKFRAME_-74011_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74011_MATRIX        = ( 0.0117851901010
                                       0.0004938596793
                                       0.9999304302785
                                       0.9999301878218
                                      -0.0008595641829
                                      -0.0117847627097
                                       0.0008536843642
                                       0.9999995086259
                                      -0.0005039553291 )
   \begintext


   The CRISM articulation frame is defined such that Z is the gimbal axis of
   rotation, X is the projection of the instrument slit center onto the plane
   normal to the gimbal axis at theta degrees, and Y completes the right-
   handed frame. At gimbal home (0 degree), the articulation frame is
   identical to the CRISM base frame, MRO_CRISM_BASE. The articulation frame
   rotates the base frame about the gimbal axis and is C-kernel based
   (see [5]).

   CRISM Articulation Frame (MRO_CRISM_ART):

   \begindata

      FRAME_MRO_CRISM_ART          = -74012
      FRAME_-74012_NAME            = 'MRO_CRISM_ART'
      FRAME_-74012_CLASS           = 3
      FRAME_-74012_CLASS_ID        = -74012
      FRAME_-74012_CENTER          = -74
      CK_-74012_SCLK               = -74999
      CK_-74012_SPK                = -74

   \begintext


   The MRO_CRISM_VNIR frame is defined such that the Z axis is the boresight
   (the instrument slit center), the -X axis is the gimbal rotation axis, and
   the Y axis completes a right-handed frame. The nominal mapping of
   CRISM VNIR coordinates to CRISM articulation frame coordinates is

      X      = -Z
       VNIR      art

      Y      =  Y
       VNIR      art

      Z      =  X
       VNIR      art,

   or as a rotation matrix:

      [   ]articulation   [  0  0  1 ]
      [ R ]             = [  0  1  0 ]
      [   ]nominal VNIR   [ -1  0  0 ]

   We will use the following measured alignments given in [3] to adjust the
   nominal frame:

      gimbal axis:         [0.0008242301, 0.9999940951, 0.0033362399]
      slit center at home: [0.013894654, -0.003154637,  0.999898488 ]

   To determine the VNIR boresight, we rotate the gimbal axis (Z   ) by
                                                                art
   theta degrees about Y    , where theta is the angle between the measured
                        art
   gimbal axis and measured slit center at home. Note that rotation of the
   slit center at home about the gimbal axis was adjusted for in the CRISM
   base frame's intermediate matrix [ R2 ]. The angular separation, theta,
   is calculated to be:

      theta = 89.9889570902 degrees

   Rotating Z    by theta about Y   , we obtain
             art                 art

      Z     = [ 0.9999999814  0.0000000000  0.0001927351 ]
       VNIR

   We obtain the VNIR Y axis by taking the cross product of the gimbal axis,
   Z  , with the boresight, Z   :
    art                      VNIR

      Y     = [ 0.0  1.0  0.0 ]
       VNIR

   The VNIR X axis completes the right-handed frame:

      X     = [ 0.0001927351  0.0000000000  -0.9999999814 ]
       VNIR

   Thus, the rotation matrix taking vectors from the VNIR frame to the
   articulation frame is

      [   ]   [  0.0001927351  0.0000000000  0.9999999814 ]
      [ R ] = [  0.0000000000  1.0000000000  0.0000000000 ]
      [   ]   [ -0.9999999814  0.0000000000  0.0001927351 ]

   CRISM VNIR Frame (MRO_CRISM_VNIR):

   \begindata

      FRAME_MRO_CRISM_VNIR         = -74017
      FRAME_-74017_NAME            = 'MRO_CRISM_VNIR'
      FRAME_-74017_CLASS           = 4
      FRAME_-74017_CLASS_ID        = -74017
      FRAME_-74017_CENTER          = -74
      TKFRAME_-74017_SPEC          = 'MATRIX'
      TKFRAME_-74017_RELATIVE      = 'MRO_CRISM_ART'
      TKFRAME_-74017_MATRIX        = ( 0.0001927351
                                       0.0000000000
                                      -0.9999999814
                                       0.0000000000
                                       1.0000000000
                                       0.0000000000
                                       0.9999999814
                                       0.0000000000
                                       0.0001927351 )

   \begintext

   The MRO_CRISM_IR frame is defined identically to the MRO_CRISM_VNIR
   frame. Any offsets between the VNIR and IR are accounted for in the
   camera model described in the MRO CRISM Instrument Kernel.

   \begindata

      FRAME_MRO_CRISM_IR           = -74018
      FRAME_-74018_NAME            = 'MRO_CRISM_IR'
      FRAME_-74018_CLASS           = 4
      FRAME_-74018_CLASS_ID        = -74018
      FRAME_-74018_CENTER          = -74
      TKFRAME_-74018_SPEC          = 'MATRIX'
      TKFRAME_-74018_RELATIVE      = 'MRO_CRISM_VNIR'
      TKFRAME_-74018_MATRIX        = ( 1.0
                                       0.0
                                       0.0
                                       0.0
                                       1.0
                                       0.0
                                       0.0
                                       0.0
                                       1.0 )

   \begintext


CTX Frames
----------

   The following frames are defined for CTX:

      -  CTX base frame (MRO_CTX_BASE) -- fixed w.r.t. and nominally
         co-aligned with the MRO_SPACECRAFT frame;

      -  CTX apparent FOV frame (MRO_CTX) -- fixed w.r.t. MRO_CTX_BASE
         and nominally co-aligned with it; it has +Z along boresight,
         +Y along the detector line, and +X completing the right hand
         frame;

   This diagram illustrates CTX frames:

                                 . ---- .
                              .'         `. HGA
                            .'             `.
                           /                 \
                          .     .-------.     .
                          |     |   o   |     |
                          .     \       /     .
                           \     \     /     /
                            `.    \   /    .'
                              `.   \ /   .'
      SAPX                      ` --o-- '                      SAMX
      ========================o_____H_____o========================
                              |   / _ \   |
                              |    ' '    |
                              |--- <----o |
                              | +Xctx*  | |
       Direction              |         | |
       of flight              |         V +Yctx*
       <-------                <----o ____.            +Zsc and +Zctx*
                           +Xsc   \_|_                   are out of the  
                                   /|\                      the page
                                    V 
                                     +Ysc

   The keyword sets below define CTX frames. Except cases were the
   source of the alignment data is specifically noted, these frame
   definitions incorporate the nominal alignment.
 
   The following CTX to HIRISE Direction Cosine Matrix (DCM) was
   provided in [14]:

         0.99999994   0.00021198  -0.00026011
        -0.00021196   0.99999998   0.00005486
         0.00026012  -0.00005481   0.99999996

   This matrix is incorporated in the MRO_CTX_BASE definition below.

   \begindata

      FRAME_MRO_CTX_BASE           = -74020
      FRAME_-74020_NAME            = 'MRO_CTX_BASE'
      FRAME_-74020_CLASS           = 4
      FRAME_-74020_CLASS_ID        = -74020
      FRAME_-74020_CENTER          = -74
      TKFRAME_-74020_SPEC          = 'MATRIX'
      TKFRAME_-74020_RELATIVE      = 'MRO_HIRISE_LOOK_DIRECTION'
      TKFRAME_-74020_MATRIX        = ( 
                                       0.99999994
                                      -0.00021196
                                       0.00026012
                                       0.00021198
                                       0.99999998
                                      -0.00005481
                                      -0.00026011
                                       0.00005486
                                       0.99999996
                                     )

      FRAME_MRO_CTX                = -74021
      FRAME_-74021_NAME            = 'MRO_CTX'
      FRAME_-74021_CLASS           = 4
      FRAME_-74021_CLASS_ID        = -74021
      FRAME_-74021_CENTER          = -74
      TKFRAME_-74021_SPEC          = 'ANGLES'
      TKFRAME_-74021_RELATIVE      = 'MRO_CTX_BASE'
      TKFRAME_-74021_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74021_AXES          = ( 1,   2,   3   )
      TKFRAME_-74021_UNITS         = 'DEGREES'

   \begintext


HIRISE Frames
-------------

   The following frames are defined for HIRISE:

      -  HIRISE ``look direction'' frame (MRO_HIRISE_LOOK_DIRECTION) --
         fixed w.r.t. and nominally co-aligned with the MRO_SPACECRAFT
         frame; it has +Z along the camera boresight, nominally defined
         as the view direction of the detector pixel 0 of CCD 5/Channel
         1 for Mars in-focus observations (*), +Y along the detector lines,
         and +X completing the right hand frame;

      -  HIRISE optical axis frame (MRO_HIRISE_OPTICAL_AXIS) -- fixed
         w.r.t. MRO_SPACECRAFT and is nominally rotated from it by
         +0.45 degrees about +Y axis; it has +Z along the camera
         optical axis, +Y along the detector lines, and +X completing
         the right hand frame;

         (*) the actual boresight direction shifts by up to 1 arcsecond
             (5 pixels) depending on the instrument focus setting.

   This diagram illustrates HIRISE frames:

                                   .o.      HGA
                                 .' | `.
                               .'   |   `.
                           -------------------
                            `.             .'
                              `-._______.-'
                                    o
                               ____/_\____
                              /           \
                             /             \
       Direction            /               \
       of flight            \ +Xsc   +Ysc (into the page)
       <-------              \ <----x      / 
                    ..........o     |     o..........
                 .-'         /|     |     |\         `-.
      SAPX    .-'                   V +Zsc| \           `-.     SAMX
           .-\\             +Xh_*  -------|  \            //-.     
        .-'   \\               <----x HiRISE  \          //   `-.
      -'       \\       ./    .___| | |___.    \.       //       `-  
     \          \\   .-'          ._|_.          `-.   //          /
      \          \\-'               V               `-//          /
       \       .-'                    +Zh_*            `-.       /
        \   .-'                                           `-.   /
         `-'            0.45 deg ->||<-                      `-'
                                   || 
                                   VV 
                            +Zh_oa    +Zh_ld (co-aligned with s/c +Z)

                                    |
                                    | Nadir
                                    V

   The keyword sets below define HIRISE frames. Except cases were the
   source of the alignment data is specifically noted, these frame
   definitions incorporate the nominal alignment.


MRO_HIRISE_LOOK_DIRECTION Frame Rotation Provided in FK Versions 0.0-0.5

   In the FK versions 0.0-0.6 this frame was named MRO_HIRISE. It was
   defined as zero offset frame relative to the MRO_SPACECRAFT frame.


MRO_HIRISE_LOOK_DIRECTION Frame Rotation Provided in FK Version 0.6

   In the FK version 0.6 this frame was named MRO_HIRISE. It was
   defined as a fixed offset frame relative to the MRO_SPACECRAFT frame
   using pre-launch, ground calibrated alignment shown below.

   Combining the following Tracker 1 cube to S/C Direction Cosine
   Matrix (DCM) (from [14]):

        -0.96621811   0.00025732  -0.25772564
        -0.08772412   0.93995995   0.32981780
         0.24233665   0.34128468  -0.90818375

   with Tracker 1 cube to HiRISE DCM (from [14]):

        -0.96607109   0.00067328  -0.25827542
        -0.08754160   0.93994921   0.32989689
         0.24298789   0.34131369  -0.90799882

   results in this HIRISE ``look direction'' to S/C DCM:

         0.99999975  -0.00019674  -0.00067690
         0.00019676   0.99999999   0.00003113
         0.00067689  -0.00003127   0.99999978

   which formatted as the frame definition keyword looks like this:

      TKFRAME_-74699_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74699_MATRIX        = ( 
                                       0.99999975
                                       0.00019676
                                       0.00067689
                                      -0.00019674
                                       0.99999999
                                      -0.00003127
                                      -0.00067690
                                       0.00003113
                                       0.99999978
                                     )


MRO_HIRISE_LOOK_DIRECTION Frame Rotation Provided in FK Version 0.7+

   The MRO_HIRISE_LOOK_DIRECTION frame definition below includes the
   following rotation with respect to the MRO_HIRISE_OPTICAL_AXIS frame
   derived by combining preliminary in-flight alignment for
   MRO_HIRISE_OPTICAL_AXIS frame with with pre-launch, ground
   calibrated alignment for the HIRISE boresight:
   
      TKFRAME_-74699_RELATIVE      = 'MRO_HIRISE_OPTICAL_AXIS'
      TKFRAME_-74699_MATRIX        = ( 
                                       0.99996484
                                       0.00020285
                                       0.00838273
                                      -0.00019649
                                       0.99999969
                                      -0.00075976
                                      -0.00838288
                                       0.00075809
                                       0.99996458
                                     )

MRO_HIRISE_LOOK_DIRECTION Frame Rotation Provided in FK Version 1.4+

   The MRO_HIRISE_LOOK_DIRECTION frame definition below includes the
   following rotation with respect to the MRO_HIRISE_OPTICAL_AXIS frame,
   updated from the previous rotation to compensate for the change in
   the orientation of the MRO_HIRISE_OPTICAL_AXIS frame in the FK 1.4
   (to keep the MRO_HIRISE_LOOK_DIRECTION frame oriented with respect to
   the spacecraft frame the same way as it was in the versions 0.7-1.3
   of the FK):

      TKFRAME_-74699_RELATIVE      = 'MRO_HIRISE_OPTICAL_AXIS'
      TKFRAME_-74699_MATRIX        = ( 
                                       0.999964843747
                                       0.000206345964
                                       0.008382642316
                                      -0.000196487983
                                       0.999999288260
                                      -0.001176805736
                                      -0.008382879179
                                       0.001175117275
                                       0.999964172576
                                     )


   This matrix is currently incorporated in the MRO_HIRISE_LOOK_DIRECTION 
   frame definition.

   \begindata

      FRAME_MRO_HIRISE_LOOK_DIRECTION = -74699
      FRAME_-74699_NAME            = 'MRO_HIRISE_LOOK_DIRECTION'
      FRAME_-74699_CLASS           = 4
      FRAME_-74699_CLASS_ID        = -74699
      FRAME_-74699_CENTER          = -74
      TKFRAME_-74699_SPEC          = 'MATRIX'
      TKFRAME_-74699_RELATIVE      = 'MRO_HIRISE_OPTICAL_AXIS'
      TKFRAME_-74699_MATRIX        = ( 
                                       0.999964843747
                                       0.000206345964
                                       0.008382642316
                                      -0.000196487983
                                       0.999999288260
                                      -0.001176805736
                                      -0.008382879179
                                       0.001175117275
                                       0.999964172576
                                     )
   \begintext


MRO_HIRISE_OPTICAL_AXIS Frame Rotation Provided in FK Versions 0.0-0.6

   In the FK versions 0.0-0.6 this frame was named MRO_HIRISE_IF. It
   was defined relative to the MRO_SPACECRAFT with single rotation of
   -0.45 degrees about Y axis. This rotation rotation angle was an
   error; it should have been by +0.45 degrees.

   
MRO_HIRISE_OPTICAL_AXIS Frame Rotation Provided in FK Version 0.7+

   The MRO_HIRISE_OPTICAL_AXIS frame definition below incorporates the
   following preliminary in-flight alignment with respect to the
   spacecraft frame provided by Jeff Anderson, USGS on September 22,
   2005:

      TKFRAME_-74690_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74690_MATRIX        = (
                                       0.999970
                                       0.000000
                                      -0.007706
                                       0.000000
                                       1.000000
                                       0.000727
                                       0.007706
                                      -0.000727
                                       0.999970
                                     )

MRO_HIRISE_OPTICAL_AXIS Frame Rotation Provided in FK Version 1.4+

   The MRO_HIRISE_OPTICAL_AXIS frame definition below incorporates the
   updated alignment provided by Laszlo Keszthelyi, USGS on February
   12, 2009, compensating for the correction to the optical distortion
   model made in 2008:

      TKFRAME_-74690_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74690_MATRIX        = (
                                       0.99997031
                                       0.00000000
                                      -0.00770600
                                       0.00000882
                                       0.99999935
                                       0.00114399
                                       0.00770599
                                      -0.00114402
                                       0.99996965
                                     )

   This matrix is currently incorporated in the MRO_HIRISE_OPTICAL_AXIS
   frame definition.

   \begindata

      FRAME_MRO_HIRISE_OPTICAL_AXIS = -74690
      FRAME_-74690_NAME            = 'MRO_HIRISE_OPTICAL_AXIS'
      FRAME_-74690_CLASS           = 4
      FRAME_-74690_CLASS_ID        = -74690
      FRAME_-74690_CENTER          = -74
      TKFRAME_-74690_SPEC          = 'MATRIX'
      TKFRAME_-74690_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74690_MATRIX        = (
                                       0.99997031
                                       0.00000000
                                      -0.00770600
                                       0.00000882
                                       0.99999935
                                       0.00114399
                                       0.00770599
                                      -0.00114402
                                       0.99996965
                                     )

   \begintext


MARCI Frames
----------

   The following frames are defined for MARCI:

      -  MARCI base frame (MRO_MARCI_BASE) -- fixed w.r.t. and rotated
         by +95 degrees about +Z axis w.r.t. the MRO_SPACECRAFT frame;

      -  MARCI apparent VIS FOV frame (MRO_MARCI_VIS) -- fixed w.r.t.
         MRO_MARCI_BASE and nominally co-aligned with it; it has +Z
         along boresight, +X along the detector lines, and +Y
         completing the right hand frame;

      -  MARCI apparent UV FOV frame (MRO_MARCI_UV) -- fixed w.r.t.
         MRO_MARCI_BASE and nominally co-aligned with it; it has +Z
         along boresight, +X along the detector lines, and +Y
         completing the right hand frame;

   This diagram illustrates MARCI frames:

                                 . ---- .
                              .'         `. HGA
                            .'             `.
                           /                 \
                          .     .-------.     .
                          |     |   o   |     |
                          .     \       /     .
                           \     \     /     /
                            `.    \   /    .'
                              `.   \ /   .'
      SAPX                      ` --o-- '                      SAMX
      ========================o_____H_____o========================
                              |   / _ \   |
                              |  | '_' | HiRISE
                              |-- \___/ --|
                              |   .->     |
           -----              |o-'   +Ymarci*
            5 deg             |`          |
             .-'          +Xsc <`---o ____.            +Zsc and +Zmarci*
            '                    V\_|_/                  are out of the  
                          +Xmarci* /|\                      the page
                                    V 
       <-------                      +Ysc
       Direction
       of flight

   Nominally the following set of rotations can be used to align the 
   MRO spacecraft frame with the MARCI base frame:

      Msc->marci = [ 95.0 ]z * [ 0.0 ]y * [ 0.0 ]x

   By co-locating pixels for several overlapping images taken on
   different orbits the MARCI team at MSSS derived the following
   updated alignment angles (from [19]):

      Msc->marci = [ 95.5 ]z * [ 0.5 ]y * [ 0.475 ]x

   These angles are used in the MARCI base frame definition below.

   The keyword sets below define MARCI frames. Except cases were the
   source of the alignment data is specifically noted, these frame
   definitions incorporate the nominal alignment.
 
   \begindata

      FRAME_MRO_MARCI_BASE         = -74400
      FRAME_-74400_NAME            = 'MRO_MARCI_BASE'
      FRAME_-74400_CLASS           = 4
      FRAME_-74400_CLASS_ID        = -74400
      FRAME_-74400_CENTER          = -74
      TKFRAME_-74400_SPEC          = 'ANGLES'
      TKFRAME_-74400_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74400_ANGLES        = ( -0.475, -0.5, -95.5 )
      TKFRAME_-74400_AXES          = (  1,      2,     3   )
      TKFRAME_-74400_UNITS         = 'DEGREES'

      FRAME_MRO_MARCI_VIS          = -74410
      FRAME_-74410_NAME            = 'MRO_MARCI_VIS'
      FRAME_-74410_CLASS           = 4
      FRAME_-74410_CLASS_ID        = -74410
      FRAME_-74410_CENTER          = -74
      TKFRAME_-74410_SPEC          = 'ANGLES'
      TKFRAME_-74410_RELATIVE      = 'MRO_MARCI_BASE'
      TKFRAME_-74410_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74410_AXES          = ( 1,   2,   3   )
      TKFRAME_-74410_UNITS         = 'DEGREES'

      FRAME_MRO_MARCI_UV           = -74420
      FRAME_-74420_NAME            = 'MRO_MARCI_UV'
      FRAME_-74420_CLASS           = 4
      FRAME_-74420_CLASS_ID        = -74420
      FRAME_-74420_CENTER          = -74
      TKFRAME_-74420_SPEC          = 'ANGLES'
      TKFRAME_-74420_RELATIVE      = 'MRO_MARCI_BASE'
      TKFRAME_-74420_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74420_AXES          = ( 1,   2,   3   )
      TKFRAME_-74420_UNITS         = 'DEGREES'

   \begintext


MCS Frames
----------

   The following frames are defined for MCS:

      -  MCS base frame (MRO_MCS_BASE) -- fixed w.r.t., and nominally
         co-aligned with the MRO_SPACECRAFT frame; the definition of
         this frame incorporates instrument misalignment determined by
         measuring the alignment cube orientation w.r.t. to the
         spacecraft at the time of instrument installation;

      -  MCS Azimuth Gimbal "Reference" position frame
         (MRO_MCS_AZ_GIMBAL_REF) -- fixed w.r.t., and nominally
         coalligned with the MCS_BASE frame, this frame is defined by
         requiring the MRO_MCS_AZ_GIMBAL_REF +Z axis be coalligned with
         the MCS Azimuth physical rotation axis, while at the same time
         minimizing the angle between the MRO_MCS_BASE +X axis and the
         MRO_MCS_AZ_GIMBAL_REF +X axis.

      -  MCS Azimuth Gimbal frame (MRO_MCS_AZ_GIMBAL) -- rotates about
         the +Z axis by AZ angle w.r.t. MCS_AZ_GIMBAL_REF frame (and,
         therefore, is defined as a CK-based frame) and is co-aligned
         with the MCS_AZ_GIMBAL_REF frame at an azimuth scan angle of
         180 degrees (2782 counts).

      -  MCS Elevation Gimbal "Reference" position frame
         (MRO_MCS_EL_GIMBAL_REF) -- fixed w.r.t., and nominally
         coaligned with the MCS_AZ_GIMBAL frame, this frame is defined
         by requiring the MRO_MCS_EL_GIMBAL_REF +Y axis be coaligned
         with the MCS Elevation physical rotation axis, while at the
         same time minimizing the angle between the MRO_MCS_AZ_GIMBAL
         +Z axis and the MRO_MCS_EL_GIMBAL_REF +Z axis.

      -  MCS Elevation Gimbal frame (MRO_MCS_EL_GIMBAL) -- rotates
         about -Y axis by EL angle w.r.t. MCS_EL_GIMBAL_REF frame (and,
         therefore, is defined as a CK-based frame) and is co-aligned
         with the MCS_EL_GIMBAL_REF frame at an elevation scan angle of
         90 degrees (1891 counts).

      -  MCS telescope boresight frame (MRO_MCS) -- fixed w.r.t, and
         nominally coaligned with the MRO_MCS_EL_GIMBAL frame, this
         frame is defined by requiring the MRO_MCS +X axis be in the
         direction of the telescope boresight, and requiring that the
         MRO_MCS Z axis be aligned with the detector arrays in such a
         sense that, when viewing the forward limb (near the +X axis),
         positive rotations about the MRO_MCS +Y axis cause Z to
         increase.

      -  MCS solar target frame (MRO_MCS_SOLAR_TARGET) -- fixed w.r.t.
         the MRO_MCS_AZ_GIMBAL frame, is defined such that its +Z axis
         is normal to the solar target plate and +Y axis is co-aligned
         with the AZ_GIMBAL frame's +Y axis. This frame is rotated from
         the AZ_GIMBAL frame by 15 degrees about +Y axis.

   Assuming that in (180,90) (AZ,EL) angle position the telescope
   boresight is pointing along the s/c +X axis (nominal), all six MCS
   frames -- BASE, AZ_GIMBAL_REF, AZ_GIMBAL, EL_GIMBAL_REF, EL_GIMBAL, and
   MRO_MCS -- will be co-aligned as shown in this diagram (SOLAR_TARGET
   frame is not shown):

                                 . ---- .
                              .'         `. HGA
                            .'             `.
                           /                 \
                          .     .-------.     .
                          |     |   o   |     |
                          .     \       /     .
                           \     \     /     /
                            `.    \   /    .'
                              `.   \ /   .'
      SAPX                      ` --o-- '                      SAMX
      ========================o_____H_____o========================
                              |   / _ \   |
                              |  | '_' | HiRISE
                              |-- \___/ --|
                              |           |
                     +Xmcs    | <----o    |
                              |      |    |
       <-------           +Xsc <----o|____.           
       Direction                  \_|V +Ymcs              +Zsc, +Zmcs
       of flight                   /|\                 and nadir are out of 
                                    V                         the page
                                     +Ysc          
                                                       Azimuth rotation is
                                                             about +Z

                                                       Elevation rotation is
                                                             about +Y

   The keyword sets below define MCS frames. Except cases were the
   source of the alignment data is specifically noted, these frame
   definitions incorporate the nominal alignment.
 
   The following MCS to HIRISE Direction Cosine Matrix (DCM) was
   provided in [14]:

         0.99996956  -0.00780193  -0.00010482
         0.00780199   0.99996939   0.00059227
         0.00010020  -0.00059307   0.99999982

   This DCM was provided in the MRO_MCS_BASE frame definition as the 
   following keyword in the FK versions 0.6 to 1.1:

      TKFRAME_-74501_MATRIX        = (  
                                       0.99996956
                                       0.00780199
                                       0.00010020
                                      -0.00780193
                                       0.99996939
                                      -0.00059307
                                      -0.00010482
                                       0.00059227
                                       0.99999982
                                     )

   Based on the analysis on the flight data the offsets for the MCS AZ
   and EL gimbals were determined with respect to the spacecraft. To
   make FK consistent with these results, starting with the FK version
   1.2 the MRO_MCS_BASE frame was re-defined to be with respect to the
   MRO_SPACECRAFT frame with zero offset rotation.

   Based on the analysis of the off-track observations carried out in
   2008 (see [20]), the MCS Team, JPL determined that the azimuth axis
   was tilted relative to the s/c +Z axis by 0.431 degrees towards the
   line 25.8 degrees off the s/c -Y axis towards the s/c +X axis. The
   following set of rotations aligning the s/c frame with the MCS base
   frame was incorporated into the MRO_MCS_BASE frame definition below
   to account for this tilt:
 
       base
      M    = [-25.8 deg]   [+0.431 deg]   [+25.8 deg]
       sc               Z              X             Z

   Incorporating the tilt into the MRO_MCS_BASE frame "raised" the
   frame's +X axis above the s/c XY plane, invalidating the previous
   zero EL angle offset included in the definition of the
   MRO_MCS_EL_GIMBAL_REF frame. To fix this, the offset was re-set from
   -0.208 degrees to -0.03 degrees.
      
   Based on the analysis of the off-track observations carried out in
   2012 (see [21]), the MCS Team, JPL determined that the MCS base was
   additionally rotated -0.118 degrees about a vector in the S/C XY
   plane oriented 18.5 degrees from the -Y axis towards the +X axis.
   The set of rotations in the MRO_MCS_BASE frame definition aligning
   the s/c frame with the MCS base was chaged to be
 
       base
      M    = [169.913679 deg]   [-0.432158 deg]   [-169.914119 deg]
       sc                    Z                 X                   Z

   to account for this additional tilt.

   In addition to tilting the base frame, the alignment of the
   MRO_MCS_EL_GIMBAL_REF frame was updated to apply a -0.132 degree
   correction to MCS elevation angle by changing the offset rotation
   from -0.03 degrees to -0.162 degrees.
      
   (The frame definition below contain the opposite of these rotations 
   because Euler angles specified in it define transformations from MCS
   base frame to the s/c frame -- see [1].)

   \begindata

      FRAME_MRO_MCS_BASE           = -74501
      FRAME_-74501_NAME            = 'MRO_MCS_BASE'
      FRAME_-74501_CLASS           = 4
      FRAME_-74501_CLASS_ID        = -74501
      FRAME_-74501_CENTER          = -74
      TKFRAME_-74501_SPEC          = 'ANGLES'
      TKFRAME_-74501_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74501_ANGLES        = (  169.914119,    0.432158, -169.913679 )
      TKFRAME_-74501_AXES          = (    3,           1,           3        )
      TKFRAME_-74501_UNITS         = 'DEGREES'

      FRAME_MRO_MCS_AZ_GIMBAL_REF  = -74502
      FRAME_-74502_NAME            = 'MRO_MCS_AZ_GIMBAL_REF'
      FRAME_-74502_CLASS           = 4
      FRAME_-74502_CLASS_ID        = -74502
      FRAME_-74502_CENTER          = -74
      TKFRAME_-74502_SPEC          = 'ANGLES'
      TKFRAME_-74502_RELATIVE      = 'MRO_MCS_BASE'
      TKFRAME_-74502_ANGLES        = ( 0.0, 0.0, -0.46 )
      TKFRAME_-74502_AXES          = ( 1,   2,    3    )
      TKFRAME_-74502_UNITS         = 'DEGREES'

      FRAME_MRO_MCS_AZ_GIMBAL      = -74503
      FRAME_-74503_NAME            = 'MRO_MCS_AZ_GIMBAL'
      FRAME_-74503_CLASS           = 3
      FRAME_-74503_CLASS_ID        = -74503
      FRAME_-74503_CENTER          = -74
      CK_-74503_SCLK               = -74
      CK_-74503_SPK                = -74

      FRAME_MRO_MCS_EL_GIMBAL_REF  = -74504
      FRAME_-74504_NAME            = 'MRO_MCS_EL_GIMBAL_REF'
      FRAME_-74504_CLASS           = 4
      FRAME_-74504_CLASS_ID        = -74504
      FRAME_-74504_CENTER          = -74
      TKFRAME_-74504_SPEC          = 'ANGLES'
      TKFRAME_-74504_RELATIVE      = 'MRO_MCS_AZ_GIMBAL'
      TKFRAME_-74504_ANGLES        = ( 0.0, -0.162,  0.0  )
      TKFRAME_-74504_AXES          = ( 1,    2,     3    )
      TKFRAME_-74504_UNITS         = 'DEGREES'

      FRAME_MRO_MCS_EL_GIMBAL      = -74505
      FRAME_-74505_NAME            = 'MRO_MCS_EL_GIMBAL'
      FRAME_-74505_CLASS           = 3
      FRAME_-74505_CLASS_ID        = -74505
      FRAME_-74505_CENTER          = -74
      CK_-74505_SCLK               = -74
      CK_-74505_SPK                = -74

      FRAME_MRO_MCS                = -74500
      FRAME_-74500_NAME            = 'MRO_MCS'
      FRAME_-74500_CLASS           = 4
      FRAME_-74500_CLASS_ID        = -74500
      FRAME_-74500_CENTER          = -74
      TKFRAME_-74500_SPEC          = 'ANGLES'
      TKFRAME_-74500_RELATIVE      = 'MRO_MCS_EL_GIMBAL'
      TKFRAME_-74500_ANGLES        = ( 0.0,  0.0,   0.0 )
      TKFRAME_-74500_AXES          = ( 1,    2,     3   )
      TKFRAME_-74500_UNITS         = 'DEGREES'

      FRAME_MRO_MCS_SOLAR_TARGET   = -74506
      FRAME_-74506_NAME            = 'MRO_MCS_SOLAR_TARGET'
      FRAME_-74506_CLASS           = 4
      FRAME_-74506_CLASS_ID        = -74506
      FRAME_-74506_CENTER          = -74
      TKFRAME_-74506_SPEC          = 'ANGLES'
      TKFRAME_-74506_RELATIVE      = 'MRO_MCS_AZ_GIMBAL'
      TKFRAME_-74506_ANGLES        = ( 0.0, -15.0,  0.0 )
      TKFRAME_-74506_AXES          = ( 1,     2,    3   )
      TKFRAME_-74506_UNITS         = 'DEGREES'

   \begintext


ONC Frames
----------

   The following frame is defined for ONC:

      -  ONC apparent FOV frame (MRO_ONC) -- fixed w.r.t.
         MRO_SPACECRAFT and has +Z along boresight, +X along the
         detector lines, and +Y completing the right hand frame;

   ONC is mounted on the -Z side of the s/c and points approximately 30
   degrees off the s/c +Y axis towards s/c -Z axis and sightly to the
   +X s/c side. These diagrams illustrate the ONC frame orientation:

   +X side view:
   -------------

                  HGA |`.
                      |  \
                    .'|   .._ 
                  ,'  |    | |                     +Xsc and +Xonc are 
                 o    |    | |                       put of the page.
                  `-. |    | |
                     `|   '|.'
                      |  / ||
                      |.'  |/
                           o                  SAPX
                            \_====================
                             \   /      \_____.             Nadir
                              \ /        \___/  HiRISE      --->
                               .----------\.
                   +Yonc <.    |           |
                           `.  |           |
                             `o| +Xsc      |
                             / .____ o----> +Zsc
                            /      \_|_/ 
                           V        /|\
                     +Zonc           V
                                      +Ysc
                        /     |
                       /<---->|
                          29.6 deg (projected on s/c Y-Z plane)


   -Z side view:
   -------------
                                 . ---- .
                              .'         `. HGA
                            .'      o      `.
                           /        |        \
                          .         |         .
                          |         o         |
                          .       .' `.        .
                           \    o'     `o    /
                            `.             .'
                              `.         .'
      SAMX                      ` --o-- '                      SAPX
      ========================o_____H_____o========================
                              |           |
                              |           |
                              |-----------|
                              |         .>|
                              |       .' +Xonc
                              |     o'    |
                              .____ x\---> +Xsc 
                                  \_|_\    
          +Zsc is into             /|\ V
            the page                V   +Zonc
                                     +Ysc
          +Yonc is out
           of the page              |     \
                                    |<---->\
                                       5.7 deg (projected on s/c X-Y plane)

   The s/c frame can be transformed into the ONC frame in nominal
   orientation by the following three rotations (derived from the
   projected angles shown above): first by -119.6 degrees about +X,
   second by +4.96 degrees about +Y and finally by "?" degrees about
   +Z. (The third rotation is not derivable from projected angles and
   is assumed to be zero.)

   Based on the cruise observations the ONC team determined the actual
   ONC alignment relative to the s/c frame. According to [15] the
   following rotations are required to align the s/c spacecraft frame
   with the ONC frame:

       onc
      M    = [1.109164]   [-61.065813]   [169.075079]
       sc              Z              X              Y

   The definition below incorporates these rotations.

   (The frame definitions below contain the opposite of these rotations 
   because Euler angles specified in them define transformations from ONC
   frames to the s/c frame -- see [1].)

   \begindata

      FRAME_MRO_ONC                = -74030
      FRAME_-74030_NAME            = 'MRO_ONC'
      FRAME_-74030_CLASS           = 4
      FRAME_-74030_CLASS_ID        = -74030
      FRAME_-74030_CENTER          = -74
      TKFRAME_-74030_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74030_SPEC          = 'ANGLES'
      TKFRAME_-74030_ANGLES        = ( -169.075079, 61.065813, -1.109164 )
      TKFRAME_-74030_AXES          = (    2,         1,         3        )
      TKFRAME_-74030_UNITS         = 'DEGREES'

   \begintext


SHARAD Frames
----------

   The following frame is defined for SHARAD:

      -  SHARAD frame (MRO_SHARAD) -- fixed w.r.t. MRO_SPACECRAFT and
         nominally co-aligned with it;

   The keyword set below defines the SHARAD frame. In this version of
   the FK it incorporates the nominal alignments.

   \begindata

      FRAME_MRO_SHARAD             = -74070
      FRAME_-74070_NAME            = 'MRO_SHARAD'
      FRAME_-74070_CLASS           = 4
      FRAME_-74070_CLASS_ID        = -74070
      FRAME_-74070_CENTER          = -74
      TKFRAME_-74070_SPEC          = 'ANGLES'
      TKFRAME_-74070_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74070_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74070_AXES          = ( 1,   2,   3   )
      TKFRAME_-74070_UNITS         = 'DEGREES'

   \begintext


MRO Antenna Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for MRO antennas -- HGA,
   LGA1, LGA2, and UHF.


High Gain Antenna Frame
-----------------------

   The HGA boresight frame -- MRO_HGA -- is defined as follows ([4],[13]):

      -  Z axis is along the HGA reflector central symmetry axis (boresight 
         axis) and points from the reflector surface towards the feed horn;
         
      -  X axis is parallel to the inner gimbal rotation axis and 
         points from the gimbal towards the antenna center;
         
      -  Y axis completes to the right hand frame;
      
      -  the origin of this frame is located at the intersection of the 
         antenna reflector symmetry axis and a plane containing HGA  
         reflector rim circle.
     
   In stowed configuration HGA boresight (+Z axis) points approximately
   along S/C -Y axis (14.5 degrees off it towards +Z.) In deployed
   configuration orientation of the HGA with respect to the s/c varies
   as the HGA moves constantly using two gimbals to track Earth.
 

HGA Baseplate and Gimbal Drive Frames
-------------------------------------

   The frame chain for HGA includes:

      -  baseplate frame that is fixed w.r.t. to the s/c frame

      -  inner gimbal frame that rotates w.r.t. to the baseplate frame

      -  outer gimbal frame rotates w.r.t. to the inner gimbal frame

      -  boresight frame (described above) that is fixed w.r.t. to the
         outer gimbal frame.

   In "0" angle position the baseplate frame, both gimbal frames, and
   the boresight frame are co-aligned.

   The MRO HGA baseplate frame is defined as follows:
   
      -  +Z axis is s/c -Z axis;
         
      -  +Y axis is s/c -Y axis;

      -  +X axis completes the right hand frame and is parallel to
         the s/c +X axis
      
      -  the origin of this frame is located at the intersection of the 
         inner gimbal rotation axis and a plane perpendicular to this 
         rotation axis and containing the outer gimbal rotation axis.

   The MRO HGA inner gimbal frame:
   
      -  Y axis is along the inner gimbal rotation axis; in deployed
         configuration with the inner and outer gimbal angles set to
         zero it points along the baseplate frame +Y axis;
         
      -  X axis is such that in deployed configuration with 
         the inner and outer gimbal angles set to zero it points along 
         the baseplate frame +X axis;

      -  Z axis completes the right hand frame and in deployed
         configuration with the inner and outer gimbal angles set to
         zero it points along the baseplate frame +Z axis;
      
      -  the origin of this frame is located at the intersection of the 
         inner gimbal rotation axis and a plane perpendicular to this 
         rotation axis and containing the outer gimbal rotation axis.
            
   The MRO HGA outer gimbal frame:
   
      -  X axis is along the outer gimbal rotation axis and points
         along the baseplate +X in deployed configuration with the
         inner and outer gimbal angles set to zero;
         
      -  Y axis is such that in deployed configuration with the inner
         and outer gimbal angles set to zero it points along the
         baseplate +Y axis;

      -  Z axis completes to the right hand frame and in deployed
         configuration with the inner and outer gimbal angles set to
         zero it points along the baseplate +Z axis;
      
      -  the origin of this frame is located at the intersection of the
         outer gimbal rotation axis and a plane perpendicular to this
         rotation axis and containing the HGA frame origin;

   When antenna is deployed and both gimbals are in zero position, the
   axes of the baseplate, inner gimbal, and outer gimbal frames are
   co-aligned while the HGA frame is rotated by +90 degrees about +Z
   axis with respect to them. The diagram below illustrates this:

                                    
                                    |  HGA Inner 
                                    . Gimbal Axis
                                    |

                                 . ---- .
                              .'  +Xhga  `. HGA (shown in "0" angle 
                            .'      ^      `.        position)
                           /        |        \
                          .     .---|---. +Yhga
                          |    |    x---->    |
                          .     \             .
                           \     \  ^ +Yhgabp/
                       +Xhgabp    \ | +Yhgaig
                       +Xhgaig.    \| +Yhgaog
      -- . -- . -      +Xhgaog <----x --'                      SAMX
      HGA Outer         ======o_____H_____o========================
     Gimbal Axis              |   / _ \   |
                              |  | '_' | HiRISE
                              |---\___/---|
                              |           |
       Direction              |           |
       of flight              |      +Zsc (out of the page)
       <-------                <----o ____.
                           +Xsc   \_|_/ 
                                   /|\
                                    V 
                                     +Ysc

                                             +Zhga, +Zhgabp, +Zhgaig, 
                                               and +Zhgaog are 
                                                into the page 
  
   The gimbal frames are defined such that rotation axis designations 
   are consistent with [4].

   [22] provided the following HGA baseplate DCM w.r.t. to the s/c frame

       9.999999804376672e-01 -1.977993966601806e-04  2.295934537796784e-07
      -1.977992524732815e-04 -9.999998719444498e-01 -4.658181490872693e-04
       3.217319645753285e-07  4.658180945613699e-04 -9.999998915066927e-01

   and the following antenna boresight direction in the ``ideal'' HGA
   frame:

       5.221214000806500e-05 -1.497243425405999e-03  9.999988777673789e-01

   The baseplate DCM and the rotations aligning the antenna's Z axis
   with the vector above are provided in the corresponding frame
   definitons below.


HGA Frame Definitions
---------------------

   The sets of keywords below contain definitions for the HGA frames.   

   \begindata

      FRAME_MRO_HGA_BASEPLATE    = -74211
      FRAME_-74211_NAME          = 'MRO_HGA_BASEPLATE'
      FRAME_-74211_CLASS         = 4
      FRAME_-74211_CLASS_ID      = -74211
      FRAME_-74211_CENTER        = -74
      TKFRAME_-74211_SPEC        = 'MATRIX'
      TKFRAME_-74211_RELATIVE    = 'MRO_SPACECRAFT'
      TKFRAME_-74211_MATRIX      = ( 

       9.999999804376672e-01 -1.977993966601806e-04  2.295934537796784e-07
      -1.977992524732815e-04 -9.999998719444498e-01 -4.658181490872693e-04
       3.217319645753285e-07  4.658180945613699e-04 -9.999998915066927e-01

                                   )

      FRAME_MRO_HGA_INNER_GIMBAL = -74212
      FRAME_-74212_NAME          = 'MRO_HGA_INNER_GIMBAL'
      FRAME_-74212_CLASS         = 3
      FRAME_-74212_CLASS_ID      = -74212
      FRAME_-74212_CENTER        = -74
      CK_-74212_SCLK             = -74
      CK_-74212_SPK              = -74

      FRAME_MRO_HGA_OUTER_GIMBAL = -74213
      FRAME_-74213_NAME          = 'MRO_HGA_OUTER_GIMBAL'
      FRAME_-74213_CLASS         = 3
      FRAME_-74213_CLASS_ID      = -74213
      FRAME_-74213_CENTER        = -74
      CK_-74213_SCLK             = -74
      CK_-74213_SPK              = -74

      FRAME_MRO_HGA              = -74214
      FRAME_-74214_NAME          = 'MRO_HGA'
      FRAME_-74214_CLASS         = 4
      FRAME_-74214_CLASS_ID      = -74214
      FRAME_-74214_CENTER        = -74
      TKFRAME_-74214_SPEC        = 'ANGLES'
      TKFRAME_-74214_RELATIVE    = 'MRO_HGA_OUTER_GIMBAL'
      TKFRAME_-74214_ANGLES      = ( -0.08578576, -0.00299154, -90.0 )
      TKFRAME_-74214_AXES        = (  1,           2,            3   )
      TKFRAME_-74214_UNITS       = 'DEGREES'

   \begintext


Low Gain Antennas
-----------------

   Both LGA boresight frames -- MRO_LGA1 and MRO_LGA2 -- are defined as
   follows:

      -  +Z axis is along the LGA boresight vector;
         
      -  +Y axis is along the HGA +Y axis;
         
      -  +X completes the right hand frame;
      
      -  the origin of the frame is located at the center of the LGA
         patch.
     
   Both LGAs are mounted on and do not move with respect to the HGA.
   Therefore their frames are specified as fixed offset frames with
   respect to the HGA boresight frame. 

   According to [4] the LGA boresights point along the following directions
   in HGA outer gimbal frame:

      LGA1 (truss-mounted LGA) --  (0.0, -0.422618,  0.906308)
      LGA2 (TWTA-mounted LGA)  --  (0.0,  0.906308, -0.422618)

   The diagram below illustrates the LGA1 and LGA2 frames:


                   ^ +Xlga1
                    \
                     \             
                     .x LGA1                       HGA is shown in  
          +Zlga1  .-' |`.                        "0" angle position.
                <'    ^ +Xhga    
                      |   .._                   +Xsc is out of the page 
           +Zhga      |    | |                    
                 <----x    | |  ^ +Zlga2       +Yhga, +Ylga1, and +Ylga2
                      |    | | /                   are into the  page.
                      |   '|.'/
                      | LGA2 x
                      |.'  |/ `.   +Xlga2
                           o    `>             SAPX
                            \_====================
                             \   /      \_____.
                              \ /        \___/  HiRISE
                               .----------\.
                               |           |
                               |           |
                               | +Xsc      |
                               .____ o----> +Zsc      ------->
                                   \_|_/               Nadir
                                    /|\
                                     V
                                      +Ysc


   As seen on the diagram the LGA1 frame is rotated from the HGA frame
   by -25 degrees about +Y while the LGA2 frame is rotated by +115
   degrees from HGA frame about +Y.

   (The frame definitions below contain the opposite of these rotations 
   because Euler angles specified in them define transformations from LGA
   frames to the HGA frame -- see [1].)

   \begindata

      FRAME_MRO_LGA1             = -74220
      FRAME_-74220_NAME          = 'MRO_LGA1'
      FRAME_-74220_CLASS         = 4
      FRAME_-74220_CLASS_ID      = -74220
      FRAME_-74220_CENTER        = -74
      TKFRAME_-74220_SPEC        = 'ANGLES'
      TKFRAME_-74220_RELATIVE    = 'MRO_HGA'
      TKFRAME_-74220_ANGLES      = ( 0.0, 0.0,   25.0 )
      TKFRAME_-74220_AXES        = ( 3,   1,      2   )
      TKFRAME_-74220_UNITS       = 'DEGREES'

      FRAME_MRO_LGA2             = -74230
      FRAME_-74230_NAME          = 'MRO_LGA2'
      FRAME_-74230_CLASS         = 4
      FRAME_-74230_CLASS_ID      = -74230
      FRAME_-74230_CENTER        = -74
      TKFRAME_-74230_SPEC        = 'ANGLES'
      TKFRAME_-74230_RELATIVE    = 'MRO_HGA'
      TKFRAME_-74230_ANGLES      = ( 0.0, 0.0, -115.0 )
      TKFRAME_-74230_AXES        = ( 3,   1,      2   )
      TKFRAME_-74230_UNITS       = 'DEGREES'

   \begintext


UHF Antenna
-----------

   The UHF frame -- MRO_UHF -- is defined as follows:

      -  +Z axis is along the antenna boresight and co-aligned with the
         s/c +Z axis;
         
      -  +Y axis is co-aligned with the s/c +Y axis;
      
      -  +X completes the right hand frame;
         
      -  the origin of this frame is located at the geometric center of 
         the antenna.
     
   Since UHF antenna is rigidly mounted on the s/c bus, it is defined as
   a fixed offset frame co-aligned with the s/c frame.

   (The frame definition below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from antenna 
   to s/c frame -- see [1].)

   \begindata

      FRAME_MRO_UHF              = -74240
      FRAME_-74240_NAME          = 'MRO_UHF'
      FRAME_-74240_CLASS         = 4
      FRAME_-74240_CLASS_ID      = -74240
      FRAME_-74240_CENTER        = -74
      TKFRAME_-74240_SPEC        = 'ANGLES'
      TKFRAME_-74240_RELATIVE    = 'MRO_SPACECRAFT'
      TKFRAME_-74240_ANGLES      = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74240_AXES        = ( 3,   2,   1   )
      TKFRAME_-74240_UNITS       = 'DEGREES'

   \begintext


MRO Solar Array Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for MRO Solar Array frames.


Solar Array Frames
------------------

   Both SA frames -- MRO_SAPX and MRO_SAMX -- are defined as follows:

      -  +Z axis is perpendicular to and points away from the array
         solar cell side (note that this is different from [4] where
         SAMX +Z axis is defined to point away from the non-cell side
         of the array);
         
      -  +X axis parallel to the long side of the array and points from
         the end of the array towards the gimbal;
      
      -  +Y axis completes the right hand frame;
         
      -  the origin of this frame is located at the intersection of the
         inner gimbal rotation axis and a plane perpendicular to this
         rotation axis and containing the outer gimbal rotation axis.

   When SAs are deployed they move constantly using two gimbals to
   track Sun.
 

Solar Array Gimbal Drive Frames
-------------------------------

   The frame chain for each of the arrays includes:

       -  baseplate frame that is fixed w.r.t. to the s/c frame

       -  inner gimbal frame that rotates w.r.t. to the baseplate frame

       -  outer gimbal frame that rotates w.r.t. to the inner gimbal
          frame

       -  boresight frame (described above) that is fixed w.r.t. to the
          outer gimbal frame.

   When SAPX is in "0" angle position its baseplate frame, both gimbal
   frames, and the boresight frame are co-aligned. When SAMX is in "0"
   angle position its baseplate frame and both gimbal frames are
   co-aligned while the boresight frame is rotated by 180 degrees about
   +X axis w.r.t. to them.

   The MRO SAPX baseplate frame is defined as follows:
   
      -  +Z axis is s/c -Y axis;
         
      -  +Y axis is along the inner gimbal rotation axis and points 
         towards the HGA side of the deck;

      -  +X axis completes the right hand frame and is along the outer
         gimbal rotation axis;
      
      -  the origin of this frame is located at the intersection of the 
         inner gimbal rotation axis and a plane perpendicular to this 
         rotation axis and containing the outer gimbal rotation axis.

   The MRO SAMX baseplate frame is defined as follows:
   
      -  +Z axis is s/c +Y axis;
         
      -  +Y axis is along the inner gimbal rotation axis and points 
         towards HGA side of the deck;

      -  +X axis completes the right hand frame and is along the outer
         gimbal rotation axis;
      
      -  the origin of this frame is located at the intersection of the 
         inner gimbal rotation axis and a plane perpendicular to this 
         rotation axis and containing the outer gimbal rotation axis.

   The MRO SAPX and SAMX inner gimbal frame:
   
      -  +Y axis is along the inner gimbal rotation axis; in deployed
         configuration with the inner and outer gimbal angles set to
         zero it points along the baseplate +Y axis;
         
      -  +X axis is such that in deployed configuration with the inner
         and outer gimbal angles set to zero it points along the
         baseplate +X axis;

      -  +Z axis completes to the right hand frame and in deployed
         configuration wit the inner and outer gimbal angles set to
         zero it points along the baseplate +Z axis;
      
      -  the origin of this frame is located at the intersection of the 
         inner gimbal rotation axis and a plane perpendicular to this 
         rotation axis and containing the outer gimbal rotation axis.
      
   The MRO SA outer gimbal frame:
   
      -  +X axis is along the outer gimbal rotation axis and points
         along the baseplate +X in deployed configuration with the
         inner and outer gimbal angles set to zero;
         
      -  +Y axis is such that in deployed configuration with 
         the inner and outer gimbal angles set to zero it points along 
         the baseplate +Y axis;

      -  Z axis completes to the right hand frame and in deployed
         configuration with the inner and outer gimbal angles set to
         zero it points along the s/c +Z axis;
      
      -  the origin of this frame is located at the intersection of the
         outer gimbal rotation axis and a plane perpendicular to this
         rotation axis and containing the solar array frame origin;

   The diagram below illustrates the solar array baseplate, gimbal and
   cell-side frames in deployed "0" angle configuration:

          
                                   .o.      HGA
                                 .' | `.
                               .'   |   `.           +Zsapx** and +Zsamx 
                           -------------------       are out of the page
                            `.             .'
                              `-._______.-'           +Zsamx** are into
                                    o                     the page
                               ____/_\____
                              /           \
                  +Ysapxbp   /   +Xsa*x**  \   +Ysamxbp
                  +Ysapxig ^                 ^ +Ysamxig
                  +Ysapxog  \     .> <.     /  +Ysamxog
                  +Ysapx     \  .'     `.  / 
                    ..........o'         `x..........
                 .-'         /|          /|\         `-.
      SAPX    .-'           / |  +Ysamx / | \           `-.     SAMX
           .-\\            /  |------- v -|  \            //-.     
        .-'   \\          /   |   |       |   \          //   `-.
      -'       \\       ./    .___|   |___.    \.       //       `-  
     \          \\   .-'          .___.          `-.   //          /
      \          \\-'                  HiRISE       `-//          /
       \       .-'                                     `-.       /
        \   .-'                                           `-.   /
          -'              +Xsc <----x +Ysc (into the page)   `-
                                    | 
                                    |
       <-------                     V 
       Direction                      +Zsc
       of flight
                                    |
                                    | Nadir
                                    V

   The gimbal frames are defined such that rotation axis designations
   are consistent with [4]. Also according to [4] the SAPX and SAMX
   baseplate frames are rotated w.r.t. to the s/c frame as follows:

      SAPX: first by +165 degrees about +Y, then by +90 deg about +X

      SAPX: first by  +15 degrees about +Y, then by -90 deg about +X


Solar Array Frames Definitions
-----------------------------

   Two sets of keywords below contain definitions for these frames.   

   \begindata

      FRAME_MRO_SAPX_BASEPLATE     = -74311
      FRAME_-74311_NAME            = 'MRO_SAPX_BASEPLATE'
      FRAME_-74311_CLASS           = 4
      FRAME_-74311_CLASS_ID        = -74311
      FRAME_-74311_CENTER          = -74
      TKFRAME_-74311_SPEC          = 'ANGLES'
      TKFRAME_-74311_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74311_ANGLES        = ( 0.0, -165.0, -90.0 )
      TKFRAME_-74311_AXES          = (   3,    2,     1   )
      TKFRAME_-74311_UNITS         = 'DEGREES'

      FRAME_MRO_SAPX_INNER_GIMBAL  = -74312
      FRAME_-74312_NAME            = 'MRO_SAPX_INNER_GIMBAL'
      FRAME_-74312_CLASS           = 3
      FRAME_-74312_CLASS_ID        = -74312
      FRAME_-74312_CENTER          = -74
      CK_-74312_SCLK               = -74
      CK_-74312_SPK                = -74

      FRAME_MRO_SAPX_OUTER_GIMBAL  = -74313
      FRAME_-74313_NAME            = 'MRO_SAPX_OUTER_GIMBAL'
      FRAME_-74313_CLASS           = 3
      FRAME_-74313_CLASS_ID        = -74313
      FRAME_-74313_CENTER          = -74
      CK_-74313_SCLK               = -74
      CK_-74313_SPK                = -74

      FRAME_MRO_SAPX               = -74314
      FRAME_-74314_NAME            = 'MRO_SAPX'
      FRAME_-74314_CLASS           = 4
      FRAME_-74314_CLASS_ID        = -74314
      FRAME_-74314_CENTER          = -74
      TKFRAME_-74314_SPEC          = 'ANGLES'
      TKFRAME_-74314_RELATIVE      = 'MRO_SAPX_OUTER_GIMBAL'
      TKFRAME_-74314_ANGLES        = ( 0.0, 0.0, 0.0 )
      TKFRAME_-74314_AXES          = (   3,   2,   1 )
      TKFRAME_-74314_UNITS         = 'DEGREES'

      FRAME_MRO_SAMX_BASEPLATE     = -74321
      FRAME_-74321_NAME            = 'MRO_SAMX_BASEPLATE'
      FRAME_-74321_CLASS           = 4
      FRAME_-74321_CLASS_ID        = -74321
      FRAME_-74321_CENTER          = -74
      TKFRAME_-74321_SPEC          = 'ANGLES'
      TKFRAME_-74321_RELATIVE      = 'MRO_SPACECRAFT'
      TKFRAME_-74321_ANGLES        = ( 0.0, -15.0, 90.0 )
      TKFRAME_-74321_AXES          = (   3,   2,    1   )
      TKFRAME_-74321_UNITS         = 'DEGREES'

      FRAME_MRO_SAMX_INNER_GIMBAL  = -74322
      FRAME_-74322_NAME            = 'MRO_SAMX_INNER_GIMBAL'
      FRAME_-74322_CLASS           = 3
      FRAME_-74322_CLASS_ID        = -74322
      FRAME_-74322_CENTER          = -74
      CK_-74322_SCLK               = -74
      CK_-74322_SPK                = -74

      FRAME_MRO_SAMX_OUTER_GIMBAL  = -74323
      FRAME_-74323_NAME            = 'MRO_SAMX_OUTER_GIMBAL'
      FRAME_-74323_CLASS           = 3
      FRAME_-74323_CLASS_ID        = -74323
      FRAME_-74323_CENTER          = -74
      CK_-74323_SCLK               = -74
      CK_-74323_SPK                = -74

      FRAME_MRO_SAMX               = -74324
      FRAME_-74324_NAME            = 'MRO_SAMX'
      FRAME_-74324_CLASS           = 4
      FRAME_-74324_CLASS_ID        = -74324
      FRAME_-74324_CENTER          = -74
      TKFRAME_-74324_SPEC          = 'ANGLES'
      TKFRAME_-74324_RELATIVE      = 'MRO_SAMX_OUTER_GIMBAL'
      TKFRAME_-74324_ANGLES        = ( 0.0, 0.0, 180.0 )
      TKFRAME_-74324_AXES          = (   3,   2,   1   )
      TKFRAME_-74324_UNITS         = 'DEGREES'

   \begintext


Mars Reconnaissance Orbiter NAIF ID Codes -- Definitions
========================================================================

   This section contains name to NAIF ID mappings for the MRO mission.
   Once the contents of this file is loaded into the KERNEL POOL, these 
   mappings become available within SPICE, making it possible to use 
   names instead of ID code in the high level SPICE routine calls. 

   Spacecraft:
   -----------

      MARS RECONNAISSANCE ORBITER   -74
      MRO                           -74
      MRO_SPACECRAFT                -74000
      MRO_SPACECRAFT_BUS            -74000
      MRO_SC_BUS                    -74000

   Science Instruments:
   --------------------      

      MRO_CRISM                     -74010
      MRO_CRISM_VNIR                -74017
      MRO_CRISM_IR                  -74018

      MRO_CTX                       -74021

      MRO_HIRISE                    -74699
      MRO_HIRISE_CCD0               -74600
      MRO_HIRISE_CCD1               -74601
      MRO_HIRISE_CCD2               -74602
      MRO_HIRISE_CCD3               -74603
      MRO_HIRISE_CCD4               -74604
      MRO_HIRISE_CCD5               -74605
      MRO_HIRISE_CCD6               -74606
      MRO_HIRISE_CCD7               -74607
      MRO_HIRISE_CCD8               -74608
      MRO_HIRISE_CCD9               -74609
      MRO_HIRISE_CCD10              -74610
      MRO_HIRISE_CCD11              -74611
      MRO_HIRISE_CCD12              -74612
      MRO_HIRISE_CCD13              -74613

      MRO_MARCI                     -74400
      MRO_MARCI_VIS                 -74410
      MRO_MARCI_VIS_BLUE            -74411
      MRO_MARCI_VIS_GREEN           -74412
      MRO_MARCI_VIS_ORANGE          -74413
      MRO_MARCI_VIS_RED             -74414
      MRO_MARCI_VIS_NIR             -74415
      MRO_MARCI_UV                  -74420
      MRO_MARCI_UV_SHORT_UV         -74421
      MRO_MARCI_UV_LONG_UV          -74422

      MRO_MCS                       -74500
      MRO_MCS_A                     -74510
      MRO_MCS_A1                    -74511
      MRO_MCS_A2                    -74512
      MRO_MCS_A3                    -74513
      MRO_MCS_A4                    -74514
      MRO_MCS_A5                    -74515
      MRO_MCS_A6                    -74516
      MRO_MCS_B                     -74520
      MRO_MCS_B1                    -74521
      MRO_MCS_B2                    -74522
      MRO_MCS_B3                    -74523

      MRO_ONC                       -74030

      MRO_SHARAD                    -74070

   Antennas:
   ---------

      MRO_HGA_BASEPLATE             -74211
      MRO_HGA_INNER_GIMBAL          -74212
      MRO_HGA_OUTER_GIMBAL          -74213
      MRO_HGA                       -74214
      MRO_LGA1                      -74220
      MRO_LGA2                      -74230
      MRO_UHF                       -74240

   Solar Arrays:
   -------------

      MRO_SAPX_BASEPLATE            -74311
      MRO_SAPX_INNER_GIMBAL         -74312
      MRO_SAPX_OUTER_GIMBAL         -74313
      MRO_SAPX                      -74314
      MRO_SAPX_C1                   -74315
      MRO_SAPX_C2                   -74316
      MRO_SAPX_C3                   -74317
      MRO_SAPX_C4                   -74318

      MRO_SAMX_BASEPLATE            -74321
      MRO_SAMX_INNER_GIMBAL         -74322
      MRO_SAMX_OUTER_GIMBAL         -74323
      MRO_SAMX                      -74324
      MRO_SAMX_C1                   -74325
      MRO_SAMX_C2                   -74326
      MRO_SAMX_C3                   -74327
      MRO_SAMX_C4                   -74328

   The mappings summarized in this table are implemented by the keywords 
   below.

   \begindata

      NAIF_BODY_NAME += ( 'MARS RECONNAISSANCE ORBITER' )
      NAIF_BODY_CODE += ( -74                           )

      NAIF_BODY_NAME += ( 'MRO'                         )
      NAIF_BODY_CODE += ( -74                           )

      NAIF_BODY_NAME += ( 'MRO_SPACECRAFT'              )
      NAIF_BODY_CODE += ( -74000                        )

      NAIF_BODY_NAME += ( 'MRO_SPACECRAFT_BUS'          )
      NAIF_BODY_CODE += ( -74000                        )

      NAIF_BODY_NAME += ( 'MRO_SC_BUS'                  )
      NAIF_BODY_CODE += ( -74000                        )

      NAIF_BODY_NAME += ( 'MRO_CRISM'                   )
      NAIF_BODY_CODE += ( -74010                        )
 
      NAIF_BODY_NAME += ( 'MRO_CRISM_VNIR'              )
      NAIF_BODY_CODE += ( -74017                        )
 
      NAIF_BODY_NAME += ( 'MRO_CRISM_IR'                )
      NAIF_BODY_CODE += ( -74018                        )
 
      NAIF_BODY_NAME += ( 'MRO_CTX'                     )
      NAIF_BODY_CODE += ( -74021                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE'                  )
      NAIF_BODY_CODE += ( -74699                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD0'             )
      NAIF_BODY_CODE += ( -74600                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD1'             )
      NAIF_BODY_CODE += ( -74601                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD2'             )
      NAIF_BODY_CODE += ( -74602                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD3'             )
      NAIF_BODY_CODE += ( -74603                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD4'             )
      NAIF_BODY_CODE += ( -74604                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD5'             )
      NAIF_BODY_CODE += ( -74605                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD6'             )
      NAIF_BODY_CODE += ( -74606                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD7'             )
      NAIF_BODY_CODE += ( -74607                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD8'             )
      NAIF_BODY_CODE += ( -74608                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD9'             )
      NAIF_BODY_CODE += ( -74609                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD10'            )
      NAIF_BODY_CODE += ( -74610                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD11'            )
      NAIF_BODY_CODE += ( -74611                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD12'            )
      NAIF_BODY_CODE += ( -74612                        )
      
      NAIF_BODY_NAME += ( 'MRO_HIRISE_CCD13'            )
      NAIF_BODY_CODE += ( -74613                        )      
      
      NAIF_BODY_NAME += ( 'MRO_MARCI'                   )
      NAIF_BODY_CODE += ( -74400                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS'               )
      NAIF_BODY_CODE += ( -74410                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS_BLUE'          )
      NAIF_BODY_CODE += ( -74411                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS_GREEN'         )
      NAIF_BODY_CODE += ( -74412                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS_ORANGE'        )
      NAIF_BODY_CODE += ( -74413                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS_RED'           )
      NAIF_BODY_CODE += ( -74414                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_VIS_NIR'           )
      NAIF_BODY_CODE += ( -74415                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_UV'                )
      NAIF_BODY_CODE += ( -74420                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_UV_SHORT_UV'       )
      NAIF_BODY_CODE += ( -74421                        )

      NAIF_BODY_NAME += ( 'MRO_MARCI_UV_LONG_UV'        )
      NAIF_BODY_CODE += ( -74422                        )
      
      NAIF_BODY_NAME += ( 'MRO_MCS'                     )
      NAIF_BODY_CODE += ( -74500                        )
      
      NAIF_BODY_NAME += ( 'MRO_MCS_A'                   )
      NAIF_BODY_CODE += ( -74510                        )
      
      NAIF_BODY_NAME += ( 'MRO_MCS_A1'                  )
      NAIF_BODY_CODE += ( -74511                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_A2'                  )
      NAIF_BODY_CODE += ( -74512                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_A3'                  )
      NAIF_BODY_CODE += ( -74513                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_A4'                  )
      NAIF_BODY_CODE += ( -74514                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_A5'                  )
      NAIF_BODY_CODE += ( -74515                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_A6'                  )
      NAIF_BODY_CODE += ( -74516                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_B'                   )
      NAIF_BODY_CODE += ( -74520                        )
      
      NAIF_BODY_NAME += ( 'MRO_MCS_B1'                  )
      NAIF_BODY_CODE += ( -74521                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_B2'                  )
      NAIF_BODY_CODE += ( -74522                        )

      NAIF_BODY_NAME += ( 'MRO_MCS_B3'                  )
      NAIF_BODY_CODE += ( -74523                        )

      NAIF_BODY_NAME += ( 'MRO_ONC'                     )
      NAIF_BODY_CODE += ( -74030                        )
      
      NAIF_BODY_NAME += ( 'MRO_SHARAD'                  )
      NAIF_BODY_CODE += ( -74070                        )
      
      NAIF_BODY_NAME += ( 'MRO_HGA_BASEPLATE'           )
      NAIF_BODY_CODE += ( -74211                        )

      NAIF_BODY_NAME += ( 'MRO_HGA_INNER_GIMBAL'        )
      NAIF_BODY_CODE += ( -74212                        )

      NAIF_BODY_NAME += ( 'MRO_HGA_OUTER_GIMBAL'        )
      NAIF_BODY_CODE += ( -74213                        )

      NAIF_BODY_NAME += ( 'MRO_HGA'                     )
      NAIF_BODY_CODE += ( -74214                        )

      NAIF_BODY_NAME += ( 'MRO_LGA1'                    )
      NAIF_BODY_CODE += ( -74220                        )

      NAIF_BODY_NAME += ( 'MRO_LGA2'                    )
      NAIF_BODY_CODE += ( -74230                        )

      NAIF_BODY_NAME += ( 'MRO_UHF'                     )
      NAIF_BODY_CODE += ( -74240                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_BASEPLATE'          )
      NAIF_BODY_CODE += ( -74311                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_INNER_GIMBAL'       )
      NAIF_BODY_CODE += ( -74312                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_OUTER_GIMBAL'       )
      NAIF_BODY_CODE += ( -74313                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX'                    )
      NAIF_BODY_CODE += ( -74314                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_C1'                 )
      NAIF_BODY_CODE += ( -74315                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_C2'                 )
      NAIF_BODY_CODE += ( -74316                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_C3'                 )
      NAIF_BODY_CODE += ( -74317                        )

      NAIF_BODY_NAME += ( 'MRO_SAPX_C4'                 )
      NAIF_BODY_CODE += ( -74318                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_BASEPLATE'          )
      NAIF_BODY_CODE += ( -74321                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_INNER_GIMBAL'       )
      NAIF_BODY_CODE += ( -74322                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_OUTER_GIMBAL'       )
      NAIF_BODY_CODE += ( -74323                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX'                    )
      NAIF_BODY_CODE += ( -74324                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_C1'                 )
      NAIF_BODY_CODE += ( -74325                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_C2'                 )
      NAIF_BODY_CODE += ( -74326                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_C3'                 )
      NAIF_BODY_CODE += ( -74327                        )

      NAIF_BODY_NAME += ( 'MRO_SAMX_C4'                 )
      NAIF_BODY_CODE += ( -74328                        )

   \begintext
