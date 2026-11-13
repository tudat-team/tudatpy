KPL/FK

\beginlabel
PDS_VERSION_ID               = PDS3
RECORD_TYPE                  = STREAM
RECORD_BYTES                 = "N/A"
^SPICE_KERNEL                = "grail_v07.tf"
MISSION_NAME                 = "GRAVITY RECOVERY AND INTERIOR LABORATORY"
SPACECRAFT_NAME              = {
                               "GRAVITY RECOVERY AND INTERIOR LABORATORY A",
                               "GRAVITY RECOVERY AND INTERIOR LABORATORY B"
                               }
DATA_SET_ID                  = "GRAIL-L-SPICE-6-V1.0"
KERNEL_TYPE_ID               = FK
PRODUCT_ID                   = "grail_v07.tf"
PRODUCT_CREATION_TIME        = 2012-11-01T14:02:02
PRODUCER_ID                  = "NAIF/JPL"
MISSION_PHASE_NAME           = {
                               LAUNCH,
                               CRUISE,
                               "ORBIT INSERTION",
                               "TRANSITION TO SCIENCE",
                               SCIENCE,
                               HIATUS,
                               "EXTENDED MISSION",
                               DECOMMISSIONING
                               }
PRODUCT_VERSION_TYPE         = ACTUAL
PLATFORM_OR_MOUNTING_NAME    = "N/A"
START_TIME                   = "N/A"
STOP_TIME                    = "N/A"
SPACECRAFT_CLOCK_START_COUNT = "N/A"
SPACECRAFT_CLOCK_STOP_COUNT  = "N/A"
TARGET_NAME                  = MOON
INSTRUMENT_NAME              = "N/A"
NAIF_INSTRUMENT_ID           = "N/A"
SOURCE_PRODUCT_ID            = "N/A"
NOTE                         = "See comments in the file for details"
OBJECT                       = SPICE_KERNEL
  INTERCHANGE_FORMAT         = ASCII
  KERNEL_TYPE                = FRAMES
  DESCRIPTION                = "GRAIL SPICE FK file providing the complete
set of frame definitions for both GRAIL spacecraft and their science
instruments, created by NAIF, JPL. "
END_OBJECT                   = SPICE_KERNEL
\endlabel


Grail Spacecraft Frame Definition FK
===============================================================================

   This frame kernel contains the frame definitions for the Grail
   spacecraft frames.


Version and Date
-------------------------------------------------------------------------------

   Version 0.7 -- October 23, 2012 -- Boris Semenov, NAIF

      Redefined GRAIL-B MK frames to rotate them 180 degrees about 
      boresights.

   Version 0.6 -- January 24, 2012 -- Boris Semenov, NAIF

      Added name-ID mappings for EBB (synonym for GRAIL-A) and FLOW
      (synonym for GRAIL-B). Corrected comments describing *_SCIENCE
      frames.

   Version 0.5 -- January 11, 2012 -- Boris Semenov, NAIF

      Added star tracker frames and name-ID mappings.

   Version 0.4 -- December 20, 2011 -- Boris Semenov, NAIF

      Corrected MK frame definitions to be based on as-built boresight 
      directions and to have correct twists about boresights.

   Version 0.3 -- October 18, 2011 -- Boris Semenov, NAIF

      Fixed LGA1 definitions to match antenna pattern clock reference
      axes based on input from Martin Schmitzer, GRAIL Telecom, LMCO
      and Chi-Wung Lau, TFP Team, JPL.

   Version 0.2 -- June 6, 2011 -- Boris Semenov, NAIF

      Added payload frames, structure frames, and dynamic frames
      emulating science pointing.

   Version 0.1 -- May 5, 2010 -- Boris Semenov, NAIF

      Changed bogus IDs used in version 0.0 (GRAIL-A/-651, GRAIL-B/-652)
      to actual flight IDs (GRAIL-A/-177, GRAIL-B/-181). Added 
      a few synonyms to the name/ID mapping section.

   Version 0.0 -- June 26, 2007 -- Boris Semenov, NAIF

      
Contact Information
-------------------------------------------------------------------------------

   Boris V. Semenov, NAIF/JPL, (818)-354-8136, Boris.Semenov@jpl.nasa.gov


References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. ``C-Kernel Required Reading''

   4. ``LIB-10: ACS Hardware Coordinate Frame Definitions and
      Transformations'', GRA-AC-09-0013, Rev 3, 7/26/2010

   5. E-mail from Ryan Olds, LMCO; 03/28/11

   6. E-mails from Martin Schmitzer, 09/30/11 and Chi-Wung Lau,
      10/17/11.
  
   7. E-mail from Ralph Roncoli, 10/29/11

   8. E-mail re. as-built MK alignments from Neil Dahya to Ralph
      Roncoli, 12/15/2011


Grail Frames
-------------------------------------------------------------------------------

   The following Grail frames are defined in this kernel file:

           Name                  Relative to           Type       NAIF ID
      ======================  =====================  ============   =======

   Spacecraft frames:
   ------------------
      GRAIL-A_SPACECRAFT      rel.to J2000           CK             -177000

      GRAIL-B_SPACECRAFT      rel.to J2000           CK             -181000

   LGRS Payload frames:
   --------------------
      GRAIL-A_RSB1            rel.to SPACECRAFT      FIXED          -177010
      GRAIL-A_RSB2            rel.to SPACECRAFT      FIXED          -177020
      GRAIL-A_TTS             rel.to SPACECRAFT      FIXED          -177050
      GRAIL-A_KBR             rel.to SPACECRAFT      FIXED          -177060

      GRAIL-B_RSB1            rel.to SPACECRAFT      FIXED          -181010
      GRAIL-B_RSB2            rel.to SPACECRAFT      FIXED          -181020
      GRAIL-B_TTS             rel.to SPACECRAFT      FIXED          -181050
      GRAIL-B_KBR             rel.to SPACECRAFT      FIXED          -181060

   MoonKAM frames:
   ---------------
      GRAIL-A_MK1             rel.to SPACECRAFT      FIXED          -177110
      GRAIL-A_MK2             rel.to SPACECRAFT      FIXED          -177120
      GRAIL-A_MK3             rel.to SPACECRAFT      FIXED          -177130
      GRAIL-A_MK4             rel.to SPACECRAFT      FIXED          -177140

      GRAIL-B_MK1             rel.to SPACECRAFT      FIXED          -181110
      GRAIL-B_MK2             rel.to SPACECRAFT      FIXED          -181120
      GRAIL-B_MK3             rel.to SPACECRAFT      FIXED          -181130
      GRAIL-B_MK4             rel.to SPACECRAFT      FIXED          -181140

   S/C Structure frames:
   ---------------------
      GRAIL-A_LGA1            rel.to SPACECRAFT      FIXED          -177510
      GRAIL-A_LGA2            rel.to SPACECRAFT      FIXED          -177520
      GRAIL-A_STA             rel.to SPACECRAFT      FIXED          -177530

      GRAIL-B_LGA1            rel.to SPACECRAFT      FIXED          -181510
      GRAIL-B_LGA2            rel.to SPACECRAFT      FIXED          -181520
      GRAIL-B_STA             rel.to SPACECRAFT      FIXED          -181530


   Nominal pointing frames:
   -----------------------
      GRAIL-A_TRACK           rel. to J2000          DYNAMIC        -177911
      GRAIL-A_SCIENCE         rel. to GRAIL-A_TRACK  FIXED          -177912 

      GRAIL-B_TRACK           rel. to J2000          DYNAMIC        -181911
      GRAIL-B_SCIENCE         rel. to GRAIL-B_TRACK  FIXED          -181912 


Grail Spacecraft Frames
-------------------------------------------------------------------------------

   The GRAIL spacecraft frame (aka "Mechanical Frames") is defined by the 
   s/c design as follows:

      -  Z is normal to bus plate where star tracker is mounted;
 
      -  X is parallel to and in opposite direction of solar array normal;

      -  Y completes the right handed frame;

      -  the origin of the frame is at the center of the launch
         vehicle interface ring.

   This diagram illustrates the s/c frames:

      +X side view:
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |     
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     o---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |     
         `----------------'|`-'| |   |`----------------'
                           `---| |---'
                               `-' KBR         +Xsc is out of the page.

      +Z side view:
                              +Zsc 
                                o---------> +Ysc
         [----------------][----|----][----------------]
                          .-----|-----.
                          |     |     |
                          |     |     |
                          |     v     |
                          |      +Xsc |
                          |           |
                          |           |
                          |           |        +Zsc is out of the page.
                          `-----------'

   Since the S/C bus attitude with respect to an inertial frame is
   provided by a CK kernel, this frame is defined as a CK-based frame.

   \begindata

      FRAME_GRAIL-A_SPACECRAFT  = -177000
      FRAME_-177000_NAME        = 'GRAIL-A_SPACECRAFT'
      FRAME_-177000_CLASS       = 3
      FRAME_-177000_CLASS_ID    = -177000
      FRAME_-177000_CENTER      = -177
      CK_-177000_SCLK           = -177
      CK_-177000_SPK            = -177

      FRAME_GRAIL-B_SPACECRAFT  = -181000
      FRAME_-181000_NAME        = 'GRAIL-B_SPACECRAFT'
      FRAME_-181000_CLASS       = 3
      FRAME_-181000_CLASS_ID    = -181000
      FRAME_-181000_CENTER      = -181
      CK_-181000_SCLK           = -181
      CK_-181000_SPK            = -181

   \begintext


LGRS Payload Frames
-------------------------------------------------------------------------------

   The frames for all Gravity Recovery and Climate Experiment (LGRS)
   antennas -- RSB, TTS, and KBR -- are defined as follows:

      -  Z is in the direction of the antenna boresight;
 
      -  X is along the antenna clock angle reference direction;

      -  Y completes the right handed frame;

      -  the origins of the frames are at the center of antenna horn outer 
         rim or patch.

   The frames for the Radio Science Beacon (RSB) antennas --
   GRAIL-A_RSB1 and GRAIL-B_RSB1 -- mounted on the -X side of the
   spacecraft bus are nominally rotated from the spacecraft frame by
   -90 degrees about +Y as shown on this diagram:

      +Z side view:
                                  ^ +Zrsb1
                                  |
                                  |
                                  |
                            +Zsc  | +Xrsb1
                                o-o------->-> +Ysc,+Yrsb1 
         [----------------][----|----][----------------]
                          .-----|-----.
                          |     |     |
                          |     |     |
                          |     v     |
                          |      +Xsc |
                          |           |
                          |           |        +Xrsb1 is out of the page.
                          |           |        +Zsc is out of the page.
                          `-----------'

   (The frame definitions below contain the opposite of this rotation
   because Euler angles specified in them define transformations from
   the antenna frames to the spacecraft frame -- see [1].)

   \begindata
 
      FRAME_GRAIL-A_RSB1        = -177010
      FRAME_-177010_NAME        = 'GRAIL-A_RSB1'
      FRAME_-177010_CLASS       = 4
      FRAME_-177010_CLASS_ID    = -177010
      FRAME_-177010_CENTER      = -177
      TKFRAME_-177010_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177010_SPEC      = 'ANGLES'
      TKFRAME_-177010_UNITS     = 'DEGREES'
      TKFRAME_-177010_ANGLES    = (   0.0,  90.0,   0.0 )
      TKFRAME_-177010_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-B_RSB1        = -181010
      FRAME_-181010_NAME        = 'GRAIL-B_RSB1'
      FRAME_-181010_CLASS       = 4
      FRAME_-181010_CLASS_ID    = -181010
      FRAME_-181010_CENTER      = -181
      TKFRAME_-181010_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181010_SPEC      = 'ANGLES'
      TKFRAME_-181010_UNITS     = 'DEGREES'
      TKFRAME_-181010_ANGLES    = (   0.0,  90.0,   0.0 )
      TKFRAME_-181010_AXES      = (   1,     2,     3   )

   \begintext

   The frames for the Radio Science Beacon (RSB) antennas --
   GRAIL-A_RSB2 and GRAIL-B_RSB2 -- mounted on the +X side of the
   spacecraft bus are nominally rotated from the spacecraft frame by
   +90 degrees about +Y as shown on this diagram:

      +Z side view:
                                  +Zsc
                                o---------> +Ysc
         [----------------][----|----][----------------]
                          .-----|-----.
                          |     |     |
                          |     |     |
                          |     v     |  
                          |      +Xsc |
                          |           |
                          |           |        +Xrsb2 is into the page.
                          |           |        +Zsc is out of the page.
                          `-----------'
                           +Xrsb2 x---------> +Yrsb2
                                  |
                                  |
                                  |
                                  |
                                  v +Zrsb2

   (The frame definitions below contain the opposite of this rotation
   because Euler angles specified in them define transformations from
   the antenna frames to the spacecraft frame -- see [1].)

   \begindata

      FRAME_GRAIL-A_RSB2        = -177020
      FRAME_-177020_NAME        = 'GRAIL-A_RSB2'
      FRAME_-177020_CLASS       = 4
      FRAME_-177020_CLASS_ID    = -177020
      FRAME_-177020_CENTER      = -177
      TKFRAME_-177020_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177020_SPEC      = 'ANGLES'
      TKFRAME_-177020_UNITS     = 'DEGREES'
      TKFRAME_-177020_ANGLES    = (   0.0, -90.0,   0.0 )
      TKFRAME_-177020_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-B_RSB2        = -181020
      FRAME_-181020_NAME        = 'GRAIL-B_RSB2'
      FRAME_-181020_CLASS       = 4
      FRAME_-181020_CLASS_ID    = -181020
      FRAME_-181020_CENTER      = -181
      TKFRAME_-181020_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181020_SPEC      = 'ANGLES'
      TKFRAME_-181020_UNITS     = 'DEGREES'
      TKFRAME_-181020_ANGLES    = (   0.0, -90.0,   0.0 )
      TKFRAME_-181020_AXES      = (   1,     2,     3   )

   \begintext

   The frames for the Time Transfer System (TTS) antennas --
   GRAIL-A_TTS and GRAIL-B_TTS -- mounted on the -Z side of the
   spacecraft bus are nominally rotated from the spacecraft frame by
   180 degrees about +Y as shown on this diagram:

      +X side view:
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |     
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     o---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |     
         `----------------'|`x---------> +Ytts  -------'
                       +Xtts-|-| |---'
                             | `-' KBR          +Xsc is out of the page.
                             |                  +Xtts is into the page.
                             |
                             v +Ztts
       
   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in them define transformations from the
   antenna frames to the spacecraft frame -- see [1].)

   \begindata
 
      FRAME_GRAIL-A_TTS         = -177050
      FRAME_-177050_NAME        = 'GRAIL-A_TTS'
      FRAME_-177050_CLASS       = 4
      FRAME_-177050_CLASS_ID    = -177050
      FRAME_-177050_CENTER      = -177
      TKFRAME_-177050_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177050_SPEC      = 'ANGLES'
      TKFRAME_-177050_UNITS     = 'DEGREES'
      TKFRAME_-177050_ANGLES    = (   0.0, 180.0,   0.0 )
      TKFRAME_-177050_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-B_TTS         = -181050
      FRAME_-181050_NAME        = 'GRAIL-B_TTS'
      FRAME_-181050_CLASS       = 4
      FRAME_-181050_CLASS_ID    = -181050
      FRAME_-181050_CENTER      = -181
      TKFRAME_-181050_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181050_SPEC      = 'ANGLES'
      TKFRAME_-181050_UNITS     = 'DEGREES'
      TKFRAME_-181050_ANGLES    = (   0.0, 180.0,   0.0 )
      TKFRAME_-181050_AXES      = (   1,     2,     3   )
 
   \begintext

   The frames for the Ka-Band Ranging (KBR) antennas -- GRAIL-A_KBR and
   GRAIL-B_KBR -- mounted on the -Z side of the spacecraft bus are
   nominally rotated from the spacecraft frame by 180 degrees about +Y
   and then canted by 2.1 deg towards the -Ysc axis for GRAIL-A and 2.1
   deg towards the +Ysc for GRAIL-B as shown on these diagrams:

      GRAIL-A +X side view:
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |     
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     o---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |     
         `----------------'|`-'| |KBR|`----------------'
                           `---| |---'
                         +Xkbr `x---------> +Ykbr
                                |
                               '|                +Xsc is out of the page.
                              . |                +Xkbr is into the page.
                                |
                       +Zkbr V  |
                              
      GRAIL-B +X side view:
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |     
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     *---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |     
         `----------------'|`-'| |KBR|`----------------'
                           `---| |---'
                         +Xkbr `x---------> +Ykbr
                                |
                                |'               +Xsc is out of the page.
                                | .              +Xkbr is into the page.
                                |  
                                |  v +Zkbr                              

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in them define transformations from the
   antenna frames to the spacecraft frame -- see [1].)

   \begindata

      FRAME_GRAIL-A_KBR         = -177060
      FRAME_-177060_NAME        = 'GRAIL-A_KBR'
      FRAME_-177060_CLASS       = 4
      FRAME_-177060_CLASS_ID    = -177060
      FRAME_-177060_CENTER      = -177
      TKFRAME_-177060_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177060_SPEC      = 'ANGLES'
      TKFRAME_-177060_UNITS     = 'DEGREES'
      TKFRAME_-177060_ANGLES    = (   0.0, 180.0,  -2.1 )
      TKFRAME_-177060_AXES      = (   3,     2,     1   )
 
      FRAME_GRAIL-B_KBR         = -181060
      FRAME_-181060_NAME        = 'GRAIL-B_KBR'
      FRAME_-181060_CLASS       = 4
      FRAME_-181060_CLASS_ID    = -181060
      FRAME_-181060_CENTER      = -181
      TKFRAME_-181060_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181060_SPEC      = 'ANGLES'
      TKFRAME_-181060_UNITS     = 'DEGREES'
      TKFRAME_-181060_ANGLES    = (   0.0, 180.0,   2.1 )
      TKFRAME_-181060_AXES      = (   3,     2,     1   ) 

   \begintext


MoonKAM Frames
-------------------------------------------------------------------------------

   All MoonKAM (MK) frames are defined as follows:

      -  Z is in the direction of the camera boresight boresight;
 
      -  X is along the longer (horizontal) dimension of the image, 
         pointing from the center to the left of the image.

      -  Y is along the shorter (vertical) dimension of the image, 
         pointing from the center to the top of the image.

      -  the origins of the frames are at the camera focal points.

   These nominal MK view directions in the spacecraft frame are given
   in [5]:

      Camera       S/c          X       Y       Z
      -----------  ---------  ------  ------  ------
      MK_Camera_1  (GRAIL-A)   0.0    -1.0     0.0
      MK_Camera_2  (GRAIL-A)   0.13   -0.491   0.856
      MK_Camera_3  (GRAIL-A)   0.13   -0.491  -0.856
      MK_Camera_4  (GRAIL-A)   0.0    -1.0     0.0

      MK_Camera_1  (GRAIL-B)   0.0     1.0     0.0
      MK_Camera_2  (GRAIL-B)   0.13    0.491   0.856
      MK_Camera_3  (GRAIL-B)   0.13    0.491  -0.856
      MK_Camera_4  (GRAIL-B)   0.0     1.0     0.0

   For MK2 and MK3 these directions are based on the incorrect
   +Y/-Y-to-boresight and ZY-plane-to-boresight angles of 60.16149823
   and 7.50470626 degrees. The frames definitions based on these angles
   were used in the FK versions 0.2-0.3.

   The correct as-built +Y/-Y-to-boresight and ZY-plane-to-boresight
   angles are 60.0 and 8.0 degrees ([8]). Given these as-built angles
   and the fact that the X axes of all camera frames are in the s/c YZ
   plane (see [7]) and point either exactly (MK1 & MK4) or
   approximately (MK2 & MK3) in the s/c -Z direction for GRAIL-A or +Z
   direction for GRAIL-B, the following three rotations are needed to
   co-align the spacecraft frame with the camera frames:
   
      Camera       S/c        Rot1 about X   Rot2 about Z  Rot3 about X
      -----------  ---------  ------------   ------------  ------------
      MK_Camera_1  (GRAIL-A)      90.0          -90.0           0.0
      MK_Camera_2  (GRAIL-A)      30.0          -90.0          -8.0
      MK_Camera_3  (GRAIL-A)     150.0          -90.0          -8.0
      MK_Camera_4  (GRAIL-A)      90.0          -90.0           0.0
      MK_Camera_1  (GRAIL-B)     -90.0          -90.0           0.0
      MK_Camera_2  (GRAIL-B)     -30.0          -90.0          -8.0
      MK_Camera_3  (GRAIL-B)    -150.0          -90.0          -8.0
      MK_Camera_4  (GRAIL-B)     -90.0          -90.0           0.0

   The following diagrams illustrate the MK frames:

      GRAIL-A +X side view:

              +Zmk2 ^
                     \
                      \    .---------.
         .-------------\--.|     +Zsc|.----------------.
         |          MK2 \ ||    ^    ||                |     
         |              .o|-----|-----.                |
         |     +Xmk2 .-'  |     |     |                |
         |         <'  MK1|     |     |                |
         |     <---------o|     |     |                |
   <---- |  +Zmk1     MK4||     o---------> +Ysc       |
   Nadir |     <---------o|      +Xsc |                |
         |  +Zmk4        ||           |                |
         |               v +Xmk1,4    |                |
         |           MK3 o.   --------'                |
         |              /  `-.                         |     
         `-------------/--'|  `> +Xmk3 ----------------'
                      /   TTS -   --- 
                     /         `-'
              +Zmk3 v        KBR           +Xsc  is out of the page.
                                           +Ymk1 is out of the page.
                                           +Ymk2 is out of the page and
                                                 tilted 8 deg towards -Zsc.
                                           +Ymk3 is out of the page and
                                                 tilted 8 deg towards +Zsc.
                                           +Ymk4 is out of the page.
                              

      GRAIL-B +X side view:

                                            ^  +Zmk2
                                           /
                           .---- +Xmk2    /
         .----------------.|+Zsc <.   .--/-------------.
         |                ||    ^  `-.| / MK2          |     
         |                .-----|-----`o               |
         |                |     |     |^ +Xmk1,4       |
         |                |     |     ||MK4      +Zmk4 |
         |                |     |     |o--------->     |
         |                |     o---------> +Ysc +Zmk1 |   -----> Nadir
         |                |   +Xsc    |o--------->     |
         |                |           | MK1 .>         |
         |                |           |  .-'  +Xmk3    |
         |                '-------     o'              |
         |                ||| || |  MK3 \              |     
         `----------------'|`-'| |     --\-------------'
                         TTS --| |        \
                               `-'         \
                            KBR             v +Zmk3

                                           +Xsc  is out of the page.
                                           +Ymk1 is out of the page.
                                           +Ymk2 is out of the page and
                                                 tilted 8 deg towards -Zsc.
                                           +Ymk3 is out of the page and
                                                 tilted 8 deg towards +Zsc.
                                           +Ymk4 is out of the page.

   (The frame definitions below contain the opposites of these
   rotations because Euler angles specified in them define
   transformations from the camera frames to the spacecraft frame --
   see [1].)

   \begindata

      FRAME_GRAIL-A_MK1         = -177110
      FRAME_-177110_NAME        = 'GRAIL-A_MK1'
      FRAME_-177110_CLASS       = 4
      FRAME_-177110_CLASS_ID    = -177110
      FRAME_-177110_CENTER      = -177
      TKFRAME_-177110_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177110_SPEC      = 'ANGLES'
      TKFRAME_-177110_UNITS     = 'DEGREES'
      TKFRAME_-177110_ANGLES    = (  -90.0,  90.0,   0.0 )
      TKFRAME_-177110_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-A_MK2         = -177120
      FRAME_-177120_NAME        = 'GRAIL-A_MK2'
      FRAME_-177120_CLASS       = 4
      FRAME_-177120_CLASS_ID    = -177120
      FRAME_-177120_CENTER      = -177
      TKFRAME_-177120_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177120_SPEC      = 'ANGLES'
      TKFRAME_-177120_UNITS     = 'DEGREES'
      TKFRAME_-177120_ANGLES    = ( -30.0    90.0   8.0 )
      TKFRAME_-177120_AXES      = (    1,     3,    1   )
 
      FRAME_GRAIL-A_MK3         = -177130
      FRAME_-177130_NAME        = 'GRAIL-A_MK3'
      FRAME_-177130_CLASS       = 4
      FRAME_-177130_CLASS_ID    = -177130
      FRAME_-177130_CENTER      = -177
      TKFRAME_-177130_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177130_SPEC      = 'ANGLES'
      TKFRAME_-177130_UNITS     = 'DEGREES'
      TKFRAME_-177130_ANGLES    = ( -150.0   90.0    8.0 )
      TKFRAME_-177130_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-A_MK4         = -177140
      FRAME_-177140_NAME        = 'GRAIL-A_MK4'
      FRAME_-177140_CLASS       = 4
      FRAME_-177140_CLASS_ID    = -177140
      FRAME_-177140_CENTER      = -177
      TKFRAME_-177140_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177140_SPEC      = 'ANGLES'
      TKFRAME_-177140_UNITS     = 'DEGREES'
      TKFRAME_-177140_ANGLES    = ( -90.0,   90.0,   0.0 )
      TKFRAME_-177140_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-B_MK1         = -181110
      FRAME_-181110_NAME        = 'GRAIL-B_MK1'
      FRAME_-181110_CLASS       = 4
      FRAME_-181110_CLASS_ID    = -181110
      FRAME_-181110_CENTER      = -181
      TKFRAME_-181110_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181110_SPEC      = 'ANGLES'
      TKFRAME_-181110_UNITS     = 'DEGREES'
      TKFRAME_-181110_ANGLES    = (  90.0    90.0    0.0 )
      TKFRAME_-181110_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-B_MK2         = -181120
      FRAME_-181120_NAME        = 'GRAIL-B_MK2'
      FRAME_-181120_CLASS       = 4
      FRAME_-181120_CLASS_ID    = -181120
      FRAME_-181120_CENTER      = -181
      TKFRAME_-181120_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181120_SPEC      = 'ANGLES'
      TKFRAME_-181120_UNITS     = 'DEGREES'
      TKFRAME_-181120_ANGLES    = (  30.0    90.0    8.0 )
      TKFRAME_-181120_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-B_MK3         = -181130
      FRAME_-181130_NAME        = 'GRAIL-B_MK3'
      FRAME_-181130_CLASS       = 4
      FRAME_-181130_CLASS_ID    = -181130
      FRAME_-181130_CENTER      = -181
      TKFRAME_-181130_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181130_SPEC      = 'ANGLES'
      TKFRAME_-181130_UNITS     = 'DEGREES'
      TKFRAME_-181130_ANGLES    = ( 150.0    90.0    8.0 )
      TKFRAME_-181130_AXES      = (    1,     3,     1   )
 
      FRAME_GRAIL-B_MK4         = -181140
      FRAME_-181140_NAME        = 'GRAIL-B_MK4'
      FRAME_-181140_CLASS       = 4
      FRAME_-181140_CLASS_ID    = -181140
      FRAME_-181140_CENTER      = -181
      TKFRAME_-181140_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181140_SPEC      = 'ANGLES'
      TKFRAME_-181140_UNITS     = 'DEGREES'
      TKFRAME_-181140_ANGLES    = (  90.0    90.0    0.0 )
      TKFRAME_-181140_AXES      = (    1,     3,     1   )
 
   \begintext


S/C Structure Frames
-------------------------------------------------------------------------------

LGA Frames

   The low gain antenna (LGA) frames -- GRAIL-A_LGA1, GRAIL-A_LGA2,
   GRAIL-B_LGA1, and GRAIL-B_LGA2 -- are defined as follows:

      -  Z is in the direction of the antenna boresight;
 
      -  X is along the antenna clock angle reference direction;

      -  Y completes the right handed frame;

      -  the origins of the frames are at the center of antenna patch.

   The frame for the LGA2 mounted on the +X side of the spacecraft bus
   is nominally rotated from the spacecraft frame by +90 degrees about
   +Y while the frame for the LGA1 mounted on the -X side of the
   spacecraft bus is nominally rotated from the spacecraft frame by
   first -90 degrees about X, then by -90 degrees about +Y as shown on
   this diagram:

      +Z side view:
                                  ^ +Zlga1
                                  |
                                  |
                                  |
                            +Zsc  | +Ylga1
                                o-x------->-> +Ysc, +Xlga1
         [----------------][----|----][----------------]
                          .-----|-----.
                          |     |     |
                          |     |     |
                          |     v     |  
                          |      +Xsc |
                          |           |        +Ylga1 is into the page.
                          |           |        +Xlga2 is into the page.
                          |           |        +Zsc is out of the page.
                          `-----------'
                       +Xlga2 x---------> +Ylga2
                              |
                              |
                              |
                              |
                              v +Zlga2

   (The frame definitions below contain the opposite of this rotation
   because Euler angles specified in them define transformations from
   the antenna frames to the spacecraft frame -- see [1].)

   LGA frames.

   \begindata

      FRAME_GRAIL-A_LGA2        = -177520
      FRAME_-177520_NAME        = 'GRAIL-A_LGA2'
      FRAME_-177520_CLASS       = 4
      FRAME_-177520_CLASS_ID    = -177520
      FRAME_-177520_CENTER      = -177
      TKFRAME_-177520_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177520_SPEC      = 'ANGLES'
      TKFRAME_-177520_UNITS     = 'DEGREES'
      TKFRAME_-177520_ANGLES    = (   0.0, -90.0,   0.0 )
      TKFRAME_-177520_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-B_LGA2        = -181520
      FRAME_-181520_NAME        = 'GRAIL-B_LGA2'
      FRAME_-181520_CLASS       = 4
      FRAME_-181520_CLASS_ID    = -181520
      FRAME_-181520_CENTER      = -181
      TKFRAME_-181520_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181520_SPEC      = 'ANGLES'
      TKFRAME_-181520_UNITS     = 'DEGREES'
      TKFRAME_-181520_ANGLES    = (   0.0, -90.0,   0.0 )
      TKFRAME_-181520_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-A_LGA1        = -177510
      FRAME_-177510_NAME        = 'GRAIL-A_LGA1'
      FRAME_-177510_CLASS       = 4
      FRAME_-177510_CLASS_ID    = -177510
      FRAME_-177510_CENTER      = -177
      TKFRAME_-177510_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177510_SPEC      = 'ANGLES'
      TKFRAME_-177510_UNITS     = 'DEGREES'
      TKFRAME_-177510_ANGLES    = (   90.0,  90.0,   0.0 )
      TKFRAME_-177510_AXES      = (   1,     2,     3   )
 
      FRAME_GRAIL-B_LGA1        = -181510
      FRAME_-181510_NAME        = 'GRAIL-B_LGA1'
      FRAME_-181510_CLASS       = 4
      FRAME_-181510_CLASS_ID    = -181510
      FRAME_-181510_CENTER      = -181
      TKFRAME_-181510_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181510_SPEC      = 'ANGLES'
      TKFRAME_-181510_UNITS     = 'DEGREES'
      TKFRAME_-181510_ANGLES    = (   90.0,  90.0,   0.0 )
      TKFRAME_-181510_AXES      = (   1,     2,     3   )
 
   \begintext


STA Frames

   The star tracker assembly (STA) frames -- GRAIL-A_STA and
   GRAIL-B_STA -- are defined as follows:

      -  Z is in the direction of the star tracker boresight; for
         GRAIL-A it is nominally canted by 30 degrees off s/c +Z towards
         s/c +Y, for GRAIL-B it is nominally canted by 30 degrees off
         s/c +Z towards s/c -Y
 
      -  X is towards the STA electronics; for GRAIL-A it nominally
         points in the direction opposite of the s/c +X, for GRAIL-B it
         nominally points in the same direction as s/c +X

      -  Y completes the right handed frame;

      -  the origin of the frame are at the STA focal point.

   The frame for the GRAIL-A STA, mounted on the +Z side of the
   spacecraft bus, is nominally rotated from the spacecraft frame by
   180 degrees about +Z, then by +30 degrees about X, while the frame
   the GRAIL-B STA, also mounted on the +Z side of the spacecraft bus,
   is nominally rotated from the spacecraft frame by +30 degrees X:

      GRAIL-A +X side view:

                                    ^ +Zsta
                    +Ysta          /
                        <-.       /
                           `-.   /
                              `-x +Xsta
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     o---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |
         `----------------'|`-'| |   |`----------------'
                           `---| |---'
                               `-' KBR         +Xsc is out of the page.
                                               +Xsta is into of the page.
                              
      GRAIL-B +X side view:

                            ^ +Zsta
                             \            +Ysta
                              \       .->
                               \   .-'
                         +Xsta  o-'
                           .---------.
         .----------------.|     +Zsc|.----------------.
         |                ||    ^    ||                |     
         |                .-----|-----.                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     |     |                |
         |                |     o---------> +Ysc       |
         |                |   +Xsc    |                |
         |                |           |                |
         |                |           |                |
         |                '-----------'                |
         |               TTS| || |   ||                |     
         `----------------'|`-'| |   |`----------------'
                           `---| |---'
                               `-' KBR         +Xsc is out of the page.
                                               +Xsta is out of the page.

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in them define transformations from the
   antenna frames to the spacecraft frame -- see [1].)

   \begindata

      FRAME_GRAIL-A_STA         = -177530
      FRAME_-177530_NAME        = 'GRAIL-A_STA'
      FRAME_-177530_CLASS       = 4
      FRAME_-177530_CLASS_ID    = -177530
      FRAME_-177530_CENTER      = -177
      TKFRAME_-177530_RELATIVE  = 'GRAIL-A_SPACECRAFT'
      TKFRAME_-177530_SPEC      = 'ANGLES'
      TKFRAME_-177530_UNITS     = 'DEGREES'
      TKFRAME_-177530_ANGLES    = (   0.0, 180.0, -30.0 )
      TKFRAME_-177530_AXES      = (   1,     3,     1   )
 
      FRAME_GRAIL-B_STA         = -181530
      FRAME_-181530_NAME        = 'GRAIL-B_STA'
      FRAME_-181530_CLASS       = 4
      FRAME_-181530_CLASS_ID    = -181530
      FRAME_-181530_CENTER      = -181
      TKFRAME_-181530_RELATIVE  = 'GRAIL-B_SPACECRAFT'
      TKFRAME_-181530_SPEC      = 'ANGLES'
      TKFRAME_-181530_UNITS     = 'DEGREES'
      TKFRAME_-181530_ANGLES    = (   0.0,   0.0, -30.0 )
      TKFRAME_-181530_AXES      = (   1,     3,     1   )
   
   \begintext


Nominal Pointing Frames
-------------------------------------------------------------------------------

   The GRAIL nominal pointing frames -- GRAIL-A_SCIENCE (ID -177912) and
   GRAIL-B_SCIENCE (ID -181912) -- are intended to emulate the nominal
   spacecraft orientations during the mapping phase of the mission.
   These frames are defined as follows:

      GRAIL-A_SCIENCE:  KA-band antenna boresight toward GRAIL-B
                        -Y towards the Moon center

      GRAIL-B_SCIENCE:  KA-band antenna boresight toward GRAIL-A, 
                        +Y towards the Moon center

   These frames are defined as fixed-offset frames rotated by -2.1 and
   2.1 degrees about X correspondingly relative to the tracking dynamic
   frames defined later in this FK.

   \begindata

      FRAME_GRAIL-A_SCIENCE     = -177912 
      FRAME_-177912_NAME        = 'GRAIL-A_SCIENCE'
      FRAME_-177912_CLASS       = 4
      FRAME_-177912_CLASS_ID    = -177912
      FRAME_-177912_CENTER      = -177
      TKFRAME_-177912_RELATIVE  = 'GRAIL-A_TRACK'
      TKFRAME_-177912_SPEC      = 'ANGLES'
      TKFRAME_-177912_UNITS     = 'DEGREES'
      TKFRAME_-177912_ANGLES    = (  -2.1,  0.0,   0.0 )
      TKFRAME_-177912_AXES      = (   1,    2,     3   )

      FRAME_GRAIL-B_SCIENCE     = -181912 
      FRAME_-181912_NAME        = 'GRAIL-B_SCIENCE'
      FRAME_-181912_CLASS       = 4
      FRAME_-181912_CLASS_ID    = -181912
      FRAME_-181912_CENTER      = -181
      TKFRAME_-181912_RELATIVE  = 'GRAIL-B_TRACK'
      TKFRAME_-181912_SPEC      = 'ANGLES'
      TKFRAME_-181912_UNITS     = 'DEGREES'
      TKFRAME_-181912_ANGLES    = (   2.1,  0.0,   0.0 )
      TKFRAME_-181912_AXES      = (   1,    2,     3   )

   \begintext

   The GRAIL tracking dynamic frames -- GRAIL-A_TRACK (ID -177911) and
   GRAIL-B_TRACK (ID -181911)-- are axillary frames needed to
   implement the nominal dynamic frames defined above. These frames are
   defined as follows:

      GRAIL-A_TRACK:  -Z toward GRAIL-B
                      -Y towards the Moon center

      GRAIL-B_TRACK:  -Z toward GRAIL-A
                      +Y towards the Moon center

   \begindata


      FRAME_GRAIL-A_TRACK          = -177911
      FRAME_-177911_NAME           = 'GRAIL-A_TRACK'
      FRAME_-177911_CLASS          = 5
      FRAME_-177911_CLASS_ID       = -177911
      FRAME_-177911_CENTER         = -177
      FRAME_-177911_RELATIVE       = 'J2000'
      FRAME_-177911_DEF_STYLE      = 'PARAMETERIZED'
      FRAME_-177911_FAMILY         = 'TWO-VECTOR'
      FRAME_-177911_PRI_AXIS       = '-Z'
      FRAME_-177911_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
      FRAME_-177911_PRI_OBSERVER   = 'GRAIL-A'
      FRAME_-177911_PRI_TARGET     = 'GRAIL-B'
      FRAME_-177911_PRI_ABCORR     = 'NONE'
      FRAME_-177911_SEC_AXIS       = '-Y'
      FRAME_-177911_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
      FRAME_-177911_SEC_OBSERVER   = 'GRAIL-A'
      FRAME_-177911_SEC_TARGET     = 'MOON'
      FRAME_-177911_SEC_ABCORR     = 'NONE'

      FRAME_GRAIL-B_TRACK          = -181911
      FRAME_-181911_NAME           = 'GRAIL-B_TRACK'
      FRAME_-181911_CLASS          = 5
      FRAME_-181911_CLASS_ID       = -181911
      FRAME_-181911_CENTER         = -181
      FRAME_-181911_RELATIVE       = 'J2000'
      FRAME_-181911_DEF_STYLE      = 'PARAMETERIZED'
      FRAME_-181911_FAMILY         = 'TWO-VECTOR'
      FRAME_-181911_PRI_AXIS       = '-Z'
      FRAME_-181911_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
      FRAME_-181911_PRI_OBSERVER   = 'GRAIL-B'
      FRAME_-181911_PRI_TARGET     = 'GRAIL-A'
      FRAME_-181911_PRI_ABCORR     = 'NONE'
      FRAME_-181911_SEC_AXIS       = '+Y'
      FRAME_-181911_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
      FRAME_-181911_SEC_OBSERVER   = 'GRAIL-B'
      FRAME_-181911_SEC_TARGET     = 'MOON'
      FRAME_-181911_SEC_ABCORR     = 'NONE'

   \begintext


GRAIL Name-ID Mappings
-------------------------------------------------------------------------------
   
   This section contains name to NAIF ID mappings for the GRAIL mission.
   Once the contents of this file is loaded into the KERNEL POOL, these
   mappings become available within SPICE, making it possible to use
   names instead of ID code in the high level SPICE routine calls.

   GRAIL-A name-ID mappings:

   \begindata

      NAIF_BODY_NAME += ( 'EBB'                     )
      NAIF_BODY_CODE += ( -177                      )

      NAIF_BODY_NAME += ( 'GRA'                     )
      NAIF_BODY_CODE += ( -177                      )

      NAIF_BODY_NAME += ( 'GRAIL_A'                 )
      NAIF_BODY_CODE += ( -177                      )

      NAIF_BODY_NAME += ( 'GRAIL-A'                 )
      NAIF_BODY_CODE += ( -177                      )

      NAIF_BODY_NAME += ( 'GRAIL-A_RSB1'            )
      NAIF_BODY_CODE += ( -177010                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_RSB2'            )
      NAIF_BODY_CODE += ( -177020                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_TTS'             )
      NAIF_BODY_CODE += ( -177050                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_KBR'             )
      NAIF_BODY_CODE += ( -177060                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_MK1'             )
      NAIF_BODY_CODE += ( -177110                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_MK2'             )
      NAIF_BODY_CODE += ( -177120                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_MK3'             )
      NAIF_BODY_CODE += ( -177130                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_MK4'             )
      NAIF_BODY_CODE += ( -177140                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_LGA1'            )
      NAIF_BODY_CODE += ( -177510                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_LGA2'            )
      NAIF_BODY_CODE += ( -177520                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_STA'             )
      NAIF_BODY_CODE += ( -177530                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_STA_SUN_ZONE'    )
      NAIF_BODY_CODE += ( -177531                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_STA_EARTH_ZONE'  )
      NAIF_BODY_CODE += ( -177532                   )

      NAIF_BODY_NAME += ( 'GRAIL-A_STA_MOON_ZONE'   )
      NAIF_BODY_CODE += ( -177533                   )

   \begintext

   GRAIL-B name-ID mappings:

   \begindata

      NAIF_BODY_NAME += ( 'FLOW'                    )
      NAIF_BODY_CODE += ( -181                      )

      NAIF_BODY_NAME += ( 'GRB'                     )
      NAIF_BODY_CODE += ( -181                      )

      NAIF_BODY_NAME += ( 'GRAIL_B'                 )
      NAIF_BODY_CODE += ( -181                      )

      NAIF_BODY_NAME += ( 'GRAIL-B'                 )
      NAIF_BODY_CODE += ( -181                      )

      NAIF_BODY_NAME += ( 'GRAIL-B_RSB1'            )
      NAIF_BODY_CODE += ( -181010                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_RSB2'            )
      NAIF_BODY_CODE += ( -181020                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_TTS'             )
      NAIF_BODY_CODE += ( -181050                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_KBR'             )
      NAIF_BODY_CODE += ( -181060                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_MK1'             )
      NAIF_BODY_CODE += ( -181110                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_MK2'             )
      NAIF_BODY_CODE += ( -181120                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_MK3'             )
      NAIF_BODY_CODE += ( -181130                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_MK4'             )
      NAIF_BODY_CODE += ( -181140                   )
  
      NAIF_BODY_NAME += ( 'GRAIL-B_LGA1'            )
      NAIF_BODY_CODE += ( -181510                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_LGA2'            )
      NAIF_BODY_CODE += ( -181520                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_STA'             )
      NAIF_BODY_CODE += ( -181530                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_STA_SUN_ZONE'    )
      NAIF_BODY_CODE += ( -181531                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_STA_EARTH_ZONE'  )
      NAIF_BODY_CODE += ( -181532                   )

      NAIF_BODY_NAME += ( 'GRAIL-B_STA_MOON_ZONE'   )
      NAIF_BODY_CODE += ( -181533                   )

   \begintext

End of FK file.
