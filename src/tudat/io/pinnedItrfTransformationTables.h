// clang-format off
// Pinned official IGN input tables; full source notices retained. Parsed only in data_input.
// SHA256 provenance accompanies readOfficialItrf20xxTransformationData.
namespace tudat { namespace input_output {
const char pinnedItrf2014Table[] = R"ITRF(Transformation parameters from ITRF2014 to past ITRFs.
----------------------------------------------------------------------------------------
SOLUTION         Tx       Ty       Tz        D        Rx        Ry        Rz      EPOCH
UNITS----------> mm       mm       mm       ppb       .001"     .001"     .001"
                 .        .        .         .        .         .         .
       RATES     Tx       Ty       Tz        D        Rx        Ry        Rz
UNITS----------> mm/y     mm/y     mm/y     ppb/y    .001"/y   .001"/y   .001"/y
-----------------------------------------------------------------------------------------
  ITRF2008        1.6      1.9      2.4     -0.02      0.00      0.00      0.00    2010.0
       rates      0.0      0.0     -0.1      0.03      0.00      0.00      0.00
  ITRF2005        2.6      1.0     -2.3      0.92      0.00      0.00      0.00    2010.0
       rates      0.3      0.0     -0.1      0.03      0.00      0.00      0.00
  ITRF2000        0.7      1.2    -26.1      2.12      0.00      0.00      0.00    2010.0
       rates      0.1      0.1     -1.9      0.11      0.00      0.00      0.00
  ITRF97          7.4     -0.5    -62.8      3.80      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF96          7.4     -0.5    -62.8      3.80      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF94          7.4     -0.5    -62.8      3.80      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF93        -50.4      3.3    -60.2      4.29     -2.81     -3.38      0.40    2010.0
       rates     -2.8     -0.1     -2.5      0.12     -0.11     -0.19      0.07
  ITRF92         15.4      1.5    -70.8      3.09      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF91         27.4     15.5    -76.8      4.49      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF90         25.4     11.5    -92.8      4.79      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF89         30.4     35.5   -130.8      8.19      0.00      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
  ITRF88         25.4     -0.5   -154.8     11.29      0.10      0.00      0.26    2010.0
       rates      0.1     -0.5     -3.3      0.12      0.00      0.00      0.02
_________________________________________________________________________________________

Note : These parameters are derived from those already published in the IERS
Technical Notes and Annual Reports. The transformation parameters should be
used with the standard model (1) given below and are valid at the indicated
epoch.


: XS :    : X :   : Tx :   :  D   -Rz   Ry : : X :
:    :    :   :   :    :   :               : :   :
: YS :  = : Y : + : Ty : + :  Rz   D   -Rx : : Y :                       (1)
:    :    :   :   :    :   :               : :   :
: ZS :    : Z :   : Tz :   : -Ry   Rx   D  : : Z :


Where X,Y,Z are the coordinates in ITRF2014 and XS,YS,ZS are the coordinates in
the other frames.

On the other hand, for a given parameter P, its value at any epoch t
is obtained by using equation (2).

                  .
P(t) = P(EPOCH) + P * (t - EPOCH)                                        (2)


where EPOCH is the epoch indicated in the above table (currently 2010.0)
    .
and P is the rate of that parameter.

)ITRF";
const char pinnedItrf2020Table[] = R"ITRF(Transformation parameters from ITRF2020 to past ITRFs.
----------------------------------------------------------------------------------------
SOLUTION         Tx       Ty       Tz        D        Rx        Ry        Rz      EPOCH
UNITS----------> mm       mm       mm       ppb       .001"     .001"     .001"
                 .        .        .         .        .         .         .
       RATES     Tx       Ty       Tz        D        Rx        Ry        Rz
UNITS----------> mm/y     mm/y     mm/y     ppb/y    .001"/y   .001"/y   .001"/y
-----------------------------------------------------------------------------------------
  ITRF2014       -1.4     -0.9      1.4     -0.42      0.00      0.00      0.00    2015.0  
       rates      0.0     -0.1      0.2      0.00      0.00      0.00      0.00     
  ITRF2008        0.2      1.0      3.3     -0.29      0.00      0.00      0.00    2015.0  
       rates      0.0     -0.1      0.1      0.03      0.00      0.00      0.00     
  ITRF2005        2.7      0.1     -1.4      0.65      0.00      0.00      0.00    2015.0  
       rates      0.3     -0.1      0.1      0.03      0.00      0.00      0.00     
  ITRF2000       -0.2      0.8    -34.2      2.25      0.00      0.00      0.00    2015.0  
       rates      0.1      0.0     -1.7      0.11      0.00      0.00      0.00     
  ITRF97          6.5     -3.9    -77.9      3.98      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF96          6.5     -3.9    -77.9      3.98      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF94          6.5     -3.9    -77.9      3.98      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF93        -65.8      1.9    -71.3      4.47     -3.36     -4.33      0.75    2015.0  
       rates     -2.8     -0.2     -2.3      0.12     -0.11     -0.19      0.07     
  ITRF92         14.5     -1.9    -85.9      3.27      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF91         26.5     12.1    -91.9      4.67      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF90         24.5      8.1   -107.9      4.97      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF89         29.5     32.1   -145.9      8.37      0.00      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
  ITRF88         24.5     -3.9   -169.9     11.47      0.10      0.00      0.36    2015.0  
       rates      0.1     -0.6     -3.1      0.12      0.00      0.00      0.02     
_________________________________________________________________________________________

Note : These parameters are derived from those already published in the IERS Technical 
Notes and Annual Reports, and from the transformation parameters between ITRF2020, ITRF2014
and ITRF2008. They supersede all values published in the past IERS/ITRF documentations.
The transformation parameters should be used with the standard model (1) given below and 
are valid at the indicated epoch.


: XS :    : X :   : Tx :   :  D   -Rz   Ry : : X :
:    :    :   :   :    :   :               : :   :
: YS :  = : Y : + : Ty : + :  Rz   D   -Rx : : Y :                       (1)
:    :    :   :   :    :   :               : :   :
: ZS :    : Z :   : Tz :   : -Ry   Rx   D  : : Z :


Where X,Y,Z are the coordinates in ITRF2020 and XS,YS,ZS are the coordinates in
the other frames.

On the other hand, for a given parameter P, its value at any epoch t
is obtained by using equation (2).

                  .
P(t) = P(EPOCH) + P * (t - EPOCH)                                        (2)


where EPOCH is the epoch indicated in the above table (currently 2015.0)
    .
and P is the rate of that parameter.

)ITRF";
}}
// clang-format on
