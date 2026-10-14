#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"

namespace tudat
{

namespace basic_astrodynamics
{

namespace
{

/*
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
 */

const std::array< double, 342 > oceanLoadingEquilibriumAmplitudes = {
    0.632208,  0.294107,  0.121046,  0.079915,  0.023818,  -0.023589, 0.022994,  0.019333,  -0.017871, 0.017192,  0.016018,  0.004671,
    -0.004662, -0.004519, 0.004470,  0.004467,  0.002589,  -0.002455, -0.002172, 0.001972,  0.001947,  0.001914,  -0.001898, 0.001802,
    0.001304,  0.001170,  0.001130,  0.001061,  -0.001022, -0.001017, 0.001014,  0.000901,  -0.000857, 0.000855,  0.000855,  0.000772,
    0.000741,  0.000741,  -0.000721, 0.000698,  0.000658,  0.000654,  -0.000653, 0.000633,  0.000626,  -0.000598, 0.000590,  0.000544,
    0.000479,  -0.000464, 0.000413,  -0.000390, 0.000373,  0.000366,  0.000366,  -0.000360, -0.000355, 0.000354,  0.000329,  0.000328,
    0.000319,  0.000302,  0.000279,  -0.000274, -0.000272, 0.000248,  -0.000225, 0.000224,  -0.000223, -0.000216, 0.000211,  0.000209,
    0.000194,  0.000185,  -0.000174, -0.000171, 0.000159,  0.000131,  0.000127,  0.000120,  0.000118,  0.000117,  0.000108,  0.000107,
    0.000105,  -0.000102, 0.000102,  0.000099,  -0.000096, 0.000095,  -0.000089, -0.000085, -0.000084, -0.000081, -0.000077, -0.000072,
    -0.000067, 0.000066,  0.000064,  0.000063,  0.000063,  0.000063,  0.000062,  0.000062,  -0.000060, 0.000056,  0.000053,  0.000051,
    0.000050,  0.368645,  -0.262232, -0.121995, -0.050208, 0.050031,  -0.049470, 0.020620,  0.020613,  0.011279,  -0.009530, -0.009469,
    -0.008012, 0.007414,  -0.007300, 0.007227,  -0.007131, -0.006644, 0.005249,  0.004137,  0.004087,  0.003944,  0.003943,  0.003420,
    0.003418,  0.002885,  0.002884,  0.002160,  -0.001936, 0.001934,  -0.001798, 0.001690,  0.001689,  0.001516,  0.001514,  -0.001511,
    0.001383,  0.001372,  0.001371,  -0.001253, -0.001075, 0.001020,  0.000901,  0.000865,  -0.000794, 0.000788,  0.000782,  -0.000747,
    -0.000745, 0.000670,  -0.000603, -0.000597, 0.000542,  0.000542,  -0.000541, -0.000469, -0.000440, 0.000438,  0.000422,  0.000410,
    -0.000374, -0.000365, 0.000345,  0.000335,  -0.000321, -0.000319, 0.000307,  0.000291,  0.000290,  -0.000289, 0.000286,  0.000275,
    0.000271,  0.000263,  -0.000245, 0.000225,  0.000225,  0.000221,  -0.000202, -0.000200, -0.000199, 0.000192,  0.000183,  0.000183,
    0.000183,  -0.000170, 0.000169,  0.000168,  0.000162,  0.000149,  -0.000147, -0.000141, 0.000138,  0.000136,  0.000136,  0.000127,
    0.000127,  -0.000126, -0.000121, -0.000121, 0.000117,  -0.000116, -0.000114, -0.000114, -0.000114, 0.000114,  0.000113,  0.000109,
    0.000108,  0.000106,  -0.000106, -0.000106, 0.000105,  0.000104,  -0.000103, -0.000100, -0.000100, -0.000100, 0.000099,  -0.000098,
    0.000093,  0.000093,  0.000090,  -0.000088, 0.000083,  -0.000083, -0.000082, -0.000081, -0.000079, -0.000077, -0.000075, -0.000075,
    -0.000075, 0.000071,  0.000071,  -0.000071, 0.000068,  0.000068,  0.000065,  0.000065,  0.000064,  0.000064,  0.000064,  -0.000064,
    -0.000060, 0.000056,  0.000056,  0.000053,  0.000053,  0.000053,  -0.000053, 0.000053,  0.000053,  0.000052,  0.000050,  -0.066607,
    -0.035184, -0.030988, 0.027929,  -0.027616, -0.012753, -0.006728, -0.005837, -0.005286, -0.004921, -0.002884, -0.002583, -0.002422,
    0.002310,  0.002283,  -0.002037, 0.001883,  -0.001811, -0.001687, -0.001004, -0.000925, -0.000844, 0.000766,  0.000766,  -0.000700,
    -0.000495, -0.000492, 0.000491,  0.000483,  0.000437,  -0.000416, -0.000384, 0.000374,  -0.000312, -0.000288, -0.000273, 0.000259,
    0.000245,  -0.000232, 0.000229,  -0.000216, 0.000206,  -0.000204, -0.000202, 0.000200,  0.000195,  -0.000190, 0.000187,  0.000180,
    -0.000179, 0.000170,  0.000153,  -0.000137, -0.000119, -0.000119, -0.000112, -0.000110, -0.000110, 0.000107,  -0.000095, -0.000095,
    -0.000091, -0.000090, -0.000081, -0.000079, -0.000079, 0.000077,  -0.000073, 0.000069,  -0.000067, -0.000066, 0.000065,  0.000064,
    -0.000062, 0.000060,  0.000059,  -0.000056, 0.000055,  -0.000051,
};
const std::array< std::array< int, 6 >, 342 > oceanLoadingDoodsonNumbers = { {
        { { 2, 0, 0, 0, 0, 0 } },    { { 2, 2, -2, 0, 0, 0 } },   { { 2, -1, 0, 1, 0, 0 } },   { { 2, 2, 0, 0, 0, 0 } },
        { { 2, 2, 0, 0, 1, 0 } },    { { 2, 0, 0, 0, -1, 0 } },   { { 2, -1, 2, -1, 0, 0 } },  { { 2, -2, 2, 0, 0, 0 } },
        { { 2, 1, 0, -1, 0, 0 } },   { { 2, 2, -3, 0, 0, 1 } },   { { 2, -2, 0, 2, 0, 0 } },   { { 2, -3, 2, 1, 0, 0 } },
        { { 2, 1, -2, 1, 0, 0 } },   { { 2, -1, 0, 1, -1, 0 } },  { { 2, 3, 0, -1, 0, 0 } },   { { 2, 1, 0, 1, 0, 0 } },
        { { 2, 2, 0, 0, 2, 0 } },    { { 2, 2, -1, 0, 0, -1 } },  { { 2, 0, -1, 0, 0, 1 } },   { { 2, 1, 0, 1, 1, 0 } },
        { { 2, 3, 0, -1, 1, 0 } },   { { 2, 0, 1, 0, 0, -1 } },   { { 2, 0, -2, 2, 0, 0 } },   { { 2, -3, 0, 3, 0, 0 } },
        { { 2, -2, 3, 0, 0, -1 } },  { { 2, 4, 0, 0, 0, 0 } },    { { 2, -1, 1, 1, 0, -1 } },  { { 2, -1, 3, -1, 0, -1 } },
        { { 2, 2, 0, 0, -1, 0 } },   { { 2, -1, -1, 1, 0, 1 } },  { { 2, 4, 0, 0, 1, 0 } },    { { 2, -3, 4, -1, 0, 0 } },
        { { 2, -1, 2, -1, -1, 0 } }, { { 2, 3, -2, 1, 0, 0 } },   { { 2, 1, 2, -1, 0, 0 } },   { { 2, -4, 2, 2, 0, 0 } },
        { { 2, 4, -2, 0, 0, 0 } },   { { 2, 0, 2, 0, 0, 0 } },    { { 2, -2, 2, 0, -1, 0 } },  { { 2, 2, -4, 0, 0, 2 } },
        { { 2, 2, -2, 0, -1, 0 } },  { { 2, 1, 0, -1, -1, 0 } },  { { 2, -1, 1, 0, 0, 0 } },   { { 2, 2, -1, 0, 0, 1 } },
        { { 2, 2, 1, 0, 0, -1 } },   { { 2, -2, 0, 2, -1, 0 } },  { { 2, -2, 4, -2, 0, 0 } },  { { 2, 2, 2, 0, 0, 0 } },
        { { 2, -4, 4, 0, 0, 0 } },   { { 2, -1, 0, -1, -2, 0 } }, { { 2, 1, 2, -1, 1, 0 } },   { { 2, -1, -2, 3, 0, 0 } },
        { { 2, 3, -2, 1, 1, 0 } },   { { 2, 4, 0, -2, 0, 0 } },   { { 2, 0, 0, 2, 0, 0 } },    { { 2, 0, 2, -2, 0, 0 } },
        { { 2, 0, 2, 0, 1, 0 } },    { { 2, -3, 3, 1, 0, -1 } },  { { 2, 0, 0, 0, -2, 0 } },   { { 2, 4, 0, 0, 2, 0 } },
        { { 2, 4, -2, 0, 1, 0 } },   { { 2, 0, 0, 0, 0, 2 } },    { { 2, 1, 0, 1, 2, 0 } },    { { 2, 0, -2, 0, -2, 0 } },
        { { 2, -2, 1, 0, 0, 1 } },   { { 2, -2, 1, 2, 0, -1 } },  { { 2, -1, 1, -1, 0, 1 } },  { { 2, 5, 0, -1, 0, 0 } },
        { { 2, 1, -3, 1, 0, 1 } },   { { 2, -2, -1, 2, 0, 1 } },  { { 2, 3, 0, -1, 2, 0 } },   { { 2, 1, -2, 1, -1, 0 } },
        { { 2, 5, 0, -1, 1, 0 } },   { { 2, -4, 0, 4, 0, 0 } },   { { 2, -3, 2, 1, -1, 0 } },  { { 2, -2, 1, 1, 0, 0 } },
        { { 2, 4, 0, -2, 1, 0 } },   { { 2, 0, 0, 2, 1, 0 } },    { { 2, -5, 4, 1, 0, 0 } },   { { 2, 0, 2, 0, 2, 0 } },
        { { 2, -1, 2, 1, 0, 0 } },   { { 2, 5, -2, -1, 0, 0 } },  { { 2, 1, -1, 0, 0, 0 } },   { { 2, 2, -2, 0, 0, 2 } },
        { { 2, -5, 2, 3, 0, 0 } },   { { 2, -1, -2, 1, -2, 0 } }, { { 2, -3, 5, -1, 0, -1 } }, { { 2, -1, 0, 0, 0, 1 } },
        { { 2, -2, 0, 0, -2, 0 } },  { { 2, 0, -1, 1, 0, 0 } },   { { 2, -3, 1, 1, 0, 1 } },   { { 2, 3, 0, -1, -1, 0 } },
        { { 2, 1, 0, 1, -1, 0 } },   { { 2, -1, 2, 1, 1, 0 } },   { { 2, 0, -3, 2, 0, 1 } },   { { 2, 1, -1, -1, 0, 1 } },
        { { 2, -3, 0, 3, -1, 0 } },  { { 2, 0, -2, 2, -1, 0 } },  { { 2, -4, 3, 2, 0, -1 } },  { { 2, -1, 0, 1, -2, 0 } },
        { { 2, 5, 0, -1, 2, 0 } },   { { 2, -4, 5, 0, 0, -1 } },  { { 2, -2, 4, 0, 0, -2 } },  { { 2, -1, 0, 1, 0, 2 } },
        { { 2, -2, -2, 4, 0, 0 } },  { { 2, 3, -2, -1, -1, 0 } }, { { 2, -2, 5, -2, 0, -1 } }, { { 2, 0, -1, 0, -1, 1 } },
        { { 2, 5, -2, -1, 1, 0 } },  { { 1, 1, 0, 0, 0, 0 } },    { { 1, -1, 0, 0, 0, 0 } },   { { 1, 1, -2, 0, 0, 0 } },
        { { 1, -2, 0, 1, 0, 0 } },   { { 1, 1, 0, 0, 1, 0 } },    { { 1, -1, 0, 0, -1, 0 } },  { { 1, 2, 0, -1, 0, 0 } },
        { { 1, 0, 0, 1, 0, 0 } },    { { 1, 3, 0, 0, 0, 0 } },    { { 1, -2, 2, -1, 0, 0 } },  { { 1, -2, 0, 1, -1, 0 } },
        { { 1, -3, 2, 0, 0, 0 } },   { { 1, 0, 0, -1, 0, 0 } },   { { 1, 1, 0, 0, -1, 0 } },   { { 1, 3, 0, 0, 1, 0 } },
        { { 1, 1, -3, 0, 0, 1 } },   { { 1, -3, 0, 2, 0, 0 } },   { { 1, 1, 2, 0, 0, 0 } },    { { 1, 0, 0, 1, 1, 0 } },
        { { 1, 2, 0, -1, 1, 0 } },   { { 1, 0, 2, -1, 0, 0 } },   { { 1, 2, -2, 1, 0, 0 } },   { { 1, 3, -2, 0, 0, 0 } },
        { { 1, -1, 2, 0, 0, 0 } },   { { 1, 1, 1, 0, 0, -1 } },   { { 1, 1, -1, 0, 0, 1 } },   { { 1, 4, 0, -1, 0, 0 } },
        { { 1, -4, 2, 1, 0, 0 } },   { { 1, 0, -2, 1, 0, 0 } },   { { 1, -2, 2, -1, -1, 0 } }, { { 1, 3, 0, -2, 0, 0 } },
        { { 1, -1, 0, 2, 0, 0 } },   { { 1, -1, 0, 0, -2, 0 } },  { { 1, 3, 0, 0, 2, 0 } },    { { 1, -3, 2, 0, -1, 0 } },
        { { 1, 4, 0, -1, 1, 0 } },   { { 1, 0, 0, -1, -1, 0 } },  { { 1, 1, -2, 0, -1, 0 } },  { { 1, -3, 0, 2, -1, 0 } },
        { { 1, 1, 0, 0, 2, 0 } },    { { 1, 1, -1, 0, 0, -1 } },  { { 1, -1, -1, 0, 0, 1 } },  { { 1, 0, 2, -1, 1, 0 } },
        { { 1, -1, 1, 0, 0, -1 } },  { { 1, -1, -2, 2, 0, 0 } },  { { 1, 2, -2, 1, 1, 0 } },   { { 1, -4, 0, 3, 0, 0 } },
        { { 1, -1, 2, 0, 1, 0 } },   { { 1, 3, -2, 0, 1, 0 } },   { { 1, 2, 0, -1, -1, 0 } },  { { 1, 0, 0, 1, -1, 0 } },
        { { 1, -2, 2, 1, 0, 0 } },   { { 1, 4, -2, -1, 0, 0 } },  { { 1, -3, 3, 0, 0, -1 } },  { { 1, -2, 1, 1, 0, -1 } },
        { { 1, -2, 3, -1, 0, -1 } }, { { 1, 0, -2, 1, -1, 0 } },  { { 1, -2, -1, 1, 0, 1 } },  { { 1, 4, -2, 1, 0, 0 } },
        { { 1, -4, 4, -1, 0, 0 } },  { { 1, -4, 2, 1, -1, 0 } },  { { 1, 5, -2, 0, 0, 0 } },   { { 1, 3, 0, -2, 1, 0 } },
        { { 1, -5, 2, 2, 0, 0 } },   { { 1, 2, 0, 1, 0, 0 } },    { { 1, 1, 3, 0, 0, -1 } },   { { 1, -2, 0, 1, -2, 0 } },
        { { 1, 4, 0, -1, 2, 0 } },   { { 1, 1, -4, 0, 0, 2 } },   { { 1, 5, 0, -2, 0, 0 } },   { { 1, -1, 0, 2, 1, 0 } },
        { { 1, -2, 1, 0, 0, 0 } },   { { 1, 4, -2, 1, 1, 0 } },   { { 1, -3, 4, -2, 0, 0 } },  { { 1, -1, 3, 0, 0, -1 } },
        { { 1, 3, -3, 0, 0, 1 } },   { { 1, 5, -2, 0, 1, 0 } },   { { 1, 1, 2, 0, 1, 0 } },    { { 1, 2, 0, 1, 1, 0 } },
        { { 1, -5, 4, 0, 0, 0 } },   { { 1, -2, 0, -1, -2, 0 } }, { { 1, 5, 0, -2, 1, 0 } },   { { 1, 1, 2, -2, 0, 0 } },
        { { 1, 1, -2, 2, 0, 0 } },   { { 1, -2, 2, 1, 1, 0 } },   { { 1, 0, 3, -1, 0, -1 } },  { { 1, 2, -3, 1, 0, 1 } },
        { { 1, -2, -2, 3, 0, 0 } },  { { 1, -1, 2, -2, 0, 0 } },  { { 1, -4, 3, 1, 0, -1 } },  { { 1, -4, 0, 3, -1, 0 } },
        { { 1, -1, -2, 2, -1, 0 } }, { { 1, -2, 0, 3, 0, 0 } },   { { 1, 4, 0, -3, 0, 0 } },   { { 1, 0, 1, 1, 0, -1 } },
        { { 1, 2, -1, -1, 0, 1 } },  { { 1, 2, -2, 1, -1, 0 } },  { { 1, 0, 0, -1, -2, 0 } },  { { 1, 2, 0, 1, 2, 0 } },
        { { 1, 2, -2, -1, -1, 0 } }, { { 1, 0, 0, 1, 2, 0 } },    { { 1, 0, 1, 0, 0, 0 } },    { { 1, 2, -1, 0, 0, 0 } },
        { { 1, 0, 2, -1, -1, 0 } },  { { 1, -1, -2, 0, -2, 0 } }, { { 1, -3, 1, 0, 0, 1 } },   { { 1, 3, -2, 0, -1, 0 } },
        { { 1, -1, -1, 0, -1, 1 } }, { { 1, 4, -2, -1, 1, 0 } },  { { 1, 2, 1, -1, 0, -1 } },  { { 1, 0, -1, 1, 0, 1 } },
        { { 1, -2, 4, -1, 0, 0 } },  { { 1, 4, -4, 1, 0, 0 } },   { { 1, -3, 1, 2, 0, -1 } },  { { 1, -3, 3, 0, -1, -1 } },
        { { 1, 1, 2, 0, 2, 0 } },    { { 1, 1, -2, 0, -2, 0 } },  { { 1, 3, 0, 0, 3, 0 } },    { { 1, -1, 2, 0, -1, 0 } },
        { { 1, -2, 1, -1, 0, 1 } },  { { 1, 0, -3, 1, 0, 1 } },   { { 1, -3, -1, 2, 0, 1 } },  { { 1, 2, 0, -1, 2, 0 } },
        { { 1, 6, -2, -1, 0, 0 } },  { { 1, 2, 2, -1, 0, 0 } },   { { 1, -1, 1, 0, -1, -1 } }, { { 1, -2, 3, -1, -1, -1 } },
        { { 1, -1, 0, 0, 0, 2 } },   { { 1, -5, 0, 4, 0, 0 } },   { { 1, 1, 0, 0, 0, -2 } },   { { 1, -2, 1, 1, -1, -1 } },
        { { 1, 1, -1, 0, 1, 1 } },   { { 1, 1, 2, 0, 0, -2 } },   { { 1, -3, 1, 1, 0, 0 } },   { { 1, -4, 4, -1, -1, 0 } },
        { { 1, 1, 0, -2, -1, 0 } },  { { 1, -2, -1, 1, -1, 1 } }, { { 1, -3, 2, 2, 0, 0 } },   { { 1, 5, -2, -2, 0, 0 } },
        { { 1, 3, -4, 2, 0, 0 } },   { { 1, 1, -2, 0, 0, 2 } },   { { 1, -1, 4, -2, 0, 0 } },  { { 1, 2, 2, -1, 1, 0 } },
        { { 1, -5, 2, 2, -1, 0 } },  { { 1, 1, -3, 0, -1, 1 } },  { { 1, 1, 1, 0, 1, -1 } },   { { 1, 6, -2, -1, 1, 0 } },
        { { 1, -2, 2, -1, -2, 0 } }, { { 1, 4, -2, 1, 2, 0 } },   { { 1, -6, 4, 1, 0, 0 } },   { { 1, 5, -4, 0, 0, 0 } },
        { { 1, -3, 4, 0, 0, 0 } },   { { 1, 1, 2, -2, 1, 0 } },   { { 1, -2, 1, 0, -1, 0 } },  { { 0, 2, 0, 0, 0, 0 } },
        { { 0, 1, 0, -1, 0, 0 } },   { { 0, 0, 2, 0, 0, 0 } },    { { 0, 0, 0, 0, 1, 0 } },    { { 0, 2, 0, 0, 1, 0 } },
        { { 0, 3, 0, -1, 0, 0 } },   { { 0, 1, -2, 1, 0, 0 } },   { { 0, 2, -2, 0, 0, 0 } },   { { 0, 3, 0, -1, 1, 0 } },
        { { 0, 0, 1, 0, 0, -1 } },   { { 0, 2, 0, -2, 0, 0 } },   { { 0, 2, 0, 0, 2, 0 } },    { { 0, 3, -2, 1, 0, 0 } },
        { { 0, 1, 0, -1, -1, 0 } },  { { 0, 1, 0, -1, 1, 0 } },   { { 0, 4, -2, 0, 0, 0 } },   { { 0, 1, 0, 1, 0, 0 } },
        { { 0, 0, 3, 0, 0, -1 } },   { { 0, 4, 0, -2, 0, 0 } },   { { 0, 3, -2, 1, 1, 0 } },   { { 0, 3, -2, -1, 0, 0 } },
        { { 0, 4, -2, 0, 1, 0 } },   { { 0, 0, 2, 0, 1, 0 } },    { { 0, 1, 0, 1, 1, 0 } },    { { 0, 4, 0, -2, 1, 0 } },
        { { 0, 3, 0, -1, 2, 0 } },   { { 0, 5, -2, -1, 0, 0 } },  { { 0, 1, 2, -1, 0, 0 } },   { { 0, 1, -2, 1, -1, 0 } },
        { { 0, 1, -2, 1, 1, 0 } },   { { 0, 2, -2, 0, -1, 0 } },  { { 0, 2, -3, 0, 0, 1 } },   { { 0, 2, -2, 0, 1, 0 } },
        { { 0, 0, 2, -2, 0, 0 } },   { { 0, 1, -3, 1, 0, 1 } },   { { 0, 0, 0, 0, 2, 0 } },    { { 0, 0, 1, 0, 0, 1 } },
        { { 0, 1, 2, -1, 1, 0 } },   { { 0, 3, 0, -3, 0, 0 } },   { { 0, 2, 1, 0, 0, -1 } },   { { 0, 1, -1, -1, 0, 1 } },
        { { 0, 1, 0, 1, 2, 0 } },    { { 0, 5, -2, -1, 1, 0 } },  { { 0, 2, -1, 0, 0, 1 } },   { { 0, 2, 2, -2, 0, 0 } },
        { { 0, 1, -1, 0, 0, 0 } },   { { 0, 5, 0, -3, 0, 0 } },   { { 0, 2, 0, -2, 1, 0 } },   { { 0, 1, 1, -1, 0, -1 } },
        { { 0, 3, -4, 1, 0, 0 } },   { { 0, 0, 2, 0, 2, 0 } },    { { 0, 2, 0, -2, -1, 0 } },  { { 0, 4, -3, 0, 0, 1 } },
        { { 0, 3, -1, -1, 0, 1 } },  { { 0, 0, 2, 0, 0, -2 } },   { { 0, 3, -3, 1, 0, 1 } },   { { 0, 2, -4, 2, 0, 0 } },
        { { 0, 4, -2, -2, 0, 0 } },  { { 0, 3, 1, -1, 0, -1 } },  { { 0, 5, -4, 1, 0, 0 } },   { { 0, 3, -2, -1, -1, 0 } },
        { { 0, 3, -2, 1, 2, 0 } },   { { 0, 4, -4, 0, 0, 0 } },   { { 0, 6, -2, -2, 0, 0 } },  { { 0, 5, 0, -3, 1, 0 } },
        { { 0, 4, -2, 0, 2, 0 } },   { { 0, 2, 2, -2, 1, 0 } },   { { 0, 0, 4, 0, 0, -2 } },   { { 0, 3, -1, 0, 0, 0 } },
        { { 0, 3, -3, -1, 0, 1 } },  { { 0, 4, 0, -2, 2, 0 } },   { { 0, 1, -2, -1, -1, 0 } }, { { 0, 2, -1, 0, 0, -1 } },
        { { 0, 4, -4, 2, 0, 0 } },   { { 0, 2, 1, 0, 1, -1 } },   { { 0, 3, -2, -1, 1, 0 } },  { { 0, 4, -3, 0, 1, 1 } },
        { { 0, 2, 0, 0, 3, 0 } },    { { 0, 6, -4, 0, 0, 0 } },
} };

// This renamed C++ implementation is derived from the IERS HARDISP package,
// corrected 19 June 2015 and updated 19 December 2016. It is neither distributed
// nor endorsed by the IERS Conventions Center. ADMINT/SPLINE/EVAL/TDFRPH have
// been translated to the functions below; single-precision admittance arithmetic
// is promoted to double, Shell sort uses std::sort, and each requested epoch is
// summed directly rather than through the equally spaced RECURS time series.
// TT and UTC come from the caller's TT epoch and SOFA conversion, replacing
// the reference ETUTC decimal-year lookup with the maintained time conversion.
// Pinned source: https://iers-conventions.obspm.fr/content/chapter7/software/hardisp/
// ADMINT.F SHA256: 478e3f0c001c09dd4d9e2e9920033d29e76eb60b808238e416e6db9ffaeae6c8
// TDFRPH.F SHA256: 24068e0cd1e2e210fab7dd8d4473bb60ff8335bfcb53a2fc12dad7e6d1cccd19

struct OceanLoadingAstronomicalArguments {
    std::array< double, 6 > phasesDegrees;
    std::array< double, 6 > frequenciesCyclesPerDay;
};

OceanLoadingAstronomicalArguments computeOceanLoadingAstronomicalArguments( const double tt )
{
    const double utc = sofa_interface::convertTTtoUTC( tt );
    const double utcDaysSinceMidnight = ( utc + 43200.0 ) / 86400.0;
    const double dayFraction = utcDaysSinceMidnight - std::floor( utcDaysSinceMidnight );
    const double t = tt / ( 86400.0 * 36525.0 );
    const double f1 = 134.9634025100 + t * ( 477198.8675605000 + t * ( 0.0088553333 + t * ( 0.0000143431 + t * -0.0000000680 ) ) );
    const double f2 = 357.5291091806 + t * ( 35999.0502911389 + t * ( -0.0001536667 + t * ( 0.0000000378 + t * -0.0000000032 ) ) );
    const double f3 = 93.2720906200 + t * ( 483202.0174577222 + t * ( -0.0035420000 + t * ( -0.0000002881 + t * 0.0000000012 ) ) );
    const double f4 = 297.8501954694 + t * ( 445267.1114469445 + t * ( -0.0017696111 + t * ( 0.0000018314 + t * -0.0000000088 ) ) );
    const double f5 = 125.0445550100 + t * ( -1934.1362619722 + t * ( 0.0020756111 + t * ( 0.0000021394 + t * -0.0000000165 ) ) );
    const double fd1 = 0.0362916471 + 0.0000000013 * t;
    const double fd2 = 0.0027377786;
    const double fd3 = 0.0367481951 - 0.0000000005 * t;
    const double fd4 = 0.0338631920 - 0.0000000003 * t;
    const double fd5 = -0.0001470938 + 0.0000000003 * t;
    OceanLoadingAstronomicalArguments result;
    result.phasesDegrees = { { 360.0 * dayFraction - f4, f3 + f5, f3 + f5 - f4, f3 + f5 - f1, -f5, f3 + f5 - f4 - f2 } };
    result.frequenciesCyclesPerDay = { { 1.0 - fd4, fd3 + fd5, fd3 + fd5 - fd4, fd3 + fd5 - fd1, -fd5, fd3 + fd5 - fd4 - fd2 } };
    return result;
}

struct OceanLoadingAdmittanceSpline {
    std::vector< double > frequencies;
    std::vector< double > values;
    std::vector< double > secondDerivatives;
};

OceanLoadingAdmittanceSpline makeOceanLoadingSpline( const std::vector< double >& x, const std::vector< double >& y )
{
    const std::size_t n = x.size( );
    OceanLoadingAdmittanceSpline spline{ x, y, std::vector< double >( n, 0.0 ) };
    if( n <= 3 )
    {
        return spline;
    }
    const auto endpointSlope = []( const double u1, const double x1, const double u2, const double x2 ) {
        return ( u1 / ( x1 * x1 ) - u2 / ( x2 * x2 ) ) / ( 1.0 / x1 - 1.0 / x2 );
    };
    const double q1 = endpointSlope( y[ 1 ] - y[ 0 ], x[ 1 ] - x[ 0 ], y[ 2 ] - y[ 0 ], x[ 2 ] - x[ 0 ] );
    const double qn = endpointSlope( y[ n - 2 ] - y[ n - 1 ], x[ n - 2 ] - x[ n - 1 ], y[ n - 3 ] - y[ n - 1 ], x[ n - 3 ] - x[ n - 1 ] );
    auto& s = spline.secondDerivatives;
    s[ 0 ] = 6.0 * ( ( y[ 1 ] - y[ 0 ] ) / ( x[ 1 ] - x[ 0 ] ) - q1 );
    for( std::size_t i = 1; i + 1 < n; ++i )
    {
        s[ i ] = 6.0 *
                ( y[ i - 1 ] / ( x[ i ] - x[ i - 1 ] ) - y[ i ] * ( 1.0 / ( x[ i ] - x[ i - 1 ] ) + 1.0 / ( x[ i + 1 ] - x[ i ] ) ) +
                  y[ i + 1 ] / ( x[ i + 1 ] - x[ i ] ) );
    }
    s[ n - 1 ] = 6.0 * ( qn + ( y[ n - 2 ] - y[ n - 1 ] ) / ( x[ n - 1 ] - x[ n - 2 ] ) );
    std::vector< double > a( n );
    a[ 0 ] = 2.0 * ( x[ 1 ] - x[ 0 ] );
    a[ 1 ] = 1.5 * ( x[ 1 ] - x[ 0 ] ) + 2.0 * ( x[ 2 ] - x[ 1 ] );
    s[ 1 ] -= 0.5 * s[ 0 ];
    for( std::size_t i = 2; i + 1 < n; ++i )
    {
        const double c = ( x[ i ] - x[ i - 1 ] ) / a[ i - 1 ];
        a[ i ] = 2.0 * ( x[ i + 1 ] - x[ i - 1 ] ) - c * ( x[ i ] - x[ i - 1 ] );
        s[ i ] -= c * s[ i - 1 ];
    }
    const double c = ( x[ n - 1 ] - x[ n - 2 ] ) / a[ n - 2 ];
    a[ n - 1 ] = ( 2.0 - c ) * ( x[ n - 1 ] - x[ n - 2 ] );
    s[ n - 1 ] = ( s[ n - 1 ] - c * s[ n - 2 ] ) / a[ n - 1 ];
    for( std::size_t i = n - 1; i-- > 0; )
    {
        s[ i ] = ( s[ i ] - ( x[ i + 1 ] - x[ i ] ) * s[ i + 1 ] ) / a[ i ];
    }
    return spline;
}

double evaluateOceanLoadingSpline( const OceanLoadingAdmittanceSpline& spline, const double frequency )
{
    const auto& x = spline.frequencies;
    const auto& y = spline.values;
    const auto& s = spline.secondDerivatives;
    if( frequency <= x.front( ) )
    {
        return y.front( );
    }
    if( frequency >= x.back( ) )
    {
        return y.back( );
    }
    const auto upper = std::lower_bound( x.begin( ), x.end( ), frequency );
    const std::size_t right = static_cast< std::size_t >( upper - x.begin( ) );
    const std::size_t left = right - 1;
    const double dy = x[ right ] - frequency;
    const double dy1 = frequency - x[ left ];
    const double dk = x[ right ] - x[ left ];
    return ( s[ left ] * dy * dy * dy + s[ right ] * dy1 * dy1 * dy1 ) / ( 6.0 * dk ) + dy1 * ( y[ right ] / dk - s[ right ] * dk / 6.0 ) +
            dy * ( y[ left ] / dk - s[ left ] * dk / 6.0 );
}

Eigen::Vector3d synthesizeExpandedOceanLoading( const double tt, const input_output::OceanLoadingStation& station )
{
    const auto astronomy = computeOceanLoadingAstronomicalArguments( tt );
    const double degreesToRadians = mathematical_constants::PI / 180.0;
    std::array< double, 342 > frequencies;
    std::array< double, 342 > phases;
    for( std::size_t tide = 0; tide < 342; ++tide )
    {
        double phase = 0.0;
        frequencies[ tide ] = 0.0;
        for( std::size_t argument = 0; argument < 6; ++argument )
        {
            phase += oceanLoadingDoodsonNumbers[ tide ][ argument ] * astronomy.phasesDegrees[ argument ];
            frequencies[ tide ] += oceanLoadingDoodsonNumbers[ tide ][ argument ] * astronomy.frequenciesCyclesPerDay[ argument ];
        }
        phase = std::fmod( phase, 360.0 );
        const int species = oceanLoadingDoodsonNumbers[ tide ][ 0 ];
        phases[ tide ] = ( phase + ( species == 0 ? 180.0 : ( species == 1 ? 90.0 : 0.0 ) ) ) * degreesToRadians;
    }
    // Canonical BLQ constituents, matched against the pinned reference Doodson table.
    const std::array< std::array< int, 6 >, 11 > mainDoodsonNumbers = { { { { 2, 0, 0, 0, 0, 0 } },
                                                                          { { 2, 2, -2, 0, 0, 0 } },
                                                                          { { 2, -1, 0, 1, 0, 0 } },
                                                                          { { 2, 2, 0, 0, 0, 0 } },
                                                                          { { 1, 1, 0, 0, 0, 0 } },
                                                                          { { 1, -1, 0, 0, 0, 0 } },
                                                                          { { 1, 1, -2, 0, 0, 0 } },
                                                                          { { 1, -2, 0, 1, 0, 0 } },
                                                                          { { 0, 2, 0, 0, 0, 0 } },
                                                                          { { 0, 1, 0, -1, 0, 0 } },
                                                                          { { 0, 0, 2, 0, 0, 0 } } } };
    std::array< std::size_t, 11 > mainIndices;
    for( std::size_t column = 0; column < 11; ++column )
    {
        const auto found =
                std::find( oceanLoadingDoodsonNumbers.begin( ), oceanLoadingDoodsonNumbers.end( ), mainDoodsonNumbers[ column ] );
        if( found == oceanLoadingDoodsonNumbers.end( ) )
        {
            throw std::runtime_error( "Internal ocean loading table is missing a main constituent" );
        }
        mainIndices[ column ] = static_cast< std::size_t >( found - oceanLoadingDoodsonNumbers.begin( ) );
    }
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );
    for( int component = 0; component < 3; ++component )
    {
        std::array< OceanLoadingAdmittanceSpline, 3 > realSplines;
        std::array< OceanLoadingAdmittanceSpline, 3 > imaginarySplines;
        for( int species = 0; species < 3; ++species )
        {
            std::vector< std::size_t > columns;
            for( std::size_t column = 0; column < 11; ++column )
            {
                if( mainDoodsonNumbers[ column ][ 0 ] == species )
                {
                    columns.push_back( column );
                }
            }
            std::sort( columns.begin( ), columns.end( ), [ & ]( const std::size_t left, const std::size_t right ) {
                return frequencies[ mainIndices[ left ] ] < frequencies[ mainIndices[ right ] ];
            } );
            std::vector< double > x, real, imaginary;
            for( const std::size_t column : columns )
            {
                const std::size_t tide = mainIndices[ column ];
                const double amplitude =
                        station.amplitudesMetres( component, column ) / std::abs( oceanLoadingEquilibriumAmplitudes[ tide ] );
                const double lag = station.phaseLagsDegrees( component, column ) * degreesToRadians;
                x.push_back( frequencies[ tide ] );
                real.push_back( amplitude * std::cos( lag ) );
                imaginary.push_back( -amplitude * std::sin( lag ) );
            }
            realSplines[ species ] = makeOceanLoadingSpline( x, real );
            imaginarySplines[ species ] = makeOceanLoadingSpline( x, imaginary );
        }
        for( std::size_t tide = 0; tide < 342; ++tide )
        {
            const int species = oceanLoadingDoodsonNumbers[ tide ][ 0 ];
            const double real = evaluateOceanLoadingSpline( realSplines[ species ], frequencies[ tide ] );
            const double imaginary = evaluateOceanLoadingSpline( imaginarySplines[ species ], frequencies[ tide ] );
            // Signed Cartwright-Edden amplitude and complex admittance retain the
            // reference phase convention, without a singular atan2 at zero amplitude.
            displacement( component ) += oceanLoadingEquilibriumAmplitudes[ tide ] *
                    ( real * std::cos( phases[ tide ] ) - imaginary * std::sin( phases[ tide ] ) );
        }
    }
    return displacement;
}

}  // namespace

//! Compatibility wrapper; neutral input records retain phases in their original degrees.
std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > readBlqFile( const std::string& fileName )
{
    const auto data = input_output::readOceanLoadingData( fileName );
    std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > result( 2 );
    for( const auto& station : data.stations )
    {
        result[ 0 ].emplace( station.first, station.second.amplitudesMetres );
        result[ 1 ].emplace( station.first, unit_conversions::convertDegreesToRadians( station.second.phaseLagsDegrees ) );
    }
    return result;
}

Eigen::Vector3d convertLocalOceanTideDisplacementToPlentocentricDisplacement(
        const Eigen::Vector3d& localDisplacement,
        const std::vector< Eigen::Vector3d >& earthNorthradialUnitVectors )
{
    if( earthNorthradialUnitVectors.size( ) != 3 )
    {
        throw std::runtime_error( "Ocean loading Cartesian conversion requires east/north/radial basis vectors" );
    }
    return -localDisplacement( 1 ) * earthNorthradialUnitVectors[ 0 ] - localDisplacement( 2 ) * earthNorthradialUnitVectors[ 1 ] +
            localDisplacement( 0 ) * earthNorthradialUnitVectors[ 2 ];
}

OceanTideEarthDeformation::OceanTideEarthDeformation( const std::vector< std::string >& blqFiles, const OceanTideSynthesis synthesis ):
    synthesis_( synthesis ),
    doodsonArgumentFunction_( []( const double time ) { return sofa_interface::calculateDoodsonFundamentalArguments( time ); } )
{
    if( synthesis != OceanTideSynthesis::iersHardisp && synthesis != OceanTideSynthesis::legacyElevenConstituents )
    {
        throw std::runtime_error( "Unknown ocean tidal loading synthesis convention" );
    }
    initializeLegacyConstituents( );
    for( const auto& file : blqFiles )
    {
        addBlqFile( file );
    }
}

OceanTideEarthDeformation::OceanTideEarthDeformation( const input_output::OceanLoadingData& data, const OceanTideSynthesis synthesis ):
    OceanTideEarthDeformation( std::vector< std::string >( ), synthesis )
{
    addOceanLoadingData( data );
}

OceanTideEarthDeformation::OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                                                      const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
    OceanTideEarthDeformation( blqFiles, OceanTideSynthesis::legacyElevenConstituents )
{
    if( !doodsonArgumentFunction )
    {
        throw std::runtime_error( "Ocean loading legacy argument callback is empty" );
    }
    doodsonArgumentFunction_ = doodsonArgumentFunction;
}

OceanTideEarthDeformation::OceanTideEarthDeformation( const input_output::OceanLoadingData& data,
                                                      const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
    OceanTideEarthDeformation( std::vector< std::string >( ), doodsonArgumentFunction )
{
    addOceanLoadingData( data );
}

void OceanTideEarthDeformation::initializeLegacyConstituents( )
{
    // Initialize characteristics of 11 main tides, from Moyer(2000), Table 5-1 and Eq. 5-88
    doodsonMultipliers_ << 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0,

            1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -2.0, 0.0, 1.0, -2.0, 0.0, 1.0,

            0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0, 0.0;

    schwiderskiFactors_ << 0.0, 0.0, 0.0, 0.0, 0.5, -.50, -.50, -.50, 0.0, 0.0, 0.0;
    schwiderskiFactors_ *= mathematical_constants::PI;
}

void OceanTideEarthDeformation::addBlqFile( const std::string& blqFile )
{
    addOceanLoadingData( input_output::readOceanLoadingData( blqFile ) );
}

void OceanTideEarthDeformation::addOceanLoadingData( const input_output::OceanLoadingData& data )
{
    input_output::validateOceanLoadingData( data );
    for( const auto& item : data.stations )
    {
        if( oceanLoadingData_.stations.count( item.first ) != 0 )
        {
            throw std::runtime_error( "Duplicate ocean loading station: " + item.first );
        }
    }
    oceanLoadingData_.source = oceanLoadingData_.source.empty( ) ? data.source : "multiple BLQ inputs";
    oceanLoadingData_.sourceHeaders.insert(
            oceanLoadingData_.sourceHeaders.end( ), data.sourceHeaders.begin( ), data.sourceHeaders.end( ) );
    for( const auto& item : data.stations )
    {
        oceanLoadingData_.stations.emplace( item.first, item.second );
        siteOceanTideAmplitudes_.emplace( item.first, item.second.amplitudesMetres );
        siteOceanTidePhases_.emplace( item.first, unit_conversions::convertDegreesToRadians( item.second.phaseLagsDegrees ) );
    }
}

Eigen::Vector3d OceanTideEarthDeformation::calculateGeocentricEnuDisplacement( const double time, const std::string& siteIdentifier )
{
    if( !std::isfinite( time ) )
    {
        throw std::runtime_error( "Ocean loading evaluation time must be finite" );
    }
    if( synthesis_ == OceanTideSynthesis::iersHardisp )
    {
        const Eigen::Vector3d uws = synthesizeExpandedOceanLoading( time, oceanLoadingData_.getStation( siteIdentifier ) );
        return Eigen::Vector3d( -uws( 1 ), -uws( 2 ), uws( 0 ) );
    }
    // Check to see if requested site is available.
    if( siteOceanTideAmplitudes_.count( siteIdentifier ) == 0 )
    {
        throw std::runtime_error( "Error, requested site: " + siteIdentifier +
                                  ", which is not available for ocean tide displacement calculation (blq file missing)" );
    }

    // Calculate doodson arguments at given time.
    Eigen::Vector6d doodsonArguments = doodsonArgumentFunction_( time );

    if( !doodsonArguments.allFinite( ) )
    {
        throw std::runtime_error( "Ocean loading legacy arguments must be finite" );
    }
    // Initialize displacement to zero.
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );

    Eigen::Matrix< double, 11, 1 > arguments;

    // Iterate over all tidal constituents and calculate current phase and resulting displacement.
    for( int i = 0; i < 11; i++ )
    {
        arguments[ i ] = 0.0;

        // Combine influence of doodson arguments (argument 5 and 6 are zero for all 11 constituent tides).
        for( int j = 0; j < 4; j++ )
        {
            arguments[ i ] += doodsonArguments[ j ] * doodsonMultipliers_( i, j );
        }

        // Add schwiderski phase factor for applicable tide.s
        arguments[ i ] += schwiderskiFactors_( i );

        // Calculate site displacement due to tide in three dimensions.
        for( int j = 0; j < 3; j++ )
        {
            displacement( j ) += siteOceanTideAmplitudes_[ siteIdentifier ]( j, i ) *
                    std::cos( arguments[ i ] - siteOceanTidePhases_[ siteIdentifier ]( j, i ) );
        }
    }

    return ( Eigen::Vector3d( ) << -displacement( 1 ), -displacement( 2 ), displacement( 0 ) ).finished( );
}

Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacementInEnuFrame( const double time, const std::string& siteIdentifier )
{
    return calculateGeocentricEnuDisplacement( time, siteIdentifier );
}

Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacementInEnuFrame(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    return stationState->getRotationFromBodyFixedToTopocentricFrame( time ) * calculateDisplacement( time, stationState );
}

//! Function to calculate site displacement at given time.
Eigen::Vector3d OceanTideEarthDeformation::calculateDisplacement(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    const Eigen::Vector3d enu = calculateGeocentricEnuDisplacement( time, stationState->getSiteId( ) );
    const auto basis = stationState->getEnuGeocentricUnitVectors( );
    return enu( 0 ) * basis[ 0 ] + enu( 1 ) * basis[ 1 ] + enu( 2 ) * basis[ 2 ];
}

}  // namespace basic_astrodynamics

}  // namespace tudat
