/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <iostream>

#include "Tudat/Basics/utilityMacros.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"

namespace tudat
{

namespace unit_tests
{

Eigen::Matrix3d convertArrayToMatrix(
        const double array[3][3] )
{
    Eigen::Matrix3d matrix;
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            matrix( i, j ) = array[ i ][ j ];
        }
    }
    return matrix;
}

Eigen::Matrix3d getSofaEarthOrientationExamples(
        const int calculationCase,
        const double dXInMas,
        const double dYInMas,
        const double xPInAs,
        const double yPInAs,
        const double ut1Correction )
{

//    std::cout << "Angles: " << calculationCase << " " << dXInMas << " " << dYInMas << " " << xPInAs << " " << yPInAs << " " << ut1Correction << std::endl;

    double AS2R = 4.848136811095359935899141E-6;

    int IY, IM, ID, IH, MIN;
    double SEC, XP, YP, DUT1,
            DDP80, DDE80, DX00, DY00, DX06, DY06,
            DJMJD0, DATE, TIME, UTC, DAT,
            TAI, TT, TUT, UT1, DP80, DE80,
            DPSI, DEPS, EPSA, RP[3][3], RN[3][3], RNPB[3][3],
            EE, GST, RC2TI[3][3], RPOM[3][3],
            RC2IT[3][3], X, Y, S,
            RC2I[3][3], ERA, DP00, DE00, RB[3][3],
            RPB[3][3], V1[3], V2[3], DDP00, DDE00;

    TUDAT_UNUSED_PARAMETER( UT1 );

    //  UTC.
    IY = 2007;
    IM = 4;
    ID = 5;
    IH = 12;
    MIN = 0;
    SEC = 0.0;

    //  Polar motion (arcsec->radians).
    XP = xPInAs * AS2R;
    YP = yPInAs * AS2R;

    //  UT1-UTC (s).
    DUT1 = ut1Correction;

    //  Nutation corrections wrt IAU 1976/1980 (mas->radians).
    DDP80 = -55.0655 * AS2R/1000.0;
    DDE80 =  -6.3580 * AS2R/1000.0;

    //  CIP offsets wrt IAU 2000A (mas->radians).
    DX00 =  0.1725 * AS2R/1000.0;
    DY00 = -0.2650 * AS2R/1000.0;

    //  CIP offsets wrt IAU 2006/2000A (mas->radians).
    DX06 =  dXInMas * AS2R/1000.0;
    DY06 =  dYInMas * AS2R/1000.0;

    //  TT (MJD).
    iauCal2jd ( IY, IM, ID, &DJMJD0, &DATE );
    TIME = ( 60.0*(60.0*double(IH) + double(MIN)) + SEC ) / 86400.0;
    UTC = DATE + TIME;
    iauDat ( IY, IM, ID, TIME, &DAT );
    TAI = UTC + DAT/86400.0;
    TT = TAI + 32.184/86400.0;

    //  UT1.
    TUT = TIME + DUT1/86400.0;
    UT1 = DATE + TUT;

    switch( calculationCase )
    {

    //  =============
    //  IAU 1976/1980
    //  =============

    case 0:
        iauPmat76 ( DJMJD0, TT, RP );

        //  IAU 1980 nutation.
        iauNut80 ( DJMJD0, TT, &DP80, &DE80 );

        //  Add adjustments: frame bias, precession-rates, geophysical.
        DPSI = DP80 + DDP80;
        DEPS = DE80 + DDE80;

        //  Mean obliquity.
        EPSA = iauObl80 ( DJMJD0, TT );

        //  Build the rotation matrix.
        iauNumat( EPSA, DPSI, DEPS, RN );

        //  Combine the matrices:  PN = N x P.
        iauRxr ( RN, RP, RNPB );

        //  Equation of the equinoxes, including nutation correction.
        EE = iauEqeq94 ( DJMJD0, TT ) + DDP80 * cos( EPSA );

        //  Greenwich apparent sidereal time (IAU 1982/1994).
        GST = iauAnp ( iauGmst82( DJMJD0+DATE, TUT ) + EE );

        //  Form celestial-terrestrial matrix (no polar motion yet).
        iauCr ( RNPB, RC2TI );
        iauRz ( GST, RC2TI );

        //  Polar motion matrix (TIRS->ITRS, IERS 1996).
        iauIr ( RPOM );
        iauRx ( -YP, RPOM );
        iauRy ( -XP, RPOM );

        //  Form celestial-terrestrial matrix (including polar motion).
        iauRxr ( RPOM, RC2TI, RC2IT );

        break;

        //  ====================
        //  IAU 2000A, CIO based
        //  ====================
    case 1:
        //  CIP and CIO, IAU 2000A.
        iauXys00a( DJMJD0, TT, &X, &Y, &S );

        //  Add CIP corrections.
        X = X + DX00;
        Y = Y + DY00;

        //  GCRS to CIRS matrix.
        iauC2ixys ( X, Y, S, RC2I );

        //  Earth rotation angle.
        ERA = iauEra00 ( DJMJD0+DATE, TUT );

        //  Form celestial-terrestrial matrix (no polar motion yet).
        iauCr ( RC2I, RC2TI );
        iauRz ( ERA, RC2TI );

        //  Polar motion matrix (TIRS->ITRS, IERS 2003).
        iauPom00 ( XP, YP, iauSp00(DJMJD0,TT), RPOM );

        //  Form celestial-terrestrial matrix (including polar motion).
        iauRxr ( RPOM, RC2TI, RC2IT );

        break;

        //  ========================
        //  IAU 2000A, equinox based
        //  ========================
    case 2:
        //  Nutation, IAU 2000A.
        iauNut00a ( DJMJD0, TT, &DP00, &DE00 );

        //  Precession-nutation quantities, IAU 2000.
        iauPn00 ( DJMJD0, TT, DP00, DE00,
                  &EPSA, RB, RP, RPB, RN, RNPB );

        //  Transform dX,dY corrections from GCRS to mean of date.
        V1[0] = DX00;
        V1[1] = DY00;
        V1[2] = 0.0;
        iauRxp ( RNPB, V1, V2 );
        DDP00 = V2[0] / sin( EPSA );
        DDE00 = V2[1];

        //  Corrected nutation.
        DPSI = DP00 + DDP00;
        DEPS = DE00 + DDE00;

        //  Build the rotation matrix.
        iauNumat ( EPSA, DPSI, DEPS, RN );

        //  Combine the matrices:  PN = N x P.
        iauRxr ( RN, RPB, RNPB );

        //  Greenwich apparent sidereal time (IAU 2000).
        GST = iauAnp ( iauGmst00 ( DJMJD0+DATE, TUT, DJMJD0, TT )
                       + iauEe00 ( DJMJD0, TT, EPSA, DPSI ) );

        //  Form celestial-terrestrial matrix (no polar motion yet).
        iauCr ( RNPB, RC2TI );
        iauRz ( GST, RC2TI );

        //  Polar motion matrix (TIRS->ITRS, IERS 2003).
        iauPom00 ( XP, YP, iauSp00(DJMJD0,TT), RPOM );

        //  Form celestial-terrestrial matrix (including polar motion).
        iauRxr ( RPOM, RC2TI, RC2IT );

        break;
        //  =========================
        //  IAU 2006/2000A, CIO based
        //  =========================

    case 3:
        //  CIP and CIO, IAU 2006/2000A.
        iauXys06a ( DJMJD0, TT, &X, &Y, &S );

        //  Add CIP corrections.
        X = X + DX06;
        Y = Y + DY06;

        //std::cout << "Corrections Nut: " << std::setprecision( 16 ) << TT << " " << X << " " << Y << " " << S << " " << DX06 << " " << DY06 << std::endl;

        //  GCRS to CIRS matrix.
        iauC2ixys ( X, Y, S, RC2I );

        //  Earth rotation angle.
        ERA = iauEra00 ( DJMJD0+DATE, TUT );

        //std::cout << "Corrections ERA: " << std::setprecision( 16 ) << ERA << " " << DJMJD0+DATE << " " << TUT << std::endl;

        //  Form celestial-terrestrial matrix (no polar motion yet).
        iauCr ( RC2I, RC2TI );
        iauRz ( ERA, RC2TI );

        //  Polar motion matrix (TIRS->ITRS, IERS 2003).
        iauPom00 ( XP, YP, iauSp00(DJMJD0,TT), RPOM );

        //  Form celestial-terrestrial matrix (including polar motion).
        iauRxr ( RPOM, RC2TI, RC2IT );

//        std::cout << "Matrix Sofa: " << std::endl << std::setprecision( 16 )
//                  << convertArrayToMatrix( RC2I ) << std::endl;
//        std::cout << "Matrix Sofa: " << std::endl << std::setprecision( 16 )
//                  << convertArrayToMatrix( RC2TI ) << std::endl;
//        std::cout << "Matrix Sofa: " << std::endl << std::setprecision( 16 )
//                  << convertArrayToMatrix( RC2IT ) << std::endl;
        break;

        //  ===========================================
        //  IAU 2006/2000A, CIO based, using X,Y series
        //  ===========================================
    case 4:
        //  CIP and CIO, IAU 2006/2000A.
        iauXy06 ( DJMJD0, TT, &X, &Y );
        S = iauS06 ( DJMJD0, TT, X, Y );

        //  Add CIP corrections.
        X = X + DX06;
        Y = Y + DY06;

        //  GCRS to CIRS matrix.
        iauC2ixys ( X, Y, S, RC2I );

        //  Earth rotation angle.
        ERA = iauEra00 ( DJMJD0+DATE, TUT );

        //  Form celestial-terrestrial matrix (no polar motion yet).
        iauCr ( RC2I, RC2TI );
        iauRz ( ERA, RC2TI );

        //  Polar motion matrix (TIRS->ITRS, IERS 2003).
        iauPom00 ( XP, YP, iauSp00(DJMJD0,TT), RPOM );

        //  Form celestial-terrestrial matrix (including polar motion).
        iauRxr ( RPOM, RC2TI, RC2IT );

        break;
    default:
        std::cerr << "Error, did not recognize case " << calculationCase << " in sofa earth orientation cookbook" << std::endl;
    }

    return convertArrayToMatrix( RC2IT );

}

}

}
