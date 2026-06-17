/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          ITU-R Recommendation P.531-15 (2023), Ionospheric propagation data and prediction methods
 *              required for the design of satellite networks and systems.
 *          B. Nava, P. Coisson, and S.M. Radicella (2008), A new version of the NeQuick ionosphere
 *              electron density model, Journal of Atmospheric and Solar-Terrestrial Physics, 70(15),
 *              pp. 1856-1862, doi:10.1016/j.jastp.2008.01.015
 *          G. Di Giovanni and S.M. Radicella (1990), An analytical model of the electron density
 *              profile in the ionosphere, Advances in Space Research, 10(11), pp. 27-30.
 *          S.M. Radicella and R. Leitinger (2001), The evolution of the DGR approach to model
 *              electron density profiles, Advances in Space Research, 27(1), pp. 35-40.
 *          R. Leitinger, M.L. Zhang, and S.M. Radicella (2005), An improved bottomside for the
 *              ionospheric electron density model NeQuick, Annals of Geophysics, 48(3).
 *          CCIR (1967), Atlas of Ionospheric Characteristics, Report 340-6, International
 *              Telecommunication Union, Geneva.
 *
 *    Original Fortran reference implementation: NeQuick2 P.531-15 (release 2.1.0, October 2011)
 *          Authors: B. Nava, P. Coisson, J. Wagner, Y. Migoya Orue (ICTP Trieste)
 *          Validation: R. Prieto-Cerdeira, R. Orus-Perez (ESA/ESTEC)
 *    C++ transliteration for Tudat.
 */

#include "tudat/astro/basic_astro/neQuick2Model.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/quadrature/gaussianQuadrature.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace tudat
{
namespace environment
{

namespace
{
const double DR = 1.74532925199433e-2;  // degrees to radians
const double RD = 57.29577951308232;    // radians to degrees
const double PI12 = 0.261799387788149;  // pi/12

// foF2 harmonic order table (QF) and M(3000)F2 harmonic order table (QM)
const int QF[ 9 ] = { 11, 11, 8, 4, 1, 0, 0, 0, 0 };
const int QM[ 7 ] = { 6, 7, 5, 2, 1, 0, 0 };
}  // namespace

// ============================================================================
// Static utility functions
// ============================================================================

double NeQuick2Model::safeExp( double a )
{
    if( a > 80.0 ) return 5.5406e34;
    if( a < -80.0 ) return 1.8049e-35;
    return std::exp( a );
}

double NeQuick2Model::djoin( double f1, double f2, double alpha, double x )
{
    double ee = safeExp( alpha * x );
    return ( f1 * ee + f2 ) / ( ee + 1.0 );
}

double NeQuick2Model::fluxToR12( double flx )
{
    return std::sqrt( 167273.0 + ( flx - 63.7 ) * 1123.6 ) - 408.99;
}

double NeQuick2Model::finter3( const double z[ 4 ], double x )
{
    double dx = x + x - 1.0;
    if( std::abs( dx + 1.0 ) < 1.0e-10 )
    {
        return z[ 1 ];
    }

    double g1 = z[ 2 ] + z[ 1 ];
    double g2 = z[ 2 ] - z[ 1 ];
    double g3 = z[ 3 ] + z[ 0 ];
    double g4 = ( z[ 3 ] - z[ 0 ] ) / 3.0;
    double a0 = 9.0 * g1 - g3;
    double a1 = 9.0 * g2 - g4;
    double a2 = g3 - g1;
    double a3 = g4 - g2;

    double zi = 0.0;
    zi = zi * dx + a3;
    zi = zi * dx + a2;
    zi = zi * dx + a1;
    zi = zi * dx + a0;
    return zi / 16.0;
}

void NeQuick2Model::solarDeclination( int month, double ut, double& sinDelta, double& cosDelta )
{
    double doy = month * 30.5 - 15.0;
    double t = doy + ( 18.0 - ut ) / 24.0;
    double amrad = ( 0.9856 * t - 3.289 ) * DR;
    double aLrad = amrad + ( 1.916 * std::sin( amrad ) + 0.020 * std::sin( 2.0 * amrad ) + 282.634 ) * DR;
    sinDelta = 0.39782 * std::sin( aLrad );
    cosDelta = std::sqrt( 1.0 - sinDelta * sinDelta );
}

double NeQuick2Model::computeHmF2( double foE, double foF2, double M3000 )
{
    double sqM = M3000 * M3000;
    double MF = M3000 * std::sqrt( ( 0.0196 * sqM + 1.0 ) / ( 1.2967 * sqM - 1.0 ) );
    double dM;
    if( foE >= 1.0e-30 )
    {
        double ratio = foF2 / foE;
        ratio = djoin( ratio, 1.75, 20.0, ratio - 1.75 );
        dM = 0.253 / ( ratio - 1.215 ) - 0.012;
    }
    else
    {
        dM = -0.012;
    }
    return 1490.0 * MF / ( M3000 + dM ) - 176.0;
}

void NeQuick2Model::timeToMonthAndUT( double time, int& month, double& ut )
{
    // Convert seconds since J2000 to Julian Day
    double julianDay = basic_astrodynamics::convertSecondsSinceEpochToJulianDay( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Extract calendar date (day, month, year)
    int day, year;
    basic_astrodynamics::convertJulianDayToCalendarDate( julianDay, day, month, year );

    // Extract UT as fractional hours from the Julian Day
    // JD 0.0 = noon, so JD fraction 0.5 = midnight (0h UT)
    double jdFraction = julianDay - std::floor( julianDay );
    ut = ( jdFraction + 0.5 ) * 24.0;
    if( ut >= 24.0 ) ut -= 24.0;
}

// ============================================================================
// Constructor
// ============================================================================

NeQuick2Model::NeQuick2Model( const input_output::CcirData& ccirData,
                              const Eigen::MatrixXd& modipGrid,
                              std::function< double( double ) > solarFluxFunction,
                              bool inputIsEffectiveIonizationLevel ):
    ccirData_( ccirData ), modipGrid_( modipGrid ), solarFluxFunction_( solarFluxFunction )
{}

// ============================================================================
// MODIP computation
// ============================================================================

double NeQuick2Model::computeModip( double latitudeDeg, double longitudeDeg ) const
{
    // Transliteration of Fortran amodip function
    const int lngp = 180;
    const double dlatp = 1.0;
    const double dlngp = 2.0;

    double along1 = std::fmod( longitudeDeg + 360.0, 360.0 );
    double dlng1 = ( along1 + 180.0 ) / dlngp;
    dlng1 = dlng1 - std::floor( dlng1 );
    int j1 = static_cast< int >( ( along1 + 180.0 ) / dlngp ) - 2;
    if( j1 < 0 ) j1 += lngp;
    if( j1 > lngp - 3 ) j1 -= lngp;

    double a = ( latitudeDeg + 90.0 ) / dlatp + 1.0;
    int iLat = static_cast< int >( a ) - 2;
    a = a - static_cast< double >( iLat + 2 );

    double z[ 4 ];
    for( int k = 0; k < 4; ++k )
    {
        double z1[ 4 ];
        for( int j = 0; j < 4; ++j )
        {
            // Fortran: pmodip(i+j, j1+k) with 0-based offset +1 for our storage
            z1[ j ] = modipGrid_( iLat + j + 1, j1 + k + 1 );
        }
        z[ k ] = finter3( z1, a );
    }
    return finter3( z, dlng1 );
}

// ============================================================================
// CCIR coefficient evaluation (gamma1 and cciri)
// ============================================================================

double NeQuick2Model::gamma1( double modip,
                              double latDeg,
                              double lonDeg,
                              double hour,
                              int iharm,
                              const int* nq,
                              int k1,
                              int m,
                              int mm,
                              int m3,
                              const double* sfe ) const
{
    double c[ 12 ], s[ 12 ];
    double coef[ 100 ];

    double hou = ( 15.0 * hour - 180.0 ) * DR;
    s[ 0 ] = std::sin( hou );
    c[ 0 ] = std::cos( hou );
    for( int i = 1; i < iharm; ++i )
    {
        c[ i ] = c[ 0 ] * c[ i - 1 ] - s[ 0 ] * s[ i - 1 ];
        s[ i ] = c[ 0 ] * s[ i - 1 ] + s[ 0 ] * c[ i - 1 ];
    }

    for( int i = 0; i < m; ++i )
    {
        int mi = i * mm;
        coef[ i ] = sfe[ mi ];
        for( int j = 0; j < iharm; ++j )
        {
            coef[ i ] += sfe[ mi + 2 * j + 1 ] * s[ j ] + sfe[ mi + 2 * j + 2 ] * c[ j ];
        }
    }

    double sum = coef[ 0 ];
    double ss = std::sin( modip * DR );
    double s3 = ss;
    double xsinx[ 14 ];
    xsinx[ 0 ] = 1.0;
    int index = nq[ 0 ];
    bool numok = true;
    for( int j = 0; j < index; ++j )
    {
        numok = std::abs( ss ) >= 1.0e-30;
        if( numok )
        {
            sum += coef[ 1 + j ] * ss;
            xsinx[ j + 1 ] = ss;
            ss *= s3;
        }
        else
        {
            xsinx[ j + 1 ] = 0.0;
        }
    }
    xsinx[ nq[ 0 ] + 1 ] = numok ? ss : 0.0;

    int np = nq[ 0 ] + 1;
    ss = std::cos( latDeg * DR );
    s3 = ss;
    for( int j = 1; j < k1; ++j )
    {
        double s0 = lonDeg * j * DR;
        double s1 = std::cos( s0 );
        double s2 = std::sin( s0 );
        index = nq[ j ] + 1;
        for( int L = 0; L < index; ++L )
        {
            sum += coef[ np ] * xsinx[ L ] * ss * s1;
            np++;
            sum += coef[ np ] * xsinx[ L ] * ss * s2;
            np++;
        }
        ss *= s3;
    }

    return sum;
}

void NeQuick2Model::computeFoF2AndM3000( double modip,
                                         int month,
                                         double ut,
                                         double R12,
                                         double latDeg,
                                         double lonDeg,
                                         double& foF2,
                                         double& M3000 ) const
{
    // Interpolate CCIR coefficients between low and high solar activity
    int monthIdx = month - 1;  // 0-based
    double RR2 = R12 / 100.0;
    double RR1 = 1.0 - RR2;

    // Build foF2 coefficient array (988 values = 13 * 76)
    std::vector< double > FF0( 988 );
    for( int i = 0; i < 76; ++i )
    {
        for( int j = 0; j < 13; ++j )
        {
            int k = j + 13 * i;
            FF0[ k ] =
                    ccirData_.foF2Coefficients[ monthIdx ][ 0 ]( i, j ) * RR1 + ccirData_.foF2Coefficients[ monthIdx ][ 1 ]( i, j ) * RR2;
        }
    }

    // Build M(3000)F2 coefficient array (441 values = 9 * 49)
    std::vector< double > xm0( 441 );
    for( int i = 0; i < 49; ++i )
    {
        for( int j = 0; j < 9; ++j )
        {
            int k = j + 9 * i;
            xm0[ k ] = ccirData_.m3000F2Coefficients[ monthIdx ][ 0 ]( i, j ) * RR1 +
                    ccirData_.m3000F2Coefficients[ monthIdx ][ 1 ]( i, j ) * RR2;
        }
    }

    foF2 = gamma1( modip, latDeg, lonDeg, ut, 6, QF, 9, 76, 13, 988, FF0.data( ) );
    M3000 = gamma1( modip, latDeg, lonDeg, ut, 4, QM, 7, 49, 9, 441, xm0.data( ) );
}

// ============================================================================
// E-layer and F1-layer critical frequencies
// ============================================================================

void NeQuick2Model::computeELayerAndF1( double latDeg, int month, double flx, double chi, double foF2, double& foE, double& foF1 )
{
    const double chi0 = 86.23292796211615;

    double chin = djoin( 90.0 - 0.24 * safeExp( 20.0 - 0.20 * chi ), chi, 12.0, chi - chi0 );
    double ee = safeExp( 0.3 * latDeg );
    double seas = 1.0 - 2.0 / ( ee + 1.0 );

    double sfac;
    if( month == 1 || month == 2 || month == 11 || month == 12 )
    {
        sfac = ( 1.112 + 0.019 * seas ) * std::pow( flx, 0.25 );
    }
    else if( month == 3 || month == 4 || month == 9 || month == 10 )
    {
        sfac = 1.112 * std::pow( flx, 0.25 );
    }
    else
    {
        sfac = ( 1.112 - 0.019 * seas ) * std::pow( flx, 0.25 );
    }

    double fa = sfac * std::pow( std::cos( chin * DR ), 0.3 );
    foE = std::sqrt( fa * fa + 0.49 );

    foF1 = djoin( 1.4 * foE, 0.0, 1000.0, foE - 2.0 );
    foF1 = djoin( 0.0, foF1, 1000.0, foE - foF1 );
    foF1 = djoin( foF1, 0.85 * foF1, 60.0, 0.85 * foF2 - foF1 );
    if( foF1 < 1.0e-6 ) foF1 = 0.0;
}

// ============================================================================
// Layer parameter preparation (prepmdgr)
// ============================================================================

void NeQuick2Model::prepareLayerParameters( int month,
                                            double R12,
                                            double foF2,
                                            double foF1,
                                            double foE,
                                            double M3000,
                                            NeQuick2LayerParameters& params )
{
    auto FNe = []( double x ) { return 0.124 * x * x; };
    auto FEpst = []( double amp, double hm, double B, double h ) {
        double ex = safeExp( ( h - hm ) / B );
        return amp * ex / ( ( 1.0 + ex ) * ( 1.0 + ex ) );
    };

    params.foF2 = foF2;
    params.foF1 = foF1;
    params.foE = foE;
    params.M3000 = M3000;

    params.NmF2 = FNe( foF2 );
    params.NmF1 = FNe( foF1 );
    params.NmE = FNe( foE );

    params.hmF2 = computeHmF2( foE, foF2, M3000 );
    params.hmE = 120.0;
    params.hmF1 = 0.5 * ( params.hmF2 + params.hmE );

    double dNdHmx = std::exp( -3.467 + 1.714 * std::log( foF2 ) + 2.02 * std::log( M3000 ) );
    double B2bot = 38.5 * params.NmF2 / dNdHmx;
    double hF2F1E = params.hmF2 - params.hmF1;
    double B1top = 0.3 * hF2F1E;
    double B1bot = 0.5 * hF2F1E;
    double Betop = B1bot;
    if( Betop < 7.0 ) Betop = 7.0;
    double Bebot = 5.0;

    params.aep[ 0 ] = 4.0 * params.NmF2;
    params.aep[ 1 ] = 4.0 * params.NmF1;
    params.aep[ 2 ] = 4.0 * params.NmE;

    if( foF1 <= 0.5 )
    {
        params.aep[ 1 ] = 0.0;
        params.aep[ 2 ] = 4.0 * ( params.NmE - FEpst( params.aep[ 0 ], params.hmF2, B2bot, params.hmE ) );
    }
    else
    {
        for( int iter = 0; iter < 5; ++iter )
        {
            params.aep[ 1 ] = 4.0 *
                    ( params.NmF1 - FEpst( params.aep[ 0 ], params.hmF2, B2bot, params.hmF1 ) -
                      FEpst( params.aep[ 2 ], params.hmE, Betop, params.hmF1 ) );
            params.aep[ 1 ] = djoin( params.aep[ 1 ], 0.8 * params.NmF1, 1.0, params.aep[ 1 ] - 0.8 * params.NmF1 );
            params.aep[ 2 ] = 4.0 *
                    ( params.NmE - FEpst( params.aep[ 1 ], params.hmF1, B1bot, params.hmE ) -
                      FEpst( params.aep[ 0 ], params.hmF2, B2bot, params.hmE ) );
        }
        if( params.aep[ 1 ] < 0.0 )
        {
            params.aep[ 1 ] = 0.0;
            params.aep[ 2 ] = 4.0 * ( params.NmE - FEpst( params.aep[ 0 ], params.hmF2, B2bot, params.hmE ) );
        }
    }
    params.aep[ 2 ] = djoin( params.aep[ 2 ], 0.005, 60.0, params.aep[ 2 ] - 0.005 );

    double b2k = 3.22 - 0.0538 * foF2 - 0.00664 * params.hmF2 + 0.113 * params.hmF2 / B2bot + 0.00257 * R12;
    b2k = djoin( b2k, 1.0, 2.0, b2k - 1.0 );
    double B2top = b2k * B2bot;

    params.bb[ 0 ] = Bebot;
    params.bb[ 1 ] = Betop;
    params.bb[ 2 ] = B1bot;
    params.bb[ 3 ] = B1top;
    params.bb[ 4 ] = B2bot;
    params.bb[ 5 ] = B2top;
}

// ============================================================================
// Electron density computation
// ============================================================================

double NeQuick2Model::bottomsideElectronDensity( const NeQuick2LayerParameters& params, double heightKm )
{
    // Transliteration of Fortran NeMdGR
    const double f1 = 10.0, f2 = 1.0;
    const double h0 = 90.0, Hd = 5.0;

    double hm[ 3 ] = { params.hmF2, params.hmF1, params.hmE };
    double B[ 3 ];
    B[ 0 ] = params.bb[ 4 ];  // B2bot
    B[ 1 ] = params.bb[ 2 ];  // B1bot
    B[ 2 ] = params.bb[ 0 ];  // Bebot

    if( heightKm > params.hmE ) B[ 2 ] = params.bb[ 1 ];   // Betop
    if( heightKm > params.hmF1 ) B[ 1 ] = params.bb[ 3 ];  // B1top

    double sum = 0.0;
    for( int jj = 0; jj < 3; ++jj )
    {
        double arg = ( heightKm - hm[ jj ] ) / B[ jj ];
        if( jj > 0 )
        {
            double d = std::abs( heightKm - hm[ 0 ] );
            arg *= std::exp( f1 / ( 1.0 + f2 * d ) );
        }
        if( heightKm < h0 )
        {
            arg *= ( Hd + h0 - heightKm ) / Hd;
        }
        if( std::abs( arg ) > 25.0 )
        {
            // s0 = 0
        }
        else
        {
            double ee = std::exp( arg );
            double s0 = params.aep[ jj ] * ee / ( ( 1.0 + ee ) * ( 1.0 + ee ) );
            sum += s0;
        }
    }
    return sum * 1.0e11;
}

double NeQuick2Model::topsideElectronDensity( double heightKm, double NmaxAtPeak, double hmaxKm, double H0 )
{
    // Transliteration of Fortran topq
    const double g = 0.125;
    const double rfac = 100.0;

    double dh = heightKm - hmaxKm;
    double gh = g * dh;
    double z = dh / ( H0 * ( 1.0 + rfac * gh / ( rfac * H0 + gh ) ) );
    double ee = safeExp( z );
    double ep;
    if( ee > 1.0e11 )
    {
        ep = 4.0 / ee;
    }
    else
    {
        ep = 4.0 * ee / ( ( 1.0 + ee ) * ( 1.0 + ee ) );
    }
    return NmaxAtPeak * ep;
}

double NeQuick2Model::computeElectronDensityFromParams( double heightKm, const NeQuick2LayerParameters& params )
{
    if( heightKm > params.hmF2 )
    {
        double NmaxAtPeak = bottomsideElectronDensity( params, params.hmF2 );
        return topsideElectronDensity( heightKm, NmaxAtPeak, params.hmF2, params.bb[ 5 ] );
    }
    return bottomsideElectronDensity( params, heightKm );
}

// ============================================================================
// Public API
// ============================================================================

NeQuick2LayerParameters NeQuick2Model::computeLayerParameters( double latitudeDeg,
                                                               double longitudeDeg,
                                                               int month,
                                                               double flx,
                                                               double ut ) const
{
    // Clamp flux
    double flx1 = std::max( 63.0, std::min( flx, 193.0 ) );
    double R12 = fluxToR12( flx1 );
    double ut1 = std::fmod( ut + 24.0, 24.0 );

    // Solar declination
    double sinDelta, cosDelta;
    solarDeclination( month, ut1, sinDelta, cosDelta );

    // Local time and solar zenith angle
    double xlt = std::fmod( ut1 + longitudeDeg / 15.0 + 24.0, 24.0 );
    double modip = computeModip( latitudeDeg, longitudeDeg );
    double cchi = std::sin( latitudeDeg * DR ) * sinDelta - std::cos( latitudeDeg * DR ) * cosDelta * std::cos( PI12 * xlt );
    double chi = std::atan2( std::sqrt( 1.0 - cchi * cchi ), cchi ) * RD;

    // CCIR lookup
    double foF2, M3000;
    computeFoF2AndM3000( modip, month, ut1, R12, latitudeDeg, longitudeDeg, foF2, M3000 );

    // E-layer and F1-layer
    double foE, foF1;
    computeELayerAndF1( latitudeDeg, month, flx1, chi, foF2, foE, foF1 );

    // Prepare layer parameters
    NeQuick2LayerParameters params;
    prepareLayerParameters( month, R12, foF2, foF1, foE, M3000, params );

    return params;
}

double NeQuick2Model::computeElectronDensity( double heightKm, double latitudeDeg, double longitudeDeg, int month, double flx, double ut )
        const
{
    NeQuick2LayerParameters params = computeLayerParameters( latitudeDeg, longitudeDeg, month, flx, ut );
    return computeElectronDensityFromParams( heightKm, params );
}

double NeQuick2Model::computeElectronDensityAtTime( double heightKm, double latitudeDeg, double longitudeDeg, double time ) const
{
    int month;
    double ut;
    timeToMonthAndUT( time, month, ut );
    double flx = solarFluxFunction_( time );
    return computeElectronDensity( heightKm, latitudeDeg, longitudeDeg, month, flx, ut );
}

double NeQuick2Model::computeElectronDensityRescaled( double heightKm,
                                                      double latitudeDeg,
                                                      double longitudeDeg,
                                                      int month,
                                                      double flx,
                                                      double ut,
                                                      double nmF2ScaleFactor ) const
{
    NeQuick2LayerParameters params = computeLayerParameters( latitudeDeg, longitudeDeg, month, flx, ut );

    // Rescale: scale foF2 by sqrt(k), which scales NmF2 by k (since NmF2 = 0.124 * foF2^2)
    double scaledFoF2 = params.foF2 * std::sqrt( nmF2ScaleFactor );

    // Recompute layer parameters with rescaled foF2
    prepareLayerParameters(
            month, fluxToR12( std::max( 63.0, std::min( flx, 193.0 ) ) ), scaledFoF2, params.foF1, params.foE, params.M3000, params );

    return computeElectronDensityFromParams( heightKm, params );
}

double NeQuick2Model::getVerticalTotalElectronContent( double latitudeDeg, double longitudeDeg, double time )
{
    int month;
    double ut;
    timeToMonthAndUT( time, month, ut );
    double flx = solarFluxFunction_( time );
    flx = std::max( 63.0, std::min( flx, 193.0 ) );

    NeQuick2LayerParameters params = computeLayerParameters( latitudeDeg, longitudeDeg, month, flx, ut );

    // Integrate electron density vertically from 0 to 20000 km
    std::function< double( double ) > electronDensityProfile = [ &params ]( double h ) {
        return computeElectronDensityFromParams( h, params );
    };

    numerical_quadrature::GaussianQuadrature< double, double > quadrature( electronDensityProfile, 0.0, 20000.0, 50 );

    // Result in el/m^2, multiply by 1e3 to convert km integration to meters, then to TECU (1e16 el/m^2 = 1 TECU)
    double vtecElPerM2 = quadrature.getQuadrature( ) * 1.0e3;  // km -> m conversion
    return vtecElPerM2 / 1.0e16;                               // Convert to TECU
}

// ============================================================================
// IonexConstrainedNeQuick2Model
// ============================================================================

IonexConstrainedNeQuick2Model::IonexConstrainedNeQuick2Model( std::shared_ptr< NeQuick2Model > neQuick2Model,
                                                              std::shared_ptr< IonosphereModel > ionexModel ):
    neQuick2Model_( neQuick2Model ), ionexModel_( ionexModel )
{}

double IonexConstrainedNeQuick2Model::computeRescaledSlantTec( const Eigen::Vector3d& txPositionEarthFixed,
                                                               const Eigen::Vector3d& rxPositionEarthFixed,
                                                               double time,
                                                               double earthRadius,
                                                               int quadratureOrder ) const
{
    int month;
    double ut;
    NeQuick2Model::timeToMonthAndUT( time, month, ut );
    double flx = neQuick2Model_->getSolarFluxFunction( )( time );
    flx = std::max( 63.0, std::min( flx, 193.0 ) );

    // Compute IPP: midpoint of the ray projected to the surface
    Eigen::Vector3d midpoint = 0.5 * ( txPositionEarthFixed + rxPositionEarthFixed );
    double midRadius = midpoint.norm( );
    double midLatRad = std::asin( midpoint.z( ) / midRadius );
    double midLonRad = std::atan2( midpoint.y( ), midpoint.x( ) );
    double ippLatDeg = midLatRad * RD;
    double ippLonDeg = midLonRad * RD;

    // Get IONEX VTEC at IPP [TECU]
    double vtecIonex = ionexModel_->getVerticalTotalElectronContent( ippLatDeg, ippLonDeg, time );

    // Get NeQuick VTEC at IPP [TECU]
    double vtecNeQuick = neQuick2Model_->getVerticalTotalElectronContent( ippLatDeg, ippLonDeg, time );

    // Compute rescaling factor
    double k = 1.0;
    if( vtecNeQuick > 0.01 )
    {
        k = vtecIonex / vtecNeQuick;
        if( k < 0.0 ) k = 0.0;  // safety clamp
    }

    // Integrate rescaled N_e along the ray using Gauss-Legendre quadrature
    Eigen::Vector3d rayDir = rxPositionEarthFixed - txPositionEarthFixed;
    double pathLength = rayDir.norm( );

    std::function< double( double ) > electronDensityAlongRay =
            [ this, &txPositionEarthFixed, &rayDir, earthRadius, month, flx, ut, k ]( double t ) {
                Eigen::Vector3d pos = txPositionEarthFixed + t * rayDir;
                double radius = pos.norm( );
                double heightKm = ( radius - earthRadius ) / 1.0e3;

                if( heightKm < 0.0 ) return 0.0;

                double latRad = std::asin( pos.z( ) / radius );
                double lonRad = std::atan2( pos.y( ), pos.x( ) );
                double latDeg = latRad * RD;
                double lonDeg = lonRad * RD;

                return neQuick2Model_->computeElectronDensityRescaled( heightKm, latDeg, lonDeg, month, flx, ut, k );
            };

    numerical_quadrature::GaussianQuadrature< double, double > quadrature( electronDensityAlongRay, 0.0, 1.0, quadratureOrder );

    // STEC = integral * path_length (path_length converts fractional to meters)
    return quadrature.getQuadrature( ) * pathLength;
}

}  // namespace environment
}  // namespace tudat
