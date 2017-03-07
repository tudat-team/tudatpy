#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "External/SofaInterface/sofaTimeConversions.h"

int main( )
{

    int latnd, latnm, lonwd, lonwm, j, iy, mo, id, ih, im, ihmsf[4];
    double slatn, slonw, hm, elon, phi, xyz[2], u, v, sec,
            utc1, utc2, dut, ut11, ut12, ut, tai1, tai2, tt1, tt2,
            tcg1, tcg2, dtr, tdb1, tdb2, tcb1, tcb2;
    /* Site terrestrial coordinates (WGS84). */
    latnd = 19;
    latnm = 28;
    slatn = 52.5;
    lonwd = 155;
    lonwm = 55;
    slonw = 59.6;
    hm = 0.0;
    /* Transform to geocentric. */
    j = iauAf2a ( '+', latnd, latnm, slatn, &phi );
    if ( j ) return 1;
    j = iauAf2a ( '-', lonwd, lonwm, slonw, &elon );
    if ( j ) return 1;
    j = iauGd2gc ( 1, elon, phi, hm, xyz );
    if ( j ) return 1;
    u = sqrt ( xyz[0]*xyz[0] + xyz[1]*xyz[1] );
    v = xyz[2];
    std::cout<<std::setprecision( 16 )<<"xyz "<<xyz[ 0 ]<<" "<<xyz[ 1 ]<<" "<<xyz[ 2 ]<<std::endl;
    /* UTC date and time. */
    iy = 2006;
    mo = 1;
    id = 15;
    ih = 21;
    im = 24;
    sec = 37.5;
    /* Transform into internal format. */
    j = iauDtf2d ( "UTC", iy, mo, id, ih, im, sec, &utc1, &utc2 );
    if ( j ) return 1;
    /* UT1-UTC (s, from IERS). */
    dut = 0.3341;
    /* UTC -> UT1. */
    j = iauUtcut1 ( utc1, utc2, dut, &ut11, &ut12 );
    if ( j ) return 1;
    /* Extract fraction for TDB-TT calculation, later. */
    ut = fmod ( fmod(ut11,1.0) + fmod(ut12,1.0), 1.0 );
            /* UTC -> TAI -> TT -> TCG. */
            j = iauUtctai ( utc1, utc2, &tai1, &tai2 );
    if ( j ) return 1;
    j = iauTaitt ( tai1, tai2, &tt1, &tt2 );
    if ( j ) return 1;
    j = iauTttcg ( tt1, tt2, &tcg1, &tcg2 );
    if ( j ) return 1;
    /* TDB-TT (using TT as a substitute for TDB). */
    dtr = iauDtdb ( tt1, tt2, ut, elon, u, v );
    std::cout<<tt1-2.45375e+006<<" "<<tt2<<" "<<ut<<" "<<elon<<" "<<u<<" "<<v<<" "<<dtr<<std::endl;

    /* TT -> TDB -> TCB. */
    j = iauTttdb ( tt1, tt2, dtr, &tdb1, &tdb2 );
    if ( j ) return 1;
    j = iauTdbtcb ( tdb1, tdb2, &tcb1, &tcb2 );
    if ( j ) return 1;
    /* Report. */
    j = iauD2dtf ( "UTC", 6, utc1, utc2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "UTC%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
    j = iauD2dtf ( "ut1", 6, ut11, ut12, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "UT1%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
    j = iauD2dtf ( "tai", 6, tai1, tai2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "TAI%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
    j = iauD2dtf ( "tt", 6, tt1, tt2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "TT %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
    j = iauD2dtf ( "tcg", 6, tcg1, tcg2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "TCG%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
    j = iauD2dtf ( "tdb", 6, tdb1, tdb2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "TDB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
            j = iauD2dtf ( "tcb", 6, tcb1, tcb2, &iy, &mo, &id, ihmsf );
    if ( j ) return 1;
    printf ( "TCB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n",
             iy, mo, id , ihmsf[0], ihmsf[1], ihmsf[2], ihmsf[3] );
}
