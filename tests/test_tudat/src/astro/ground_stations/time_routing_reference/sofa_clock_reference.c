/* Original Tudat test driver; calls the unmodified SOFA 2023-10-11 library. */
#include <math.h>
#include <stdio.h>
#include <sofa/sofa.h>
static double delta( double day )
{
    return iauDtdb( 2451545.0 + floor( day ), day - floor( day ), 0, 0, 0, 0 );
}
static double peak( double center, int sign )
{
    double a = center - .25, b = center + .25;
    for( int i = 0; i < 100; i++ )
    {
        double c = a + ( b - a ) / 3, d = b - ( b - a ) / 3;
        if( sign * delta( c ) < sign * delta( d ) )
            a = c;
        else
            b = d;
    }
    return ( a + b ) / 2;
}
int main( void )
{
    const double start = 7304.5;
    for( int sign = -1; sign <= 1; sign += 2 )
    {
        double best = start;
        for( int i = 0; i <= 366 * 4; i++ )
        {
            double day = start + i * .25;
            if( sign * delta( day ) > sign * delta( best ) ) best = day;
        }
        double day = peak( best, sign ), dtr = delta( day ), tai1, tai2, utc1, utc2;
        iauTttai( 2451545.0 + floor( day ), day - floor( day ), &tai1, &tai2 );
        iauTaiutc( tai1, tai2, &utc1, &utc2 );
        int y, m, d, clock[ 4 ];
        iauD2dtf( "UTC", 9, utc1, utc2, &y, &m, &d, clock );
        printf( "peak sign=%d TT_seconds=%.17g TT_JD1=%.17g TT_JD2=%.17g TDB_minus_TT=%.17g UTC=%04d-%02d-%02dT%02d:%02d:%02d.%09d\n",
                sign,
                day * 86400.0,
                2451545.0 + floor( day ),
                day - floor( day ),
                dtr,
                y,
                m,
                d,
                clock[ 0 ],
                clock[ 1 ],
                clock[ 2 ],
                clock[ 3 ] );
    }
    double utc1, utc2, tai1, tai2, tt1, tt2;
    for( int i = 0; i < 3; i++ )
    {
        iauDtf2d( "UTC", 2016, 12, 31, 23, 59, 59.5 + i * .5, &utc1, &utc2 );
        iauUtctai( utc1, utc2, &tai1, &tai2 );
        iauTaitt( tai1, tai2, &tt1, &tt2 );
        double dtr = iauDtdb( tt1, tt2, 0, 0, 0, 0 );
        printf( "leap UTC_second=%.1f TT_JD1=%.17g TT_JD2=%.17g TDB_minus_TT=%.17g\n", 59.5 + i * .5, tt1, tt2, dtr );
    }
    return 0;
}
