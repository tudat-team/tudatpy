#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

#include "Tudat/External/SofaInterface/fundamentalArguments.h"

namespace tudat
{

namespace sofa_interface
{

//! Function to calculate mean elongation of the Moon from the Sun
double calculateMeanElongationOfMoonFromSun( const double julianCenturiesSinceJ2000 )
{
    return iauFad03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean argument of latitude of the Moon
double calculateMeanArgumentOfLatitudeOfMoon( const double julianCenturiesSinceJ2000 )
{
    return iauFaf03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean anomaly of the Moon
double calculateMeanAnomalyOfMoon( const double julianCenturiesSinceJ2000 )
{
    return iauFal03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean anomaly of the Sun
double calculateMeanAnomalyOfSun( const double julianCenturiesSinceJ2000 )
{
    return iauFalp03( julianCenturiesSinceJ2000 );
}

//! Function to calculate longitude of the Moon's ascending node
double calculateLongitudeOfMoonsAscendingNode( const double julianCenturiesSinceJ2000 )
{
    return iauFaom03( julianCenturiesSinceJ2000 );
}

//Eigen::Vector3d calculateOceanTideArgumets( const double ephemerisTime )
//{
//    double julianCenturiesSinceJ2000 = basic_astrodynamics::convertSecondsSinceEpochToJulianCenturiesSinceEpoch( ephemerisTime );

//    double meanElongationOfMoonFromSun = calculateMeanElongationOfMoonFromSun( julianCenturiesSinceJ2000 );
//    double meanArgumentOfLatitudeOfMoon = calculateMeanArgumentOfLatitudeOfMoon( julianCenturiesSinceJ2000 );
//    double meanAnomalyOfMoon = calculateMeanAnomalyOfMoon( julianCenturiesSinceJ2000 );
//    double longitudeOfMoonsAscendingNode = calculateLongitudeOfMoonsAscendingNode( julianCenturiesSinceJ2000 );

//    Eigen::Vector3d oceanTideArguments;
//    oceanTideArguments( 0 ) = meanArgumentOfLatitudeOfMoon + longitudeOfMoonsAscendingNode;
//    oceanTideArguments( 1 ) = oceanTideArguments( 0 ) - meanElongationOfMoonFromSun;
//    oceanTideArguments( 2 ) = oceanTideArguments( 0 ) - meanAnomalyOfMoon;

//    return oceanTideArguments;
//}

Eigen::Matrix< double, 5, 1 > calculateDelaunayFundamentalArguments(
        const double ephemerisTime )
{
    double julianCenturiesSinceJ2000 = basic_astrodynamics::convertSecondsSinceEpochToJulianCenturiesSinceEpoch( ephemerisTime );

    return ( Eigen::Matrix< double, 5, 1 >( ) <<
             calculateMeanAnomalyOfMoon( julianCenturiesSinceJ2000 ),
             calculateMeanAnomalyOfSun( julianCenturiesSinceJ2000 ),
             calculateMeanArgumentOfLatitudeOfMoon( julianCenturiesSinceJ2000 ),
             calculateMeanElongationOfMoonFromSun( julianCenturiesSinceJ2000 ),
             calculateLongitudeOfMoonsAscendingNode( julianCenturiesSinceJ2000 ) ).finished( );
}


basic_mathematics::Vector6d  calculateFundamentalArguments(
        const double ephemerisTime )
{
    double utc = convertTTtoUTC( ephemerisTime );

    double greenwichMeanSiderealTime = iauGmst06( basic_astrodynamics::JULIAN_DAY_ON_J2000, utc / physical_constants::JULIAN_DAY,
                                                  basic_astrodynamics::JULIAN_DAY_ON_J2000, ephemerisTime / physical_constants::JULIAN_DAY );

    basic_mathematics::Vector6d fundamentalArguments;
    fundamentalArguments << greenwichMeanSiderealTime + mathematical_constants::PI, calculateDelaunayFundamentalArguments( ephemerisTime );

    return fundamentalArguments;

}

basic_mathematics::Vector6d calculateDoodsonFundamentalArguments(
        const double ephemerisTime, const double terrestrialTime, const double universalTime1, const double ttAndUt1EpochShift )
{
    using mathematical_constants::PI;

    double julianCenturiesSinceJ2000 = basic_astrodynamics::convertSecondsSinceEpochToJulianCenturiesSinceEpoch( ephemerisTime );

    //std::cout<<"Centuries since J2000 "<<std::setprecision( 26 )<<julianCenturiesSinceJ2000 - 0.07995893223819302<<std::endl;
    double greenwichMeanSiderealTime = iauGmst06( ttAndUt1EpochShift, universalTime1 / physical_constants::JULIAN_DAY,
                                                  ttAndUt1EpochShift, terrestrialTime / physical_constants::JULIAN_DAY );

    double meanElongationOfMoonFromSun = calculateMeanElongationOfMoonFromSun( julianCenturiesSinceJ2000 );
    double meanArgumentOfLatitudeOfMoon = calculateMeanArgumentOfLatitudeOfMoon( julianCenturiesSinceJ2000 );
    double meanAnomalyOfMoon = calculateMeanAnomalyOfMoon( julianCenturiesSinceJ2000 );
    double meanAnomalyOfSun = calculateMeanAnomalyOfSun( julianCenturiesSinceJ2000 );
    double longitudeOfMoonsAscendingNode = calculateLongitudeOfMoonsAscendingNode( julianCenturiesSinceJ2000 );

    basic_mathematics::Vector6d doodsonArguments;
    doodsonArguments( 1 ) = meanArgumentOfLatitudeOfMoon + longitudeOfMoonsAscendingNode;
    doodsonArguments( 0 ) = greenwichMeanSiderealTime + PI - doodsonArguments( 1 );
    doodsonArguments( 2 ) = doodsonArguments( 1 ) - meanElongationOfMoonFromSun;
    doodsonArguments( 3 ) = doodsonArguments( 1 ) - meanAnomalyOfMoon;
    doodsonArguments( 5 ) = doodsonArguments( 2 ) - meanAnomalyOfSun;
    doodsonArguments( 4 ) = -longitudeOfMoonsAscendingNode;

    return doodsonArguments;
}

// Aproximations TT = TDB and UT1 = UTC used!
basic_mathematics::Vector6d calculateDoodsonFundamentalArguments( const double ephemerisTime )
{
    double utc = convertTTtoUTC( ephemerisTime );

    return calculateDoodsonFundamentalArguments( ephemerisTime, ephemerisTime,  utc, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

} // namespace sofa_interfaces

} // namespace tudat
