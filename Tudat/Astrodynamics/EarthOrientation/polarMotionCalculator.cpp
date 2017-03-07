#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to calculate the mean Celestial Intermediate Pole (CIP) position in the ITRS.
Eigen::Vector2d calculateMeanCipPositionInItrs( const double julianYearsSinceJ2000 )
{
    using namespace tudat::unit_conversions;

    Eigen::Vector2d meanCipPositionInItrs = Eigen::Vector2d::Zero( );
    static double milliAsToRadians = convertArcSecondsToRadians< double >( 1.0E-3 );

    // If date is before 2010, use measured values.
    if( julianYearsSinceJ2000 < 10.0 )
    {
        meanCipPositionInItrs.x( ) = 55.974 +  ( 1.8243 + julianYearsSinceJ2000 *
                                                 ( 0.18413 + julianYearsSinceJ2000 * 0.007024 ) ) * julianYearsSinceJ2000;
        meanCipPositionInItrs.y( ) = 346.346 +  ( 1.7896 + julianYearsSinceJ2000 *
                                                  ( -0.10729 + julianYearsSinceJ2000 * -0.000908 ) ) * julianYearsSinceJ2000;

    }
    // If date is after beginning of 2010, use linear extrapolation model.
    else
    {
        meanCipPositionInItrs.x( ) = 23.513 + 7.6141 * julianYearsSinceJ2000;
        meanCipPositionInItrs.y( ) = 358.891 + -0.6287 * julianYearsSinceJ2000;
    }

    // Scale values to radians.
    meanCipPositionInItrs *= milliAsToRadians;

    return meanCipPositionInItrs;
}

//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs(
        const double ttSinceEpoch,
        const double utcSinceEpoch,
        const double epochShift )
{
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Check if input is consistent with (required) interpolator settings.
    if( epochShift != basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        std::cerr<<"Warning, polar motion calculation currently only possible for time arguments referenced to J2000."<<std::endl;
    }

    // Add interpolated measured offsets.
    poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );

    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrections(
                ttSinceEpoch );

    return poleOffsetsInItrs;
}

//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs(
        Eigen::Vector6d doodsonArguments,
        const double utcSinceEpoch,
        const double epochShift )
{
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Check if input is consistent with (required) interpolator settings.
    if( epochShift != basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        std::cerr<<"Warning, polar motion calculation currently only possible for time arguments referenced to J2000."<<std::endl;
    }

    // Add interpolated measured offsets.
    poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );

    std::cout<<"Offset: "<<poleOffsetsInItrs<<std::endl;
    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrections(
                doodsonArguments );

    return poleOffsetsInItrs;
}

Eigen::Vector2d PolarMotionCalculator::getPolarMotionWobbleVariables(
        Eigen::Vector6d doodsonArguments,
        const double utcSinceEpoch,
        const double epochShift )
{
    double julianYearsSinceJ2000 = ( basic_astrodynamics::JULIAN_DAY_ON_J2000 - epochShift ) /
            physical_constants::JULIAN_YEAR_IN_DAYS + utcSinceEpoch / physical_constants::JULIAN_YEAR;
    return ( getPositionOfCipInItrs( doodsonArguments, utcSinceEpoch, epochShift ) -
            calculateMeanCipPositionInItrs( julianYearsSinceJ2000 ) );
}

Eigen::Vector2d PolarMotionCalculator::getPolarMotionWobbleVariables(
        const double ttSinceEpoch,
        const double utcSinceEpoch,
        const double epochShift )
{
    double julianYearsSinceJ2000 = ( basic_astrodynamics::JULIAN_DAY_ON_J2000 - epochShift ) /
            physical_constants::JULIAN_YEAR_IN_DAYS + utcSinceEpoch / physical_constants::JULIAN_YEAR;
    return ( getPositionOfCipInItrs( ttSinceEpoch, utcSinceEpoch, epochShift ) -
            calculateMeanCipPositionInItrs( julianYearsSinceJ2000 ) );
}

}

}
