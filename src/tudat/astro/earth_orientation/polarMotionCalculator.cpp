/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{

namespace earth_orientation
{

Eigen::Vector2d getSecularPolePositionInMas( const double time )
{
    double timeInYears = time / physical_constants::JULIAN_YEAR;
    return ( Eigen::Vector2d( ) << 55.0 + 1.677 * timeInYears, 320.5 + 3.460 * timeInYears ).finished( );
}

Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrsAtTt( const double ttSinceEpoch )
{
    if( !evaluator_ ) return getPositionOfCipInItrs( ttSinceEpoch, sofa_interface::convertTTtoUTC( ttSinceEpoch ) );
    const Time tai = evaluator_->convertTtToTai( ttSinceEpoch );
    const Time ut1 = ut1FromTai_ ? ut1FromTai_( tai ) : tai + evaluator_->evaluate( tai ).ut1MinusTai;
    const Eigen::Vector6d arguments =
            sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( ttSinceEpoch, ttSinceEpoch, static_cast< double >( ut1 ) );
    return getPositionOfCipInItrsFromTai( tai, arguments );
}

Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrsAtTt( const Time& tt )
{
    if( !evaluator_ )
        return getPositionOfCipInItrs( static_cast< double >( tt ), sofa_interface::convertTTtoUTC( static_cast< double >( tt ) ) );
    const Time tai = evaluator_->convertTtToTai( tt );
    const Time ut1 = ut1FromTai_ ? ut1FromTai_( tai ) : tai + evaluator_->evaluate( tai ).ut1MinusTai;
    const double ttDouble = static_cast< double >( tt );
    const Eigen::Vector6d arguments =
            sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( ttDouble, ttDouble, static_cast< double >( ut1 ) );
    return getPositionOfCipInItrsFromTai( tai, arguments );
}

Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrsFromTai( const Time& tai, const Eigen::Vector6d& arguments )
{
    if( !evaluator_ ) throw std::runtime_error( "TAI polar motion requires a continuous EOP model" );
    Eigen::Vector2d value = evaluator_->evaluate( tai ).polarMotion;
    if( shortPeriodPolarMotionCalculator_ ) value += shortPeriodPolarMotionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
    return value;
}

//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs( const double ttSinceEpoch, const double utcSinceEpoch )
{
    if( evaluator_ ) return getPositionOfCipInItrsAtTt( ttSinceEpoch );
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Add interpolated measured offsets.
    try
    {
        poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in pole offset interpolation.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }

    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrections( ttSinceEpoch );

    return poleOffsetsInItrs;
}

//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs( Eigen::Vector6d fundamentalArguments, const double utcSinceEpoch )
{
    if( evaluator_ ) return getPositionOfCipInItrsFromTai( sofa_interface::convertUTCtoTAI( Time( utcSinceEpoch ) ), fundamentalArguments );
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Add interpolated measured offsets.
    try
    {
        poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in pole offset interpolation.\nOriginal error: " + std::string( caughtException.what( ) ) );
    }

    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrectionsFromFundamentalArgument( fundamentalArguments );

    return poleOffsetsInItrs;
}

}  // namespace earth_orientation

}  // namespace tudat
