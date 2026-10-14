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

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/precessionNutationCalculator.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include <cmath>

namespace tudat
{

namespace earth_orientation
{

//! Constructor for (CIO-based) precession-nutation calculation object.
PrecessionNutationCalculator::PrecessionNutationCalculator(
        basic_astrodynamics::IAUConventions precessionNutationTheory,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > dailyCorrectionInterpolator,
        const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > angleInterpolatorSettings,
        const std::shared_ptr< EarthOrientationDataEvaluator > evaluator ):
    evaluator_( evaluator ), precessionNutationTheory_( precessionNutationTheory ),
    dailyCorrectionInterpolator_( dailyCorrectionInterpolator )
{
    if( angleInterpolatorSettings != nullptr )
    {
        std::function< Eigen::Vector3d( const double ) > angleFunction = std::bind( &sofa_interface::getPositionOfCipInGcrs,
                                                                                    std::placeholders::_1,
                                                                                    basic_astrodynamics::JULIAN_DAY_ON_J2000,
                                                                                    precessionNutationTheory );
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > angleInterpolator =
                interpolators::createOneDimensionalInterpolator< double, Eigen::Vector3d >( angleFunction, angleInterpolatorSettings );
        typedef interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > LocalInterpolator;
        nominalCipPositionFunction_ =
                std::bind( static_cast< Eigen::Vector3d ( LocalInterpolator::* )( const double ) >( &LocalInterpolator::interpolate ),
                           angleInterpolator,
                           std::placeholders::_1 );
    }
    else
    {
        // Link selected SOFA function wrapper for direct calculation of precession-nutation.
        nominalCipPositionFunction_ = std::bind( sofa_interface::getPositionOfCipInGcrs,
                                                 std::placeholders::_1,
                                                 basic_astrodynamics::JULIAN_DAY_ON_J2000,
                                                 precessionNutationTheory );
    }
}

// Differentiate the direct SOFA model, also when nominal positions are interpolated.
std::pair< Eigen::Vector3d, Eigen::Vector3d > PrecessionNutationCalculator::getNominalPositionAndRate( const Time& tt )
{
    // Nominal X, Y and s in radians at a TT offset in seconds; the epoch offset is formed in Time.
    const auto directPosition = [ this, &tt ]( const double offset ) -> Eigen::Vector3d {
        return sofa_interface::getPositionOfCipInGcrs(
                static_cast< double >( tt + offset ), basic_astrodynamics::JULIAN_DAY_ON_J2000, precessionNutationTheory_ );
    };

    // Fourth-order central difference over 32 s, pairing opposite samples before weighting.
    const double step = 32.0;

    const Eigen::Vector3d rate = ( 8.0 * ( directPosition( step ) - directPosition( -step ) ) -
                                   ( directPosition( 2.0 * step ) - directPosition( -2.0 * step ) ) ) /
            ( 12.0 * step );

    if( !rate.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite nominal CIP rate" );
    }
    return { nominalCipPositionFunction_( static_cast< double >( tt ) ), rate };
}

// Add the selected measured celestial-pole corrections to the nominal TT model.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrsFromTai( const double terrestrialTime, const Time& tai )
{
    if( !evaluator_ )
    {
        throw std::runtime_error( "TAI CPO evaluation requires a continuous EOP model" );
    }

    // Nominal X, Y and s in radians; measured dX/dY affect the first two components.
    Eigen::Vector3d value = nominalCipPositionFunction_( terrestrialTime );

    value.segment< 2 >( 0 ) += evaluator_->getCelestialPoleOffsets( tai );
    return value;
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrs( const double terrestrialTime )
{
    if( evaluator_ )
    {
        return getPositionOfCipInGcrsFromTai( terrestrialTime, basic_astrodynamics::convertTTtoTAI( Time( terrestrialTime ) ) );
    }

    // Calculate current UTC from SOFA.
    double utc = sofa_interface::convertTTtoUTC( terrestrialTime );

    // Call function to compte precession-nutation from UTC and TT.
    return getPositionOfCipInGcrs( terrestrialTime, utc );
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrs( const double terrestrialTime, const double utc )
{
    if( evaluator_ )
    {
        return getPositionOfCipInGcrsFromTai( terrestrialTime, basic_astrodynamics::convertTTtoTAI( Time( terrestrialTime ) ) );
    }

    // Calculate nominal precession-nutation values.
    Eigen::Vector3d nominalCipPosition = nominalCipPositionFunction_( terrestrialTime );

    // Retrieve measured corrections to model.
    try
    {
        nominalCipPosition.segment( 0, 2 ) += dailyCorrectionInterpolator_->interpolate( utc );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error in precession-nutation interpolation.\nOriginal error: " +
                                  std::string( caughtException.what( ) ) );
    }

    // Add nominal values and corrections and return.
    return nominalCipPosition;
}

}  // namespace earth_orientation

}  // namespace tudat
