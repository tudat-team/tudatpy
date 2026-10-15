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
#include "tudat/math/basic/numericalDerivative.h"
#include <cmath>

namespace tudat
{

namespace earth_orientation
{

//! Constructor for (CIO-based) precession-nutation calculation object.
PrecessionNutationCalculator::PrecessionNutationCalculator(
        const basic_astrodynamics::IAUConventions precessionNutationTheory,
        const std::shared_ptr< EarthOrientationDataEvaluator > eopData,
        const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > angleInterpolatorSettings ):
    eopData_( eopData ), precessionNutationTheory_( precessionNutationTheory )
{
    if( eopData_ == nullptr )
    {
        throw std::invalid_argument( "Precession-nutation calculator requires Earth orientation data" );
    }
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

// X, Y and s from the SOFA implementation of the IAU model.
Eigen::Vector3d PrecessionNutationCalculator::getNominalPositionFromSofa( const Time& tt )
{
    return sofa_interface::getPositionOfCipInGcrs(
            static_cast< double >( tt ), basic_astrodynamics::JULIAN_DAY_ON_J2000, precessionNutationTheory_ );
}

// Differentiate the SOFA model numerically, also when the nominal positions are interpolated.
std::pair< Eigen::Vector3d, Eigen::Vector3d > PrecessionNutationCalculator::getNominalPositionAndRate( const Time& tt )
{
    // Fourth-order central difference over 32 s.
    const std::function< Eigen::Vector3d( const Time& ) > nominalPositionFunction =
            std::bind( &PrecessionNutationCalculator::getNominalPositionFromSofa, this, std::placeholders::_1 );
    const Eigen::Vector3d rate = numerical_derivatives::computeFourthOrderCentralDifference( nominalPositionFunction, tt, 32.0 );

    if( !rate.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite nominal CIP rate" );
    }
    return { nominalCipPositionFunction_( static_cast< double >( tt ) ), rate };
}

// IAU model at TT, plus the measured celestial pole offsets at TAI.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrsFromTai( const double terrestrialTime, const Time& tai )
{
    // Nominal X, Y and s in radians; measured dX/dY affect the first two components.
    Eigen::Vector3d value = nominalCipPositionFunction_( terrestrialTime );

    value.segment< 2 >( 0 ) += eopData_->getCelestialPoleOffsets( tai );
    return value;
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrs( const double terrestrialTime )
{
    return getPositionOfCipInGcrsFromTai( terrestrialTime, basic_astrodynamics::convertTTtoTAI( Time( terrestrialTime ) ) );
}

}  // namespace earth_orientation

}  // namespace tudat
