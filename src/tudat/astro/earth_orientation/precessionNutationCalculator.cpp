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
#include <exception>
#include <limits>
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
        nominalCipInterpolator_ = angleInterpolator;
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

void PrecessionNutationCalculator::setNominalDerivativeTimeStep( const double step )
{
    if( !std::isfinite( step ) || step <= 0.0 ) throw std::invalid_argument( "Nominal CIP derivative step must be positive and finite" );
    nominalDerivativeTimeStep_ = step;
}

std::pair< Eigen::Vector3d, Eigen::Vector3d > PrecessionNutationCalculator::getNominalPositionAndRate( const Time& tt )
{
    const double centre = static_cast< double >( tt );
    if( !std::isfinite( centre ) ) throw std::invalid_argument( "Non-finite nominal CIP derivative epoch" );
    if( nominalCipInterpolator_ )
    {
        const auto nodes = nominalCipInterpolator_->getIndependentValues( );
        if( nodes.size( ) < 2 || centre < nodes.front( ) || centre > nodes.back( ) )
            throw std::runtime_error( "Nominal CIP interpolation rate outside its coverage; position policy remains unchanged" );
    }
    Eigen::Vector3d central = nominalCipPositionFunction_( centre );
    if( !central.allFinite( ) ) throw std::runtime_error( "Non-finite nominal CIP position" );
    Eigen::Vector3d rate = Eigen::Vector3d::Zero( );
    try
    {
        if( !nominalCipInterpolator_ )
        {
            const double h = nominalDerivativeTimeStep_;
            const auto sample = [ & ]( double delta ) -> Eigen::Vector3d {
                const Time shifted = tt + delta;
                const long double separation = static_cast< long double >( shifted - tt );
                if( !std::isfinite( separation ) ||
                    std::abs( separation - delta ) > 16.0L * std::numeric_limits< double >::epsilon( ) * std::abs( delta ) )
                    throw std::invalid_argument( "Unrepresentable nominal CIP derivative stencil" );
                const double shiftedDouble = static_cast< double >( shifted );
                if( std::abs( ( shiftedDouble - centre ) - delta ) > 1e-6 * std::abs( delta ) )
                    throw std::invalid_argument( "Nominal CIP derivative step is too small for the scalar SOFA epoch" );
                const Eigen::Vector3d value = nominalCipPositionFunction_( shiftedDouble );
                if( !value.allFinite( ) ) throw std::runtime_error( "Non-finite nominal CIP derivative sample" );
                return value;
            };
            rate = ( 8.0 * ( sample( h ) - sample( -h ) ) - ( sample( 2 * h ) - sample( -2 * h ) ) ) / ( 12 * h );
        }
        else
        {
            const auto x = nominalCipInterpolator_->getIndependentValues( );
            if( x.size( ) < 2 || centre < x.front( ) || centre > x.back( ) )
                throw std::runtime_error( "Nominal CIP interpolation rate outside its coverage; position policy remains unchanged" );
            // The derivative at a knot is the branch on its right, except at the final endpoint.
            const std::size_t left =
                    std::min< std::size_t >( std::upper_bound( x.begin( ), x.end( ), centre ) - x.begin( ) - 1, x.size( ) - 2 );
            const double width = x[ left + 1 ] - x[ left ];
            if( !std::isfinite( width ) || width <= 0.0 ) throw std::runtime_error( "Invalid nominal CIP interpolation interval" );
            const auto type = nominalCipInterpolator_->getInterpolatorType( );
            std::size_t count = 4;
            if( type == interpolators::linear_interpolator )
            {
                const auto y = nominalCipInterpolator_->getDependentValues( );
                rate = ( y[ left + 1 ] - y[ left ] ) / width;
            }
            else
            {
                if( type == interpolators::lagrange_interpolator )
                {
                    const auto lagrange = std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< double, Eigen::Vector3d > >(
                            nominalCipInterpolator_ );
                    if( !lagrange ) throw std::runtime_error( "Unsupported nominal CIP Lagrange scalar type" );
                    count = lagrange->getNumberOfStages( );
                    const std::size_t offset = count / 2 - 1;
                    if( left < offset || left >= x.size( ) - offset - 1 )
                    {
                        const auto boundary = lagrange->getLagrangeBoundaryHandling( );
                        if( boundary != interpolators::lagrange_cubic_spline_boundary_interpolation &&
                            boundary != interpolators::lagrange_cubic_spline_boundary_interpolation_with_warning )
                            throw std::runtime_error( "Nominal CIP boundary interpolation rate unavailable" );
                        count = 4;
                    }
                }
                else if( type != interpolators::cubic_spline_interpolator && type != interpolators::hermite_spline_interpolator )
                    throw std::runtime_error( "Selected nominal CIP interpolation rate unavailable" );
                // Recover the selected polynomial branch from interior samples only.
                // Actual representable abscissae give nonuniform weights when needed.
                std::vector< double > nodes;
                std::vector< Eigen::Vector3d > values;
                for( std::size_t i = 0; i < count; ++i )
                {
                    const double sampleTime = x[ left ] + width * ( static_cast< double >( i ) + 0.5 ) / count;
                    if( sampleTime <= x[ left ] || sampleTime >= x[ left + 1 ] ||
                        ( !nodes.empty( ) && ( sampleTime - x[ left ] ) / width <= nodes.back( ) ) )
                        throw std::runtime_error( "Unrepresentable nominal CIP polynomial rate samples" );
                    nodes.push_back( ( sampleTime - x[ left ] ) / width );
                    const Eigen::Vector3d value = nominalCipPositionFunction_( sampleTime );
                    if( !value.allFinite( ) ) throw std::runtime_error( "Non-finite nominal CIP polynomial sample" );
                    values.push_back( value );
                }
                const double u = ( centre - x[ left ] ) / width;
                for( std::size_t j = 1; j < count; ++j )
                {
                    double derivative = 0.0;
                    for( std::size_t k = 0; k < count; ++k )
                    {
                        if( k == j ) continue;
                        double term = 1.0 / ( nodes[ j ] - nodes[ k ] );
                        for( std::size_t m = 0; m < count; ++m )
                            if( m != j && m != k ) term *= ( u - nodes[ m ] ) / ( nodes[ j ] - nodes[ m ] );
                        derivative += term;
                    }
                    rate += ( derivative / width ) * ( values[ j ] - values.front( ) );
                }
            }
        }
        if( !rate.allFinite( ) ) throw std::runtime_error( "Non-finite nominal CIP derivative" );
        central = nominalCipPositionFunction_( centre );
        return { central, rate };
    }
    catch( ... )
    {
        const auto original = std::current_exception( );
        try
        {
            nominalCipPositionFunction_( centre );
        }
        catch( ... )
        {}
        std::rethrow_exception( original );
    }
}

Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrsFromTai( const double terrestrialTime, const Time& tai )
{
    if( !evaluator_ ) throw std::runtime_error( "TAI CPO evaluation requires a continuous EOP model" );
    Eigen::Vector3d value = nominalCipPositionFunction_( terrestrialTime );
    value.segment< 2 >( 0 ) += evaluator_->getCelestialPoleOffsets( tai );
    return value;
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrs( const double terrestrialTime )
{
    if( evaluator_ )
        return getPositionOfCipInGcrsFromTai( terrestrialTime, basic_astrodynamics::convertTTtoTAI( Time( terrestrialTime ) ) );
    // Calculate current UTC from SOFA.
    double utc = sofa_interface::convertTTtoUTC( terrestrialTime );

    // Call function to compte precession-nutation from UTC and TT.
    return getPositionOfCipInGcrs( terrestrialTime, utc );
}

//! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
Eigen::Vector3d PrecessionNutationCalculator::getPositionOfCipInGcrs( const double terrestrialTime, const double utc )
{
    if( evaluator_ )
        return getPositionOfCipInGcrsFromTai( terrestrialTime, basic_astrodynamics::convertTTtoTAI( Time( terrestrialTime ) ) );
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
