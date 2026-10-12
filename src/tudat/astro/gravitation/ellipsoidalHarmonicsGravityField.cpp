/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityField.h"

#include <cmath>
#include <stdexcept>

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace gravitation
{

using namespace tudat::basic_mathematics;

EllipsoidalHarmonicsGravityCache::EllipsoidalHarmonicsGravityCache( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                                    const double firstFocalParameter,
                                                                    const double secondFocalParameter,
                                                                    const int quadratureOrder,
                                                                    const bool computeLegendreDerivatives ):
    lameFunctionCoefficients_( lameFunctionCoefficients ),
    legendreCache_( lameFunctionCoefficients.getMaximumDegree( ), lameFunctionCoefficients.getMaximumDegree( ), true ),
    firstFocalParameter_( firstFocalParameter ), secondFocalParameter_( secondFocalParameter )
{
    radialFunctions_ = createEllipsoidalHarmonicCoefficients( lameFunctionCoefficients.getMaximumDegree( ) );

    // The derivatives are singular on the long axis of the body, so they are only computed when
    // something actually needs them.
    legendreCache_.setComputeFirstDerivatives( computeLegendreDerivatives );
    legendreCache_.setComputeSecondDerivatives( false );

    computeGaussLegendreNodesAndWeights( quadratureOrder, quadratureNodes_, quadratureWeights_ );

    // Starting from NaN makes the first update always run.
    currentBodyFixedPosition_ = ( Eigen::Vector3d( ) << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    ellipsoidalCoordinates_ = Eigen::Vector3d::Zero( );
}

void EllipsoidalHarmonicsGravityCache::update( const Eigen::Vector3d& currentBodyFixedPosition )
{
    if( currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        ellipsoidalCoordinates_ =
                convertCartesianToEllipsoidalCoordinates( currentBodyFixedPosition, firstFocalParameter_, secondFocalParameter_ );

        computeRadialFunctions( lameFunctionCoefficients_,
                                ellipsoidalCoordinates_( 0 ),
                                firstFocalParameter_,
                                secondFocalParameter_,
                                quadratureNodes_,
                                quadratureWeights_,
                                radialFunctions_ );

        legendreCache_.update( std::sin( ellipsoidalCoordinates_( 1 ) ) );
    }
}

EllipsoidalHarmonicsGravityField::EllipsoidalHarmonicsGravityField( const double gravitationalParameter,
                                                                    const Eigen::Vector3d& referenceEllipsoidSemiAxes,
                                                                    const EllipsoidalHarmonicCoefficients& coefficients,
                                                                    const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                                    const EllipsoidalHarmonicMatrices& transformFactors,
                                                                    const int quadratureOrder,
                                                                    const std::string& fixedReferenceFrame,
                                                                    const std::function< void( ) > updateInertiaTensor ):
    GravityFieldModel( gravitationalParameter, updateInertiaTensor ), referenceEllipsoidSemiAxes_( referenceEllipsoidSemiAxes ),
    firstFocalParameter_( std::sqrt( referenceEllipsoidSemiAxes( 0 ) * referenceEllipsoidSemiAxes( 0 ) -
                                     referenceEllipsoidSemiAxes( 1 ) * referenceEllipsoidSemiAxes( 1 ) ) ),
    secondFocalParameter_( std::sqrt( referenceEllipsoidSemiAxes( 0 ) * referenceEllipsoidSemiAxes( 0 ) -
                                      referenceEllipsoidSemiAxes( 2 ) * referenceEllipsoidSemiAxes( 2 ) ) ),
    maximumDegree_( coefficients.getMaximumDegree( ) ), quadratureOrder_( quadratureOrder ),
    lameFunctionCoefficients_( lameFunctionCoefficients ), transformFactors_( transformFactors ),
    fixedReferenceFrame_( fixedReferenceFrame )
{
    // The semi-axes must be strictly ordered, otherwise the focal parameters are not real and
    // the coordinates this model uses are not defined.
    if( !( referenceEllipsoidSemiAxes( 0 ) > referenceEllipsoidSemiAxes( 1 ) &&
           referenceEllipsoidSemiAxes( 1 ) > referenceEllipsoidSemiAxes( 2 ) && referenceEllipsoidSemiAxes( 2 ) > 0.0 ) )
    {
        throw std::runtime_error(
                "Error when creating ellipsoidal harmonics gravity field, the reference ellipsoid "
                "semi-axes must be given longest first and strictly decreasing, but they are " +
                std::to_string( referenceEllipsoidSemiAxes( 0 ) ) + ", " + std::to_string( referenceEllipsoidSemiAxes( 1 ) ) + ", " +
                std::to_string( referenceEllipsoidSemiAxes( 2 ) ) + "." );
    }

    if( lameFunctionCoefficients.getMaximumDegree( ) != maximumDegree_ || transformFactors.getMaximumDegree( ) != maximumDegree_ )
    {
        throw std::runtime_error(
                "Error when creating ellipsoidal harmonics gravity field, the coefficients, the Lame "
                "polynomials and the transform factors must all have the same maximum degree, but "
                "they have " +
                std::to_string( maximumDegree_ ) + ", " + std::to_string( lameFunctionCoefficients.getMaximumDegree( ) ) + " and " +
                std::to_string( transformFactors.getMaximumDegree( ) ) + "." );
    }

    gravityCache_ = std::make_shared< EllipsoidalHarmonicsGravityCache >(
            lameFunctionCoefficients_, firstFocalParameter_, secondFocalParameter_, quadratureOrder_, false );

    // The coefficients are supplied normalised to the surface of the reference ellipsoid, so the
    // radial parts are evaluated there once and divided out.
    Eigen::VectorXd quadratureNodes;
    Eigen::VectorXd quadratureWeights;
    computeGaussLegendreNodesAndWeights( quadratureOrder_, quadratureNodes, quadratureWeights );

    EllipsoidalHarmonicCoefficients radialFunctionsAtSurface = createEllipsoidalHarmonicCoefficients( maximumDegree_ );
    computeRadialFunctions( lameFunctionCoefficients_,
                            referenceEllipsoidSemiAxes_( 0 ),
                            firstFocalParameter_,
                            secondFocalParameter_,
                            quadratureNodes,
                            quadratureWeights,
                            radialFunctionsAtSurface );

    coefficients_ = normaliseCoefficientsToReferenceSurface( coefficients, radialFunctionsAtSurface );
}

void EllipsoidalHarmonicsGravityField::checkPositionIsOutsideReferenceEllipsoid( const Eigen::Vector3d& bodyFixedPosition,
                                                                                 const double ellipsoidalRadius ) const
{
    if( ellipsoidalRadius < referenceEllipsoidSemiAxes_( 0 ) )
    {
        throw std::runtime_error( "Error when evaluating ellipsoidal harmonics gravity field, the position (" +
                                  std::to_string( bodyFixedPosition( 0 ) ) + ", " + std::to_string( bodyFixedPosition( 1 ) ) + ", " +
                                  std::to_string( bodyFixedPosition( 2 ) ) +
                                  ") lies inside the reference ellipsoid, where the expansion does not converge. Its "
                                  "ellipsoidal radius is " +
                                  std::to_string( ellipsoidalRadius ) + " and the reference ellipsoid reaches " +
                                  std::to_string( referenceEllipsoidSemiAxes_( 0 ) ) + "." );
    }
}

double EllipsoidalHarmonicsGravityField::getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    gravityCache_->update( bodyFixedPosition );
    checkPositionIsOutsideReferenceEllipsoid( bodyFixedPosition, gravityCache_->getEllipsoidalRadius( ) );

    return sumEllipsoidalHarmonicTerms( coefficients_,
                                        gravityCache_->getRadialFunctions( ),
                                        transformFactors_,
                                        gravityCache_->getLegendreCache( ),
                                        gravityCache_->getAngleAboutLongAxis( ) );
}

Eigen::Vector3d EllipsoidalHarmonicsGravityField::getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    throw std::runtime_error(
            "Error, the gradient of an ellipsoidal harmonics gravity field is not implemented yet, so "
            "this field cannot be used to compute an acceleration." );
}

}  // namespace gravitation

}  // namespace tudat
