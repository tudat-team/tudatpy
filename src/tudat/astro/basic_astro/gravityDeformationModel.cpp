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

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{

namespace basic_astrodynamics
{

MaxwellGravityDeformationModel::MaxwellGravityDeformationModel(
        const StateFunction stateOfDeformingBodyFunction,
        const std::vector< std::string > perturbingBody,
        const double maxwellRelaxationTime,
        const double globalRelaxationTime,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldModel,
        const std::vector< std::function< double( ) > > gravitationalParameterPerturbingBody,
        const std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody,
        const std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody,
        const double k2,
        const std::vector< StateFunction > stateOfPerturbingBodyFunction,
        const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction,
        const std::function< Eigen::Matrix3d( ) > rotationToLocalFrameDerivativeFunction,
        const bool includeOrder1,
        const bool includeCentrifugalPotential ):
    GravityDeformationModel( ), perturbingBody_( perturbingBody ), maxwellRelaxationTime_( maxwellRelaxationTime ),
    globalRelaxationTime_( globalRelaxationTime ), gravityFieldModel_( gravityFieldModel ),
    gravitationalParameterPerturbingBody_( gravitationalParameterPerturbingBody ), k2_( k2 ),
    equilibriumCoefficients_( Eigen::VectorXd::Zero( 5 ) ), derivativeEquilibriumCoefficients_( Eigen::VectorXd::Zero( 5 ) ),
    currentCoefficientVariation_( Eigen::VectorXd::Zero( 5 ) ),
    rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
    rotationToBodyFixedDerivativeFunction_( rotationToLocalFrameDerivativeFunction ),
    stateOfDeformingBodyFunction_( stateOfDeformingBodyFunction ), stateOfPerturbingBodyFunction_( stateOfPerturbingBodyFunction ),
    includeOrder1_( includeOrder1 ), includeCentrifugalPotential_( includeCentrifugalPotential ),
    angularVelocityDeformingBody_( angularVelocityDeformingBody ),
    angularVelocityDerivativeDeformingBody_( angularVelocityDerivativeDeformingBody )
{
    if( gravityFieldModel_ == nullptr )
    {
        throw std::runtime_error( "Error when creating Maxwell gravity deformation model: gravity field is null." );
    }
    if( gravitationalParameterPerturbingBody_.size( ) != perturbingBody_.size( ) ||
        stateOfPerturbingBodyFunction_.size( ) != perturbingBody_.size( ) )
    {
        throw std::runtime_error( "Error when creating Maxwell gravity deformation model: inconsistent perturber input sizes." );
    }
    const Eigen::MatrixXd cosineCoefficients = gravityFieldModel_->getCosineCoefficientsWithoutVariations( );
    const Eigen::MatrixXd sineCoefficients = gravityFieldModel_->getSineCoefficientsWithoutVariations( );
    if( cosineCoefficients.rows( ) < 3 || cosineCoefficients.cols( ) < 3 || sineCoefficients.rows( ) < 3 || sineCoefficients.cols( ) < 3 )
    {
        throw std::runtime_error( "Error when creating Maxwell gravity deformation model: degree-two gravity coefficients are absent." );
    }

    const unsigned int numberPerturbingBodies = perturbingBody_.size( );
    stateOfPerturbingBody_.resize( numberPerturbingBodies );
    currentRelativePosition_.resize( numberPerturbingBodies );
    currentRelativeVelocity_.resize( numberPerturbingBodies );
    currentInertialRelativeState_.resize( numberPerturbingBodies );
    currentLongitude_.resize( numberPerturbingBodies );
    currentLatitude_.resize( numberPerturbingBodies );
    currentLongitudeDerivative_.resize( numberPerturbingBodies );
    currentLatitudeDerivative_.resize( numberPerturbingBodies );
}

void MaxwellGravityDeformationModel::updateMembers( const double currentTime )
{
    // The integrated state is a variation, so remove the gravity field's own variation-free
    // baseline from its current total coefficients. The Maxwell model does not retain a second
    // copy of the static field.
    const Eigen::MatrixXd cosineHarmonicCoefficients =
            gravityFieldModel_->getCosineCoefficients( ) - gravityFieldModel_->getCosineCoefficientsWithoutVariations( );
    const Eigen::MatrixXd sineHarmonicCoefficients =
            gravityFieldModel_->getSineCoefficients( ) - gravityFieldModel_->getSineCoefficientsWithoutVariations( );

    currentCoefficientVariation_[ 0 ] = cosineHarmonicCoefficients( 2, 0 );
    currentCoefficientVariation_[ 1 ] = cosineHarmonicCoefficients( 2, 1 );
    currentCoefficientVariation_[ 2 ] = cosineHarmonicCoefficients( 2, 2 );
    currentCoefficientVariation_[ 3 ] = sineHarmonicCoefficients( 2, 1 );
    currentCoefficientVariation_[ 4 ] = sineHarmonicCoefficients( 2, 2 );

    // The environment stores geodesy-normalised coefficients; the propagated Maxwell state is unnormalised.
    currentCoefficientVariation_[ 0 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    currentCoefficientVariation_[ 1 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    currentCoefficientVariation_[ 2 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    currentCoefficientVariation_[ 3 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    currentCoefficientVariation_[ 4 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

    // Update rotation and positions
    rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
    stateOfDeformingBodyFunction_( stateOfDeformingBody_ );

    for( unsigned int k = 0; k < perturbingBody_.size( ); k++ )
    {
        stateOfPerturbingBodyFunction_.at( k )( stateOfPerturbingBody_.at( k ) );

        // Compute relative inertial state
        currentInertialRelativeState_.at( k ) = stateOfPerturbingBody_.at( k ) - stateOfDeformingBody_;
    }

    Eigen::Matrix3d currentRotationToLocalFrameDerivative = rotationToBodyFixedDerivativeFunction_( );

    for( unsigned int k = 0; k < perturbingBody_.size( ); k++ )
    {
        // Compute current relative state in body-fixed frame
        currentRelativePosition_.at( k ) = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativeState_.at( k ).segment( 0, 3 );
        currentRelativeVelocity_.at( k ) = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativeState_.at( k ).segment( 3, 3 ) +
                currentRotationToLocalFrameDerivative * currentInertialRelativeState_.at( k ).segment( 0, 3 );

        // Compute spherical coordinates of perturbing body in body-fixed frame
        Eigen::Vector3d currentSphericalPositionPerturbingBody =
                coordinate_conversions::convertCartesianToSpherical( currentRelativePosition_.at( k ) );
        currentLongitude_.at( k ) = currentSphericalPositionPerturbingBody[ 2 ];
        currentLatitude_.at( k ) = mathematical_constants::PI / 2.0 - currentSphericalPositionPerturbingBody.y( );

        // Compute current derivative of the perturbing body's body-fixed longitude
        currentLongitudeDerivative_.at( k ) = ( ( currentRelativeVelocity_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 0 ] -
                                                  currentRelativeVelocity_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 1 ] ) /
                                                ( currentRelativePosition_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 0 ] +
                                                  currentRelativePosition_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 1 ] ) );

        // Compute current derivative of the perturbing body's body-fixed latitude
        double currentDistance = currentRelativePosition_.at( k ).segment( 0, 3 ).norm( );
        double currentDistanceDerivative = ( currentRelativePosition_.at( k )[ 0 ] * currentRelativeVelocity_.at( k )[ 0 ] +
                                             currentRelativePosition_.at( k )[ 1 ] * currentRelativeVelocity_.at( k )[ 1 ] +
                                             currentRelativePosition_.at( k )[ 2 ] * currentRelativeVelocity_.at( k )[ 2 ] ) /
                currentDistance;
        currentLatitudeDerivative_.at( k ) = ( currentRelativeVelocity_.at( k )[ 2 ] * currentDistance -
                                               currentRelativePosition_.at( k )[ 2 ] * currentDistanceDerivative ) /
                ( currentDistance *
                  std::sqrt( currentRelativePosition_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 0 ] +
                             currentRelativePosition_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 1 ] ) );
    }

    updateEquilibriumDeformation( );

    currentDeformation_ = ( 1.0 / globalRelaxationTime_ ) *
            ( equilibriumCoefficients_ - currentCoefficientVariation_ + maxwellRelaxationTime_ * derivativeEquilibriumCoefficients_ );
}

void MaxwellGravityDeformationModel::updateEquilibriumDeformation( const double currentTime )
{
    // Reset equilibrium coefficients and derivatives to zero
    equilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );
    derivativeEquilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );

    const Eigen::Vector3d angularVelocity = angularVelocityDeformingBody_( );
    const Eigen::Vector3d angularVelocityDerivative = angularVelocityDerivativeDeformingBody_( );
    const double rotationRate = angularVelocity.norm( );
    const double rotationRateDerivative = ( rotationRate > 0.0 ) ? angularVelocity.dot( angularVelocityDerivative ) / rotationRate : 0.0;
    const double gravitationalParameterDeformingBody = gravityFieldModel_->getGravitationalParameter( );
    const double referenceRadius = gravityFieldModel_->getReferenceRadius( );
    if( includeCentrifugalPotential_ )
    {
        const double centrifugalCoefficientFactor =
                -k2_ * referenceRadius * referenceRadius * referenceRadius / ( 3.0 * gravitationalParameterDeformingBody );
        equilibriumCoefficients_[ 0 ] += centrifugalCoefficientFactor * rotationRate * rotationRate;
        derivativeEquilibriumCoefficients_[ 0 ] += 2.0 * centrifugalCoefficientFactor * rotationRate * rotationRateDerivative;
    }

    for( unsigned int k = 0; k < perturbingBody_.size( ); k++ )
    {
        double relativeDistance = currentRelativePosition_.at( k ).segment( 0, 3 ).norm( );
        double radiusRatioPowerThree =
                referenceRadius * referenceRadius * referenceRadius / ( relativeDistance * relativeDistance * relativeDistance );

        double gravitationalParametersRatio = gravitationalParameterPerturbingBody_.at( k )( ) / gravitationalParameterDeformingBody;

        equilibriumCoefficients_[ 0 ] += k2_ *
                ( 0.5 * gravitationalParametersRatio * radiusRatioPowerThree *
                  ( 3.0 * std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) - 1.0 ) );
        equilibriumCoefficients_[ 2 ] += k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree *
                ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                std::cos( 2.0 * currentLongitude_.at( k ) );
        equilibriumCoefficients_[ 4 ] += k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree *
                ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                std::sin( 2.0 * currentLongitude_.at( k ) );
        if( includeOrder1_ )
        {
            equilibriumCoefficients_[ 1 ] += -k2_ * gravitationalParametersRatio * radiusRatioPowerThree *
                    ( -std::cos( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                    std::cos( currentLongitude_.at( k ) );
            equilibriumCoefficients_[ 3 ] += -k2_ * gravitationalParametersRatio * radiusRatioPowerThree *
                    ( -std::cos( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                    std::sin( currentLongitude_.at( k ) );
        }
        double currentDistanceDerivative = ( currentRelativePosition_.at( k )[ 0 ] * currentRelativeVelocity_.at( k )[ 0 ] +
                                             currentRelativePosition_.at( k )[ 1 ] * currentRelativeVelocity_.at( k )[ 1 ] +
                                             currentRelativePosition_.at( k )[ 2 ] * currentRelativeVelocity_.at( k )[ 2 ] ) /
                relativeDistance;

        derivativeEquilibriumCoefficients_[ 0 ] += -k2_ *
                ( 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree * 3.0 * currentDistanceDerivative / relativeDistance *
                          ( 3.0 * std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) - 1.0 ) -
                  1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree *
                          ( 6.0 * currentLatitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) *
                            std::cos( currentLatitude_.at( k ) ) ) );

        derivativeEquilibriumCoefficients_[ 2 ] += -k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree *
                ( 3.0 * currentDistanceDerivative / relativeDistance *
                          ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                          std::cos( 2.0 * currentLongitude_.at( k ) ) +
                  2.0 * currentLongitudeDerivative_.at( k ) * std::sin( 2.0 * currentLongitude_.at( k ) ) *
                          ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) +
                  2.0 * std::cos( 2.0 * currentLongitude_.at( k ) ) * currentLatitudeDerivative_.at( k ) *
                          std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) );

        derivativeEquilibriumCoefficients_[ 4 ] += -k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree *
                ( 3.0 * currentDistanceDerivative / relativeDistance *
                          ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                          std::sin( 2.0 * currentLongitude_.at( k ) ) -
                  2.0 * currentLongitudeDerivative_.at( k ) * std::cos( 2.0 * currentLongitude_.at( k ) ) *
                          ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) +
                  2.0 * std::sin( 2.0 * currentLongitude_.at( k ) ) * currentLatitudeDerivative_.at( k ) *
                          std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) );

        if( includeOrder1_ )
        {
            derivativeEquilibriumCoefficients_[ 1 ] += -k2_ * gravitationalParametersRatio * radiusRatioPowerThree *
                    ( 3.0 * currentDistanceDerivative / relativeDistance * std::sin( currentLatitude_.at( k ) ) *
                              std::cos( currentLatitude_.at( k ) ) * std::cos( currentLongitude_.at( k ) ) +
                      currentLongitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) *
                              std::sin( currentLongitude_.at( k ) ) -
                      currentLatitudeDerivative_.at( k ) *
                              ( std::cos( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) -
                                std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                              std::cos( currentLongitude_.at( k ) ) );

            derivativeEquilibriumCoefficients_[ 3 ] += -k2_ * gravitationalParametersRatio * radiusRatioPowerThree *
                    ( 3.0 * currentDistanceDerivative / relativeDistance * std::sin( currentLatitude_.at( k ) ) *
                              std::cos( currentLatitude_.at( k ) ) * std::sin( currentLongitude_.at( k ) ) -
                      currentLongitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) *
                              std::cos( currentLongitude_.at( k ) ) -
                      currentLatitudeDerivative_.at( k ) *
                              ( std::cos( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) -
                                std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) *
                              std::sin( currentLongitude_.at( k ) ) );
        }
    }
}

}  // namespace basic_astrodynamics

}  // namespace tudat
