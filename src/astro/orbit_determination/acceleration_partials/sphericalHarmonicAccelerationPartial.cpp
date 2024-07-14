/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicPartialFunctions.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravityFieldVariationParameters.h"

namespace tudat
{

namespace acceleration_partials
{

//! Contructor.
SphericalHarmonicsGravityPartial::SphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials,
        const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >&
        tidalLoveNumberPartialInterfaces ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::spherical_harmonic_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    bodyReferenceRadius_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getReferenceRadius,
                                     accelerationModel ) ),
    cosineCoefficients_( accelerationModel->getCosineHarmonicCoefficientsFunction( ) ),
    sineCoefficients_( accelerationModel->getSineHarmonicCoefficientsFunction( ) ),
    sphericalHarmonicCache_( accelerationModel->getSphericalHarmonicsCache( ) ),
    positionFunctionOfAcceleratedBody_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                   getCurrentPositionOfBodySubjectToAcceleration, accelerationModel ) ),
    positionFunctionOfAcceleratingBody_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                    getCurrentPositionOfBodyExertingAcceleration, accelerationModel ) ),
    fromBodyFixedToIntegrationFrameRotation_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                         getCurrentRotationToIntegrationFrameMatrix, accelerationModel ) ),
    accelerationFunction_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAcceleration,
                                      accelerationModel ) ),
    updateFunction_( std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::updateMembers,
                                accelerationModel, std::placeholders::_1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials ),
    tidalLoveNumberPartialInterfaces_( tidalLoveNumberPartialInterfaces ),
    accelerationUsesMutualAttraction_( accelerationModel->getIsMutualAttractionUsed( ) )
{
    sphericalHarmonicCache_->getLegendreCache( )->setComputeSecondDerivatives( 1 );

    // Update number of degrees and orders in legendre cache for calculation of position partials

    maximumDegree_ = cosineCoefficients_( ).rows( ) - 1;
    maximumOrder_ = sineCoefficients_( ).cols( ) - 1;

    if( sphericalHarmonicCache_->getMaximumDegree( ) < maximumDegree_ ||
            sphericalHarmonicCache_->getMaximumOrder( ) < maximumOrder_ + 2 )
    {
        sphericalHarmonicCache_->resetMaximumDegreeAndOrder( maximumDegree_, maximumOrder_ + 2 );
    }
}

//! Function to create a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > SphericalHarmonicsGravityPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;


    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair =
                getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
        partialFunction = partialFunctionPair.first;
        numberOfRows = partialFunctionPair.second;
    }
    else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                               parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = std::bind(
                            &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                            this, std::placeholders::_1, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = 1;
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of sh acceleration wrt rotational parameter" +
                        std::to_string( parameter->getParameterName( ).first ) + " of " +
                        parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
            }
        }

        // Check if partial is a tidal property of body exerting acceleration.
        else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::TidalLoveNumber< double > > tidalLoveNumber =
                    std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< double >  >( parameter );
            if( tidalLoveNumber == nullptr )
            {
                throw std::runtime_error( "Error when getting tidal Love number vector parameter, object is nullptr" );
            }

            // Get degree and order(s) of tidal variations
            int degree = tidalLoveNumber->getDegree( );
            std::vector< int > orders = tidalLoveNumber->getOrders( );
            int sumOrders = tidalLoveNumber->getSumOrders( );

            std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
            for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
            {
                // Check dependency on current partial object
                currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction(
                            parameter, maximumDegree_, maximumOrder_ );

                // Check consistency
                if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
                {
                    throw std::runtime_error( "Error when getting double tidal parameter partial, multiple dependencies found " +
                                              std::to_string( numberOfRows ) + ", " +
                                              std::to_string( currentTidalPartialOutput.first ) );
                }
                else
                {
                    // If tidal dependency esists, set partial function
                    if( currentTidalPartialOutput.first > 0 )
                    {
                        std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
                                std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentDoubleParameterPartial,
                                           tidalLoveNumberPartialInterfaces_.at( i ),
                                           parameter, currentTidalPartialOutput.second );
                        partialFunction = std::bind(
                                    &SphericalHarmonicsGravityPartial::wrtTidalModelParameter, this, coefficientPartialFunction,
                                    degree, orders, sumOrders, parameter->getParameterSize( ), std::placeholders::_1 );
                        numberOfRows = currentTidalPartialOutput.first;
                    }
                }
            }
        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > SphericalHarmonicsGravityPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                               parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = std::bind(
                            &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                            this, std::placeholders::_1, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = parameter->getParameterSize( );
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of sh acceleration wrt rotational parameter" +
                        std::to_string( parameter->getParameterName( ).first ) + " of " +
                        parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
            }
        }
        // Check if partial is a tidal property of body exerting acceleration.
        else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > tidalLoveNumber =
                    std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd >  >( parameter );
            if( tidalLoveNumber == nullptr )
            {
                throw std::runtime_error( "Error when getting tidal Love number vector parameter, object is nullptr" );
            }

            // Get degree and order(s) of tidal variations
            int degree = tidalLoveNumber->getDegree( );
            std::vector< int > orders = tidalLoveNumber->getOrders( );
            int sumOrders = tidalLoveNumber->getSumOrders( );

            std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
            for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
            {
                // Check dependency on current partial object
                currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction(
                            parameter, maximumDegree_, maximumOrder_ );
                if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
                {
                    throw std::runtime_error( "Error when getting vector tidal parameter partial, inconsistent output" +
                                              std::to_string( numberOfRows ) + ", " +
                                              std::to_string( currentTidalPartialOutput.first ) );
                }
                else
                {
                    // If tidal dependency esists, set partial function
                    if( currentTidalPartialOutput.first > 0 )
                    {
                        std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
                                std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentVectorParameterPartial,
                                           tidalLoveNumberPartialInterfaces_.at( i ), parameter, currentTidalPartialOutput.second );
                        partialFunction = std::bind(
                                    &SphericalHarmonicsGravityPartial::wrtTidalModelParameter, this, coefficientPartialFunction,
                                    degree, orders, sumOrders, parameter->getParameterSize( ), std::placeholders::_1 );

                        numberOfRows = currentTidalPartialOutput.first;

                    }
                }
            }
        }
        else if( estimatable_parameters::isParameterNonTidalGravityFieldVariationProperty( parameter->getParameterName( ).first ) )
        {
            switch( parameter->getParameterName( ).first )
            {
            case polynomial_gravity_field_variation_amplitudes:
            {
                std::shared_ptr< PolynomialGravityFieldVariationsParameters > polynomialVariationParameter =
                    std::dynamic_pointer_cast< PolynomialGravityFieldVariationsParameters >( parameter );
                std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex =
                    polynomialVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
                std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex =
                    polynomialVariationParameter->getIndexAndPowerPerSineBlockIndex( );

                partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtPolynomialGravityFieldVariations, this,
                                             utilities::createVectorFromMapKeys( indexAndPowerPerCosineBlockIndex ),
                                             utilities::createVectorFromMapKeys( indexAndPowerPerSineBlockIndex ),
                                             utilities::createVectorFromMapValues( indexAndPowerPerCosineBlockIndex ),
                                             utilities::createVectorFromMapValues( indexAndPowerPerSineBlockIndex ),
                                             polynomialVariationParameter->getPolynomialVariationModel( )->getReferenceEpoch( ),
                                             std::placeholders::_1 );

                numberOfRows = parameter->getParameterSize( );
                break;
            }
            case periodic_gravity_field_variation_amplitudes:
            {
                std::shared_ptr< PeriodicGravityFieldVariationsParameters > periodicVariationParameter =
                    std::dynamic_pointer_cast< PeriodicGravityFieldVariationsParameters >( parameter );
                std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex =
                    periodicVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
                std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex =
                    periodicVariationParameter->getIndexAndPowerPerSineBlockIndex( );

                partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtPeriodicGravityFieldVariations, this,
                                             utilities::createVectorFromMapKeys( indexAndPowerPerCosineBlockIndex ),
                                             utilities::createVectorFromMapKeys( indexAndPowerPerSineBlockIndex ),
                                             utilities::createVectorFromMapValues( indexAndPowerPerCosineBlockIndex ),
                                             utilities::createVectorFromMapValues( indexAndPowerPerSineBlockIndex ),
                                             periodicVariationParameter->getPeriodicVariationModel( )->getFrequencies( ),
                                             periodicVariationParameter->getPeriodicVariationModel( )->getReferenceEpoch( ),
                                             std::placeholders::_1 );

                numberOfRows = parameter->getParameterSize( );
                break;
            }
            default:
                break;
            }
        }
        // Check non-rotational parameters.
        else
        {
            switch( parameter->getParameterName( ).first )
            {
            case spherical_harmonics_cosine_coefficient_block:
            {
                // Cast parameter object to required type.
                std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtCosineCoefficientBlock, this,
                                             coefficientsParameter->getBlockIndices( ), std::placeholders::_1 );
                numberOfRows = coefficientsParameter->getParameterSize( );

                break;
            }
            case spherical_harmonics_sine_coefficient_block:
            {
                // Cast parameter object to required type.

                std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtSineCoefficientBlock, this,
                                             coefficientsParameter->getBlockIndices( ), std::placeholders::_1 );
                numberOfRows = coefficientsParameter->getParameterSize( );

                break;
            }
            default:
                break;
            }
        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning the current partial w.r.t. a gravitational parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
SphericalHarmonicsGravityPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check for dependency
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                         this, std::placeholders::_1, 0 );
            numberOfColumns = 1;

        }

        if( parameterId.second.first == acceleratedBody_ )
        {
            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = std::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                             this, std::placeholders::_1, 0 );
                numberOfColumns = 1;
            }
        }
    }
    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function for updating the partial object to current state and time.
void SphericalHarmonicsGravityPartial::update( const double currentTime )
{
    using namespace tudat::coordinate_conversions;

    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        Eigen::Matrix3d currentRotationToBodyFixedFrame_ = fromBodyFixedToIntegrationFrameRotation_( ).inverse( );
        bodyFixedPosition_ = currentRotationToBodyFixedFrame_ *
                ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) );

        // Calculate spherical position in frame fixed to body exerting acceleration
        bodyFixedSphericalPosition_ = convertCartesianToSpherical( bodyFixedPosition_ );
        bodyFixedSphericalPosition_( 1 ) = mathematical_constants::PI / 2.0 - bodyFixedSphericalPosition_( 1 );

        // Get spherical harmonic coefficients
        currentCosineCoefficients_ = cosineCoefficients_( );
        currentSineCoefficients_ = sineCoefficients_( );

        // Update trogonometric functions of multiples of longitude.
        sphericalHarmonicCache_->update(
                    bodyFixedSphericalPosition_( 0 ), std::sin( bodyFixedSphericalPosition_( 1 ) ),
                    bodyFixedSphericalPosition_( 2 ), bodyReferenceRadius_( ) );

        // Calculate partial of acceleration wrt position of body undergoing acceleration.
        currentBodyFixedPartialWrtPosition_ = computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                    bodyFixedPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                    currentCosineCoefficients_, currentSineCoefficients_, sphericalHarmonicCache_ );

        currentPartialWrtVelocity_.setZero( );
        currentPartialWrtPosition_.setZero( );

        // Compute partial w.r.t. position in inertial frame
        currentPartialWrtPosition_ +=
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        // If rotation matrix depends on translational state, add correction partials
        if( rotationMatrixPartials_.count(
                    std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        {
            // Compute the acceleration and body-fixed position partial, without the central term (to avoid numerical errors)
            Eigen::Vector3d nonCentralAcceleration = accelerationFunction_( );
            Eigen::Matrix3d nonCentralBodyFixedPartial = currentBodyFixedPartialWrtPosition_;
            if( currentCosineCoefficients_( 0 ) > 0.0 )
            {
                nonCentralAcceleration -= gravitation::computeGravitationalAcceleration(
                            positionFunctionOfAcceleratedBody_( ), gravitationalParameterFunction_( ), positionFunctionOfAcceleratingBody_( ) );

                nonCentralBodyFixedPartial -=
                        currentRotationToBodyFixedFrame_ * calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                            positionFunctionOfAcceleratedBody_( ), positionFunctionOfAcceleratingBody_( ), gravitationalParameterFunction_( ) ) *
                        currentRotationToBodyFixedFrame_.inverse( );
            }

            // Compute rotation matrix partials
            std::vector< Eigen::Matrix3d > rotationPositionPartials =
                    rotationMatrixPartials_.at(
                        std::make_pair( estimatable_parameters::initial_body_state, "" ) )->
                    calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

            // Add correction terms to position and velocity partials
            for( unsigned int i = 0; i < 3; i++ )
            {
                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );

                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        currentRotationToBodyFixedFrame_.inverse( ) * nonCentralBodyFixedPartial *
                        ( rotationPositionPartials.at( i ).transpose( ) *
                          ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) ) );
                currentPartialWrtVelocity_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i + 3 ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );
//                        currentRotationToBodyFixedFrame_.inverse( ) * nonCentralBodyFixedPartial *
//                        ( rotationPositionPartials.at( i + 3 ).transpose( ) *
//                          ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) ) );
            }
        }

        currentTime_ = currentTime;

        // Update tidal interfaces
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->update( currentTime );
        }
    }

}

void SphericalHarmonicsGravityPartial::wrtPolynomialGravityFieldVariations(
    const std::vector< std::pair< int, int > >& cosineBlockIndices,
    const std::vector< std::pair< int, int > >& sineBlockIndices,
    const std::vector< std::vector< std::pair< int, int > > > powersPerCosineBlockIndex,
    const std::vector< std::vector< std::pair< int, int > > > powersPerSineBlockIndex,
    const double referenceEpoch,
    Eigen::MatrixXd& partialDerivatives )
{
    Eigen::MatrixXd staticCosinePartialsMatrix = Eigen::MatrixXd::Zero( 3, cosineBlockIndices.size( ) );
    Eigen::MatrixXd staticSinePartialsMatrix = Eigen::MatrixXd::Zero( 3, sineBlockIndices.size( ) );

    calculateSphericalHarmonicGravityWrtCCoefficients(
        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
        sphericalHarmonicCache_,
        cosineBlockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), staticCosinePartialsMatrix,
        maximumDegree_, maximumOrder_ );

    calculateSphericalHarmonicGravityWrtSCoefficients(
        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
        sphericalHarmonicCache_,
        sineBlockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), staticSinePartialsMatrix,
        maximumDegree_, maximumOrder_ );

    partialDerivatives.setZero( );

    int counter = 0;
    for( unsigned int i = 0; i < powersPerCosineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerCosineBlockIndex.at( i ).size( ); j++ )
        {
            partialDerivatives.block( 0, powersPerCosineBlockIndex.at( i ).at( j ).first, 3, 1 ) +=
                staticCosinePartialsMatrix.block( 0, i, 3, 1 ) * std::pow( ( currentTime_ - referenceEpoch), powersPerCosineBlockIndex.at( i ).at( j ).second );
            counter++;
        }
    }

    for( unsigned int i = 0; i < powersPerSineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerSineBlockIndex.at( i ).size( ); j++ )
        {
            partialDerivatives.block( 0, powersPerSineBlockIndex.at( i ).at( j ).first + counter, 3, 1 ) +=
                staticSinePartialsMatrix.block( 0, i, 3, 1 ) * std::pow( ( currentTime_ - referenceEpoch), powersPerSineBlockIndex.at( i ).at( j ).second );
        }
    }
}


void SphericalHarmonicsGravityPartial::wrtPeriodicGravityFieldVariations(
    const std::vector< std::pair< int, int > >& cosineBlockIndices,
    const std::vector< std::pair< int, int > >& sineBlockIndices,
    const std::vector< std::vector< std::pair< int, int > > > powersPerCosineBlockIndex,
    const std::vector< std::vector< std::pair< int, int > > > powersPerSineBlockIndex,
    const std::vector< double >& frequencies,
    const double referenceEpoch,
    Eigen::MatrixXd& partialDerivatives )
{
    Eigen::MatrixXd staticCosinePartialsMatrix = Eigen::MatrixXd::Zero( 3, cosineBlockIndices.size( ) );
    Eigen::MatrixXd staticSinePartialsMatrix = Eigen::MatrixXd::Zero( 3, sineBlockIndices.size( ) );

    calculateSphericalHarmonicGravityWrtCCoefficients(
        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
        sphericalHarmonicCache_,
        cosineBlockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), staticCosinePartialsMatrix,
        maximumDegree_, maximumOrder_ );

    calculateSphericalHarmonicGravityWrtSCoefficients(
        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
        sphericalHarmonicCache_,
        sineBlockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), staticSinePartialsMatrix,
        maximumDegree_, maximumOrder_ );

    partialDerivatives.setZero( );

    int counter = 0;
    for( unsigned int i = 0; i < powersPerCosineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerCosineBlockIndex.at( i ).size( ); j++ )
        {
            double frequency = frequencies.at( powersPerCosineBlockIndex.at( i ).at( j ).second );
            partialDerivatives.block( 0, 2 * powersPerCosineBlockIndex.at( i ).at( j ).first, 3, 1 ) +=
                staticCosinePartialsMatrix.block( 0, i, 3, 1 ) * std::cos( frequency * ( currentTime_ - referenceEpoch ) );
            partialDerivatives.block( 0, 2 * powersPerCosineBlockIndex.at( i ).at( j ).first + 1, 3, 1 ) +=
                staticCosinePartialsMatrix.block( 0, i, 3, 1 ) * std::sin( frequency * ( currentTime_ - referenceEpoch ) );
            counter++;
        }
    }

    for( unsigned int i = 0; i < powersPerSineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerSineBlockIndex.at( i ).size( ); j++ )
        {
            double frequency = frequencies.at( powersPerSineBlockIndex.at( i ).at( j ).second );
            partialDerivatives.block( 0, 2 * ( powersPerSineBlockIndex.at( i ).at( j ).first + counter ), 3, 1 ) +=
                staticSinePartialsMatrix.block( 0, i, 3, 1 ) * std::cos( frequency * ( currentTime_ - referenceEpoch ) );
            partialDerivatives.block( 0, 2 * ( powersPerSineBlockIndex.at( i ).at( j ).first + counter ) + 1, 3, 1 ) +=
                staticSinePartialsMatrix.block( 0, i, 3, 1 ) * std::sin( frequency * ( currentTime_ - referenceEpoch ) );
        }
    }
}


//! Function to calculate the partial of the acceleration wrt a set of cosine coefficients.
void SphericalHarmonicsGravityPartial::wrtCosineCoefficientBlock(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialDerivatives )
{
    calculateSphericalHarmonicGravityWrtCCoefficients(
                bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                sphericalHarmonicCache_,
                blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                    bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), partialDerivatives,
                maximumDegree_, maximumOrder_ );
}


//! Function to calculate the partial of the acceleration wrt a set of sine coefficients.
void SphericalHarmonicsGravityPartial::wrtSineCoefficientBlock(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialDerivatives )
{
    calculateSphericalHarmonicGravityWrtSCoefficients(
                bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                sphericalHarmonicCache_,
                blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                    bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), partialDerivatives,
                maximumDegree_, maximumOrder_ );
}

//! Function to calculate an acceleration partial wrt a rotational parameter.
void SphericalHarmonicsGravityPartial::wrtRotationModelParameter(
        Eigen::MatrixXd& accelerationPartial,
        const estimatable_parameters::EstimatebleParametersEnum parameterType,
        const std::string& secondaryIdentifier )
{
    // Calculate distance vector between bodies.
    Eigen::Vector3d distanceVector = positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( );

    // Get rotation matrix partial(s) wrt requested parameter
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )->
            calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );

    // Iterate for each single parameter entry partial.
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        // Calculate acceleration partial for current parameter entry.
        accelerationPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] *
                ( fromBodyFixedToIntegrationFrameRotation_( ).inverse( ) ) * accelerationFunction_( ) +
                fromBodyFixedToIntegrationFrameRotation_( ) * currentBodyFixedPartialWrtPosition_*
                rotationMatrixPartials[ i ].transpose( ) * distanceVector;
    }
}

//! Function to calculate an acceleration partial wrt a tidal parameter.
void SphericalHarmonicsGravityPartial::wrtTidalModelParameter(
        const std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunctions,
        const int degree,
        const std::vector< int >& orders,
        const bool sumOrders,
        const int parameterSize,
        Eigen::MatrixXd& partialMatrix )
{
    // Initialize partial matrix to zero values.
    partialMatrix = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, parameterSize );

    // Calculate multiplicative term found in all partial terms (partial of C,S coefficients w.r.t. parameter).
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > coefficientPartialsPerOrder_ = coefficientPartialFunctions( );
    int singleOrderPartialSize = coefficientPartialsPerOrder_.at( 0 ).cols( );

    Eigen::MatrixXd currentPartialContribution = Eigen::MatrixXd::Zero( 3, 1 );
    std::vector< std::pair< int, int > > blockIndices;
    blockIndices.resize( 1 );

    // Iterate over all required orders in current degree.
    for( unsigned int i = 0; i < orders.size( ); i++ )
    {
        // Set coefficient degree/order for current partials
        blockIndices[ 0 ] = std::make_pair( degree, orders.at( i ) );

        // Compute acceleration w.r.t. C and S coefficients, and multiply with partials of C,S coefficients w.r.t. parameter
        if( sumOrders )
        {
            calculateSphericalHarmonicGravityWrtCCoefficients(
                        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                        sphericalHarmonicCache_,
                        blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), currentPartialContribution,
                        maximumDegree_, maximumOrder_  );

            partialMatrix.block( 0, 0, 3, singleOrderPartialSize ) +=
                    currentPartialContribution * coefficientPartialsPerOrder_.at( i ).block( 0, 0, 1, singleOrderPartialSize );


            blockIndices[ 0 ] = std::make_pair( degree, orders.at( i ) );
            calculateSphericalHarmonicGravityWrtSCoefficients(
                        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                        sphericalHarmonicCache_,
                        blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), currentPartialContribution,
                        maximumDegree_, maximumOrder_  );

            partialMatrix.block( 0, 0, 3, singleOrderPartialSize ) +=
                    currentPartialContribution * coefficientPartialsPerOrder_.at( i ).block( 1, 0, 1, singleOrderPartialSize );
        }
        else
        {
            calculateSphericalHarmonicGravityWrtCCoefficients(
                        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                        sphericalHarmonicCache_,
                        blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), currentPartialContribution,
                        maximumDegree_, maximumOrder_  );

            partialMatrix.block( 0, i * singleOrderPartialSize, 3, singleOrderPartialSize ) +=
                    currentPartialContribution * coefficientPartialsPerOrder_.at( i ).block( 0, 0, 1, singleOrderPartialSize );

            calculateSphericalHarmonicGravityWrtSCoefficients(
                        bodyFixedSphericalPosition_, bodyReferenceRadius_( ), gravitationalParameterFunction_( ),
                        sphericalHarmonicCache_,
                        blockIndices, coordinate_conversions::getSphericalToCartesianGradientMatrix(
                            bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), currentPartialContribution,
                        maximumDegree_, maximumOrder_  );

            partialMatrix.block( 0, i * singleOrderPartialSize, 3, singleOrderPartialSize ) +=
                    currentPartialContribution * coefficientPartialsPerOrder_.at( i ).block( 1, 0, 1, singleOrderPartialSize );

        }
    }
}

}

}

