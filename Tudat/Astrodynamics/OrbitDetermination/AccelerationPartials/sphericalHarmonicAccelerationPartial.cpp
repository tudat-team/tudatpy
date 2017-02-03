/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicSineCoefficients.h"

namespace tudat
{

namespace acceleration_partials
{

//! Contructor.
SphericalHarmonicsGravityPartial::SphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::spherical_harmonic_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    bodyReferenceRadius_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getReferenceRadius,
                                       accelerationModel ) ),
    cosineCoefficients_( accelerationModel->getCosineHarmonicCoefficientsFunction( ) ),
    sineCoefficients_( accelerationModel->getSineHarmonicCoefficientsFunction( ) ),
    sphericalHarmonicCache_( accelerationModel->getSphericalHarmonicsCache( ) ),
    positionFunctionOfAcceleratedBody_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                     getCurrentPositionOfBodySubjectToAcceleration, accelerationModel ) ),
    positionFunctionOfAcceleratingBody_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                      getCurrentPositionOfBodyExertingAcceleration, accelerationModel ) ),
    fromBodyFixedToIntegrationFrameRotation_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                                      getCurrentRotationToIntegrationFrameMatrix, accelerationModel ) ),
    accelerationFunction_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAcceleration,
                                        accelerationModel ) ),
    updateFunction_( boost::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::updateMembers,
                                  accelerationModel, _1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials ),
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
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > SphericalHarmonicsGravityPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;


    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair =
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
                partialFunction = boost::bind(
                            &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                                               this, _1, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = 1;
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of sh acceleration wrt rotational parameter" +
                        boost::lexical_cast< std::string >( parameter->getParameterName( ).first ) + " of " +
                        parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
            }

        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning a partial w.r.t. a vector parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > SphericalHarmonicsGravityPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
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
                partialFunction = boost::bind(
                            &SphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                            this,_1, parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) );
                numberOfRows = parameter->getParameterSize( );
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of sh acceleration wrt rotational parameter" +
                        boost::lexical_cast< std::string >( parameter->getParameterName( ).first ) + " of " +
                        parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
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
                boost::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        boost::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtCosineCoefficientBlock, this,
                                               coefficientsParameter->getBlockIndices( ), _1 );
                numberOfRows = coefficientsParameter->getParameterSize( );

                break;
            }
            case spherical_harmonics_sine_coefficient_block:
            {
                // Cast parameter object to required type.

                boost::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        boost::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtSineCoefficientBlock, this,
                                               coefficientsParameter->getBlockIndices( ), _1 );
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
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
SphericalHarmonicsGravityPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check for dependency
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                           this, _1 );
            numberOfColumns = 1;

        }

        if( parameterId.second.first == acceleratedBody_ )
        {
            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = boost::bind( &SphericalHarmonicsGravityPartial::wrtGravitationalParameterOfCentralBody,
                                               this, _1 );
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
        currentPartialWrtPosition_ =
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ *
                currentRotationToBodyFixedFrame_;

        currentTime_ = currentTime;
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
                    bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), partialDerivatives );
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
                    bodyFixedPosition_ ), fromBodyFixedToIntegrationFrameRotation_( ), partialDerivatives );
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

}

}

