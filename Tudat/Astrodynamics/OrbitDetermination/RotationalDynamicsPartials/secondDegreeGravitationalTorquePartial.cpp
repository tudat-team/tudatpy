/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/secondDegreeGravitationalTorquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"


namespace tudat
{

namespace acceleration_partials
{

//! Function that computes the partial of degree 2 torque w.r.t. the current quaternion elements
Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
        const double premultiplier,
        const Eigen::Matrix3d& inertiaTensor,
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const Eigen::Vector3d& inertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions )
{
    Eigen::Matrix3d scalingMatrix = linear_algebra::getCrossProductMatrix(
                bodyFixedRelativePosition ) *  inertiaTensor - linear_algebra::getCrossProductMatrix(
                inertiaTensor * bodyFixedRelativePosition );
    Eigen::Matrix< double, 3, 4 > partialOfBodyFixedPositionWrtQuaternion = Eigen::Matrix< double, 3, 4 >::Zero( );
    for( unsigned int i = 0; i < derivativeOfRotationMatrixWrtQuaternions.size( ); i++ )
    {
        partialOfBodyFixedPositionWrtQuaternion.block( 0, i, 3, 1 ) =
                derivativeOfRotationMatrixWrtQuaternions.at( i ).transpose( ) * inertialRelativePosition;
    }
    return premultiplier * scalingMatrix * partialOfBodyFixedPositionWrtQuaternion;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
SecondDegreeGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter && (
                parameter->getParameterName( ).second.first == bodyExertingTorque_ ) )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair = std::make_pair(
                    std::bind( &SecondDegreeGravitationalTorquePartial::wrtGravitationalParameterOfCentralBody,
                                 this, std::placeholders::_1 ), 1 );
    }
    else
    {
        partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > SecondDegreeGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >  partialFunction = std::make_pair(
                std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case spherical_harmonics_cosine_coefficient_block:
        {
            // Cast parameter object to required type.
            std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

            int c20Index, c21Index, c22Index;
            coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );

            if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
            {
                if( ( getInertiaTensorNormalizationFactor_ == nullptr ) )
                {
                    throw std::runtime_error( "Error when getting partial of 2nd degree grac torque w.r.t. cosine sh parameters, inertia tensor normalization function not found." );
                }
                partialFunction = std::make_pair(
                            std::bind( &SecondDegreeGravitationalTorquePartial::
                                         wrtCosineSphericalHarmonicCoefficientsOfCentralBody, this,
                                         std::placeholders::_1, c20Index, c21Index, c22Index ), coefficientsParameter->getParameterSize( ) );
            }

            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            // Cast parameter object to required type.

            std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

            int s21Index, s22Index;
            coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );

            if( s21Index >= 0 || s22Index >= 0 )
            {
                if( getInertiaTensorNormalizationFactor_ == nullptr )
                {
                    throw std::runtime_error( "Error when getting partial of 2nd degree grac torque w.r.t. sine sh parameters, inertia tensor normalization function not found." );
                }
                partialFunction = std::make_pair(
                            std::bind( &SecondDegreeGravitationalTorquePartial::
                                         wrtSineSphericalHarmonicCoefficientsOfCentralBody, this,
                                         std::placeholders::_1, s21Index, s22Index ), coefficientsParameter->getParameterSize( ) );
            }


            break;
        }
        default:
            break;
        }
    }
    return partialFunction;
}

//! Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
void SecondDegreeGravitationalTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( ( stateReferencePoint.first == bodyExertingTorque_ ||
          stateReferencePoint.first == bodyUndergoingTorque_ ) && integratedStateType == propagators::translational_state )
    {
        partialMatrix.block( 0, 0, 3, 3 ) +=
                ( ( stateReferencePoint.first == bodyExertingTorque_ ) ? 1.0 : -1.0 ) *
                ( torqueModel_->getCurrentTorqueMagnitudePremultiplier( ) *
                  ( linear_algebra::getCrossProductMatrix( currentBodyFixedRelativePosition_ ) *
                      torqueModel_->getCurrentInertiaTensorOfRotatingBody( ) -
                      linear_algebra::getCrossProductMatrix( torqueModel_->getCurrentInertiaTensorTimesRelativePositionOfBody( ) ) ) *
                  ( torqueModel_->getCurrentRotationToBodyFixedFrame( ) ).toRotationMatrix( ) -
                  5.0 * torqueModel_->getTorque( ) *
                  torqueModel_->getCurrentRelativePositionOfBodySubjectToTorque( ).normalized( ).transpose( ) /
                  currentBodyFixedRelativePosition_.norm( ) );
    }
}

//! Update partial model to current time
void SecondDegreeGravitationalTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );

        currentQuaternionVector_ = linear_algebra::convertQuaternionToVectorFormat(
                    ( torqueModel_->getCurrentRotationToBodyFixedFrame( ) ).inverse( ) );
        linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
                    currentQuaternionVector_,  currentRotationMatrixDerivativesWrtQuaternion_ );
        currentBodyFixedRelativePosition_ = torqueModel_->getCurrentRelativeBodyFixedPositionOfBodySubjectToTorque( );
        currentCoefficientPartialPremultiplier_ = torqueModel_->getCurrentTorqueMagnitudePremultiplier( ) *
                linear_algebra::getCrossProductMatrix( currentBodyFixedRelativePosition_ );

        currentPartialDerivativeWrtQuaternion_ = getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
                    torqueModel_->getCurrentTorqueMagnitudePremultiplier( ),
                    torqueModel_->getCurrentInertiaTensorOfRotatingBody( ),
                    torqueModel_->getCurrentRelativeBodyFixedPositionOfBodySubjectToTorque(),
                    torqueModel_->getCurrentRelativePositionOfBodySubjectToTorque( ),
                    currentRotationMatrixDerivativesWrtQuaternion_ );
    }
}


//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void SecondDegreeGravitationalTorquePartial::wrtGravitationalParameterOfCentralBody(
        Eigen::MatrixXd& gravitationalParameterPartial )
{
    if( torqueModel_->getCurrentGravitationalParameterOfAttractingBody( ) == 0.0 )
    {
        throw std::runtime_error( "Error when calculating partial of SecondDegreeGravitationalTorquePartial w.r.t. mu: mu=0." );
    }
    gravitationalParameterPartial =
            torqueModel_->getTorque( ) / torqueModel_->getCurrentGravitationalParameterOfAttractingBody( );
}

//! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
void SecondDegreeGravitationalTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody(
        Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
        const int c20Index, const int c21Index, const int c22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( c20Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c20Index, 3, 1 ) =
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) *
                getInertiaTensorNormalizationFactor_( ) * currentCoefficientPartialPremultiplier_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 * currentBodyFixedRelativePosition_;
    }

    if( c21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c21Index, 3, 1 ) =
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) *
                getInertiaTensorNormalizationFactor_( ) * currentCoefficientPartialPremultiplier_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 * currentBodyFixedRelativePosition_;
    }

    if( c22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, c22Index, 3, 1 ) =
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) *
                getInertiaTensorNormalizationFactor_( ) * currentCoefficientPartialPremultiplier_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 * currentBodyFixedRelativePosition_;
    }
}

//! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
void SecondDegreeGravitationalTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody(
        Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
        const int s21Index, const int s22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( s21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, s21Index, 3, 1 ) =
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) *
                getInertiaTensorNormalizationFactor_( ) * currentCoefficientPartialPremultiplier_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 * currentBodyFixedRelativePosition_;
    }

    if( s22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial .block( 0, s22Index, 3, 1 ) =
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) *
                getInertiaTensorNormalizationFactor_( ) * currentCoefficientPartialPremultiplier_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 * currentBodyFixedRelativePosition_;
    }
}


}

}
