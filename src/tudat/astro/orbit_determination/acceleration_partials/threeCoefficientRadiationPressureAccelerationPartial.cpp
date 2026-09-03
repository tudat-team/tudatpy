/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/threeCoefficientRadiationPressureAccelerationPartial.h"

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/threeCoefficientRadiationPressureCoefficients.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{
namespace acceleration_partials
{

using linear_algebra::getCrossProductMatrix;

ThreeCoefficientRadiationPressureAccelerationPartial::ThreeCoefficientRadiationPressureAccelerationPartial(
        const std::shared_ptr< electromagnetism::ThreeCoefficientRadiationPressureAcceleration >& accelerationModel,
        const std::string& acceleratedBody,
        const std::string& acceleratingBody ):
    AccelerationPartial( acceleratedBody, acceleratingBody, accelerationModel, basic_astrodynamics::three_coefficient_radiation_pressure ),
    accelerationModel_( accelerationModel ), referenceBody_( accelerationModel->getReferenceBodyName( ) ),
    currentPartialWrtTargetPosition_( Eigen::Matrix3d::Zero( ) ), currentPartialWrtSourcePosition_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtSourceVelocity_( Eigen::Matrix3d::Zero( ) ), currentPartialWrtReferencePosition_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtReferenceVelocity_( Eigen::Matrix3d::Zero( ) )
{}

void ThreeCoefficientRadiationPressureAccelerationPartial::addPartial( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                       const Eigen::Matrix3d& partial,
                                                                       const bool addContribution,
                                                                       const int startRow,
                                                                       const int startColumn ) const
{
    partialMatrix.block( startRow, startColumn, 3, 3 ) += ( addContribution ? 1.0 : -1.0 ) * partial;
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                         const bool addContribution,
                                                                                         const int startRow,
                                                                                         const int startColumn )
{
    Eigen::Matrix3d partial = currentPartialWrtTargetPosition_;
    if( referenceBody_ == acceleratedBody_ )
    {
        partial += currentPartialWrtReferencePosition_;
    }
    addPartial( partialMatrix, partial, addContribution, startRow, startColumn );
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                         const bool addContribution,
                                                                                         const int startRow,
                                                                                         const int startColumn )
{
    if( referenceBody_ == acceleratedBody_ )
    {
        addPartial( partialMatrix, currentPartialWrtReferenceVelocity_, addContribution, startRow, startColumn );
    }
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                          const bool addContribution,
                                                                                          const int startRow,
                                                                                          const int startColumn )
{
    addPartial( partialMatrix, currentPartialWrtSourcePosition_, addContribution, startRow, startColumn );
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                          const bool addContribution,
                                                                                          const int startRow,
                                                                                          const int startColumn )
{
    addPartial( partialMatrix, currentPartialWrtSourceVelocity_, addContribution, startRow, startColumn );
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtPositionOfAdditionalBody( const std::string& bodyName,
                                                                                        Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                        const bool addContribution,
                                                                                        const int startRow,
                                                                                        const int startColumn )
{
    if( bodyName == referenceBody_ )
    {
        addPartial( partialMatrix, currentPartialWrtReferencePosition_, addContribution, startRow, startColumn );
    }
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtVelocityOfAdditionalBody( const std::string& bodyName,
                                                                                        Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                                                        const bool addContribution,
                                                                                        const int startRow,
                                                                                        const int startColumn )
{
    if( bodyName == referenceBody_ )
    {
        addPartial( partialMatrix, currentPartialWrtReferenceVelocity_, addContribution, startRow, startColumn );
    }
}

bool ThreeCoefficientRadiationPressureAccelerationPartial::isAccelerationPartialWrtAdditionalBodyNonnullptr( const std::string& bodyName )
{
    return bodyName == referenceBody_ && bodyName != acceleratedBody_ && bodyName != acceleratingBody_;
}

bool ThreeCoefficientRadiationPressureAccelerationPartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    return stateReferencePoint.first == acceleratedBody_ && stateReferencePoint.second.empty( ) &&
            integratedStateType == propagators::body_mass_state;
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtNonTranslationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType,
        const bool addContribution )
{
    if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
    {
        partialMatrix.block( 0, 0, 3, 1 ) +=
                ( addContribution ? -1.0 : 1.0 ) * accelerationModel_->getAcceleration( ) / accelerationModel_->getCurrentTargetMass( );
    }
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
ThreeCoefficientRadiationPressureAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    if( parameter->getParameterName( ).first == estimatable_parameters::three_coefficient_radiation_pressure_coefficients &&
        parameter->getParameterName( ).second.first == acceleratedBody_ &&
        parameter->getParameterName( ).second.second == acceleratingBody_ )
    {
        return std::make_pair(
                std::bind( &ThreeCoefficientRadiationPressureAccelerationPartial::wrtThreeCoefficientRadiationPressureCoefficients,
                           this,
                           std::placeholders::_1 ),
                3 );
    }
    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

void ThreeCoefficientRadiationPressureAccelerationPartial::wrtThreeCoefficientRadiationPressureCoefficients( Eigen::MatrixXd& partial )
{
    partial = accelerationModel_->getAccelerationScalingFactor( ) * accelerationModel_->getCurrentAccelerationScalingFactor( ) *
            accelerationModel_->getCurrentBasis( );
}

void ThreeCoefficientRadiationPressureAccelerationPartial::update( const double currentTime )
{
    accelerationModel_->updateMembers( currentTime );

    const Eigen::Vector3d sourceToTarget = accelerationModel_->getCurrentSourceToTargetPosition( );
    currentPartialWrtTargetPosition_ =
            -2.0 * accelerationModel_->getAcceleration( ) * sourceToTarget.transpose( ) / sourceToTarget.squaredNorm( );

    const Eigen::Vector3d sourceToReference = accelerationModel_->getCurrentSourceToReferencePosition( );
    const Eigen::Vector3d sourceToReferenceVelocity = accelerationModel_->getCurrentSourceToReferenceVelocity( );
    const Eigen::Matrix3d basis = accelerationModel_->getCurrentBasis( );
    const Eigen::Vector3d uDirection = basis.col( 0 );
    const Eigen::Vector3d wDirection = basis.col( 2 );
    const Eigen::Vector3d orbitAngularMomentum = sourceToReference.cross( sourceToReferenceVelocity );
    const Eigen::Vector3d orbitNormal = orbitAngularMomentum.normalized( );

    const Eigen::Matrix3d uPartialWrtPosition =
            -( Eigen::Matrix3d::Identity( ) - uDirection * uDirection.transpose( ) ) / sourceToReference.norm( );
    const Eigen::Matrix3d orbitNormalPartialWrtAngularMomentum =
            ( Eigen::Matrix3d::Identity( ) - orbitNormal * orbitNormal.transpose( ) ) / orbitAngularMomentum.norm( );
    const Eigen::Matrix3d orbitNormalPartialWrtPosition =
            orbitNormalPartialWrtAngularMomentum * ( -getCrossProductMatrix( sourceToReferenceVelocity ) );
    const Eigen::Matrix3d orbitNormalPartialWrtVelocity = orbitNormalPartialWrtAngularMomentum * getCrossProductMatrix( sourceToReference );

    const double obliquity = unit_conversions::convertDegreesToRadians( 23.4 );
    const Eigen::Matrix3d wPartialWrtOrbitNormal =
            std::cos( obliquity ) * Eigen::Matrix3d::Identity( ) + std::sin( obliquity ) * getCrossProductMatrix( uDirection );
    const Eigen::Matrix3d wPartialWrtU = -std::sin( obliquity ) * getCrossProductMatrix( orbitNormal );
    const Eigen::Matrix3d wPartialWrtPosition = wPartialWrtOrbitNormal * orbitNormalPartialWrtPosition + wPartialWrtU * uPartialWrtPosition;
    const Eigen::Matrix3d wPartialWrtVelocity = wPartialWrtOrbitNormal * orbitNormalPartialWrtVelocity;

    const Eigen::Matrix3d vPartialWrtPosition =
            -getCrossProductMatrix( uDirection ) * wPartialWrtPosition + getCrossProductMatrix( wDirection ) * uPartialWrtPosition;
    const Eigen::Matrix3d vPartialWrtVelocity = -getCrossProductMatrix( uDirection ) * wPartialWrtVelocity;

    const Eigen::Vector3d coefficients = accelerationModel_->getCoefficients( );
    const Eigen::Matrix3d resolvedAccelerationPartialWrtPosition =
            coefficients( 0 ) * uPartialWrtPosition + coefficients( 1 ) * vPartialWrtPosition + coefficients( 2 ) * wPartialWrtPosition;
    const Eigen::Matrix3d resolvedAccelerationPartialWrtVelocity =
            coefficients( 1 ) * vPartialWrtVelocity + coefficients( 2 ) * wPartialWrtVelocity;
    const double scale = accelerationModel_->getAccelerationScalingFactor( ) * accelerationModel_->getCurrentAccelerationScalingFactor( );
    currentPartialWrtReferencePosition_ = scale * resolvedAccelerationPartialWrtPosition;
    currentPartialWrtReferenceVelocity_ = scale * resolvedAccelerationPartialWrtVelocity;

    currentPartialWrtSourcePosition_ = -currentPartialWrtTargetPosition_ - currentPartialWrtReferencePosition_;
    currentPartialWrtSourceVelocity_ = -currentPartialWrtReferenceVelocity_;
    currentTime_ = currentTime;
}

}  // namespace acceleration_partials
}  // namespace tudat
