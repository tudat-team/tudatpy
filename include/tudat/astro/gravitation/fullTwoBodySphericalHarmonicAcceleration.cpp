#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicAcceleration.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

namespace tudat
{

namespace gravitation
{

namespace
{

bool isFullTwoBodyTorqueDebugEnabled( )
{
    static const bool isEnabled = [ ]( )
    {
        const char* flag = std::getenv( "TUDAT_DEBUG_FULL_TWO_BODY_TORQUE" );
        return ( flag != nullptr && std::string( flag ) != "0" );
    }( );
    return isEnabled;
}

std::string debugStatusFromDifference( const double difference, const double tolerance )
{
    return ( std::fabs( difference ) <= tolerance ? "OK" : "MISMATCH" );
}

bool isApproximatelyIdentityQuaternion( const Eigen::Quaterniond& quaternion, const double tolerance )
{
    const Eigen::Quaterniond identityQuaternion = Eigen::Quaterniond::Identity( );
    const double positiveDifferenceNorm = ( quaternion.coeffs( ) - identityQuaternion.coeffs( ) ).norm( );
    const double negativeDifferenceNorm = ( quaternion.coeffs( ) + identityQuaternion.coeffs( ) ).norm( );
    return ( std::min( positiveDifferenceNorm, negativeDifferenceNorm ) <= tolerance );
}

}


FullTwoBodySphericalHarmonicAcceleration::FullTwoBodySphericalHarmonicAcceleration(
        const std::function< Eigen::Vector3d( ) > positionOfBody1Function,
        const std::function< Eigen::Vector3d( ) > positionOfBody2Function,
        const std::function< double( ) > gravitationalParameterFunction,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const std::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody1Function,
        const std::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody1Function,
        const std::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody2Function,
        const std::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody2Function,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation,
        const std::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation,
        const bool useCentraBodyFrame,
        const bool areCoefficientsNormalized ):
    positionOfBody1Function_( positionOfBody1Function ), positionOfBody2Function_( positionOfBody2Function ),
    gravitationalParameterFunction_( gravitationalParameterFunction ),
    equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ), equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ),
    coefficientCombinationsToUse_( coefficientCombinationsToUse ),
    toLocalFrameOfBody1Transformation_( toLocalFrameOfBody1Transformation ),
    toLocalFrameOfBody2Transformation_( toLocalFrameOfBody2Transformation ),
    useCentraBodyFrame_( useCentraBodyFrame ),
    areCoefficientsNormalized_( areCoefficientsNormalized )
{
    maximumDegree_ = 0;
    maximumOrder_ = 0;

    unsigned int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));

        if( degreeOfBody1 + degreeOfBody2 > maximumDegree_ )
        {
            maximumDegree_ = degreeOfBody1 + degreeOfBody2;
        }

        if( orderOfBody1 + orderOfBody2 > maximumOrder_ )
        {
            maximumOrder_ = orderOfBody1 + orderOfBody2;
        }
    }

    // The geodesy-normalized acceleration evaluation precomputes P_lm up to m=l+1 for all l<=maximumDegree_.
    // Ensure the Legendre cache order budget covers this, even when selected interaction orders are low (e.g. C20-only).
    maximumOrder_ = std::max( maximumOrder_, maximumDegree_ );

    sphericalHarmonicsCache_ = std::make_shared< basic_mathematics::SphericalHarmonicsCache >(
                maximumDegree_ + 1, maximumOrder_ + 1 );
    effectiveMutualPotentialField_ =  std::make_shared< EffectiveMutualSphericalHarmonicsField >(
                coefficientCombinationsToUse_,
                cosineHarmonicCoefficientsOfBody1Function, sineHarmonicCoefficientsOfBody1Function,
                cosineHarmonicCoefficientsOfBody2Function, sineHarmonicCoefficientsOfBody2Function,
                gravitationalParameterFunction, equatorialRadiusOfBody1_, equatorialRadiusOfBody2_, areCoefficientsNormalized );
    effectiveCosineCoefficientFunction_ = std::bind(
                &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                effectiveMutualPotentialField_,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3,
                std::placeholders::_4 );
    effectiveSineCoefficientFunction_ = std::bind(
                &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                effectiveMutualPotentialField_,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3,
                std::placeholders::_4 );
    radius1Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree1( ) + 1 );
    radius2Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree2( ) + 1 );

}

void FullTwoBodySphericalHarmonicAcceleration::updateMembers( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        currentRotationFromInertialToBody1_ = toLocalFrameOfBody1Transformation_( );
        currentRotationFromBody2ToBody1_ =
                currentRotationFromInertialToBody1_ * toLocalFrameOfBody2Transformation_( ).inverse( );

        currentRelativePosition_ = positionOfBody1Function_( ) - positionOfBody2Function_( );
        currentBodyFixedRelativePosition_ =
                currentRotationFromInertialToBody1_ * ( currentRelativePosition_ );

        if( isFullTwoBodyTorqueDebugEnabled( ) )
        {
            std::cout << "[FTB-DBG][STEP acceleration_state]"
                      << " t=" << currentTime
                      << " mu=" << gravitationalParameterFunction_( )
                      << " r_inertial=" << currentRelativePosition_.transpose( )
                      << " r_body1=" << currentBodyFixedRelativePosition_.transpose( )
                      << " q_I_to_B1=(" << currentRotationFromInertialToBody1_.w( ) << ","
                      << currentRotationFromInertialToBody1_.x( ) << ","
                      << currentRotationFromInertialToBody1_.y( ) << ","
                      << currentRotationFromInertialToBody1_.z( ) << ")"
                      << " q_B2_to_B1=(" << currentRotationFromBody2ToBody1_.w( ) << ","
                      << currentRotationFromBody2ToBody1_.x( ) << ","
                      << currentRotationFromBody2ToBody1_.y( ) << ","
                      << currentRotationFromBody2ToBody1_.z( ) << ")"
                      << std::endl;

            if( isApproximatelyIdentityQuaternion( currentRotationFromInertialToBody1_, 1.0E-15 ) )
            {
                const Eigen::Vector3d frameDifference = currentBodyFixedRelativePosition_ - currentRelativePosition_;
                std::cout << "[FTB-DBG][STEP acceleration_state_expected]"
                          << " identity_rotation_expected_r_body1=r_inertial"
                          << " diff=" << frameDifference.transpose( )
                          << " status=" << debugStatusFromDifference( frameDifference.norm( ), 1.0E-15 )
                          << std::endl;
            }
        }

        effectiveMutualPotentialField_->computeCurrentEffectiveCoefficients( currentRotationFromBody2ToBody1_ );

        double currentDistance = currentRelativePosition_.norm( );
        for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree1( ); i++ )
        {
            radius1Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody1_ / currentDistance, i );

        }

        for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree2( ); i++ )
        {
            radius2Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody2_ / currentDistance, i );

        }

        if( isFullTwoBodyTorqueDebugEnabled( ) )
        {
            std::cout << "[FTB-DBG][STEP radius_powers]"
                      << " distance=" << currentDistance
                      << " (R1/r)^0..2="
                      << radius1Powers_.at( 0 ) << " "
                      << ( radius1Powers_.size( ) > 1 ? radius1Powers_.at( 1 ) : TUDAT_NAN ) << " "
                      << ( radius1Powers_.size( ) > 2 ? radius1Powers_.at( 2 ) : TUDAT_NAN )
                      << " (R2/r)^0..2="
                      << radius2Powers_.at( 0 ) << " "
                      << ( radius2Powers_.size( ) > 1 ? radius2Powers_.at( 1 ) : TUDAT_NAN ) << " "
                      << ( radius2Powers_.size( ) > 2 ? radius2Powers_.at( 2 ) : TUDAT_NAN )
                      << std::endl;
        }


        if( areCoefficientsNormalized_ )
        {
            mutualPotentialGradient_ = computeGeodesyNormalizedMutualGravitationalAccelerationSum(
                        currentBodyFixedRelativePosition_, gravitationalParameterFunction_( ),
                        equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                        effectiveCosineCoefficientFunction_, effectiveSineCoefficientFunction_,
                        coefficientCombinationsToUse_,
                        effectiveMutualPotentialField_->getMaximumDegree1( ),
                        effectiveMutualPotentialField_->getMaximumDegree2( ),
                        maximumDegree_,
                        radius1Powers_,
                        radius2Powers_,
                        sphericalHarmonicsCache_ );
        }
        else
        {            
            mutualPotentialGradient_ = computeUnnormalizedMutualGravitationalAccelerationSum(
                        currentBodyFixedRelativePosition_, gravitationalParameterFunction_( ),
                        equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                        effectiveCosineCoefficientFunction_, effectiveSineCoefficientFunction_,
                        coefficientCombinationsToUse_, sphericalHarmonicsCache_ );
        }

        if( useCentraBodyFrame_ )
        {
            currentAcceleration_ = mutualPotentialGradient_;
        }
        else
        {
            currentAcceleration_ = currentRotationFromInertialToBody1_.inverse( ) * mutualPotentialGradient_;
        }

        if( isFullTwoBodyTorqueDebugEnabled( ) )
        {
            const Eigen::Vector3d totalSpecificTorqueFromGradient =
                    currentBodyFixedRelativePosition_.cross( mutualPotentialGradient_ );
            std::cout << "[FTB-DBG][STEP acceleration_gradient]"
                      << " gradU=" << mutualPotentialGradient_.transpose( )
                      << " total_specific_torque=rXgradU=" << totalSpecificTorqueFromGradient.transpose( )
                      << " acceleration_output=" << currentAcceleration_.transpose( )
                      << std::endl;
        }
        currentTime_ = currentTime;
    }
}

}

}
