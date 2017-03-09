#include <Eigen/Geometry>

#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"

namespace tudat
{

namespace relativity
{

double calculateRelativisticAccelerationCorrectionsCommonterm(
        double centralBodyGravitationalParameter,
        double relativeDistance )
{
    return centralBodyGravitationalParameter /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT *
              relativeDistance * relativeDistance * relativeDistance );
}

Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        double centralBodyGravitationalParameter,
        Eigen::Vector3d relativePosition,
        Eigen::Vector3d relativeVelocity,
        double relativeDistance,
        double commonCorrectionTerm,
        double ppnParameterGamma,
        double ppnParameterBeta )
{
    Eigen::Vector3d acceleration = ( 2.0 * ( ppnParameterGamma + ppnParameterBeta ) *
                      centralBodyGravitationalParameter / relativeDistance - ppnParameterGamma *
                      relativeVelocity.dot( relativeVelocity ) ) * relativePosition +
            2.0 * ( 1.0 + ppnParameterGamma ) * ( relativePosition.dot( relativeVelocity ) ) * relativeVelocity;
    return commonCorrectionTerm * acceleration;
}

Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        double centralBodyGravitationalParameter,
        Eigen::Vector6d relativeState,
        double ppnParameterGamma,
        double ppnParameterBeta )
{
    return calculateScharzschildGravitationalAccelerationCorrection(
                centralBodyGravitationalParameter, relativeState.segment( 0, 3 ),
                relativeState.segment( 3, 3 ), relativeState.segment( 0, 3 ).norm( ),
                calculateRelativisticAccelerationCorrectionsCommonterm(
                    centralBodyGravitationalParameter, relativeState.segment( 0, 3 ).norm( ) ),
                ppnParameterGamma, ppnParameterBeta );
}

Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        Eigen::Vector3d relativePosition,
        Eigen::Vector3d relativeVelocity,
        double relativeDistance,
        double commonCorrectionTerm,
        Eigen::Vector3d centralBodyAngularMomentum,
        double ppnParameterGamma )
{
    Eigen::Vector3d acceleration = 3.0 / (
                relativeDistance * relativeDistance ) *
            relativePosition.cross( relativeVelocity ) *
            ( relativePosition.dot( centralBodyAngularMomentum ) ) +
            relativeVelocity.cross( centralBodyAngularMomentum );
    return acceleration * ( 1.0 + ppnParameterGamma ) * commonCorrectionTerm;
}

Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        double centralBodyGravitationalParameter,
        Eigen::Vector6d relativeState,
        Eigen::Vector3d centralBodyAngularMomentum,
        double ppnParameterGamma )
{
    return calculateLenseThirringCorrectionAcceleration(
                relativeState.segment( 0, 3 ), relativeState.segment( 3, 3 ), relativeState.segment( 0, 3 ).norm( ),
                calculateRelativisticAccelerationCorrectionsCommonterm(
                    centralBodyGravitationalParameter, relativeState.segment( 0, 3 ).norm( ) ),
                centralBodyAngularMomentum, ppnParameterGamma );

}

Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        Eigen::Vector3d orbiterRelativeVelocity,
        Eigen::Vector3d orbitedBodyPositionWrtLargerBody,
        Eigen::Vector3d orbitedBodyVelocityWrtLargerBody,
        double commonCorrectionTermOfLargerBody,
        double ppnParameterGamma )
{
    return - commonCorrectionTermOfLargerBody * ( 1.0 + 2.0 * ppnParameterGamma ) *
            ( orbitedBodyVelocityWrtLargerBody.cross( orbitedBodyPositionWrtLargerBody ) ).cross( orbiterRelativeVelocity );
}

Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        double largerBodyGravitationalParameter,
        Eigen::Vector6d orbiterRelativeState,
        Eigen::Vector6d orbitedBodyStateWrtLargerBody,
        double ppnParameterGamma )
{
    return calculateDeSitterCorrectionAcceleration(
                orbiterRelativeState.segment( 3, 3 ),
                orbitedBodyStateWrtLargerBody.segment( 0, 3 ),
                orbitedBodyStateWrtLargerBody.segment( 3, 3 ),
                calculateRelativisticAccelerationCorrectionsCommonterm(
                    largerBodyGravitationalParameter,
                    orbitedBodyStateWrtLargerBody.segment( 0, 3 ).norm( ) ),
                ppnParameterGamma );
}

void RelativisticAccelerationCorrection::updateMembers( const double currentTime )
{
    if( !( this->currentTime_ == currentTime ) )
    {
        this->currentTime_ = currentTime;

        stateOfAcceleratedBodyWrtCentralBody_ = stateFunctionOfAcceleratedBody_( ) - stateFunctionOfCentralBody_( );
        stateOfCentralBodyWrtPrimaryBody_ = stateFunctionOfCentralBody_( ) - stateFunctionOfPrimaryBody_( );
        gravitationalParameterOfCentralBody_ = gravitationalParameterFunctionOfCentralBody_( );
        gravitationalParameterOfPrimaryBody_ = gravitationalParameterFunctionOfPrimaryBody_( );

        ppnParameterGamma_ = ppnParameterGammaFunction_( );
        ppnParameterBeta_ = ppnParameterBetaFunction_( );

        commonCorrectionTerm_ = calculateRelativisticAccelerationCorrectionsCommonterm(
                    gravitationalParameterOfCentralBody_,
                    stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ).norm( ) );

        currentAcceleration_.setZero( );

        double relativeDistance = stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ).norm( );

        if( calculateSchwarzschildCorrection_ )
        {
            currentAcceleration_ = calculateScharzschildGravitationalAccelerationCorrection(
                        gravitationalParameterOfCentralBody_,
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ),
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                        relativeDistance, commonCorrectionTerm_, ppnParameterGamma_,
                        ppnParameterBeta_ );
        }


        if( calculateLenseThirringCorrection_ )
        {
            centalBodyAngularMomentum_ = centalBodyAngularMomentumFunction_( );
            currentAcceleration_ +=  calculateLenseThirringCorrectionAcceleration(
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ),
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                        relativeDistance, commonCorrectionTerm_, centalBodyAngularMomentum_,
                        ppnParameterGamma_ );
        }

        if( calculateDeSitterCorrection_ )
        {
            double primaryDistance = stateOfCentralBodyWrtPrimaryBody_.segment( 0, 3 ).norm( );

            stateOfCentralBodyWrtPrimaryBody_ = stateFunctionOfCentralBody_( ) - stateFunctionOfPrimaryBody_( );
            gravitationalParameterOfPrimaryBody_ = gravitationalParameterFunctionOfPrimaryBody_( );

            largerBodyCommonCorrectionTerm_ =  gravitationalParameterOfPrimaryBody_ / (
                    primaryDistance * primaryDistance * primaryDistance *
                        physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );

            currentAcceleration_ += calculateDeSitterCorrectionAcceleration(
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                        stateOfCentralBodyWrtPrimaryBody_.segment( 0, 3 ),
                        stateOfCentralBodyWrtPrimaryBody_.segment( 3, 3 ),
                        largerBodyCommonCorrectionTerm_,
                        ppnParameterGamma_ );
//            std::cout<<stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ).transpose( )<<" "<<
//                       stateOfCentralBodyWrtPrimaryBody_.transpose( )<<" "<<largerBodyCommonCorrectionTerm_<<" "<<
//                       ppnParameterGamma_<<" "<<gravitationalParameterOfPrimaryBody_<<" "<<primaryDistance<<" "<<(
//                           primaryDistance * primaryDistance * primaryDistance *
//                               physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT )<<std::endl;
//            std::cout<<currentAcceleration_.transpose( )<<std::endl<<std::endl;
        }


    }
}

}

}
