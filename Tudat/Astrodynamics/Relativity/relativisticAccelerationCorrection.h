#ifndef RELATIVISTICACCELERATIONCORRECTION_H
#define RELATIVISTICACCELERATIONCORRECTION_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace relativity
{

double calculateRelativisticAccelerationCorrectionsCommonterm(
        double centralBodyGravitationalParameter,
        double relativeDistance );

Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        double centralBodyGravitationalParameter,
        Eigen::Vector3d relativePosition,
        Eigen::Vector3d relativeVelocity,
        double relativeDistance,
        double commonCorrectionTerm,
        double ppnParameterGamma = 1.0,
        double ppnParameterBeta = 1.0 );

Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        double centralBodyGravitationalParameter,
        Eigen::Vector6d relativeState,
        double ppnParameterGamma = 1.0,
        double ppnParameterBeta = 1.0 );

Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        Eigen::Vector3d relativePosition,
        Eigen::Vector3d relativeVelocity,
        double relativeDistance,
        double commonCorrectionTerm,
        Eigen::Vector3d centralBodyAngularMomentum,
        double ppnParameterGamma = 1.0 );

Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        double centralBodyGravitationalParameter,
        Eigen::Vector6d relativeState,
        Eigen::Vector3d centralBodyAngularMomentum,
        double ppnParameterGamma = 1.0 );

Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        Eigen::Vector3d orbiterRelativeVelocity,
        Eigen::Vector3d orbitedBodyPositionWrtLargerBody,
        Eigen::Vector3d orbitedBodyVelocityWrtLargerBody,
        double commonCorrectionTermOfLargerBody,
        double ppnParameterGamma = 1.0 );

Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        double largerBodyGravitationalParameter,
        Eigen::Vector6d orbiterRelativeState,
        Eigen::Vector6d orbitedBodyStateWrtLargerBody,
        double ppnParameterGamma = 1.0 );

class RelativisticAccelerationCorrection: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    RelativisticAccelerationCorrection(
            boost::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            boost::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            boost::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody,
            boost::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            boost::function< double( ) > gravitationalParameterFunctionOfPrimaryBody,
            std::string primaryBodyName,
            boost::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction = boost::function< Eigen::Vector3d( ) >( ),
            boost::function< double( ) > ppnParameterGammaFunction = boost::lambda::constant( 1.0 ),
            boost::function< double( ) > ppnParameterBetaFunction = boost::lambda::constant( 1.0 ),
            const bool calculateSchwarzschildCorrection = true ):
        stateFunctionOfAcceleratedBody_( stateFunctionOfAcceleratedBody ),
        stateFunctionOfCentralBody_( stateFunctionOfCentralBody ),
        stateFunctionOfPrimaryBody_( stateFunctionOfPrimaryBody ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        gravitationalParameterFunctionOfPrimaryBody_( gravitationalParameterFunctionOfPrimaryBody ),
        primaryBodyName_( primaryBodyName ),
        centalBodyAngularMomentumFunction_( centalBodyAngularMomentumFunction ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction ),
        ppnParameterBetaFunction_( ppnParameterBetaFunction ),
        calculateSchwarzschildCorrection_( calculateSchwarzschildCorrection ),
        calculateDeSitterCorrection_( true ),
        calculateLenseThirringCorrection_( !centalBodyAngularMomentumFunction.empty( ) )
    { }

    RelativisticAccelerationCorrection(
            boost::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            boost::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            boost::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            boost::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction,
            boost::function< double( ) > ppnParameterGammaFunction = boost::lambda::constant( 1.0 ),
            boost::function< double( ) > ppnParameterBetaFunction = boost::lambda::constant( 1.0 ),
            const bool calculateSchwarzschildCorrection = true ):
        stateFunctionOfAcceleratedBody_( stateFunctionOfAcceleratedBody ),
        stateFunctionOfCentralBody_( stateFunctionOfCentralBody ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        centalBodyAngularMomentumFunction_( centalBodyAngularMomentumFunction ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction ),
        ppnParameterBetaFunction_( ppnParameterBetaFunction ),
        calculateSchwarzschildCorrection_( calculateSchwarzschildCorrection ),
        calculateDeSitterCorrection_( false ),
        calculateLenseThirringCorrection_( true )
    { }


    RelativisticAccelerationCorrection(
            boost::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            boost::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            boost::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            boost::function< double( ) > ppnParameterGammaFunction = boost::lambda::constant( 1.0 ),
            boost::function< double( ) > ppnParameterBetaFunction = boost::lambda::constant( 1.0 ) ):
        stateFunctionOfAcceleratedBody_( stateFunctionOfAcceleratedBody ),
        stateFunctionOfCentralBody_( stateFunctionOfCentralBody ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction ),
        ppnParameterBetaFunction_( ppnParameterBetaFunction ),
        calculateSchwarzschildCorrection_( true ),
        calculateDeSitterCorrection_( false ),
        calculateLenseThirringCorrection_( false )
    { }

    ~RelativisticAccelerationCorrection( ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return  currentAcceleration_;
    }

    void updateMembers( const double currentTime = TUDAT_NAN );

    boost::function< Eigen::Vector6d( ) > getStateFunctionOfAcceleratedBody( )
    { return stateFunctionOfAcceleratedBody_; }

    boost::function< Eigen::Vector6d( ) > getStateFunctionOfCentralBody( )
    { return stateFunctionOfCentralBody_; }

    boost::function< double( ) > getGravitationalParameterFunctionOfCentralBody( )
    { return gravitationalParameterFunctionOfCentralBody_; }

    boost::function< double( ) > getPpnParameterGammaFunction_( )
    { return ppnParameterGammaFunction_; }

    boost::function< double( ) > getPpnParameterBetaFunction_( )
    { return ppnParameterBetaFunction_; }

    bool getCalculateSchwarzschildCorrection( )
    { return calculateSchwarzschildCorrection_; }

    bool getCalculateDeSitterCorrection( )
    { return calculateDeSitterCorrection_; }

    bool getCalculateLenseThirringCorrection( )
    { return calculateLenseThirringCorrection_; }

    std::string getPrimaryBodyName( )
    { return primaryBodyName_; }

private:

    boost::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody_;
    boost::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody_;
    boost::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody_;

    boost::function< double( ) > gravitationalParameterFunctionOfCentralBody_;
    boost::function< double( ) > gravitationalParameterFunctionOfPrimaryBody_;

    std::string primaryBodyName_;

    boost::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction_;

    boost::function< double( ) > ppnParameterGammaFunction_;
    boost::function< double( ) > ppnParameterBetaFunction_;


    Eigen::Vector6d stateOfAcceleratedBody_;
    Eigen::Vector6d stateOfCentralBodyWrtPrimaryBody_;
    Eigen::Vector6d stateOfAcceleratedBodyWrtCentralBody_;

    double gravitationalParameterOfCentralBody_;
    double gravitationalParameterOfPrimaryBody_;

    Eigen::Vector3d centalBodyAngularMomentum_;

    double ppnParameterGamma_;
    double ppnParameterBeta_;

    double commonCorrectionTerm_;
    double largerBodyCommonCorrectionTerm_;

    bool calculateSchwarzschildCorrection_;
    bool calculateDeSitterCorrection_;
    bool calculateLenseThirringCorrection_;

    Eigen::Vector3d currentAcceleration_;

};

}

}

#endif // RELATIVISTICACCELERATIONCORRECTION_H
