/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTICACCELERATIONCORRECTION_H
#define TUDAT_RELATIVISTICACCELERATIONCORRECTION_H

#include <functional>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace relativity
{

//! Function to compute a term common to several relativistic acceleration terms
/*!
 * Function to compute a term common to several relativistic acceleration terms, defined by mu / (c^2 * r^3), with mu the
 * gravitational parameter of the central body, c the speed of light and r the distance between the bodies
 * \param centralBodyGravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param relativeDistance Distance between bodies undergoing and exerting acceleration
 * \return Common term in relativistic accelerations.
 */
double calculateRelativisticAccelerationCorrectionsCommonterm(
        const double centralBodyGravitationalParameter,
        const double relativeDistance );

//! Function to compute the Schwarzschild term of the relativistic acceleration correction.
/*!
 *  Function to compute the Schwarzschild term of the relativistic acceleration correction.
 * \param centralBodyGravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param relativePosition Position of body undergoing, w.r.t. body exerting, acceleration.
 * \param relativeVelocity Velocity of body undergoing, w.r.t. body exerting, acceleration.
 * \param relativeDistance Distance between bodies undergoing and exerting acceleration (norm of relativePosition)
 * \param commonCorrectionTerm Common term in relativistic accelerations, as computed by
 * calculateRelativisticAccelerationCorrectionsCommonterm function
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter beta
 * \return Schwarzschild term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& relativeVelocity,
        const double relativeDistance,
        const double commonCorrectionTerm,
        const double ppnParameterGamma = 1.0,
        const double ppnParameterBeta = 1.0 );

//! Function to compute the Schwarzschild term of the relativistic acceleration correction.
/*!
 *  Function to compute the Schwarzschild term of the relativistic acceleration correction.
 * \param centralBodyGravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param ppnParameterGamma PPN parameter gamma
 * \param ppnParameterBeta PPN parameter beta
 * \return Schwarzschild term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateScharzschildGravitationalAccelerationCorrection(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector6d& relativeState,
        const double ppnParameterGamma = 1.0,
        const double ppnParameterBeta = 1.0 );

//! Function to compute the Lense-Thirring term of the relativistic acceleration correction.
/*!
 *  Function to compute the Lense-Thirring term of the relativistic acceleration correction.
 * \param relativePosition Position of body undergoing, w.r.t. body exerting, acceleration.
 * \param relativeVelocity Velocity of body undergoing, w.r.t. body exerting, acceleration.
 * \param relativeDistance Distance between bodies undergoing and exerting acceleration (norm of relativePosition)
 * \param commonCorrectionTerm Common term in relativistic accelerations, as computed by
 * calculateRelativisticAccelerationCorrectionsCommonterm function
 * \param centralBodyAngularMomentum Angular momentum vector of central body.
 * \param ppnParameterGamma PPN parameter gamma
 * \return Lense-Thirring term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        const Eigen::Vector3d& relativePosition,
        const Eigen::Vector3d& relativeVelocity,
        const double relativeDistance,
        const double commonCorrectionTerm,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        const double ppnParameterGamma = 1.0 );

//! Function to compute the Lense-Thirring term of the relativistic acceleration correction.
/*!
 *  Function to compute the Lense-Thirring term of the relativistic acceleration correction.
 * \param centralBodyGravitationalParameter Gravitational parameter of body exerting acceleration.
 * \param relativeState Cartesian state of body undergoing, w.r.t. body exerting, acceleration.
 * \param centralBodyAngularMomentum Angular momentum vector of central body.
 * \param ppnParameterGamma PPN parameter gamma
 * \return Schwarzschild term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateLenseThirringCorrectionAcceleration(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector6d& relativeState,
        const Eigen::Vector3d& centralBodyAngularMomentum,
        const double ppnParameterGamma = 1.0 );

//! Function to compute the de Sitter term of the relativistic acceleration correction.
/*!
 *  Function to compute the  de Sitter term of the relativistic acceleration correction.
 * \param orbiterRelativeVelocity Velocity of body undergoing, w.r.t. body exerting, acceleration.
 * \param orbitedBodyPositionWrtLargerBody Position of body undergoing acceleration w.r.t. its central body. For an
 * acceleration on a satellite orbiting the Earth, this would be the position of the Earth w.r.t. the Sun.
 * \param orbitedBodyVelocityWrtLargerBody Velocity of body undergoing acceleration w.r.t. its central body. For an
 * acceleration on a satellite orbiting the Earth, this would be the velocity of the Earth w.r.t. the Sun.
 * \param commonCorrectionTermOfLargerBody Common term in relativistic accelerations, with properties of e.g. Earth w.r.t
 * the Sun for an Earth-orbitign satellite, as computed by calculateRelativisticAccelerationCorrectionsCommonterm function
 * \param ppnParameterGamma PPN parameter gamma
 * \return De Sitter term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        const Eigen::Vector3d& orbiterRelativeVelocity,
        const Eigen::Vector3d& orbitedBodyPositionWrtLargerBody,
        const Eigen::Vector3d& orbitedBodyVelocityWrtLargerBody,
        const double commonCorrectionTermOfLargerBody,
        const double ppnParameterGamma = 1.0 );

//! Function to compute the de Sitter term of the relativistic acceleration correction.
/*!
 *  Function to compute the  de Sitter term of the relativistic acceleration correction.
 * \param largerBodyGravitationalParameter Gravitational parameter of body primarilly responsible for the motion of the
 * body exerting the acceleration.
 * \param orbiterRelativeState State of body undergoing, w.r.t. body exerting, acceleration.
 * \param orbitedBodyStateWrtLargerBody Cartesian state of body undergoing acceleration w.r.t. its central body. For an
 * acceleration on a satellite orbiting the Earth, this would be the Cartesian state of the Earth w.r.t. the Sun.
 * \param ppnParameterGamma PPN parameter gamma
 * \return De Sitter term of the relativistic acceleration correction.
 */
Eigen::Vector3d calculateDeSitterCorrectionAcceleration(
        const double largerBodyGravitationalParameter,
        const Eigen::Vector6d& orbiterRelativeState,
        const Eigen::Vector6d& orbitedBodyStateWrtLargerBody,
        const double ppnParameterGamma = 1.0 );

//! Class to compute typical relativistic corrections to the dynamics of an orbiter.
/*!
 *  Class to compute typical relativistic corrections to the dynamics of an orbiter. The model allows the inclusion of the
 *  Schwarzschild, Lense-Thirring and de Sitter terms (any one of them may be turned on and off). An excellent introduction to
 *  these models is given in 'General Relativity and Space Geodesy' by L. Combrinck (2012).
 */
class RelativisticAccelerationCorrection: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor, used when including deSitter acceleration
    /*!
     * Constructor, used when including deSitter acceleration
     * \param stateFunctionOfAcceleratedBody State function of vehicle undergoing acceleration
     * \param stateFunctionOfCentralBody State function of main body exerting acceleration (e.g. Earth for an Earth-orbiting
     * satellite).
     * \param stateFunctionOfPrimaryBody State function of large body primarily responsible for motion of central body
     * (e.g. Sun for acceleration acting on an Earth-orbiting satellite).
     * \param gravitationalParameterFunctionOfCentralBody Function returning the gravitational parameter of the central body
     * \param gravitationalParameterFunctionOfPrimaryBody Function returning the gravitational parameter of the primary body
     * \param primaryBodyName Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite)
     * \param centalBodyAngularMomentumFunction Function returning the angular momenum of the central body, expressed in the
     * propagation frame (default empty; no Lense-Thirring acceleration if empty).
     * \param ppnParameterGammaFunction Function returning the PPN parameter gamma (default 1)
     * \param ppnParameterBetaFunction Function returning the PPN parameter beta (default 1)
     * \param calculateSchwarzschildCorrection Boolean denoting whether the schwarzschild term is to be used.
     */
    RelativisticAccelerationCorrection(
            std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            std::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody,
            std::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            std::function< double( ) > gravitationalParameterFunctionOfPrimaryBody,
            std::string primaryBodyName,
            std::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction = std::function< Eigen::Vector3d( ) >( ),
            std::function< double( ) > ppnParameterGammaFunction = [ ]( ){ return 1.0; },
            std::function< double( ) > ppnParameterBetaFunction = [ ]( ){ return 1.0; },
            const bool calculateSchwarzschildCorrection = true ):
        AccelerationModel< Eigen::Vector3d >( ),
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
        calculateLenseThirringCorrection_( !( centalBodyAngularMomentumFunction == nullptr ) )
    { }

    //! Constructor, used when including Lense-Thirring, but not de Sitter, acceleration
    /*!
     * Constructor, used when including Lense-Thirring, but not de Sitter, acceleration
     * \param stateFunctionOfAcceleratedBody State function of vehicle undergoing acceleration
     * \param stateFunctionOfCentralBody State function of main body exerting acceleration (e.g. Earth for an Earth-orbiting
     * satellite).
     * \param gravitationalParameterFunctionOfCentralBody Function returning the gravitational parameter of the central body
     * \param centalBodyAngularMomentumFunction Function returning the angular momenum of the central body, expressed in the
     * propagation frame (default empty; no Lense-Thirring acceleration if empty).
     * \param ppnParameterGammaFunction Function returning the PPN parameter gamma (default 1)
     * \param ppnParameterBetaFunction Function returning the PPN parameter beta (default 1)
     * \param calculateSchwarzschildCorrection Boolean denoting whether the schwarzschild term is to be used.
     */
    RelativisticAccelerationCorrection(
            std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            std::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            std::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction,
            std::function< double( ) > ppnParameterGammaFunction = [ ]( ){ return 1.0; },
            std::function< double( ) > ppnParameterBetaFunction = [ ]( ){ return 1.0; },
            const bool calculateSchwarzschildCorrection = true ):
        AccelerationModel< Eigen::Vector3d >( ),
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

    //! Constructor, used for Schwarzschild term only
    /*!
     * Constructor, used for Schwarzschild term only
     * \param stateFunctionOfAcceleratedBody State function of vehicle undergoing acceleration
     * \param stateFunctionOfCentralBody State function of main body exerting acceleration (e.g. Earth for an Earth-orbiting
     * satellite).
     * \param gravitationalParameterFunctionOfCentralBody Function returning the gravitational parameter of the central body
     * \param ppnParameterGammaFunction Function returning the PPN parameter gamma (default 1)
     * \param ppnParameterBetaFunction Function returning the PPN parameter beta (default 1)
     */
    RelativisticAccelerationCorrection(
            std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            std::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            std::function< double( ) > ppnParameterGammaFunction = [ ]( ){ return 1.0; },
            std::function< double( ) > ppnParameterBetaFunction = [ ]( ){ return 1.0; } ):
        AccelerationModel< Eigen::Vector3d >( ),
        stateFunctionOfAcceleratedBody_( stateFunctionOfAcceleratedBody ),
        stateFunctionOfCentralBody_( stateFunctionOfCentralBody ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        ppnParameterGammaFunction_( ppnParameterGammaFunction ),
        ppnParameterBetaFunction_( ppnParameterBetaFunction ),
        calculateSchwarzschildCorrection_( true ),
        calculateDeSitterCorrection_( false ),
        calculateLenseThirringCorrection_( false )
    { }

    //! Destructor
    ~RelativisticAccelerationCorrection( ){ }

    //! Function to return the current acceleration
    /*!
     * Returns the relativistic correction acceleration. Value is computed by updateMembers function
     * \return Acceleration.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return  currentAcceleration_;
    }

    //! Update member variables used by the relativistic correction acceleration model.
    /*!
     * Updates member variables used by the relativistic correction acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->currentTime_ = currentTime;

            // Update common variables
            stateOfAcceleratedBodyWrtCentralBody_ = stateFunctionOfAcceleratedBody_( ) - stateFunctionOfCentralBody_( );
            gravitationalParameterOfCentralBody_ = gravitationalParameterFunctionOfCentralBody_( );

            ppnParameterGamma_ = ppnParameterGammaFunction_( );
            ppnParameterBeta_ = ppnParameterBetaFunction_( );

            commonCorrectionTerm_ = calculateRelativisticAccelerationCorrectionsCommonterm(
                        gravitationalParameterOfCentralBody_,
                        stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ).norm( ) );

            currentAcceleration_.setZero( );

            double relativeDistance = stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ).norm( );

            // Compute Schwarzschild term (if requested)
            if( calculateSchwarzschildCorrection_ )
            {
                currentAcceleration_ = calculateScharzschildGravitationalAccelerationCorrection(
                            gravitationalParameterOfCentralBody_,
                            stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ),
                            stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                            relativeDistance, commonCorrectionTerm_, ppnParameterGamma_,
                            ppnParameterBeta_ );
            }

            // Compute Lense-Thirring term (if requested)
            if( calculateLenseThirringCorrection_ )
            {
                centalBodyAngularMomentum_ = centalBodyAngularMomentumFunction_( );
                currentAcceleration_ +=  calculateLenseThirringCorrectionAcceleration(
                            stateOfAcceleratedBodyWrtCentralBody_.segment( 0, 3 ),
                            stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                            relativeDistance, commonCorrectionTerm_, centalBodyAngularMomentum_,
                            ppnParameterGamma_ );

            }

            // Compute de Sitter term (if requested)
            if( calculateDeSitterCorrection_ )
            {
                stateOfCentralBodyWrtPrimaryBody_ = stateFunctionOfCentralBody_( ) - stateFunctionOfPrimaryBody_( );
                gravitationalParameterOfPrimaryBody_ = gravitationalParameterFunctionOfPrimaryBody_( );

                double primaryDistance = stateOfCentralBodyWrtPrimaryBody_.segment( 0, 3 ).norm( );

                stateOfCentralBodyWrtPrimaryBody_ = stateFunctionOfCentralBody_( ) - stateFunctionOfPrimaryBody_( );
                gravitationalParameterOfPrimaryBody_ = gravitationalParameterFunctionOfPrimaryBody_( );

                double largerBodyCommonCorrectionTerm =  gravitationalParameterOfPrimaryBody_ / (
                        primaryDistance * primaryDistance * primaryDistance *
                            physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );

                currentAcceleration_ += calculateDeSitterCorrectionAcceleration(
                            stateOfAcceleratedBodyWrtCentralBody_.segment( 3, 3 ),
                            stateOfCentralBodyWrtPrimaryBody_.segment( 0, 3 ),
                            stateOfCentralBodyWrtPrimaryBody_.segment( 3, 3 ),
                            largerBodyCommonCorrectionTerm,
                            ppnParameterGamma_ );

            }


        }
    }


    //! Function to return the current state of the body undergoing acceleration
    /*!
     * Function to return the current state of the body undergoing acceleration
     * \return Current state of the body undergoing acceleration
     */
    std::function< Eigen::Vector6d( ) > getStateFunctionOfAcceleratedBody( )
    { return stateFunctionOfAcceleratedBody_; }

    //! Function to return the current state of the main body exerting acceleration
    /*!
     * Function to return the current state of the main body exerting acceleration
     * \return Current state of the main body exerting acceleration
     */
    std::function< Eigen::Vector6d( ) > getStateFunctionOfCentralBody( )
    { return stateFunctionOfCentralBody_; }

    //! Function to return the current gravitational parameter of central body
    /*!
     * Function to return the current gravitational parameter of central body
     * \return Current gravitational parameter of central body
     */
    std::function< double( ) > getGravitationalParameterFunctionOfCentralBody( )
    { return gravitationalParameterFunctionOfCentralBody_; }

    //! Function to return the current PPN parameter gamma
    /*!
     * Function to return the current PPN parameter gamma
     * \return Current PPN parameter gamma
     */
    std::function< double( ) > getPpnParameterGammaFunction_( )
    { return ppnParameterGammaFunction_; }

    //! Function to return the current PPN parameter beta
    /*!
     * Function to return the current PPN parameter beta
     * \return Current PPN parameter beta
     */
    std::function< double( ) > getPpnParameterBetaFunction_( )
    { return ppnParameterBetaFunction_; }

    //! Function to return the boolean denoting wheter the Schwarzschild term is used.
    /*!
     * Function to return the boolean denoting wheter the Schwarzschild term is used.
     * \return Boolean denoting wheter the Schwarzschild term is used.
     */
    bool getCalculateSchwarzschildCorrection( )
    { return calculateSchwarzschildCorrection_; }

    //! Function to return the boolean denoting wheter the de Sitter term is used.
    /*!
     * Function to return the boolean denoting wheter the de Sitter term is used.
     * \return Boolean denoting wheter the de Sitter term is used.
     */
    bool getCalculateDeSitterCorrection( )
    { return calculateDeSitterCorrection_; }

    //! Function to return the boolean denoting wheter the Lense-Thirring term is used.
    /*!
     * Function to return the boolean denoting wheter the Lense-Thirring term is used.
     * \return Boolean denoting wheter the Lense-Thirring term is used.
     */
    bool getCalculateLenseThirringCorrection( )
    { return calculateLenseThirringCorrection_; }

    std::string getPrimaryBodyName( )
    { return primaryBodyName_; }

private:


    //! State function of vehicle undergoing acceleration
    std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody_;

    //! State function of main body exerting acceleration (e.g. Earth for an Earth-orbiting satellite).
    std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody_;

    //! State function of large body primarily responsible for motion of central body (e.g. Sun for acceleration acting on
    //! an Earth-orbiting satellite).
    std::function< Eigen::Vector6d( ) > stateFunctionOfPrimaryBody_;

    //! Function returning the gravitational parameter of the central body
    std::function< double( ) > gravitationalParameterFunctionOfCentralBody_;

    //! Function returning the gravitational parameter of the primary body
    std::function< double( ) > gravitationalParameterFunctionOfPrimaryBody_;

    //! Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite)
    std::string primaryBodyName_;

    //! Function returning the angular momenum of the central body (expressed in the propagation frame)
    std::function< Eigen::Vector3d( ) > centalBodyAngularMomentumFunction_;

    //! Function returning the PPN parameter gamma
    std::function< double( ) > ppnParameterGammaFunction_;

    //! Function returning the PPN parameter beta
    std::function< double( ) > ppnParameterBetaFunction_;



    //! Current state of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfAcceleratedBody_;

    //! Current state of the main body exerting acceleration, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfCentralBodyWrtPrimaryBody_;

    //! Current state of the primary body, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfAcceleratedBodyWrtCentralBody_;

    //! Current gravitational parameter of central body
    double gravitationalParameterOfCentralBody_;

    //! Current gravitational parameter of primary body
    double gravitationalParameterOfPrimaryBody_;

    //! Current angulat momentum vector of central body
    Eigen::Vector3d centalBodyAngularMomentum_;

    //! Current PPN parameter gamma
    double ppnParameterGamma_;

    //! Current PPN parameter beta
    double ppnParameterBeta_;

    //! Pre-computed common term for corrections (computed by calculateRelativisticAccelerationCorrectionsCommonterm)
    double commonCorrectionTerm_;



    //! Boolean denoting wheter the Schwarzschild term is used.
    bool calculateSchwarzschildCorrection_;

    //! Boolean denoting wheter the de Sitter term is used.
    bool calculateDeSitterCorrection_;

    //! Boolean denoting wheter the Lense-Thirring term is used.
    bool calculateLenseThirringCorrection_;



    //! Relativistic acceleration correction, as computed by last call to updateMembers function
    Eigen::Vector3d currentAcceleration_;

};

}

}

#endif // TUDAT_RELATIVISTICACCELERATIONCORRECTION_H
