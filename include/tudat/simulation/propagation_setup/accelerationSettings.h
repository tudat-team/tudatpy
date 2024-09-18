/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ACCELERATIONSETTINGS_H
#define TUDAT_ACCELERATIONSETTINGS_H

#include <functional>
#include <memory>
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/gravitation/thirdBodyPerturbation.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/deprecationWarnings.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"

// #include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

namespace simulation_setup
{


// Class for providing settings for acceleration model.
/*
 *  Class for providing settings for acceleration model. This class is a functional (base) class for
 *  settings of acceleration models that  require no information in addition to their type.
 *  Classes defining settings for acceleration models requiring additional information must be
 *  derived from this class.
 *  Bodies exerting and undergong acceleration are set externally from this class.
 *  This class can be used for the easy setup of acceleration models
 *  (see createAccelerationModels.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */
//! @get_docstring(AccelerationSettings.__docstring__)
class AccelerationSettings
{
public:

    // Constructor, sets type of acceleration.
    /*
     *  Constructor, sets type of acceleration.
     *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
     */
    AccelerationSettings( const basic_astrodynamics::AvailableAcceleration accelerationType ):
        accelerationType_( accelerationType ){ }

    // Destructor.
    virtual ~AccelerationSettings( ){ }

    // Type of acceleration from AvailableAcceleration enum.
    basic_astrodynamics::AvailableAcceleration accelerationType_;

};

class RadiationPressureAccelerationSettings: public AccelerationSettings
{
public:

    // Constructor, sets type of acceleration.
    /*
     *  Constructor, sets type of acceleration.
     *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
     */
    RadiationPressureAccelerationSettings( const RadiationPressureTargetModelType targetModelType = undefined_target ):
        AccelerationSettings( basic_astrodynamics::radiation_pressure ), targetModelType_( targetModelType )
        {
            if( targetModelType_ == multi_type_target )
            {
                throw std::runtime_error( "Error when creating radiation pressure acceleration settings, cannot select multi-type target" );
            }
        }

    // Destructor.
    virtual ~RadiationPressureAccelerationSettings( ){ }

    RadiationPressureTargetModelType targetModelType_;

};

inline std::shared_ptr< AccelerationSettings > acceleration( basic_astrodynamics::AvailableAcceleration accelerationType  )
{
    return std::make_shared< AccelerationSettings >( accelerationType );
}

//! @get_docstring(pointMassGravityAcceleration)
inline std::shared_ptr< AccelerationSettings > pointMassGravityAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity );
}

inline std::shared_ptr< AccelerationSettings > einsteinInfledHoffmannGravityAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::einstein_infeld_hoffmann_acceleration );
}


//! @get_docstring(aerodynamicAcceleration)
inline std::shared_ptr< AccelerationSettings > aerodynamicAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic );
}

//! @get_docstring(cannonBallRadiationPressureAcceleration)
inline std::shared_ptr< AccelerationSettings > cannonBallRadiationPressureAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure );
}

inline std::shared_ptr< AccelerationSettings > radiationPressureAcceleration( const RadiationPressureTargetModelType targetModelType = undefined_target )
{
    return std::make_shared< RadiationPressureAccelerationSettings >( targetModelType );
}


// Class for providing settings for spherical harmonics acceleration model.
/*
 *  Class for providing settings for spherical harmonics acceleration model,
 *  specifically the maximum degree and order up to which the field is to be expanded. Note that
 *  the minimum degree and order are currently always set to zero.
 */
//! @get_docstring(SphericalHarmonicAccelerationSettings.__docstring__)
class SphericalHarmonicAccelerationSettings: public AccelerationSettings
{
public:
    // Constructor to set maximum degree and order that is to be taken into account.
    /*
     *  Constructor to set maximum degree and order that is to be taken into account.
     *  \param maximumDegree Maximum degree
     *  \param maximumOrder Maximum order
     */
    SphericalHarmonicAccelerationSettings( const int maximumDegree,
                                           const int maximumOrder,
                                           const bool removePointMass = false ):
        AccelerationSettings( basic_astrodynamics::spherical_harmonic_gravity ),
        maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder ),
        removePointMass_( removePointMass ){ }


    // Maximum degree that is to be used for spherical harmonic acceleration
    int maximumDegree_;

    // Maximum order that is to be used for spherical harmonic acceleration
    int maximumOrder_;

    bool removePointMass_;
};

//! @get_docstring(sphericalHarmonicAcceleration)
inline std::shared_ptr< AccelerationSettings > sphericalHarmonicAcceleration(
        const int maximumDegree, const int maximumOrder )
{
    return std::make_shared< SphericalHarmonicAccelerationSettings >( maximumDegree, maximumOrder );
}

// Class for providing acceleration settings for mutual spherical harmonics acceleration model.
/*
 *  Class for providing acceleration settings for mutual spherical harmonics acceleration model,
 *  specifically the maximum degree and order up to which the fields of the bodies are be expanded.
 *  Please note that the minimum degrees and orders are currently always set to zero.
 */
//! @get_docstring(MutualSphericalHarmonicAccelerationSettings.__docstring__)
class MutualSphericalHarmonicAccelerationSettings: public AccelerationSettings
{
public:

    // Constructor to set maximum degrees and orders that are to be taken into account.
    /*
     * Constructor to set maximum degrees and orders that are to be taken into account.
     * \param maximumDegreeOfBodyExertingAcceleration Maximum degree of body exerting acceleration.
     * \param maximumOrderOfBodyExertingAcceleration Maximum order of body exerting acceleration.
     * \param maximumDegreeOfBodyUndergoingAcceleration Maximum degree of body undergoing acceleration.
     * \param maximumOrderOfBodyUndergoingAcceleration Maximum order of body undergoing acceleration.
     * \param maximumDegreeOfCentralBody Maximum degree of central body (only relevant for 3rd body acceleration).
     * \param maximumOrderOfCentralBody Maximum order of central body (only relevant for 3rd body acceleration).
     */
    MutualSphericalHarmonicAccelerationSettings( const int maximumDegreeOfBodyExertingAcceleration,
                                                 const int maximumOrderOfBodyExertingAcceleration,
                                                 const int maximumDegreeOfBodyUndergoingAcceleration,
                                                 const int maximumOrderOfBodyUndergoingAcceleration,
                                                 const int maximumDegreeOfCentralBody = 0,
                                                 const int maximumOrderOfCentralBody = 0 ):
        AccelerationSettings( basic_astrodynamics::mutual_spherical_harmonic_gravity ),
        maximumDegreeOfBodyExertingAcceleration_( maximumDegreeOfBodyExertingAcceleration ),
        maximumOrderOfBodyExertingAcceleration_( maximumOrderOfBodyExertingAcceleration ),
        maximumDegreeOfBodyUndergoingAcceleration_( maximumDegreeOfBodyUndergoingAcceleration ),
        maximumOrderOfBodyUndergoingAcceleration_( maximumOrderOfBodyUndergoingAcceleration ),
        maximumDegreeOfCentralBody_( maximumDegreeOfCentralBody ), maximumOrderOfCentralBody_( maximumOrderOfCentralBody ){ }

    // Maximum degree of body exerting acceleration.
    int maximumDegreeOfBodyExertingAcceleration_;

    // Maximum order of body exerting acceleration.
    int maximumOrderOfBodyExertingAcceleration_;

    // Maximum degree of body undergoing acceleration.
    int maximumDegreeOfBodyUndergoingAcceleration_;

    // Maximum order of body undergoing acceleration.
    int maximumOrderOfBodyUndergoingAcceleration_;

    // Maximum degree of central body (only releveant for 3rd body acceleration).
    int maximumDegreeOfCentralBody_;

    // Maximum order of central body (only releveant for 3rd body acceleration).
    int maximumOrderOfCentralBody_;

};

//! @get_docstring(mutualSphericalHarmonicAcceleration)
inline std::shared_ptr< AccelerationSettings > mutualSphericalHarmonicAcceleration(
		const int maximumDegreeOfBodyExertingAcceleration,
		const int maximumOrderOfBodyExertingAcceleration,
		const int maximumDegreeOfBodyUndergoingAcceleration,
		const int maximumOrderOfBodyUndergoingAcceleration,
		const int maximumDegreeOfCentralBody = 0,
		const int maximumOrderOfCentralBody = 0
		)
{
	return std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
			maximumDegreeOfBodyExertingAcceleration,
			maximumOrderOfBodyExertingAcceleration,
			maximumDegreeOfBodyUndergoingAcceleration,
			maximumOrderOfBodyUndergoingAcceleration,
			maximumDegreeOfCentralBody,
			maximumOrderOfCentralBody
			);
}

inline std::shared_ptr< AccelerationSettings > polyhedronAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::polyhedron_gravity );
}

inline std::shared_ptr< AccelerationSettings > ringAcceleration( )
{
    return std::make_shared< AccelerationSettings >( basic_astrodynamics::ring_gravity );
}

// Class to provide settings for typical relativistic corrections to the dynamics of an orbiter.
/*
 *  Class to provide settings for typical relativistic corrections to the dynamics of an orbiter: the
 *  Schwarzschild, Lense-Thirring and de Sitter terms. An excellent introduction to
 *  these models is given in 'General relativity and Space Geodesy' by L. Combrinck (2012).
 */
//! @get_docstring(RelativisticAccelerationCorrectionSettings.__docstring__)
class RelativisticAccelerationCorrectionSettings: public AccelerationSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param calculateSchwarzschildCorrection Boolean denoting whether the Schwarzschild term is used.
     * \param calculateLenseThirringCorrection Boolean denoting whether the Lense-Thirring term is used.
     * \param calculateDeSitterCorrection Boolean denoting whether the de Sitter term is used.
     * \param primaryBody Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite)
     * \param centralBodyAngularMomentum Constant angular momentum of central body. NOTE: Passing angular momentum through this
     * function is temporary: in the future this will be done consistently with rotation/gravity field.
     */
    RelativisticAccelerationCorrectionSettings(
            const bool calculateSchwarzschildCorrection = true,
            const bool calculateLenseThirringCorrection = false,
            const bool calculateDeSitterCorrection = false,
            const std::string primaryBody = "",
            const Eigen::Vector3d centralBodyAngularMomentum = Eigen::Vector3d::Zero( ) ):
        AccelerationSettings(  basic_astrodynamics::relativistic_correction_acceleration ),
        calculateSchwarzschildCorrection_( calculateSchwarzschildCorrection ),
        calculateLenseThirringCorrection_( calculateLenseThirringCorrection ),
        calculateDeSitterCorrection_( calculateDeSitterCorrection ),
        primaryBody_( primaryBody ),
        centralBodyAngularMomentum_( centralBodyAngularMomentum )
    {
        if( calculateDeSitterCorrection_ && primaryBody_ == "" )
        {
            throw std::runtime_error(
                        "Error when making relativistic acceleration correction, deSitter acceleration requested without primary body" );
        }
    }

    // Boolean denoting wheter the Schwarzschild term is used.
    bool calculateSchwarzschildCorrection_;

    // Boolean denoting wheter the Lense-Thirring term is used.
    bool calculateLenseThirringCorrection_;

    // Boolean denoting wheter the de Sitter term is used.
    bool calculateDeSitterCorrection_;

    // Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite)
    std::string primaryBody_;

    // Constant angular momentum of central body
    Eigen::Vector3d centralBodyAngularMomentum_;

};

//! @get_docstring(relativisticAccelerationCorrection)
inline std::shared_ptr< AccelerationSettings > relativisticAccelerationCorrection(
		const bool calculateSchwarzschildCorrection = true,
		const bool calculateLenseThirringCorrection = false,
		const bool calculateDeSitterCorrection = false,
		const std::string primaryBody = "",
        const Eigen::Vector3d centralBodyAngularMomentum = Eigen::Vector3d::Zero( ) )
{
	return std::make_shared< RelativisticAccelerationCorrectionSettings >(
			calculateSchwarzschildCorrection,
			calculateLenseThirringCorrection,
			calculateDeSitterCorrection,
			primaryBody,
			centralBodyAngularMomentum
			);
}

// Class to define settings for empirical accelerations
//! @get_docstring(EmpiricalAccelerationSettings.__docstring__)
class EmpiricalAccelerationSettings: public AccelerationSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param constantAcceleration Acceleration (in RSW frame) that is constant
     * \param sineAcceleration Acceleration (in RSW frame) that scales with sine of true anomaly
     * \param cosineAcceleration Acceleration (in RSW frame) that scales with cosine of true anomaly
     */
    EmpiricalAccelerationSettings(
            const Eigen::Vector3d& constantAcceleration = Eigen::Vector3d::Zero( ),
            const Eigen::Vector3d& sineAcceleration = Eigen::Vector3d::Zero( ),
            const Eigen::Vector3d& cosineAcceleration = Eigen::Vector3d::Zero( ) ):
        AccelerationSettings( basic_astrodynamics::empirical_acceleration ),
        constantAcceleration_( constantAcceleration ),
        sineAcceleration_( sineAcceleration ),
        cosineAcceleration_( cosineAcceleration ){ }

    // Acceleration (in RSW frame) that is constant
    Eigen::Vector3d constantAcceleration_;

    // Acceleration (in RSW frame) that scales with sine of true anomaly
    Eigen::Vector3d sineAcceleration_;

    // Acceleration (in RSW frame) that scales with cosine of true anomaly
    Eigen::Vector3d cosineAcceleration_;

};

//! @get_docstring(empiricalAcceleration)
inline std::shared_ptr< AccelerationSettings > empiricalAcceleration(
		const Eigen::Vector3d& constantAcceleration = Eigen::Vector3d::Zero( ),
		const Eigen::Vector3d& sineAcceleration = Eigen::Vector3d::Zero( ),
		const Eigen::Vector3d& cosineAcceleration = Eigen::Vector3d::Zero( )
		)
{
	return std::make_shared< EmpiricalAccelerationSettings >( constantAcceleration, sineAcceleration, cosineAcceleration );
}

// Class to define settings for yarkovsky accelerations
//! @get_docstring(YarkovskyAccelerationSettings.__docstring__)
class YarkovskyAccelerationSettings : public AccelerationSettings
{
    // Constructor
    /*
     * Constructor
     * \param yarkovskyParameter (A2) au d^{-1}
     */
    public:
    YarkovskyAccelerationSettings( const double yarkovskyParameter = 0.0 ):
    AccelerationSettings(basic_astrodynamics::yarkovsky_acceleration),
    yarkovskyParameter_(yarkovskyParameter) { }

    // Yarkovsky parameter (A2) au d^{-1}
    double yarkovskyParameter_;
};

//! @get_docstring(yarkovskyAcceleration)
inline std::shared_ptr< AccelerationSettings> yarkovskyAcceleration(
    const double yarkovskyParameter )
{
    return std::make_shared< YarkovskyAccelerationSettings >( yarkovskyParameter );
}

// Interface class that allows single interpolator to be used for thrust direction and magnitude (which are separated in
// thrust implementation)
// TODO: not exposed
class FullThrustInterpolationInterface
{
public:

    // Constructor
    /*
     * Constructor
     * \param thrustInterpolator Object that returns the total thrust vector, expressed in some reference frame B
     * \param rotationFunction Function that returns the rotation matrix from the frame B to the frame in which the
     * propagation is performed.
     */
    FullThrustInterpolationInterface(
            const std::function< Eigen::Vector3d( const double ) > thrustForceFunction,
            const std::function< Eigen::Matrix3d( ) > rotationFunction =
            [ ]( ){ return Eigen::Matrix3d::Identity( ); } ):
        thrustForceFunction_( thrustForceFunction ), rotationFunction_( rotationFunction ),
        currentThrust_( Eigen::Vector3d::Constant( TUDAT_NAN ) ), currentTime_( TUDAT_NAN ){ }

    // Function to retrieve the current thrust magnitude
    /*
     * Function to retrieve the current thrust magnitude, updates thrust to current time if needed.
     * \param time Time at which thrust must be evaluated.
     * \return  Current thrust magnitude.
     */
    double getThrustMagnitude( const double time )
    {
        updateThrust( time );

        return currentThrust_.norm( );
    }

    // Function to retrieve the current thrust direction (in the propagation frame).
    /*
     * Function to retrieve the current thrust direction (in the propagation frame)., updates thrust to current time if
     * needed.
     * \param time Time at which thrust must be evaluated.
     * \return Current thrust direction in propagation frame..
     */
    Eigen::Vector3d getThrustDirection( const double time )
    {
        updateThrust( time );
        return currentThrust_.normalized( );
    }

    // Function to reset the function to rotate to propation frame
    /*
     *  Function to reset the function to rotate to propation frame
     *  \param rotationFunction New function that returns the rotation matrix from the frame B to the frame in which the
     *  propagation is performed.
     */
    void resetRotationFunction( const std::function< Eigen::Matrix3d( ) > rotationFunction )
    {
        rotationFunction_ = rotationFunction;
    }

    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

    // Function to retrieve the thrust interpolator.
    /*
     * Function to retrieve the thrust interpolator.
     */
    std::function< Eigen::Vector3d( const double ) >  getThrustForceFunction( )
    {
        return thrustForceFunction_;
    }

private:

    // Function to update the thrust vector to the current time
    /*
     * Function to update the thrust vector to the current time
     * \param time Time at which thrust must be evaluated.
     */
    void updateThrust( const double time )
    {
        if( !( time == currentTime_ ) )
        {
            currentThrust_ = rotationFunction_( ) * thrustForceFunction_( time );
            currentTime_ = time;
        }
    }

    // Object that returns the total thrust vector, expressed in some reference frame B
    std::function< Eigen::Vector3d( const double ) > thrustForceFunction_;

    // Function that returns the rotation matrix from the frame B to the frame in which the propagation is performed.
    std::function< Eigen::Matrix3d( ) > rotationFunction_;

    // Total thrust vector (in propagation frame) computed by last call to updateThrust function.
    Eigen::Vector3d currentThrust_;

    // Time at which the last call to updateThrust was made (e.g. time associated with current thrust).
    double currentTime_;

};


// Class for providing acceleration settings for a thrust acceleration model
/*
 *  Class for providing acceleration settings for a thrust acceleration model. Settings for the direction and magnitude
 *  guidance of the thrust are provided/
 */
//! @get_docstring(ThrustAccelerationSettings.__docstring__)
class ThrustAccelerationSettings: public AccelerationSettings
{
public:


    ThrustAccelerationSettings( const std::string& engineId ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration )
    {
        engineIds_.push_back( engineId );
        useAllEngines_ = false;
    }

   ThrustAccelerationSettings( const std::vector< std::string >& engineIds ):
            AccelerationSettings( basic_astrodynamics::thrust_acceleration ), engineIds_( engineIds )
   {
        useAllEngines_ = false;
   }

   ThrustAccelerationSettings( ):
       AccelerationSettings( basic_astrodynamics::thrust_acceleration )
   {
        engineIds_ = std::vector< std::string >( );
        useAllEngines_ = true;
   }


    // Destructor.
    ~ThrustAccelerationSettings( ){ }

    std::vector< std::string > engineIds_;

    bool useAllEngines_;    


    template< typename ReturnType >
    ReturnType printDeprecationError( )
    {
        utilities::printDeprecationError(
                    "tudatpy.numerical_simulation.propagation_setup.acceleration.direction_settings/magnitude_settings",
                    "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
        return nullptr;
    }

};

inline Eigen::Vector3d applyAccelerationScalingFunction(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction,
        const std::function< double( const double) > scalingFunction,
        const double time )
{
    return accelerationFunction( time ) * scalingFunction( time );
}

//! @get_docstring(thrustAcceleration, 1)
inline std::shared_ptr< AccelerationSettings > thrustAcceleration(
        const std::vector< std::string >& engineIds )
{
    return std::make_shared< ThrustAccelerationSettings >( engineIds );
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationFromSingleEngine(
        const std::string& engineId )
{
    return std::make_shared< ThrustAccelerationSettings >( std::vector< std::string >( { engineId } ) );
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationFromAllEngines( )
{
    return std::make_shared< ThrustAccelerationSettings >( );
}


//// TODO: not exposed
//// Retrieve acceleration model (thrust).
//inline std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getLowThrustLegAccelerationSettings(
//        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::string& bodyToPropagate,
//        const std::function< double( const double ) > specificImpulseFunction,
//        const double lowThrustLegInitialTime )
//{
//    using namespace low_thrust_trajectories;

//    std::shared_ptr< simulation_setup::Body > vehicle = bodies.at( bodyToPropagate );

//    // Define thrust magnitude function from the shaped trajectory.
//    std::function< double( const double ) > thrustForceMagnitudeFunction;
//    if( lowThrustLeg->getLegModelIsForceBased( ) )
//    {
//        thrustForceMagnitudeFunction =
//                std::bind( &LowThrustLeg::getForceBasedThrustMagnitude, lowThrustLeg,
//                           std::placeholders::_1, lowThrustLegInitialTime );
//    }
//    else
//    {
//        std::function< double( const double, const double ) > thrustAccelerationMagnitudeFunction =
//                std::bind( &LowThrustLeg::getAccelerationBasedThrustMagnitude, lowThrustLeg,
//                           std::placeholders::_1, lowThrustLegInitialTime, std::placeholders::_2 );
//        std::function< double( ) > bodyMassFunction =
//                std::bind( &Body::getBodyMass, vehicle );
//        thrustForceMagnitudeFunction = [=](const double currentTime ){
//            return thrustAccelerationMagnitudeFunction( currentTime, bodyMassFunction( ) ); };
//    }

//    // Define thrust magnitude settings from thrust magnitude function.
//    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings =
//            std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
//                thrustForceMagnitudeFunction, specificImpulseFunction );


//    // Define thrust direction function from the shaped trajectory.
//    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
//            std::bind( &LowThrustLeg::getThrustDirection, lowThrustLeg,
//                       std::placeholders::_1, lowThrustLegInitialTime );

//    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
//    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
//            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

//    // Define thrust acceleration settings.
//    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
//            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
//                thrustDirectionSettings, thrustMagnitudeSettings );

//    return thrustAccelerationSettings;
//}

//! @get_docstring(CustomAccelerationSettings.__docstring__)
class CustomAccelerationSettings: public AccelerationSettings
{
public:

    CustomAccelerationSettings(
            const std::function< Eigen::Vector3d( const double ) > accelerationFunction  ):
        AccelerationSettings( basic_astrodynamics::custom_acceleration ),
        accelerationFunction_( accelerationFunction ){ }

    CustomAccelerationSettings(
            const std::function< Eigen::Vector3d( const double ) > accelerationFunction,
            const std::function< double( const double) > scalingFunction ):
        AccelerationSettings( basic_astrodynamics::custom_acceleration ),
        accelerationFunction_(
            std::bind( &applyAccelerationScalingFunction, accelerationFunction, scalingFunction,
                       std::placeholders::_1 ) ){ }

    std::function< Eigen::Vector3d( const double ) > accelerationFunction_;
};

//! @get_docstring(customAccelerationSettings)
inline std::shared_ptr< AccelerationSettings > customAccelerationSettings(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction,
        const std::function< double( const double ) > scalingFunction = nullptr )
{
    if( scalingFunction == nullptr )
    {
        return std::make_shared< CustomAccelerationSettings >(
                    accelerationFunction );
    }
    else
    {
        return std::make_shared< CustomAccelerationSettings >(
                    accelerationFunction, scalingFunction );
    }
}

// Class for providing settings for a direct tidal acceleration model, with approach of Lainey et al. (2007, 2009, ..)
/*
 *  Class for providing settings for a direct tidal acceleration model, with approach of Lainey et al. (2007, 2009, ..).
 *  Using this approach does includes the effect of tides raised by/on a planetary satelltie on the orbit of the satellite by
 *  a dedicated acceleration model, instead of modifying the gravity field coefficients of the satellite/host planet/
 */
//! @get_docstring(DirectTidalDissipationAccelerationSettings.__docstring__)
class DirectTidalDissipationAccelerationSettings: public AccelerationSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param k2LoveNumber Static k2 Love number of the satellite
     * \param timeLag Time lag of tidal bulge on satellite
     * \param includeDirectRadialComponent  True if term independent of time lag is to be included, false otherwise
     * \param useTideRaisedOnPlanet True if acceleration model is to model tide raised on planet by satellite, false if vice
     * versa
     */
    DirectTidalDissipationAccelerationSettings( const double k2LoveNumber, const double timeLag,
                                                const bool includeDirectRadialComponent = true,
                                                const bool useTideRaisedOnPlanet = true,
                                                const bool explicitLibraionalTideOnSatellite = false ):
        AccelerationSettings(
            ( useTideRaisedOnPlanet ? basic_astrodynamics::direct_tidal_dissipation_in_central_body_acceleration :
                                      basic_astrodynamics::direct_tidal_dissipation_in_orbiting_body_acceleration ) ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ),
        inverseTidalQualityFactor_( TUDAT_NAN ), tidalPeriod_( TUDAT_NAN ),
        includeDirectRadialComponent_( includeDirectRadialComponent ),
        useTideRaisedOnPlanet_( useTideRaisedOnPlanet ), explicitLibraionalTideOnSatellite_( explicitLibraionalTideOnSatellite )
    {
        if( explicitLibraionalTideOnSatellite_ && useTideRaisedOnPlanet_ )
        {
            throw std::runtime_error( "Error when creating tidal dissipation acceleration model, cannot use tide on planet and librational tide on satellite in same model" );
        }
    }

    // Constructor
    /*
     * Constructor
     * \param k2LoveNumber Static k2 Love number of the satellite
     * \param inverseTidalQualityFactor Inverse of tidal quality factor Q
     * \param period Tidal period to be considered to compute the tidal time lag
     * \param includeDirectRadialComponent  True if term independent of time lag is to be included, false otherwise
     * \param useTideRaisedOnPlanet True if acceleration model is to model tide raised on planet by satellite, false if vice
     * versa
     */
    DirectTidalDissipationAccelerationSettings( const double k2LoveNumber, const double inverseTidalQualityFactor,
                                                const double period,
                                                const bool includeDirectRadialComponent = true,
                                                const bool useTideRaisedOnPlanet = true,
                                                const bool explicitLibraionalTideOnSatellite = false ):
            AccelerationSettings(
                    ( useTideRaisedOnPlanet ? basic_astrodynamics::direct_tidal_dissipation_in_central_body_acceleration :
                      basic_astrodynamics::direct_tidal_dissipation_in_orbiting_body_acceleration ) ),
            k2LoveNumber_( k2LoveNumber ),
            timeLag_( period * std::atan( inverseTidalQualityFactor ) / ( 2.0 * mathematical_constants::PI ) ),
            inverseTidalQualityFactor_( inverseTidalQualityFactor ),
            tidalPeriod_( period ),
            includeDirectRadialComponent_( includeDirectRadialComponent ),
            useTideRaisedOnPlanet_( useTideRaisedOnPlanet ),
            explicitLibraionalTideOnSatellite_( explicitLibraionalTideOnSatellite ){ }

    // Static k2 Love number of the satellite
    double k2LoveNumber_;

    // Time lag of tidal bulge on satellite
    double timeLag_;

    // Inverse of tidal quality factor of the satellite (set to NaN if tidal lag is a direct input of the model)
    double inverseTidalQualityFactor_;

    // Period to be consider for the tides (set to Nan if tidal lag is a direct input of the model)
    double tidalPeriod_;

    // True if term independent of time lag is to be included, false otherwise
    bool includeDirectRadialComponent_;

    // True if acceleration model is to model tide raised on planet by satellite, false if vice versa
    bool useTideRaisedOnPlanet_;

    bool explicitLibraionalTideOnSatellite_;

};

//! @get_docstring(directTidalDissipationAcceleration)
inline std::shared_ptr< AccelerationSettings > directTidalDissipationAcceleration(
		const double k2LoveNumber, const double timeLag,
		const bool includeDirectRadialComponent = true,
        const bool useTideRaisedOnPlanet = true,
        const bool explicitLibraionalTideOnSatellite = false
		)
{
	return std::make_shared< DirectTidalDissipationAccelerationSettings >( k2LoveNumber, timeLag,
																		includeDirectRadialComponent,
                                                                        useTideRaisedOnPlanet,
                                                                           explicitLibraionalTideOnSatellite);
}

//! @get_docstring(directTidalDissipationAccelerationFromInvQ)
inline std::shared_ptr< AccelerationSettings > directTidalDissipationAccelerationFromInvQ(
        const double k2LoveNumber, const double inverseTidalQualityFactor,
        const double tidalPeriod,
        const bool includeDirectRadialComponent = true,
        const bool useTideRaisedOnPlanet = true,
        const bool explicitLibraionalTideOnSatellite = false
)
{
    return std::make_shared< DirectTidalDissipationAccelerationSettings >(
            k2LoveNumber, inverseTidalQualityFactor, tidalPeriod, includeDirectRadialComponent, useTideRaisedOnPlanet, explicitLibraionalTideOnSatellite);
}

// Class for providing acceleration settings for a momentum wheel desaturation acceleration model.
/*
 *  Class for providing acceleration settings for a momentum wheel desaturation acceleration model.
 *  The deltaV values for each of the desaturation maneuvers are provided by the user.
 */
//! @get_docstring(MomentumWheelDesaturationAccelerationSettings.__docstring__)
class MomentumWheelDesaturationAccelerationSettings: public AccelerationSettings
{
public:

    // Constructor.
    /*
     * Constructor.
     * \param thrustMidTimes Vector containing the midtime of each desaturation maneuver.
     * \param deltaVValues Vector containing the deltaV values of the desaturation maneuvers.
     * \param totalManeuverTime Total duration of the desaturation maneuvers.
     * \param maneuverRiseTime Rise time of the desaturation maneuvers.
     */
    MomentumWheelDesaturationAccelerationSettings(
            const std::vector< double > thrustMidTimes,
            const std::vector< Eigen::Vector3d > deltaVValues,
            const double totalManeuverTime,
            const double maneuverRiseTime ): AccelerationSettings( basic_astrodynamics::momentum_wheel_desaturation_acceleration ),
        thrustMidTimes_( thrustMidTimes ), deltaVValues_( deltaVValues ),
        totalManeuverTime_( totalManeuverTime ), maneuverRiseTime_( maneuverRiseTime ){ }

    // Vector containing the midtime of each desaturation maneuver.
    std::vector< double > thrustMidTimes_;

    // Vector containing the deltaV values of the momentum wheel desaturation maneuvers.
    std::vector< Eigen::Vector3d > deltaVValues_;

    // Total desaturation maneuver time.
    double totalManeuverTime_;

    // Desaturation maneuvers rise time.
    double maneuverRiseTime_;

};

//! @get_docstring(momentumWheelDesaturationAcceleration)
inline std::shared_ptr< AccelerationSettings > momentumWheelDesaturationAcceleration(
		const std::vector< double > thrustMidTimes,
		const std::vector< Eigen::Vector3d > deltaVValues,
		const double totalManeuverTime,
		const double maneuverRiseTime
		)
{
	return std::make_shared< MomentumWheelDesaturationAccelerationSettings >( thrustMidTimes, deltaVValues,
																		   totalManeuverTime, maneuverRiseTime);
}

// Typedef defining a list of acceleration settings, set up in the same manner as the
// AccelerationMap typedef.
typedef std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > >
SelectedAccelerationMap;

typedef std::map< std::string, std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > > >
SelectedAccelerationList;


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_ACCELERATIONSETTINGS_H
