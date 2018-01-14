/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/SimulationSetup/PropagationSetup/createThrustModelGuidance.h"
// #include "Tudat/Mathematics/Interpolators/createInterpolator.h"


namespace tudat
{

namespace simulation_setup
{


//! Class for providing settings for acceleration model.
/*!
 *  Class for providing settings for acceleration model. This class is a functional (base) class for
 *  settings of acceleration models that  require no information in addition to their type.
 *  Classes defining settings for acceleration models requiring additional information must be
 *  derived from this class.
 *  Bodies exerting and undergong acceleration are set externally from this class.
 *  This class can be used for the easy setup of acceleration models
 *  (see createAccelerationModels.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */
class AccelerationSettings
{
public:

    //! Constructor, sets type of acceleration.
    /*!
     *  Constructor, sets type of acceleration.
     *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
     */
    AccelerationSettings( const basic_astrodynamics::AvailableAcceleration accelerationType ):
        accelerationType_( accelerationType ){ }

    //! Destructor.
    virtual ~AccelerationSettings( ){ }

    //! Type of acceleration from AvailableAcceleration enum.
    basic_astrodynamics::AvailableAcceleration accelerationType_;

};

//! Class for providing settings for spherical harmonics acceleration model.
/*!
 *  Class for providing settings for spherical harmonics acceleration model,
 *  specifically the maximum degree and order up to which the field is to be expanded. Note that
 *  the minimum degree and order are currently always set to zero.
 */
class SphericalHarmonicAccelerationSettings: public AccelerationSettings
{
public:
    //! Constructor to set maximum degree and order that is to be taken into account.
    /*!
     *  Constructor to set maximum degree and order that is to be taken into account.
     *  \param maximumDegree Maximum degree
     *  \param maximumOrder Maximum order
     */
    SphericalHarmonicAccelerationSettings( const int maximumDegree,
                                           const int maximumOrder ):
        AccelerationSettings( basic_astrodynamics::spherical_harmonic_gravity ),
        maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder ){ }

    //! Maximum degree that is to be used for spherical harmonic acceleration
    int maximumDegree_;

    //! Maximum order that is to be used for spherical harmonic acceleration
    int maximumOrder_;
};

//! Class for providing acceleration settings for mutual spherical harmonics acceleration model.
/*!
 *  Class for providing acceleration settings for mutual spherical harmonics acceleration model,
 *  specifically the maximum degree and order up to which the fields of the bodies are be expanded.
 *  Please note that the minimum degrees and orders are currently always set to zero.
 */
class MutualSphericalHarmonicAccelerationSettings: public AccelerationSettings
{
public:

    //! Constructor to set maximum degrees and orders that are to be taken into account.
    /*!
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

    //! Maximum degree of body exerting acceleration.
    int maximumDegreeOfBodyExertingAcceleration_;

    //! Maximum order of body exerting acceleration.
    int maximumOrderOfBodyExertingAcceleration_;

    //! Maximum degree of body undergoing acceleration.
    int maximumDegreeOfBodyUndergoingAcceleration_;

    //! Maximum order of body undergoing acceleration.
    int maximumOrderOfBodyUndergoingAcceleration_;

    //! Maximum degree of central body (only releveant for 3rd body acceleration).
    int maximumDegreeOfCentralBody_;

    //! Maximum order of central body (only releveant for 3rd body acceleration).
    int maximumOrderOfCentralBody_;
};

//! Class to proivide settings for typical relativistic corrections to the dynamics of an orbiter.
/*!
 *  Class to proivide settings for typical relativistic corrections to the dynamics of an orbiter: the
 *  Schwarzschild, Lense-Thirring and de Sitter terms. An excellent introduction to
 *  these models is given in 'General Relativity and Space Geodesy' by L. Combrinck (2012).
 */
class RelativisticAccelerationCorrectionSettings: public AccelerationSettings
{
public:

    //! Constructor
    /*!
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

    //! Boolean denoting wheter the Schwarzschild term is used.
    bool calculateSchwarzschildCorrection_;

    //! Boolean denoting wheter the Lense-Thirring term is used.
    bool calculateLenseThirringCorrection_;

    //! Boolean denoting wheter the de Sitter term is used.
    bool calculateDeSitterCorrection_;

    //! Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite)
    std::string primaryBody_;

    //! Constant angular momentum of central body
    Eigen::Vector3d centralBodyAngularMomentum_;
};

//! Class to define settings for empirical accelerations
class EmpiricalAccelerationSettings: public AccelerationSettings
{
public:

    //! Constructor
    /*!
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

    //! Acceleration (in RSW frame) that is constant
    Eigen::Vector3d constantAcceleration_;

    //! Acceleration (in RSW frame) that scales with sine of true anomaly
    Eigen::Vector3d sineAcceleration_;

    //! Acceleration (in RSW frame) that scales with cosine of true anomaly
    Eigen::Vector3d cosineAcceleration_;
};

//! Interface class that allows single interpolator to be used for thrust direction and magnitude (which are separated in
//! thrust implementation)
class FullThrustInterpolationInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustInterpolator Object that returns the total thrust vector, expressed in some reference frame B
     * \param rotationFunction Function that returns the rotation matrix from the frame B to teh frame in which the
     * propagation is performed.
     */
    FullThrustInterpolationInterface(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator<
            double, Eigen::Vector3d > > thrustInterpolator,
            const boost::function< Eigen::Matrix3d( ) > rotationFunction =
            boost::lambda::constant( Eigen::Matrix3d::Identity( ) ) ):
        thrustInterpolator_( thrustInterpolator ), rotationFunction_( rotationFunction ),
        currentThrust_( Eigen::Vector3d::Constant( TUDAT_NAN ) ), currentTime_( TUDAT_NAN ){ }

    //! Function to retrieve the current thrust magnitude
    /*!
     * Function to retrieve the current thrust magnitude, updates thrust to current time if needed.
     * \param time Time at which thrust must be evaluated.
     * \return  Current thrust magnitude.
     */
    double getThrustMagnitude( const double time )
    {
        updateThrust( time );

        return currentThrust_.norm( );
    }

    //! Function to retrieve the current thrust direction (in the propagation frame).
    /*!
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

    //! Function to reset the function to rotate to propation frame
    /*!
     *  Function to reset the function to rotate to propation frame
     *  \param rotationFunction New function that returns the rotation matrix from the frame B to teh frame in which the
     *  propagation is performed.
     */
    void resetRotationFunction( const boost::function< Eigen::Matrix3d( ) > rotationFunction )
    {
        rotationFunction_ = rotationFunction;
    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

    //! Function to retrieve the thrust interpolator.
    /*!
     * Function to retrieve the thrust interpolator.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > getThrustInterpolator( )
    {
        return thrustInterpolator_;
    }

private:

    //! Function to update the thrust vector to the current time
    /*!
     * Function to update the thrust vector to the current time
     * \param time Time at which thrust must be evaluated.
     */
    void updateThrust( const double time )
    {
        if( !( time == currentTime_ ) )
        {
            currentThrust_ = rotationFunction_( ) * thrustInterpolator_->interpolate( time );
            currentTime_ = time;
        }
    }

    //! Object that returns the total thrust vector, expressed in some reference frame B
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustInterpolator_;

    //! Function that returns the rotation matrix from the frame B to teh frame in which the propagation is performed.
    boost::function< Eigen::Matrix3d( ) > rotationFunction_;

    //! Total thrust vector (in propagation frame) computed by last call to updateThrust function.
    Eigen::Vector3d currentThrust_;

    //! Time at which the last call to updateThrust was made (e.g. time associated with current thrust).
    double currentTime_;
};


//! Enum defining identifiers of frames in which a user-specifief thrust is defined.
enum ThrustFrames
{
    unspecified_thurst_frame = -1,
    inertial_thurst_frame = 0,
    lvlh_thrust_frame = 1
};

//! Class for providing acceleration settings for a thrust acceleration model
/*!
 *  Class for providing acceleration settings for a thrust acceleration model. Settings for the direction and magnitude
 *  guidance of the thrust are provided/
 */
class ThrustAccelerationSettings: public AccelerationSettings
{
public:

    //! Constructor from separate magnitude and diretion settings.
    /*!
     * Constructor from separate magnitude and diretion settings.
     * \param thrustDirectionGuidanceSettings Settings for the direction of the thrust
     * \param thrustMagnitudeSettings Settings for the magnitude of the thrust
     */
    ThrustAccelerationSettings(
            const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
            const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration ),
        thrustDirectionGuidanceSettings_( thrustDirectionGuidanceSettings ),
        thrustMagnitudeSettings_( thrustMagnitudeSettings ),
        thrustFrame_( unspecified_thurst_frame ){ }

    //! Constructor used for defining total thrust vector (in local or inertial frame) from interpolator using
    //! variable specific impulse
    /*!
     * Constructor used for defining total thrust vector (in local or inertial frame) from interpolator using
     * variable specific impulse
     * \param dataInterpolationSettings Settings to create the interpolator that returns the thrust as a function of
     * time in frame defined by thrustFrame
     * \param specificImpulseFunction Function returning the specific impulse as a function of time
     * \param thrustFrame Identifier of frame in which thrust returned by fullThrustInterpolator is expressed
     * \param centralBody Central body identifier for thrustFrame (if needed; empty by default).
     */
    ThrustAccelerationSettings(
            const boost::shared_ptr< interpolators::DataInterpolationSettings< double, Eigen::Vector3d > >&
            dataInterpolationSettings,
            const boost::function< double( const double ) > specificImpulseFunction,
            const ThrustFrames thrustFrame = unspecified_thurst_frame,
            const std::string centralBody = "" ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration ),
        constantSpecificImpulse_( TUDAT_NAN ), thrustFrame_( thrustFrame ),
        centralBody_( centralBody ), dataInterpolationSettings_( dataInterpolationSettings )
    {
        interpolatorInterface_ = boost::make_shared< FullThrustInterpolationInterface >(
                    interpolators::createOneDimensionalInterpolator( dataInterpolationSettings ) );
        thrustDirectionGuidanceSettings_ = boost::make_shared< CustomThrustDirectionSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustDirection, interpolatorInterface_, _1 ) );
        thrustMagnitudeSettings_ =  boost::make_shared< FromFunctionThrustEngineSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustMagnitude, interpolatorInterface_, _1 ),
                    specificImpulseFunction, boost::lambda::constant( true ),
                    boost::lambda::constant( Eigen::Vector3d::UnitX( ) ),
                    boost::bind( &FullThrustInterpolationInterface::resetTime, interpolatorInterface_, _1 ) );
    }

    //! Constructor used for defining total thrust vector (in local or inertial frame) from interpolator using constant
    //! specific impulse
    /*!
     * Constructor used for defining total thrust vector (in local or inertial frame) from interpolator using constant
     * specific impulse
     * \param dataInterpolationSettings Settings to create the interpolator that returns the thrust as a function of
     * time in frame defined by thrustFrame
     * \param constantSpecificImpulse Constant specific impulse
     * \param thrustFrame Identifier of frame in which thrust returned by fullThrustInterpolator is expressed
     * \param centralBody Central body identifier for thrustFrame (if needed; empty by default).
     */
    ThrustAccelerationSettings(
            const boost::shared_ptr< interpolators::DataInterpolationSettings< double, Eigen::Vector3d > >&
            dataInterpolationSettings,
            const double constantSpecificImpulse,
            const ThrustFrames thrustFrame = unspecified_thurst_frame,
            const std::string centralBody = "" ):
        ThrustAccelerationSettings( dataInterpolationSettings,
                                    boost::lambda::constant( constantSpecificImpulse ),
                                    thrustFrame,
                                    centralBody )
    {
        constantSpecificImpulse_ = constantSpecificImpulse;
    }


    //! Destructor.
    ~ThrustAccelerationSettings( ){ }


    //! Settings for the direction of the thrust
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings_;

    //! Settings for the magnitude of the thrust
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings_;

    //! Constant specific impulse used when determining the direction and magnitude of thrust from an interpolator.
    //! NaN if the specific impulse is not constant (i.e. is defined using a boost::function).
    double constantSpecificImpulse_ = TUDAT_NAN;

    //! Identifier of frame in which thrust returned by fullThrustInterpolator is expressed.
    /*!
     *  Identifier of frame in which thrust returned by fullThrustInterpolator is expressed. Unspecifief by default,
     *  only used if interpolatorInterface_ is set
     */
    ThrustFrames thrustFrame_;

    //! Central body identifier for thrustFrame.
    /*!
     *  Central body identifier for thrustFrame. Empty by default,
     *  only used if interpolatorInterface_ is set
     */
    std::string centralBody_;

    //! Settings to create the interpolator interface
    boost::shared_ptr< interpolators::DataInterpolationSettings< double, Eigen::Vector3d > > dataInterpolationSettings_;

    //! Interface object used when full thrust (direction and magnitude) are defined by a single user-supplied interpolation.
    boost::shared_ptr< FullThrustInterpolationInterface > interpolatorInterface_;

};

//! Class for providing settings for a direct tidal acceleration model, with approach of Lainey et al. (2007, 2009, ..)
/*!
 *  Class for providing settings for a direct tidal acceleration model, with approach of Lainey et al. (2007, 2009, ..).
 *  Using this approach does includes the effect of tides raised by/on a planetary satelltie on the orbit of the satellite by
 *  a dedicated acceleration model, instead of modifying the gravity field coefficients of the satellite/host planet/
 */
class DirectTidalDissipationAccelerationSettings: public AccelerationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param k2LoveNumber Static k2 Love number of the satellite
     * \param timeLag Time lag of tidal bulge on satellite
     * \param includeDirectRadialComponent  True if term independent of time lag is to be included, false otherwise
     * \param useTideRaisedOnPlanet True if acceleration model is to model tide raised on planet by satellite, false if vice
     * versa
     */
    DirectTidalDissipationAccelerationSettings( const double k2LoveNumber, const double timeLag,
                                                const bool includeDirectRadialComponent = true,
                                                const bool useTideRaisedOnPlanet = true ):
        AccelerationSettings( basic_astrodynamics::direct_tidal_dissipation_acceleration ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), includeDirectRadialComponent_( includeDirectRadialComponent ),
        useTideRaisedOnPlanet_( useTideRaisedOnPlanet ){ }

    //! Static k2 Love number of the satellite
    double k2LoveNumber_;

    //! Time lag of tidal bulge on satellite
    double timeLag_;

    //! True if term independent of time lag is to be included, false otherwise
    bool includeDirectRadialComponent_;

    //! True if acceleration model is to model tide raised on planet by satellite, false if vice versa
    bool useTideRaisedOnPlanet_;
};

//! Typedef defining a list of acceleration settings, set up in the same manner as the
//! AccelerationMap typedef.
typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > >
SelectedAccelerationMap;

typedef std::map< std::string, std::vector< std::pair< std::string, boost::shared_ptr< AccelerationSettings > > > >
SelectedAccelerationList;

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_ACCELERATIONSETTINGS_H
