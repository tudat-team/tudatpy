/*    Copyright (c) 2010-2016, Delft University of Technology
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
     * \param maximumDegreeOfCentralBody Maximum degree of central body (only releveant for 3rd body acceleration).
     * \param maximumOrderOfCentralBody Maximum order of central body (only releveant for 3rd body acceleration).
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

    //! Constructor used for defining total thrust vector (in local or inertial frame) from interpolator
    /*!
     * Constructor used for defining total thrust vector (in local or inertial frame) from interpolator
     * \param fullThrustInterpolator Interpolator that returns the thrust as a function of time in
     * frame defined by thrustFrame
     * \param specificImpulseFunction Function returning the specific impulse as a function of time
     * \param thrustFrame Identifier of frame in which thrust returned by fullThrustInterpolator is expressed
     * \param centralBody Central body identifier for thrustFrame (if needed; empty by default).
     */
    ThrustAccelerationSettings(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > >
            fullThrustInterpolator,
            const boost::function< double( const double ) > specificImpulseFunction,
            const ThrustFrames thrustFrame = unspecified_thurst_frame,
            const std::string centralBody = "" ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration ), thrustFrame_( thrustFrame ),
        centralBody_( centralBody )
    {
        interpolatorInterface_ =
                boost::make_shared< FullThrustInterpolationInterface >( fullThrustInterpolator );
        thrustDirectionGuidanceSettings_ = boost::make_shared< CustomThrustDirectionSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustDirection, interpolatorInterface_, _1 ) );
        thrustMagnitudeSettings_ =  boost::make_shared< FromFunctionThrustEngineSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustMagnitude, interpolatorInterface_, _1 ),
                    specificImpulseFunction, boost::lambda::constant( true ),
                    Eigen::Vector3d::UnitX( ),
                    boost::bind( &FullThrustInterpolationInterface::resetTime, interpolatorInterface_, _1 ) );
    }

    //! Destructor.
    ~ThrustAccelerationSettings( ){ }


    //! Settings for the direction of the thrust
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings_;

    //! Settings for the magnitude of the thrust
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings_;

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

    //! Interface object used when full thrust (direction and magnitude) are defined by a single user-supplied interpolation.
    boost::shared_ptr< FullThrustInterpolationInterface > interpolatorInterface_;

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
