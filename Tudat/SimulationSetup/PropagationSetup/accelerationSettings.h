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


class FullThrustInterpolationInterface
{
public:
    FullThrustInterpolationInterface(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustInterpolator ):
        thrustInterpolator_( thrustInterpolator ){ }

    double getThrustMagnitude( const double time )
    {
        return thrustInterpolator_->interpolate( time ).norm( );
    }

    Eigen::Vector3d getThrustDirection( const double time )
    {
        return thrustInterpolator_->interpolate( time ).normalized( );
    }

private:

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustInterpolator_;
};

//! Class for providing acceleration settings for a thrust acceleration model
/*!
 *  Class for providing acceleration settings for a thrust acceleration model. Settings for the direction and magnitude
 *  guidance of the thrust are provided/
 */
class ThrustAccelerationSettings: public AccelerationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustDirectionGuidanceSettings Settings for the direction of the thrust
     * \param thrustMagnitudeSettings Settings for the magnitude of the thrust
     */
    ThrustAccelerationSettings(
            const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
            const boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration ),
        thrustDirectionGuidanceSettings_( thrustDirectionGuidanceSettings ),
        thrustMagnitudeSettings_( thrustMagnitudeSettings ){ }

    ThrustAccelerationSettings(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > thrustDirectionGuidanceSettings,
            const boost::function< double( const double ) > specificImpulseFunction ):
        AccelerationSettings( basic_astrodynamics::thrust_acceleration )
    {
        boost::shared_ptr< FullThrustInterpolationInterface > interpolatorInterface =
                boost::make_shared< FullThrustInterpolationInterface >( thrustDirectionGuidanceSettings );
        thrustDirectionGuidanceSettings_ = boost::make_shared< CustomThrustDirectionSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustDirection, interpolatorInterface, _1 ) );
        thrustMagnitudeSettings_ =  boost::make_shared< FromFunctionThrustEngineSettings >(
                    boost::bind( &FullThrustInterpolationInterface::getThrustMagnitude, interpolatorInterface, _1 ),
                    specificImpulseFunction );
    }


    //! Destructor.
    ~ThrustAccelerationSettings( ){ }

    //! Settings for the direction of the thrust
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings_;

    //! Settings for the magnitude of the thrust
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings_;
};

//! Typedef defining a list of acceleration settings, set up in the same manner as the
//! AccelerationMap typedef.
typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr<
AccelerationSettings > > > > SelectedAccelerationMap;

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_ACCELERATIONSETTINGS_H
