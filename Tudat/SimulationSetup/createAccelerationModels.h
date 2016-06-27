/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEACCELERATIONMODELS_H
#define TUDAT_CREATEACCELERATIONMODELS_H

#include <vector>
#include <string>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/SimulationSetup/accelerationSettings.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create central gravity acceleration model.
/*!
 *  Function to create central gravity acceleration model from bodies exerting and undergoing
 *  acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the central gravity
 *  acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the central gravity
 *   acceleration.
 *  \param useCentralBodyFixedFrame Boolean setting whether the central attraction of body
 *  undergoing acceleration on body exerting acceleration is to be included in acceleration model.
 *  Should be set to true in case the body undergoing acceleration is a celestial body
 *  (with gravity field) and integration is performed in the frame centered at the body exerting
 *  acceleration.
 *  \return Central gravity acceleration model pointer.
 */
boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d >
createCentralGravityAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame );

//! Function to create spherical harmonic gravity acceleration model.
/*!
 *  Function to create spherical harmonic gravity acceleration model from bodies exerting and
 *  undergoing acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the spherical
 *  harmonic gravity acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the spherical harmonic
 *  gravity acceleration.
 *  \param accelerationSettings Settings for acceleration model that is to be created (should
 *  be of derived type associated with spherical harmonic acceleration.
 *  \param useCentralBodyFixedFrame Boolean setting whether the central attraction of body
 *  undergoing acceleration on body exerting acceleration is to be included in acceleration model.
 *  Should be set to true in case the body undergoing acceleration is a celestial body
 *  (with gravity field) and integration is performed in the frame centered at the body exerting
 *  acceleration.
 *  \return Spherical harmonic gravity acceleration model pointer.
 */
boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >
createSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame );

//! Function to create a third body central gravity acceleration model.
/*!
 *  Function to create a third body central gravity acceleration model from bodies exerting and
 *  undergoing acceleration, as well as the central body, w.r.t. which the integration is to be
 *  performed.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be
 *  calculated.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the central
 *  gravity acceleration.
 *  \param nameOfCentralBody Name of central body in frame cenetered at which acceleration is to
 *  be calculated.
 *  \return Pointer to object for calculating central gravity acceleration between bodies.
 */
boost::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody );

//! Function to create an aerodynamic acceleration model.
/*!
 *  Function to create an aerodynamic acceleration model, automatically creates all required
 *  links to environment models, vehicle properies and frame conversions
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration,
 *  i.e. body with the atmosphere through which the accelerated body is flying.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \return Pointer to object for calculating aerodynamic acceleration.
 */
boost::shared_ptr< aerodynamics::AerodynamicAcceleration >
createAerodynamicAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration );

//! Function to create a cannonball radiation pressure acceleration model.
/*!
 *  Function to create a cannonball radiation pressure automatically creates all required
 *  links to environment models, vehicle properies and frame conversions
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration,
 *  i.e. body emitting the radiation.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \return Pointer to object for calculating cannonball radiation pressures acceleration.
 */
boost::shared_ptr< electro_magnetism::CannonBallRadiationPressureAcceleration >
createCannonballRadiationPressureAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration );

//! Function to create acceleration model object.
/*!
 *  Function to create acceleration model object.
 *  Type of requested model is checked and corresponding factory function is called.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting acceleration,
 *  \param accelerationSettings Settings for acceleration model that is to be created.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be
 *  calculated (optional, only relevant for third body accelerations).
 *  \param nameOfCentralBody Name of central body in frame cenetered at which acceleration is to
 *  be calculated (optional, only relevant for third body accelerations).
 *  \return Acceleration model pointer.
 */
boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
createAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody = boost::shared_ptr< Body >( ),
        const std::string& nameOfCentralBody = "" );

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodyMap List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param centralBodies Map of central bodies for each body undergoing acceleration.
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies );

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodyMap List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param propagatedBodies List of bodies that are to be propagated
 *  \param centralBodies List of central bodies for each body undergoing acceleration (in same order as propagatedBodies).
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies );

} // namespace simulation_setup

} // namespace tudat
#endif // TUDAT_CREATEACCELERATIONMODELS_H
