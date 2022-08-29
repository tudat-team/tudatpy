/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/astro/electromagnetism/cannonBallRadiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/solarSailAcceleration.h"
#include "tudat/astro/gravitation/thirdBodyPerturbation.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace simulation_setup
{



//! Function to create a direct (i.e. not third-body) gravitational acceleration (of any type)
/*!
 * Function to create a direct (i.e. not third-body) gravitational acceleration of any type (i.e. point mass,
 * spherical harmonic, mutual spherical harmonic).
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the gravitational acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the gravitational acceleration.
 *  \param accelerationSettings Settings object for the gravitational acceleration.
 *  \param nameOfCentralBody Name of central body in frame centered at which acceleration is to
 *  be calculated.
 *  \param isCentralBody Boolean defining whether the body undergoing the acceleration is the central body for a
 *  third-body acceleration, of which the return object of this funciton is one of the sub-parts. Boolean is
 *  only used when creating mutual spherical harmonic acceleration, to ensure the correct usage of the acceleration
 *  settings.
 *  \return Direct gravitational acceleration model of requested settings.
 */
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createDirectGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfCentralBody = "",
        const bool isCentralBody = 0 );

//! Function to create a third-body gravitational acceleration (of any type)
/*!
 * Function to create a direct third-body gravitational acceleration of any type (i.e. point mass,
 * spherical harmonic, mutual spherical harmonic).
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the gravitational acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be calculated.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the gravitational acceleration.
 *  \param nameOfCentralBody Name of central body in frame centered at which acceleration is to
 *  be calculated.
 *  \param accelerationSettings Settings object for the gravitational acceleration.
 *  \return Third-body gravitational acceleration model of requested settings.
 */
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createThirdBodyGravitationalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings );

//! Function to create gravitational acceleration (of any type)
/*!
 * Function to create a third-body or direct gravitational acceleration of any type (i.e. point mass,
 * spherical harmonic, mutual spherical harmonic).
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the gravitational acceleration.
 *  \param accelerationSettings Settings object for the gravitational acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be calculated.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the gravitational acceleration.
 *  \param nameOfCentralBody Name of central body in frame centered at which acceleration is to
 *  be calculated.
 *  \return Gravitational acceleration model of requested settings.
 */
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > createGravitationalAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody );

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
std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d >
createCentralGravityAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
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
 *  \param useDegreeZeroTerm Boolean denoting whether to manually remove the C(0,0) term from
 *  cosine coefficients (default true).
 *  \return Spherical harmonic gravity acceleration model pointer.
 */
std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
createSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool useDegreeZeroTerm = true );

//! Function to create mutual spherical harmonic gravity acceleration model.
/*!
 *  Function to create mutual spherical harmonic gravity acceleration model from bodies exerting and
 *  undergoing acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the mutual spherical
 *  harmonic gravity acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the mutual spherical harmonic
 *  gravity acceleration.
 *  \param accelerationSettings Settings for acceleration model that is to be created (should
 *  be of derived type associated with mutual spherical harmonic acceleration).
 *  \param useCentralBodyFixedFrame Boolean setting whether the central attraction of body
 *  undergoing acceleration on body exerting acceleration is to be included in acceleration model.
 *  Should be set to true in case the body undergoing acceleration is a celestial body
 *  (with gravity field) and integration is performed in the frame centered at the body exerting
 *  acceleration.
 *  \param acceleratedBodyIsCentralBody Boolean defining whether the body undergoing the acceleration is the central body
 *  for a third-body acceleration, of which the return object of this funciton is one of the sub-parts.
 *  \return Mutual spherical harmonic gravity acceleration model pointer.
 */
std::shared_ptr< gravitation::MutualSphericalHarmonicsGravitationalAccelerationModel >
createMutualSphericalHarmonicsGravityAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame,
        const bool acceleratedBodyIsCentralBody );

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
std::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody );

//! Function to create a third body spheric harmonic gravity acceleration model.
/*!
 *  Function to create a third body spheric harmonic gravity acceleration model from bodies exerting and
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
 *  \param accelerationSettings Settings for acceleration model that is to be created (should
 *  be of derived type associated with spherical harmonic acceleration).
 *  \return Pointer to object for calculating third-body spheric harmonic gravity acceleration between bodies.
 */
std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
createThirdBodySphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings );

//! Function to create a third body mutual spheric harmonic gravity acceleration model.
/*!
 *  Function to create a third body mutual spheric harmonic gravity acceleration model from bodies exerting and
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
 *  \param accelerationSettings Settings for acceleration model that is to be created (should
 *  be of derived type associated with mutual spherical harmonic acceleration).
 *  \return Pointer to object for calculating third-body mutual spheric harmonic gravity acceleration between bodies.
 */
std::shared_ptr< gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >
createThirdBodyMutualSphericalHarmonicGravityAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody,
        const std::shared_ptr< AccelerationSettings > accelerationSettings );

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
std::shared_ptr< aerodynamics::AerodynamicAcceleration >
createAerodynamicAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
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
std::shared_ptr< electromagnetism::CannonBallRadiationPressureAcceleration >
createCannonballRadiationPressureAcceleratioModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration );

//! Function to create a panelled radiation pressure acceleration model.
/*!
 *  Function to create a panelled radiation pressure automatically creates all required
 *  links to environment models, vehicle properies and frame conversions
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration,
 *  i.e. body emitting the radiation.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \return Pointer to object for calculating cannonball radiation pressures acceleration.
 */
std::shared_ptr< electromagnetism::PanelledRadiationPressureAcceleration > createPanelledRadiationPressureAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration );

//! Function to create a non-ideal solar pressure radiation acceleration.
/*!
 * Function to create a non-ideal solar pressure radiation acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration.
 *  \param centralBody Pointer to object of central body.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the acceleration.
 *  \return Solar sail acceleration model pointer.
 */
std::shared_ptr< electromagnetism::SolarSailAcceleration >
createSolarSailAccelerationModel(const std::shared_ptr< Body > bodyUndergoingAcceleration,
                                 const std::shared_ptr< Body > bodyExertingAcceleration,
                                 const std::shared_ptr< Body > centralBody,
                                 const std::string& nameOfBodyUndergoingAcceleration,
                                 const std::string& nameOfBodyExertingAcceleration);

//! Function to create a thrust acceleration model.
/*!
 *  Function to create a thrust acceleration model. Creates all required
 *  links to environment models, vehicle properies and frame conversions.
 *  \param accelerationSettings Settings of thrust acceleration model.
 *  \param bodies List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param nameOfBodyUndergoingThrust Name of body that is undergoing the thrust acceleration
 *  \return Pointer to object for calculating thrust acceleration.
 */
std::shared_ptr< propulsion::ThrustAcceleration >
createThrustAcceleratioModel(
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const SystemOfBodies& bodies,
        const std::string& nameOfBodyUndergoingThrust );

//! Function to create a direct tidal acceleration model, according to approach of Lainey et al. (2007, 2009, ...)
/*!
 *  Function to create a direct tidal acceleration model, according to approach of Lainey et al. (2007, 2009, ...).
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of main body that is exerting the acceleration
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration
 *  \param accelerationSettings Settings for the acceleration model
 *  \return Pointer to object for calculating acceleration.
 */
std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > createDirectTidalDissipationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings );

//! Function to create a momentum wheel desaturation acceleration model.
/*!
 *  Function to create a momentum wheel desaturation acceleration model.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of main body that is exerting the acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \param accelerationSettings Settings for the acceleration model.
 *  \return Pointer to object for calculating acceleration.
 */
std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > createMomentumWheelDesaturationAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const  std::shared_ptr< AccelerationSettings > accelerationSettings );

//! Function to create an orbiter relativistic correction acceleration model
/*!
 *  Function to create an orbiter relativistic correction acceleration model (Schwarzschild, Lense-Thirring and/or
 *  de Sitter terms).
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of main body that is exerting the acceleration (e.g. Earth for an orbiter
 *  around the Earth).
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration
 *  \param accelerationSettings Settings for the acceleration model.
 *  \param bodies List of pointers to bodies that comprise the full environment.
 *  \return Pointer to object for calculating relativistic correction acceleration.
 */
std::shared_ptr< relativity::RelativisticAccelerationCorrection > createRelativisticCorrectionAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const SystemOfBodies& bodies );

//! Function to create empirical acceleration model.
/*!
 *  Function to create empirical acceleration model from body undergoing acceleration and body wrt which its orbit is determined
 *  for setting determining phase of once per orbit empirical accelerations.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body wrt whic orbit of body undergoing the acceleration is calculated.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body wrt whic orbit of body undergoing the acceleration is calculated.
 *  \param accelerationSettings Object containing additional settings for acceleration model, must be of type EmpiricalAccelerationSettings
 *  \return Pointer to object for calculating empiricalt acceleration.
 */
std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > createEmpiricalAcceleration(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings );


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
 *  \param bodies List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \return Acceleration model pointer.
 */
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
createAccelerationModel(
        const std::shared_ptr< Body > bodyUndergoingAcceleration,
        const std::shared_ptr< Body > bodyExertingAcceleration,
        const std::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::shared_ptr< Body > centralBody = std::shared_ptr< Body >( ),
        const std::string& nameOfCentralBody = "",
        const SystemOfBodies& bodies = SystemOfBodies( ) );

//! Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
/*!
 * Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
 * \param selectedAccelerationPerBody List of acceleration settings per body.
 * \return selectedAccelerationPerBody, put in order to ensure correct model creation.
 */
SelectedAccelerationList orderSelectedAccelerationMap( const SelectedAccelerationMap& selectedAccelerationPerBody );

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodies List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param centralBodies Map of central bodies for each body undergoing acceleration.
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
//basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
//        const SystemOfBodies& bodies,
//        const SelectedAccelerationMap& selectedAccelerationPerBody,
//        const std::map< std::string, std::string >& centralBodies );
//! Function to create a set of acceleration models from a map of bodies and acceleration model types.
inline basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    basic_astrodynamics::AccelerationMap accelerationModelMap;

    // Put selectedAccelerationPerBody in correct order
    SelectedAccelerationList orderedAccelerationPerBody =
            orderSelectedAccelerationMap( selectedAccelerationPerBody );

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationList::const_iterator bodyIterator =
         orderedAccelerationPerBody.begin( ); bodyIterator != orderedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        std::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !ephemerides::isFrameInertial( currentCentralBodyName ) )
        {
            if( bodies.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                            std::string( "Error, could not find non-inertial central body ") +
                            currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                            " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodies.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodies
        if( bodies.count( bodyUndergoingAcceleration ) ==  0 )
        {
            throw std::runtime_error(
                        std::string( "Error when making acceleration models, requested forces" ) +
                        "acting on body " + bodyUndergoingAcceleration  +
                        ", but no such body found in map of bodies" );
        }

        // Declare map of acceleration models acting on current body.
        basic_astrodynamics::SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        std::vector< std::pair< std::string, std::shared_ptr< AccelerationSettings > > > thrustAccelerationSettings;

        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > currentAcceleration;

        // Iterate over all bodies exerting an acceleration
        for( unsigned int i = 0; i < accelerationsForBody.size( ); i++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = accelerationsForBody.at( i ).first;

            // Check if body exerting acceleration is included in bodies
            if( bodies.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            if( !( accelerationsForBody.at( i ).second->accelerationType_ == basic_astrodynamics::thrust_acceleration ) )
            {
                currentAcceleration = createAccelerationModel( bodies.at( bodyUndergoingAcceleration ),
                                                               bodies.at( bodyExertingAcceleration ),
                                                               accelerationsForBody.at( i ).second,
                                                               bodyUndergoingAcceleration,
                                                               bodyExertingAcceleration,
                                                               currentCentralBody,
                                                               currentCentralBodyName,
                                                               bodies );


                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            currentAcceleration );
            }
            else
            {
                thrustAccelerationSettings.push_back( accelerationsForBody.at( i ) );
            }

        }

        // Create thrust accelerations last
        for( unsigned int i = 0; i < thrustAccelerationSettings.size( ); i++ )
        {
            currentAcceleration = createAccelerationModel( bodies.at( bodyUndergoingAcceleration ),
                                                           bodies.at( thrustAccelerationSettings.at( i ).first ),
                                                           thrustAccelerationSettings.at( i ).second,
                                                           bodyUndergoingAcceleration,
                                                           thrustAccelerationSettings.at( i ).first,
                                                           currentCentralBody,
                                                           currentCentralBodyName,
                                                           bodies );


            // Create acceleration model.
            mapOfAccelerationsForBody[ thrustAccelerationSettings.at( i ).first  ].push_back(
                        currentAcceleration );
        }


        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    return accelerationModelMap;
}

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodies List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param propagatedBodies List of bodies that are to be propagated
 *  \param centralBodies List of central bodies for each body undergoing acceleration (in same order as propagatedBodies).
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONMODELS_H
