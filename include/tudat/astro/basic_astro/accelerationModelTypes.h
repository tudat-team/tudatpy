/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ACCELERATIONMODELTYPES_H
#define TUDAT_ACCELERATIONMODELTYPES_H

#include "tudat/astro/basic_astro/customAccelerationModel.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"
#include "tudat/astro/gravitation/thirdBodyPerturbation.h"
#include "tudat/astro/gravitation/ringGravityModel.h"
#include "tudat/astro/gravitation/polyhedronGravityModel.h"
#include "tudat/astro/gravitation/directTidalDissipationAcceleration.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/propulsion/thrustAccelerationModel.h"
#include "tudat/astro/propulsion/massRateFromThrust.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"
#include "tudat/astro/relativity/einsteinInfeldHoffmannAcceleration.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/astro/propulsion/massRateFromThrust.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"

namespace tudat
{

namespace basic_astrodynamics
{

// List of accelerations available in simulations
/*
 *  List of accelerations available in simulations. Acceleration models not defined by this
 *  given enum cannot be used for automatic acceleration model setup.
 *
 */
//! @get_docstring(AvailableAcceleration.__docstring__)
enum AvailableAcceleration
{
    undefined_acceleration,
    point_mass_gravity,
    aerodynamic,
    cannon_ball_radiation_pressure,
    spherical_harmonic_gravity,
    mutual_spherical_harmonic_gravity,
    polyhedron_gravity,
    ring_gravity,
    third_body_point_mass_gravity,
    third_body_spherical_harmonic_gravity,
    third_body_mutual_spherical_harmonic_gravity,
    third_body_polyhedron_gravity,
    third_body_ring_gravity,
    thrust_acceleration,
    relativistic_correction_acceleration,
    empirical_acceleration,
    direct_tidal_dissipation_in_central_body_acceleration,
    direct_tidal_dissipation_in_orbiting_body_acceleration,
    radiation_pressure,
    momentum_wheel_desaturation_acceleration,
    custom_acceleration,
    einstein_infeld_hoffmann_acceleration,
    yarkovsky_acceleration
};

// Function to get a string representing a 'named identification' of an acceleration type
/*
 * Function to get a string representing a 'named identification' of an acceleration type
 * \param accelerationType Type of acceleration model.
 * \return String with acceleration id.
 */
std::string getAccelerationModelName( const AvailableAcceleration accelerationType );

// List of model types for body mass rates.
/*
*  List of model types for body mass rates available in simulations. Mass rate models not defined by this
*  given enum cannot be used for automatic mass rate model setup.
*/
//! @get_docstring(AvailableMassRateModels.__docstring__)
enum AvailableMassRateModels
{
    undefined_mass_rate_model,
    custom_mass_rate_model,
    from_thrust_mass_rate_model
};

// Function to identify the derived class type of an acceleration model.
/*
 *  Function to identify the derived class type of an acceleration model. The type must be defined
 *  in the AvailableAcceleration enum to be recognized by this function.
 *  \param accelerationModel Acceleration model of which the type is to be identified.
 *  \return Type of the accelerationModel, as identified by AvailableAcceleration enum.
 */
AvailableAcceleration getAccelerationModelType(
        const std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
        accelerationModel );

// Function to identify the type of a mass rate model.
/*
 *  Function to identify the type of a mass rate model. The type must be defined
 *  in the AvailableMassRateModels enum to be recognized by this function.
 *  \param massRateModel Mass rate model of which the type is to be identified.
 *  \return Type of the massRateModel, as identified by AvailableMassRateModels enum.
 */
AvailableMassRateModels getMassRateModelType(
        const std::shared_ptr< MassRateModel > massRateModel );

// Function to get all acceleration models of a given type from a list of models
/*
 * Function to get all acceleration models of a given type from a list of models
 * \param fullList List of acceleration models
 * \param modelType Type for which all models are to be retrieved
 * \return Subset of fullList for which the acceleration model type is modelType
 */
std::vector< std::shared_ptr< AccelerationModel3d > > getAccelerationModelsOfType(
        const std::vector< std::shared_ptr< AccelerationModel3d > >& fullList,
        const AvailableAcceleration modelType );



// Function to check whether an acceleration type is a direct gravitational acceleration
/*
 * Function to check whether an acceleration type is a direct gravitational acceleration, e.g. a gravitational
 * acceleration that is not from a third-body.
 * \param accelerationType Acceleration type for which it is to be checked whether it is direct gravitational.
 * \return True if acceleration type is direct gravitational, false otherwise.
 */
bool isAccelerationDirectGravitational( const AvailableAcceleration accelerationType );

// Function to check whether an acceleration type is a third-body gravitational acceleration
/*
 * Function to check whether an acceleration type is a third-body gravitational acceleration
 * \param accelerationType Acceleration type for which it is to be checked whether it is third-body gravitational.
 * \return True if acceleration type is third-body gravitational, false otherwise.
 */
bool isAccelerationFromThirdBody( const AvailableAcceleration accelerationType );

// Function to get the third-body counterpart of a direct gravitational acceleration type
/*
 * Function to get the third-body counterpart of a direct gravitational acceleration type, e.g. a third_body_point_mass_gravity
 * for a point_mass_gravity input. Function throws an exception is input is not direct gravitational
 * \param accelerationType Acceleration type for which the third-body counterpart is to be determined.
 * \return Third-body counterpart of accelerationType.
 */
AvailableAcceleration getAssociatedThirdBodyAcceleration( const AvailableAcceleration accelerationType );

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_ACCELERATIONMODELTYPES_H
