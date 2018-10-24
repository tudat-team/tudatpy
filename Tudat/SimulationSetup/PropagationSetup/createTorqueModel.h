/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATETORQUEMODEL_H
#define TUDAT_CREATETORQUEMODEL_H


#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/torqueSettings.h"
#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicGravitationalTorque.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicTorque.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/customTorque.h"

namespace tudat
{

namespace simulation_setup
{

typedef std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< TorqueSettings > > > > SelectedTorqueMap;

//! Function to create an aerodynamic torque model.
std::shared_ptr< basic_astrodynamics::InertialTorqueModel > createInertialTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::string& nameOfBodyUndergoingTorque );

//! Function to create an aerodynamic torque model.
/*!
 *  Function to create an aerodynamic torque model, automatically creates all required
 *  links to environment models, vehicle properies and frame conversions
 *  \param bodyUndergoingTorque Pointer to object of body that is being accelerated.
 *  \param bodyExertingTorque Pointer to object of body that is exerting the torque,
 *  i.e. body with the atmosphere through which the accelerated body is flying.
 *  \param nameOfBodyUndergoingTorque Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingTorque Name of object of body that is exerting the torque.
 *  \return Pointer to object for calculating aerodynamic torque.
 */
std::shared_ptr< aerodynamics::AerodynamicTorque > createAerodynamicTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

//! Function to create a second-degree gravitational torque.
/*!
 * Function to create a second-degree gravitational torque, exerted by a point mass
 *  \param bodyUndergoingTorque Pointer to object of body that is being accelerated.
 *  \param bodyExertingTorque Pointer to object of body that is exerting the gravitational torque.
 *  \param nameOfBodyUndergoingTorque Name of body that is being accelerated.
 *  \param nameOfBodyExertingTorque Name of body that is exerting the gravitational torque.
 *  \return Pointer to object for calculating gravitational torque.
 */
std::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > createSecondDegreeGravitationalTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

//! Function to create a spherical harmonic gravitational torque
/*!
 * Function to create spherical harmonic gravitational torque, exerted by a point mass
*  \param bodyUndergoingTorque Pointer to object of body that is being accelerated.
*  \param bodyExertingTorque Pointer to object of body that is exerting the gravitational torque.
*  \param torqueSettings Settings for the torque that is to be created
*  \param nameOfBodyUndergoingTorque Name of body that is being accelerated.
*  \param nameOfBodyExertingTorque Name of body that is exerting the gravitational torque.
*  \return Direct gravitational torque model of requested settings.
*/
std::shared_ptr< gravitation::SphericalHarmonicGravitationalTorqueModel > createSphericalHarmonicGravitationalTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

//! Function to create torque model object.
/*!
 *  Function to create torque model object.
 *  Type of requested model is checked and corresponding factory function is called.
 *  \param bodyUndergoingTorque Pointer to object of body that is being accelerated.
 *  \param bodyExertingTorque Pointer to object of body that is exerting torque,
 *  \param torqueSettings Settings for torque model that is to be created.
 *  \param nameOfBodyUndergoingTorque Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingTorque Name of object of body that is exerting the torque.
 *  \return Torque model pointer.
 */
std::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

//! Function to create torque models from a map of bodies and torque model settings.
/*!
 *  Function to create torque models from a map of bodies and torque model settings.
 *  The return type can be used to identify both the body undergoing and exerting torque.
 *  \param bodyMap List of pointers to bodies required for the creation of the torque model
 *  objects.
 *  \param selectedTorquePerBody List identifying which bodies exert which type of
 *  torque(s) on which bodies.
 *  \param propagatedBodies List of bodies that are to be propagated.
 *  \return Torque models for the input map of bodies, based on the torque model settings.
 */
basic_astrodynamics::TorqueModelMap createTorqueModelsMap(
        const NamedBodyMap& bodyMap,
        SelectedTorqueMap selectedTorquePerBody,
        const std::vector< std::string >& propagatedBodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif // TUDAT_CREATETORQUEMODEL_H
