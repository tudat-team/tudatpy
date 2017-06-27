/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicTorque.h"

namespace tudat
{

namespace simulation_setup
{

typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< TorqueSettings > > > > SelectedTorqueMap;


boost::shared_ptr< aerodynamics::AerodynamicTorque > createAerodynamicTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

boost::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > createSecondDegreeGravitationalTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
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
boost::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const boost::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

}

}

#endif // TUDAT_CREATETORQUEMODEL_H
