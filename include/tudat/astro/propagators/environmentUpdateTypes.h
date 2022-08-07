/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_ENVIRONMENTUPDATETYPES_H
#define TUDAT_ENVIRONMENTUPDATETYPES_H

#include <string>
#include <vector>
#include <map>

namespace tudat
{

namespace propagators
{

//! Enum defining types of environment model updates that can be done.
enum EnvironmentModelsToUpdate
{
    body_translational_state_update = 0,
    body_rotational_state_update = 1,
    body_mass_update = 2,
    spherical_harmonic_gravity_field_update = 3,
    vehicle_flight_conditions_update = 4,
    radiation_pressure_interface_update = 5,
    radiation_source_model_update = 6,
    radiation_pressure_target_model_update = 7
};

//! Function to extend existing list of required environment update types
/*!
 * Function to extend existing list of required environment update types
 * \param environmentUpdateList List of environment updates to extend
 * (passed by reference and modified by function)
 * \param updatesToAdd List of environment updates that are to be added to environmentUpdateList
 */
void addEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        environmentUpdateList,
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
        updatesToAdd );

} // namespace propagators

} // namespace tudat
#endif // TUDAT_ENVIRONMENTUPDATETYPES_H
