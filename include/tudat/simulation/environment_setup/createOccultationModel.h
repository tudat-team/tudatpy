/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_CREATEOCCULTATIONMODEL_H
#define TUDATBUNDLE_CREATEOCCULTATIONMODEL_H


#include "tudat/astro/electromagnetism/occultationModel.h"
#include "tudat/simulation/environment_setup/body.h"


namespace tudat
{
namespace simulation_setup
{

//! Function to create occultation model from a list of occulting bodies.
/*!
 * Function to create occultation model from a list of occulting bodies.
 *
 * \param occultingBodies Names of bodies to occult a source as seen from an observer
 * \param bodies System of bodies
 * \return Shared pointer to radiation source model
 */
std::shared_ptr<electromagnetism::OccultationModel> createOccultationModel(
        const std::vector<std::string>& occultingBodies,
        const SystemOfBodies& bodies);

} // tudat
} // electromagnetism

#endif //TUDATBUNDLE_CREATEOCCULTATIONMODEL_H
