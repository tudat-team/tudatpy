/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H
#define TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H

#include <vector>
#include <string>
#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/gravityDeformationSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create Maxwell gravity deformation model.
/*!
 *  Function to create Maxwell gravity deformation model from perturbing and deforming bodies.
 *  \param deformingBody Pointer to object of deforming body.
 *  \param perturbingBody Pointer to object of perturbing body.
 *  \param nameOfDeformingBody Name of deforming body.
 *  \param nameOfPerturbingBody Name of perturbing body.
 *  \param deformationSettings Settings for gravity deformation model that is to be created.
 *  \return Maxwell gravity deformation model pointer.
 */
std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > createMaxwellGravityFieldDeformationModel(
        const std::shared_ptr< simulation_setup::Body > deformingBody,
        const std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBody,
        const std::string& nameOfDeformingBody,
        const std::vector< std::string >& nameOfPerturbingBody,
        const std::shared_ptr< GravityDeformationSettings > deformationSettings );

basic_astrodynamics::GravityDeformationModelMap createGravityDeformationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedGravityDeformationModelMap& gravityDeformationSettings );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H
