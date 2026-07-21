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

#ifndef TUDAT_GRAVITY_DEFORMATION_MODEL_TYPES_H
#define TUDAT_GRAVITY_DEFORMATION_MODEL_TYPES_H

#include <vector>
#include <map>
#include <unordered_map>

#include <memory>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"

namespace tudat
{
namespace basic_astrodynamics
{

enum GravityDeformationType { maxwell_deformation = 0 };

GravityDeformationType getGravityDeformationModelType(
        std::shared_ptr< basic_astrodynamics::GravityDeformationModel > gravityDeformationModel );

//! Get the descriptive name of a gravity deformation model type.
std::string getGravityDeformationModelName( const GravityDeformationType gravityDeformationType );

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif  // TUDAT_GRAVITY_DEFORMATION_MODEL_TYPES_H
