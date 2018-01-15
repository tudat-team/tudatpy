/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <Eigen/Dense>

#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute Lorentz Force due to static magnetic field
Eigen::Vector3d computeLorentzForceDueToStaticMagneticField(
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration )
{
    //Return Lorentz force due to static magnetic field
    return chargeOfBodySubjectToAcceleration * velocityOfBodySubjectToAcceleration.cross( localMagneticField );

}

} // namespace electro_magnetism
} // namespace tudat
