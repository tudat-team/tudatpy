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

#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticForce.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticAcceleration.h"


namespace tudat
{
namespace electro_magnetism
{

//! Compute acceleration due to static magnetic field
Eigen::Vector3d computeLorentzAccelerationDueToStaticMagneticField(
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration,
        const double massOfBodySubjectToAcceleration )
{
    return computeLorentzForceDueToStaticMagneticField(
                velocityOfBodySubjectToAcceleration, localMagneticField,
                chargeOfBodySubjectToAcceleration ) / massOfBodySubjectToAcceleration;
}

} // namespace electro_magnetism
} // namespace tudat
