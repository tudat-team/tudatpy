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

#include "tudat/astro/electromagnetism/lorentzStaticMagneticForce.h"
#include "tudat/astro/electromagnetism/lorentzStaticMagneticAcceleration.h"


namespace tudat
{
namespace electromagnetism
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

} // namespace electromagnetism
} // namespace tudat
