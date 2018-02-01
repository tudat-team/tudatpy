/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K.F. Astrodynamics I, Delft University of Technology, 2010.
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *
 */

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

namespace tudat
{
namespace gravitation
{

//! Compute perturbing acceleration by third body
Eigen::Vector3d computeThirdBodyPerturbingAcceleration(
        const double gravitationalParameterOfPerturbingBody,
        const Eigen::Vector3d& positionOfPerturbingBody,
        const Eigen::Vector3d& positionOfAffectedBody,
        const Eigen::Vector3d& positionOfCentralBody )
// Using chapter 4 of (Wakker, 2010).
{
    // Return acceleration.
    return computeGravitationalAcceleration( positionOfAffectedBody,
                                             gravitationalParameterOfPerturbingBody,
                                             positionOfPerturbingBody ) -
            computeGravitationalAcceleration( positionOfCentralBody,
                                              gravitationalParameterOfPerturbingBody,
                                              positionOfPerturbingBody );
}

} // namespace gravitation
} // namespace tudat
