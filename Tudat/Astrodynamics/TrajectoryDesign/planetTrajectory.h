/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120606    P. Musegaas       First creation of code.
 *      120827    P. Musegaas       Adaptation to own ephemeris type.
 *
 *    References
 *
 */

// This contains some functions to facilitate making plots of planet trajectories.

#include <vector>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Return a vector of positions and times from ephemeris data for a certain epoch and duration.
void returnPlanetTrajectory( const ephemerides::EphemerisPointer& ephemerisPtr,
                             const double centralBodyGravitationalParameter,
                             const double startingEpochMJD2000,
                             const double duration,
                             const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime );

//! Return a propagated vector of positions and times from ephemeris data for one revolution.
void returnSingleRevolutionPlanetTrajectory(
        const ephemerides::EphemerisPointer& ephemerisPtr,
        const double centralBodyGravitationalParameter,
        const double startingEpochMJD2000,
        const double maximumTimeStep,
        std::vector < Eigen::Vector3d >& positionVector,
        std::vector < double >& timeVector,
        const double startingTime = 0. );

} // namespace spaceTrajectories
} // namespace tudat

