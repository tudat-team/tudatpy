/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      121025    K. Kumar         File created.
 *      140107    J. Geul          Added MEE and common accelerations.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_STATE_INDICES_H
#define TUDAT_STATE_INDICES_H

namespace tudat
{
namespace orbital_element_conversions
{

//! Cartesian elements indices.
enum CartesianElementIndices
{
    xCartesianPositionIndex = 0,
    yCartesianPositionIndex = 1,
    zCartesianPositionIndex = 2,
    xCartesianVelocityIndex = 3,
    yCartesianVelocityIndex = 4,
    zCartesianVelocityIndex = 5
};

//! Keplerian elements indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex,
    eccentricityIndex,
    inclinationIndex,
    argumentOfPeriapsisIndex,
    longitudeOfAscendingNodeIndex,
    trueAnomalyIndex,
    semiLatusRectumIndex = 0
};

//! Modified equinoctial element vector indices.
enum ModifiedEquinoctialElementVectorIndices
{
    semiParameterIndex,
    fElementIndex,
    gElementIndex,
    hElementIndex,
    kElementIndex,
    trueLongitudeIndex
};

//! Unified State Model indices.
enum UnifiedStateModelElementIndices
{
    CHodographIndex = 0,
    Rf1HodographIndex = 1,
    Rf2HodographIndex = 2,
    epsilon1QuaternionIndex = 3,
    epsilon2QuaternionIndex = 4,
    epsilon3QuaternionIndex = 5,
    etaQuaternionIndex = 6
};

//! Cartesian acceleration indices.
enum CartesianAccelerationElementIndices
{
    xCartesianAccelerationIndex,
    yCartesianAccelerationIndex,
    zCartesianAccelerationIndex
};

//! Acceleration indices in CSN frame for orbital elements.
enum CSNAccelerationElementIndices
{
    cAccelerationIndex,
    sAccelerationIndex,
    nAccelerationIndex
};

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_STATE_INDICES_H
