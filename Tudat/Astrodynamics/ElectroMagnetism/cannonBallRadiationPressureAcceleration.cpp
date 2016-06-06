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
 *      121123    D. Dirkx          File created.
 *      130124    K. Kumar          Added missing file header; updated layout; migrated force
 *                                  free function to separate file; added acceleration free
 *                                  function.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure acceleration using a cannon-ball model.
Eigen::Vector3d computeCannonBallRadiationPressureAcceleration(
        const double radiationPressure,
        const Eigen::Vector3d& vectorToSource,
        const double area,
        const double radiationPressureCoefficient,
        const double mass )
{
    return computeCannonBallRadiationPressureForce(
                radiationPressure, vectorToSource, area, radiationPressureCoefficient ) / mass;
}

//! Get radiation pressure acceleration.
Eigen::Vector3d CannonBallRadiationPressureAcceleration::getAcceleration( )
{
    return computeCannonBallRadiationPressureAcceleration(
                currentRadiationPressure_, currentVectorToSource_, currentArea_,
                currentRadiationPressureCoefficient_, currentMass_ );
}

//! Update member variables used by the radiation pressure acceleration model.
void CannonBallRadiationPressureAcceleration::updateMembers( const double currentTime )
{
    if( !( this->currentTime_ == currentTime ) )
    {
        currentVectorToSource_ = ( sourcePositionFunction_( )
                                   - acceleratedBodyPositionFunction_( ) ).normalized( );
        currentRadiationPressure_ = radiationPressureFunction_( );
        currentRadiationPressureCoefficient_ = radiationPressureCoefficientFunction_( );
        currentArea_ = areaFunction_( );
        currentMass_ = massFunction_( );
    }
}

} // namespace electro_magnetism
} // namespace tudat
