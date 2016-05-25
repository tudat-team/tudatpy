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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to identify the derived class type of an acceleration model.
AvailableAcceleration getAccelerationModelType(
        const boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
        accelerationModel )
{
    using namespace tudat::aerodynamics;
    using namespace tudat::electro_magnetism;
    using namespace tudat::gravitation;

    // Nominal type is undefined
    AvailableAcceleration accelerationType = undefined_acceleration;

    // Check for each accelerarion mdoel type implemented as AvailableAcceleration.
    if( boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                accelerationModel ) != NULL )
    {
        accelerationType = central_gravity;
    }
    else if( boost::dynamic_pointer_cast< CannonBallRadiationPressureAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = cannon_ball_radiation_pressure;
    }
    else if( boost::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = third_body_central_gravity;
    }
    else if( boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != NULL  )
    {
        accelerationType = spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel ) != NULL )
    {
        accelerationType = mutual_spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< AerodynamicAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = aerodynamic;
    }
    else
    {
        throw std::runtime_error(
                    "Error, acceleration model not identified when getting acceleration type." );
    }

    // Return identified type.
    return accelerationType;

}


} // namespace simulation_setup

} // namespace tudat
