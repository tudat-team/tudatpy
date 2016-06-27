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

#ifndef TUDAT_ACCELERATIONMODELTYPES_H
#define TUDAT_ACCELERATIONMODELTYPES_H

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"


namespace tudat
{

namespace basic_astrodynamics
{


//! List of accelerations available in simulations
/*!
 *  List of accelerations available in simulations. Acceleration models not defined by this
 *  given enum cannot be used for automatic acceleration model setup.
 */
enum AvailableAcceleration
{
    undefined_acceleration,
    central_gravity,
    aerodynamic,
    cannon_ball_radiation_pressure,
    spherical_harmonic_gravity,
    third_body_central_gravity,
    third_body_spherical_harmonic_gravity
};

//! Function to identify the derived class type of an acceleration model.
/*!
 *  Function to identify the derived class type of an acceleration model. The type must be defined
 *  in the AvailableAcceleration enum to be recognized by this function.
 *  \param accelerationModel Acceleration model of which the type is to be identified.
 *  \return Type of the accelerationModel, as identified by AvailableAcceleration enum.
 */
AvailableAcceleration getAccelerationModelType(
        const boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
        accelerationModel );

//! Function to get all acceleration models of a given type from a list of models
/*!
 * Function to get all acceleration models of a given type from a list of models
 * \param fullList List of acceleration models
 * \param modelType Type for which all models are to be retrieved
 * \return Subset of fullList for which the acceleration model type is modelType
 */
std::vector< boost::shared_ptr< AccelerationModel3d > > getAccelerationModelsOfType(
        const std::vector< boost::shared_ptr< AccelerationModel3d > >& fullList,
        const AvailableAcceleration modelType );

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_ACCELERATIONMODELTYPES_H
