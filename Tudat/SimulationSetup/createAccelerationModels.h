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

#ifndef TUIDAT_CREATEACCELERATIONMODELS_H
#define TUDAT_CREATEACCELERATIONMODELS_H

#include <vector>
#include <string>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/SimulationSetup/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to determine if a given frame is an inertial frame.
/*!
 *  Function to determine if a given frame is an inertial frame. Currently a frame identified as
 *  "SSB" (solar system barycenter), "inertial" or "" (empty) is recognized as inertial.
 *  \param frame Name of frame for which it is to be determined whether it is inertial.
 *  \return True if inertial, false if not.
 */
bool isFrameInertial( const std::string& frame );

//! Function to create central gravity acceleration model.
/*!
 *  Function to create central gravity acceleration model from bodies exerting and undergoing
 *  acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the central gravity
 *  acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the central gravity
 *   acceleration.
 *  \param useCentralBodyFixedFrame Boolean setting whether the central attraction of body
 *  undergoing acceleration on body exerting acceleration is to be included in acceleration model.
 *  Should be set to true in case the body undergoing acceleration is a celestial body
 *  (with gravity field) and integration is performed in the frame centered at the body exerting
 *  acceleration.
 */
boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d >
createCentralGravityAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame );

//! Function to create spherical harmonic gravity acceleration model.
/*!
 *  Function to create spherical harmonic gravity acceleration model from bodies exerting and
 *  undergoing acceleration.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the spherical
 *  harmonic gravity acceleration.
 *  \param nameOfBodyUndergoingAcceleration Name of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of body that is exerting the spherical harmonic
 *  gravity acceleration.
 *  \param accelerationSettings Settings for acceleration model that is to be created (should
 *  be of derived type associated with spherical harmonic acceleration.
 *  \param useCentralBodyFixedFrame Boolean setting whether the central attraction of body
 *  undergoing acceleration on body exerting acceleration is to be included in acceleration model.
 *  Should be set to true in case the body undergoing acceleration is a celestial body
 *  (with gravity field) and integration is performed in the frame centered at the body exerting
 *  acceleration.
 */
boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >
createSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame );

//! Function to create a third body central gravity acceleration model.
/*!
 *  Function to create a third body central gravity acceleration model from bodies exerting and
 *  undergoing acceleration, as well as the central body, w.r.t. which the integration is to be
 *  performed.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting the acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be
 *  calculated.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the central
 *  gravity acceleration.
 *  \param nameOfCentralBody Name of central body in frame cenetered at which acceleration is to
 *  be calculated.
 *  \return Pointer to object for calculating central gravity acceleration between bodies.
 */
boost::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody );

//! Function to create acceleration model object.
/*!
 *  Function to create acceleration model object.
 *  Type of requested model is checked and corresponding factory function is called.
 *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
 *  \param bodyExertingAcceleration Pointer to object of body that is exerting acceleration,
 *  \param accelerationSettings Settings for acceleration model that is to be created.
 *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
 *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
 *  \param centralBody Pointer to central body in frame centered at which acceleration is to be
 *  calculated (only relevant for third body accelerations).
 *  \param nameOfCentralBody Name of central body in frame cenetered at which acceleration is to
 *  be calculated (only relevant for third body accelerations).
 *  \return Acceleration model pointer.
 */
boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
createAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody = boost::shared_ptr< Body >( ),
        const std::string& nameOfCentralBody = "" );

//! Function to create acceleration models from a map of bodies and acceleration model types.
/*!
 *  Function to create acceleration models from a map of bodies and acceleration model types.
 *  The return type can be used to identify both the body undergoing and exerting acceleration.
 *  \param bodyMap List of pointers to bodies required for the creation of the acceleration model
 *  objects.
 *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
 *  acceleration(s) on which bodies.
 *  \param centralBodies Map of central bodies for each body undergoing acceleration.
 *  \return List of acceleration model objects, in form of AccelerationMap.
 */
AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies );

}

}
#endif // TUDAT_CREATEACCELERATIONMODELS_H
