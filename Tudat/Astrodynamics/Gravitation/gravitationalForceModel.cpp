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
 *      100916    K. Kumar          File created.
 *      100916    K. Kumar          Filename modified.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor corrections to include statements and comments.
 *      110113    K. Kumar          Changed setBody( ) argument to pointer; added pointer to
 *                                  GravityFieldModel.
 *      110119    K. Kumar          Changed computeStateDerivatives( ) to computeForce( ).
 *      110202    K. Kumar          Updated code to make use of the State and
 *                                  CartesianPositionElements classes.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110815    K. Kumar          Changed filename and class name; changed computeForce( )
 *                                  function and added setMass( ) function.
 *      120209    K. Kumar          Changed class into free functions.
 *      120316    D. Dirkx          Added old GravitationalForceModel class, with note that this
 *                                  will be removed shortly.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 */

#include <cmath>

#include <TudatCore/Basics/utilityMacros.h>

#include "Tudat/Astrodynamics/Gravitation/gravitationalAccelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/gravitationalForceModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace force_models
{

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double universalGravitationalParameter, const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double massOfBodyExertingForce, const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce *
            acceleration_models::computeGravitationalAcceleration(
                universalGravitationalParameter, positionOfBodySubjectToForce,
                massOfBodyExertingForce, positionOfBodyExertingForce );
}

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double massOfBodySubjectToForce, const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double gravitationalParameterOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce *
            acceleration_models::computeGravitationalAcceleration(
                positionOfBodySubjectToForce, gravitationalParameterOfBodyExertingForce,
                positionOfBodyExertingForce );
}

//! Compute force due to gravity field.
void GravitationalForceModel::computeForce( StatePointer state, const double time )
{
    TUDAT_UNUSED_PARAMETER( time );

    force_ = gravityFieldModel_->getGradientOfPotential( state->state.segment( 0, 3 ) )
            * bodySubjectToForce_->getMass( );
}

} // namespace force_models
} // namespace astrodynamics
} // namespace tudat
