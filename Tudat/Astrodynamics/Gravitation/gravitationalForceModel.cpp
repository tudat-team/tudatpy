/*    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
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
 */

#include <cmath>
#include "Tudat/Astrodynamics/Gravitation/gravitationalAccelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/gravitationalForceModel.h"

namespace tudat
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
            tudat::acceleration_models::computeGravitationalAcceleration(
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
            tudat::acceleration_models::computeGravitationalAcceleration(
                positionOfBodySubjectToForce, gravitationalParameterOfBodyExertingForce,
                positionOfBodyExertingForce );
}

//! Compute force due to gravity field.
void GravitationalForceModel::computeForce( State* pointerToState, double time = 0.0 )
{
    CartesianPositionElements cartesianPositionElements_;
    cartesianPositionElements_.state = pointerToState->state.segment( 0, 3 );
    force_ = pointerToGravityFieldModel_->getGradientOfPotential( &cartesianPositionElements_ )
            * pointerToBodySubjectToForce_->getMass( );
}

} // namespace force_models
} // namespace tudat
