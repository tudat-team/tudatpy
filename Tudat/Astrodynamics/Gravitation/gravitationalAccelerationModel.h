/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120209    K. Kumar          File created.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.

 *
 *    References
 */

#ifndef TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
#define TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H

#include <Eigen/Core>

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Compute gravitational acceleration.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with body2.
 * The basis for this gravitational acceleration is that body2 is a point mass,
 * generating acceleration due to Newton's gravitational force:
 * \f[
 *      \bar{a}_{gravity} = -\frac{G * m_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$G\f$ is the universal gravitational constant, \f$m_{2}\f$ is the mass of body2,
 * and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to body 2, with respect to an
 * inertial (barycentric) reference frame.
 * \param universalGravitationalParameter Universal gravitational constant [m^3 kg^-1 s^-2].
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to acceleration
 *          (body1) [m].
 * \param massOfBodyExertingAcceleration Mass of body exerting acceleration (body2) [kg].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 [m s^-2].
 * \sa computeGravitationalForce.
 */
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalParameter,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2). The basis for this gravitational acceleration is that the body exerting the
 * acceleration is a point mass, generating acceleration due to Newton's gravitational force:
 * \f[
 *      \bar{a}_{gravity} = -\frac{\mu_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to body 2. with respect to an
 * inertial (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 [m s^-2].
 * \sa computeGravitationalForce.
 */
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
