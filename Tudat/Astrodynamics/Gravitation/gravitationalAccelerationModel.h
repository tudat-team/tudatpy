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
 *      120209    K. Kumar          File created.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      Wakker.
 *      Melman, J. PhD Thesis, 2012.
 */

#ifndef TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
#define TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H

#include <map>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Indices for Cartesian position components.
/*!
 * Vector indices for Cartesian position components.
 */
enum cartesianPositionIndices
{
    cartesianXPositionIndex, cartesianYPositionIndex, cartesianZPositionIndex
};

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
 * \param universalGravitationalConstant Universal gravitational constant [m^3 kg^-1 s^-2].
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to acceleration
 *          (body1) [m].
 * \param massOfBodyExertingAcceleration Mass of body exerting acceleration (body2) [kg].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 [m s^-2].
 * \sa computeGravitationalForce.
 */
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalConstant,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with body2.
 * The basis for this gravitational acceleration is that body2 is a point mass, generating
 * acceleration due to Newton's gravitational force:
 * \f[
 *      \bar{a}_{gravity} = -\frac{\mu_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to body 2, with respect to an
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

//! Compute gravitational acceleration due to J2.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2) with an irregular gravity field. The body exerting the acceleration has an
 * irregular gravity field. The acceleration due to the J2-coefficient is given by (Melman, 2012):
 * \f{eqnarray*}{
 *      {a}_{gravity,x} &=& -\mu_{2}*\frac{x_{1}-x_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,y} &=& -\mu_{2}*\frac{y_{1}-y_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,z} &=& \frac{-\mu_{2}}{r^{2}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (3-5*\hat{z}^2) * \hat{z}  \\
 * \f}
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * \f$\hat{z} = \frac{z}{r}\f$, \f$x\f$, \f$y\f$ and \f$z\f$ are the Cartesian position components,
 * \f$r\f$ is the radial position, and \f$J_{2}\f$ is the second zonal coefficient of the gravity
 * field of body2. The positions and accelerations are given with respect to an inertial
 * (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to J2-acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param j2CoefficientOfGravityField J2-coefficient, describing irregularity of the gravity field
 *          of body2 [-].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 due to J2-effect [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationDueToJ2(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j2CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration due to J3.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2) with an irregular gravity field. The body exerting the acceleration has an
 * irregular gravity field. The acceleration due to the J3-coefficient is given by (Melman, 2012):
 * \f{eqnarray*}{
 *      {a}_{gravity,x} &=& -\mu_{2}*\frac{x_{1}-x_{2}}{r^{3}} * \frac{5}{2}*J_{3}
 *                          * \left(\frac{R}{r}\right)^{3} * (3-7*\hat{z}^2) * \hat{z}  \\
 *      {a}_{gravity,y} &=& -\mu_{2}*\frac{y_{1}-y_{2}}{r^{3}} * \frac{5}{2}*J_{3}
 *                          * \left(\frac{R}{r}\right)^{3} * (3-7*\hat{z}^2) * \hat{z}  \\
 *      {a}_{gravity,z} &=& \frac{-\mu_{2}}{r^{2}} * \frac{5}{2}*J_{3}
 *                          * \left(\frac{R}{r}\right)^{3}
 *                          * (-\frac{3}{5}+6*\hat{z}^2-7*\hat{z}^4)  \\
 * \f}
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * \f$\hat{z} = \frac{z}{r}\f$, \f$x\f$, \f$y\f$ and \f$z\f$ are the Cartesian position components,
 * \f$r\f$ is the radial position, and \f$J_{3}\f$ is the third zonal coefficient of the gravity
 * field of body2. The positions and accelerations are given with respect to an inertial
 * (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to J3-acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param j3CoefficientOfGravityField J3-coefficient, describing irregularity of the gravity field
 *          of body2 [-].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 due to J3-effect [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationDueToJ3(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j3CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration due to J4.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2) with an irregular gravity field. The body exerting the acceleration has an
 * irregular gravity field. The acceleration due to the J4-coefficient is given by (Melman, 2012):
 * \f{eqnarray*}{
 *      {a}_{gravity,x} &=& \mu_{2}*\frac{x_{1}-x_{2}}{r^{3}} * \frac{35}{8}*J_{4}
 *                          * \left(\frac{R}{r}\right)^{4}
 *                          * (\frac{3}{7}-6*\hat{z}^2+9*\hat{z}^{4}) \\
 *      {a}_{gravity,y} &=& \mu_{2}*\frac{y_{1}-y_{2}}{r^{3}} * \frac{35}{8}*J_{4}
 *                          * \left(\frac{R}{r}\right)^{4}
 *                          * (\frac{3}{7}-6*\hat{z}^2+9*\hat{z}^{4}) \\
 *      {a}_{gravity,z} &=& \frac{\mu_{2}}{r^{2}} * \frac{35}{8}*J_{4}
 *                          * \left(\frac{R}{r}\right)^{4}
 *                          * (\frac{15}{7}-10*\hat{z}^2+9*\hat{z}^4) * \hat{z} \\
 * \f}
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * \f$\hat{z} = \frac{z}{r}\f$, \f$x\f$, \f$y\f$ and \f$z\f$ are the Cartesian position components,
 * \f$r\f$ is the radial position, and \f$J_{4}\f$ is the fourth zonal coefficient of the gravity
 * field of body2. The positions and accelerations are given with respect to an inertial
 * (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to J4-acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param j4CoefficientOfGravityField J4-coefficient, describing irregularity of the gravity field
 *          of body2 [-].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 due to J4-effect [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationDueToJ4(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j4CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration zonal sum.
/*!
 * Computes sum of gravitational acceleration for given zonal terms experienced by body1, due to
 * its interaction with another body (body2) with an irregular gravity field. The body exerting
 * the acceleration has an irregular gravity field. This function essentially serves as a
 * wrapper to compute the total gravitational acceleration up to a certain zonal terms (specified
 * as input parameter). So given \f$n\f$ as the maximum order of the zonal terms, the total
 * acceleration is computed as:
 * \f{eqnarray*}{
 *      {a}_{gravity sum,x} &=& {a}_{gravity central,x}
 *                              + \sum\limits_{i=2}^{n} {a}_{gravity J_{n},x}  \\
 *      {a}_{gravity sum,y} &=& {a}_{gravity central,y}
 *                              + \sum\limits_{i=2}^{n} {a}_{gravity J_{n},y}  \\
 *      {a}_{gravity sum,z} &=& {a}_{gravity central,z}
 *                              + \sum\limits_{i=2}^{n} {a}_{gravity J_{n},z}  \\
 * \f}
 * where \f$J_{n}\f$ is the coefficient describing the irregularity of the gravity field of body2
 * due to the nth zonal term. All accelerations are given with respect to an inertial (barycentric)
 * reference frame. NOTE: Currently, this algorithm only computes accelerations up to \f$J_{4}\f$.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to J4-acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param coefficientsOfGravityField Coefficients of zonal terms, describing irregularity of the
 *          gravity field of body2 [-]. The map key contains the degree, and the map value contains
 *          the coefficient value. At the moment, the only degrees supported are 2 (\f$J_{2}\f$),
 *          3 (\f$J_{3}\f$), and 4 (\f$J_{4}\f$) etc. An error is thrown for any other values.
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Total gravitational acceleration exerted on body1 due to given zonal terms [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationZonalSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const std::map< int, double > zonalCoefficientsOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
