/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CENTRAL_J2_J3_J4_GRAVITY_MODEL_H
#define TUDAT_CENTRAL_J2_J3_J4_GRAVITY_MODEL_H

#include <map>

#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h"

namespace tudat
{
namespace gravitation
{

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
 * \param equatorialRadiusOfBodyExertingAcceleration Equatorial radius of body exerting
 *          acceleration (body2), in formulation of spherical harmonics expansion [m].
 * \param j4CoefficientOfGravityField J4-coefficient, describing irregularity of the gravity field
 *          of body2 [-].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 due to J4-effect [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationDueToJ4(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const double j4CoefficientOfGravityField,
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
 * \param equatorialRadiusOfBodyExertingAcceleration Equatorial radius of body exerting
 *          acceleration (body2), in formulation of spherical harmonics expansion [m].
 * \param zonalCoefficientsOfGravityField Coefficients of zonal terms, describing irregularity of the
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
        const double equatorialRadiusOfBodyExertingAcceleration,
        const std::map< int, double > zonalCoefficientsOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Central + J2 + J3 + J4 gravitational acceleration model class.
/*!
 * This class implements a gravitational acceleration model that includes the central, J2, J3, and
 * J4 (unnormalized coefficient of general spherical harmonics expansion) terms.
 */
class CentralJ2J3J4GravitationalAccelerationModel
        : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >,
        public SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d >
{
private:

    //! Typedef for base class.
    typedef SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > Base;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter, J2-, J3-, and J4-coefficients,
     * and equatorial radius of the body exerting the acceleration, and a pointer to a function
     * returning the position of the body exerting the gravitational acceleration (typically the
     * central body). This constructor uses the Boost::lambda library to create a function
     * on-the-fly that returns the constant gravitational parameter, J2-, J3-, and J4-coefficients,
     * and equatorial radius provided. The constructor also updates all the internal members. The
     * position of the body exerting the gravitational acceleration is an optional parameter; the
     * default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aJ2GravityCoefficient A (constant) J2-coefficient.
     * \param aJ3GravityCoefficient A (constant) J3-coefficient.
     * \param aJ4GravityCoefficient A (constant) J4-coefficient.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     */
    CentralJ2J3J4GravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double anEquatorialRadius,
            const double aJ2GravityCoefficient,
            const double aJ3GravityCoefficient,
            const double aJ4GravityCoefficient,
            const StateFunction positionOfBodyExertingAccelerationFunction
            = boost::lambda::constant( Eigen::Vector3d::Zero( ) ) )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameter,
                positionOfBodyExertingAccelerationFunction, false ),
          equatorialRadius( anEquatorialRadius ),
          j2GravityCoefficient( aJ2GravityCoefficient ),
          j3GravityCoefficient( aJ3GravityCoefficient ),
          j4GravityCoefficient( aJ4GravityCoefficient )
    {
        this->updateMembers( );
    }

    //! Get gravitational acceleration.
    /*!
     * Returns the gravitational acceleration computed using the input parameters provided to the
     * class. This function serves as a wrapper for the computeGravitationalAcceleration(),
     * computeGravitationalAccelerationDueToJ2(), and computeGravitationalAccelerationDueToJ3()
     * functions.
     * \return Computed gravitational acceleration vector.
     */
    Eigen::Vector3d getAcceleration( );

    //! Update members.
    /*!
     * Updates class members relevant for computing the central gravitational acceleration. In this
     * case the function simply updates the members in the base class.
     * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->updateBaseMembers( );
        }
    }

protected:

    //! Equatorial radius [m].
    /*!
     * Equatorial radius of unnormalized spherical harmonics gravity field representation [m].
     */
    const double equatorialRadius;

    //! J2 gravity coefficient.
    /*!
     * J2 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    const double j2GravityCoefficient;

    //! J3 gravity coefficient.
    /*!
     * J3 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    const double j3GravityCoefficient;

    //! J4 gravity coefficient.
    /*!
     * J4 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    const double j4GravityCoefficient;

private:
};

//! Typedef for shared-pointer to CentralJ2J3J4GravitationalAccelerationModel.
typedef boost::shared_ptr< CentralJ2J3J4GravitationalAccelerationModel >
CentralJ2J3J4GravitationalAccelerationModelPointer;

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_CENTRAL_J2_J3_J4_GRAVITY_MODEL_H
