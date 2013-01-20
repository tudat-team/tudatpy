/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      121023    K. Kumar          File created.
 *      121105    K. Kumar          Class simplified, file renamed, and contents merged from other
 *                                  files.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_CENTRAL_J2_J3_GRAVITY_H
#define TUDAT_CENTRAL_J2_J3_GRAVITY_H

#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h"

namespace tudat
{
namespace gravitation
{

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

//! Template class for central + J2 + J3 gravitational acceleration model.
/*!
 * This template class implements a gravitational acceleration model that includes the central, J2,
 * and J3 (unnormalized coefficient of general spherical harmonics expansion) terms.
 * \tparam DataType Data type (default = double).
 * \tparam PositionType Data type for position vector (default = Eigen::Vector3DataType).
 */
template< typename DataType = double, typename PositionType = Eigen::Matrix< DataType, 3, 1 > >
class CentralJ2J3GravitationalAccelerationModel
        : public basic_astrodynamics::AccelerationModel< PositionType >,
        SphericalHarmonicsGravitationalAccelerationModelBase< DataType, PositionType >
{
private:

    //! Typedef for base class.
    typedef SphericalHarmonicsGravitationalAccelerationModelBase< DataType, PositionType > Base;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter, J2- and J3-coefficients, and
     * equatorial radius of the body exerting the acceleration, and a pointer to a function
     * returning the position of the body exerting the gravitational acceleration (typically the
     * central body). This constructor uses the Boost::lambda library to create a function
     * on-the-fly that returns the constant gravitational parameter, J2- and J3-coefficients, and
     * equatorial radius provided. The constructor also updates all the internal members. The
     * position of the body exerting the gravitational acceleration is an optional parameter; the
     * default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aJ2GravityCoefficient A (constant) J2-coefficient.
     * \param aJ3GravityCoefficient A (constant) J3-coefficient.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     */
    CentralJ2J3GravitationalAccelerationModel(
            const typename Base::PositionReturningFunction
            positionOfBodySubjectToAccelerationFunction,
            const DataType aGravitationalParameter,
            const DataType anEquatorialRadius,
            const DataType aJ2GravityCoefficient,
            const DataType aJ3GravityCoefficient,
            const typename Base::PositionReturningFunction
            positionOfBodyExertingAccelerationFunction
            = boost::lambda::constant( PositionType::Zero( ) ) );

    //! Constructor taking functions for position of bodies, and parameters of spherical harmonics
    //! expansion.
    /*!
     * Constructor taking pointer to functions returning the position of the body subject to
     * gravitational acceleration, the gravitational parameter of the body exerting the
     * acceleration (central body), the equatorial radius of the central body, the J2- and
     * J3-coefficients of the spherical harmonics expansion, and the position of the central body.
     * The constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param gravitationalParameterFunction Pointer to function returning gravitational parameter.
     * \param equatorialRadiusFunction Pointer to function returning equatorial radius.
     * \param j2GravityCoefficientFunction Pointer to function returning J2-coefficient of
     *          spherical harmonics expansion.
     * \param j3GravityCoefficientFunction Pointer to function returning J3-coefficient of
     *          spherical harmonics expansion.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     */
    CentralJ2J3GravitationalAccelerationModel(
            const typename Base::PositionReturningFunction
            positionOfBodySubjectToAccelerationFunction,
            const typename Base::DataTypeReturningFunction gravitationalParameterFunction,
            const typename Base::DataTypeReturningFunction equatorialRadiusFunction,
            const typename Base::DataTypeReturningFunction j2GravityCoefficientFunction,
            const typename Base::DataTypeReturningFunction j3GravityCoefficientFunction,
            const typename Base::PositionReturningFunction
            positionOfBodyExertingAccelerationFunction
            = boost::lambda::constant( PositionType::Zero( ) ) );

    //! Get gravitational acceleration.
    /*!
     * Returns the gravitational acceleration computed using the input parameters provided to the
     * class. This function serves as a wrapper for the computeGravitationalAcceleration(),
     * computeGravitationalAccelerationDueToJ2(), and computeGravitationalAccelerationDueToJ3()
     * functions.
     * \return Computed gravitational acceleration vector.
     */
    typename Base::AccelerationType getAcceleration( );

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class
     * members of this class.
     * \return Flag indicating if update was successful or not.
     */
    bool updateMembers( );

protected:

    //! Equatorial radius.
    /*!
     * Equatorial radius of unnormalized spherical harmonics gravity field representation.
     */
    DataType equatorialRadius;

    //! J2 gravity coefficient.
    /*!
     * J2 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    DataType j2GravityCoefficient;

    //! J3 gravity coefficient.
    /*!
     * J3 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    DataType j3GravityCoefficient;

    //! Pointer to function returning equatorial radius.
    /*!
     * Pointer to function that returns the current value of the equatorial radius.
     */
    const typename Base::DataTypeReturningFunction getEquatorialRadius;

    //! Pointer to function returning J2 gravity coefficient.
    /*!
     * Pointer to function that returns the current value of the J2 gravity coefficient.
     */
    const typename Base::DataTypeReturningFunction getJ2GravityCoefficient;

    //! Pointer to function returning J3 gravity coefficient.
    /*!
     * Pointer to function that returns the current value of the J3 gravity coefficient.
     */
    const typename Base::DataTypeReturningFunction getJ3GravityCoefficient;

private:
};

//! Typedef for CentralJ2J3GravitationalAccelerationModel3d.
typedef CentralJ2J3GravitationalAccelerationModel< > CentralJ2J3GravitationalAccelerationModel3d;

//! Typedef for shared-pointer to CentralJ2J3GravitationalAccelerationModel3d.
typedef boost::shared_ptr< CentralJ2J3GravitationalAccelerationModel3d >
CentralJ2J3GravitationalAccelerationModel3dPointer;

// Template class source.
// The code given below is effectively the ".cpp file" for the template class definition, so you
// only need to look at the code below if you are interested in the source implementation.

//! Constructor taking position-functions for bodies, and constant parameters of spherical
//! harmonics expansion.
template< typename DataType, typename PositionType >
CentralJ2J3GravitationalAccelerationModel< DataType, PositionType >::
CentralJ2J3GravitationalAccelerationModel(
        const typename Base::PositionReturningFunction positionOfBodySubjectToAccelerationFunction,
        const DataType aGravitationalParameter,
        const DataType anEquatorialRadius,
        const DataType aJ2GravityCoefficient,
        const DataType aJ3GravityCoefficient,
        const typename Base::PositionReturningFunction positionOfBodyExertingAccelerationFunction )
    : Base( positionOfBodySubjectToAccelerationFunction,
            boost::lambda::constant( aGravitationalParameter ),
            positionOfBodyExertingAccelerationFunction ),
      getEquatorialRadius( boost::lambda::constant( anEquatorialRadius ) ),
      getJ2GravityCoefficient( boost::lambda::constant( aJ2GravityCoefficient ) ),
      getJ3GravityCoefficient( boost::lambda::constant( aJ3GravityCoefficient ) )
{
    this->updateMembers( );
}

//! Constructor taking functions for position of bodies, and parameters of spherical harmonics
//! expansion.
template< typename DataType, typename PositionType >
CentralJ2J3GravitationalAccelerationModel< DataType, PositionType >::
CentralJ2J3GravitationalAccelerationModel(
        const typename Base::PositionReturningFunction positionOfBodySubjectToAccelerationFunction,
        const typename Base::DataTypeReturningFunction gravitationalParameterFunction,
        const typename Base::DataTypeReturningFunction equatorialRadiusFunction,
        const typename Base::DataTypeReturningFunction j2GravityCoefficientFunction,
        const typename Base::DataTypeReturningFunction j3GravityCoefficientFunction,
        const typename Base::PositionReturningFunction positionOfBodyExertingAccelerationFunction )
    : Base( positionOfBodySubjectToAccelerationFunction,
            gravitationalParameterFunction,
            positionOfBodyExertingAccelerationFunction ),
      getEquatorialRadius( equatorialRadiusFunction ),
      getJ2GravityCoefficient( j2GravityCoefficientFunction ),
      getJ3GravityCoefficient( j3GravityCoefficientFunction )
{
    this->updateMembers( );
}

//! Get gravitational acceleration.
template< typename DataType, typename PositionType >
typename CentralJ2J3GravitationalAccelerationModel<
DataType, PositionType >::Base::AccelerationType
CentralJ2J3GravitationalAccelerationModel< DataType, PositionType >::getAcceleration( )
{
    // Sum and return constituent acceleration terms.
    return computeGravitationalAcceleration(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ2(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->j2GravityCoefficient,
                this->equatorialRadius,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ3(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->j3GravityCoefficient,
                this->equatorialRadius,
                this->positionOfBodyExertingAcceleration );
}

//! Update class members.
template< typename DataType, typename PositionType >
bool CentralJ2J3GravitationalAccelerationModel< DataType, PositionType >::updateMembers( )
{
    this->updateBaseMembers( );
    this->equatorialRadius = this->getEquatorialRadius( );
    this->j2GravityCoefficient = this->getJ2GravityCoefficient( );
    this->j3GravityCoefficient = this->getJ3GravityCoefficient( );
    return true;
}

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_CENTRAL_J2_J3_GRAVITY_H
