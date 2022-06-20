/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CENTRAL_GRAVITY_MODEL_H
#define TUDAT_CENTRAL_GRAVITY_MODEL_H

#include <boost/lambda/lambda.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModelBase.h"

namespace tudat
{
namespace gravitation
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
        const Eigen::Vector3d& positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( )  );

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
        const Eigen::Vector3d& positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( ) );

//! Compute gravitational force.
/*!
 * Computes gravitational force experienced by body1, due to its interaction with body2.
 * The basis for this gravitational force is that both body1 and body2 are point masses,
 * generating Newton's gravitational force:
 * \f[
 *      \bar{F}_{gravity} = -\frac{G * m_{1} * m_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$G\f$ is the universal gravitational constant, \f$m_{1}\f$ is the mass of body1,
 * \f$m_{2}\f$ is the mass of body2, and \f$\bar{r}_{12}\f$ is the relative position vector
 * from body1 to body 2, with respect to an inertial (barycentric) reference frame.
 * \param universalGravitationalParameter Universal gravitational constant [m^3 kg^-1 s^-2].
 * \param massOfBodySubjectToForce Mass of body subject to force (body1) [kg].
 * \param positionOfBodySubjectToForce Position vector of body subject to force (body1) [m].
 * \param massOfBodyExertingForce Mass of body exerting force (body2) [kg].
 * \param positionOfBodyExertingForce Position vector of body exerting force (body2) [m].
 * \return Gravitational force exerted on body1 [N].
 * \sa computeGravitationalAcceleration.
 */
Eigen::Vector3d computeGravitationalForce(
        const double universalGravitationalParameter,
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double massOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce = Eigen::Vector3d::Zero( )  );

//! Compute gravitational force.
/*!
 * Computes gravitational force experienced by body1, due to its interaction with another body
 * (body2) The basis for this gravitational force is that both body1 and body2 are point masses,
 * generating Newton's gravitational force:
 * \f[
 *      \bar{F}_{gravity} = -\frac{m_{1} * \mu_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$m_{1}\f$ is the mass of body1, \f$\mu_{2}\f$ is the gravitational parameter of the
 * central body (body2), and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to
 * the body exerting the force (body2), with respect to an inertial (barycentric) reference frame.
 * \param massOfBodySubjectToForce Mass of body subject to force (body1) [kg].
 * \param positionOfBodySubjectToForce Position vector of body subject to force (body1) [m].
 * \param gravitationalParameterOfBodyExertingForce Gravitational parameter of body exerting
 *          force (body2) [m^3 s^-2].
 * \param positionOfBodyExertingForce Position vector of body exerting force (body2) [m].
 * \return Gravitational force exerted on body1 [N].
 * \sa computeGravitationalAcceleration.
 */
Eigen::Vector3d computeGravitationalForce(
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double gravitationalParameterOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce = Eigen::Vector3d::Zero( )  );

//! Template class for central gravitational acceleration model.
/*!
 * This template class implements a central gravitational acceleration model, i.e., only the
 * central term of the general spherical harmonics expansion.
 * \tparam StateMatrix Data type for state matrix (default = Eigen::Vector3d).
 */
template< typename StateMatrix = Eigen::Vector3d >
class CentralGravitationalAccelerationModel
        : public basic_astrodynamics::AccelerationModel< StateMatrix >,
        public SphericalHarmonicsGravitationalAccelerationModelBase< StateMatrix >
{
private:

    //! Typedef for base class.
    typedef SphericalHarmonicsGravitationalAccelerationModelBase< StateMatrix > Base;

public:

    //! Constructor taking position-functions for bodies, and constant gravitational parameter.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, a constant gravitational parameter, and a pointer to a function
     * returning the position of the body exerting the gravitational acceleration (typically the
     * central body). This constructor uses the Boost::lambda library to create a function
     * on-the-fly that returns the constant gravitational parameter provided. The constructor also
     * updates all the internal members. The position of the body exerting the gravitational
     * acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    CentralGravitationalAccelerationModel(
            const typename Base::StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const typename Base::StateFunction positionOfBodyExertingAccelerationFunction =
            [ ]( StateMatrix& input ){ input = StateMatrix::Zero( ); },
            const bool isMutualAttractionUsed = false )
        : Base( positionOfBodySubjectToAccelerationFunction,
                [ = ]( ){ return aGravitationalParameter; },
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed )
    {
    }

    //! Constructor taking position-functions for bodies, and constant gravitational parameter.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, a constant gravitational parameter, and a pointer to a function
     * returning the position of the body exerting the gravitational acceleration (typically the
     * central body). This constructor uses the Boost::lambda library to create a function
     * on-the-fly that returns the constant gravitational parameter provided. The constructor also
     * updates all the internal members. The position of the body exerting the gravitational
     * acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameterFunction Functioning returning a (constant) gravitational
     *          parameter [m^2 s^-3].
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    CentralGravitationalAccelerationModel(
            const typename Base::StateFunction positionOfBodySubjectToAccelerationFunction,
            const std::function< double( ) > aGravitationalParameterFunction,
            const typename Base::StateFunction positionOfBodyExertingAccelerationFunction =
            [ ]( StateMatrix& input ){ input = StateMatrix::Zero( ); },
            const bool isMutualAttractionUsed = false )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameterFunction,
                positionOfBodyExertingAccelerationFunction,
                isMutualAttractionUsed )
    {
    }

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
            this->currentAcceleration_ = computeGravitationalAcceleration(
                        this->positionOfBodySubjectToAcceleration,
                        this->gravitationalParameter,
                        this->positionOfBodyExertingAcceleration );
        }
    }


protected:
private:
};

//! Typedef for CentralGravitationalAccelerationModel3d.
typedef CentralGravitationalAccelerationModel< > CentralGravitationalAccelerationModel3d;

//! Typedef for shared-pointer to CentralGravitationalAccelerationModel3d.
typedef std::shared_ptr< CentralGravitationalAccelerationModel3d >
CentralGravitationalAccelerationModel3dPointer;

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_CENTRAL_GRAVITY_MODEL_H
