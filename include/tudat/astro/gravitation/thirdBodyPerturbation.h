/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K.F. astro I, Delft University of Technology, 2010.
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *
 */

#ifndef TUDAT_THIRD_BODY_PERTURBATION_H
#define TUDAT_THIRD_BODY_PERTURBATION_H

#include <Eigen/Core>

#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"

namespace tudat
{
namespace gravitation
{

//! Compute perturbing acceleration by third body.
/*!
 * Computes the perturbing acceleration on a point mass in orbit about a central body (point mass),
 * caused by a third body (point mass). This acceleration is expressed with respect to a
 * pseudo-inertial reference frame centered at the central body body, i.e, a non-rotating frame
 * aligned with an inertial reference frame. This acceleration may be summed with other third body
 * accelerations in order to get the total gravitational perturbation acceleration of N bodies
 * (Wakker, 2010; Montenbruck & Gill, 2005).
 * \param gravitationalParameterOfPerturbingBody The gravitational parameter of the perturbing body,
 *      i.e., the third body.                                                             [m^3/s^2]
 * \param positionOfPerturbingBody The position of the third body, the body that
 *      causes the perturbations, in Cartesian coordinates.
 *          positionOfPerturbingBody( 0 ) = x-position                                          [m]
 *          positionOfPerturbingBody( 1 ) = y-position                                          [m]
 *          positionOfPerturbingBody( 2 ) = z-position                                          [m]
 * \param positionOfAffectedBody The position of the second body, the body that
 *      experiences the perturbation, in Cartesian coordinates.
 *          positionOfAffectedBody( 0 ) = x-position                                            [m]
 *          positionOfAffectedBody( 1 ) = y-position                                            [m]
 *          positionOfAffectedBody( 2 ) = z-position                                            [m]
 * \param positionOfCentralBody The position of the central body, the body w.r.t. which the
 *      acceleration is calculated the perturbation, in Cartesian coordinates (default=origin).
 *          positionOfAffectedBody( 0 ) = x-position                                            [m]
 *          positionOfAffectedBody( 1 ) = y-position                                            [m]
 *          positionOfAffectedBody( 2 ) = z-position                                            [m]
 * \return perturbingAcceleration The perturbing acceleration in Cartesian components, that gives
 *      the difference of the total acceleration with the central two-body acceleration.
 *          perturbingAcceleration( 0 ) = x-acceleration                                    [m/s^2]
 *          perturbingAcceleration( 1 ) = y-acceleration                                    [m/s^2]
 *          perturbingAcceleration( 2 ) = z-acceleration                                    [m/s^2]
 */
Eigen::Vector3d computeThirdBodyPerturbingAcceleration(
        const double gravitationalParameterOfPerturbingBody,
        const Eigen::Vector3d& positionOfPerturbingBody,
        const Eigen::Vector3d& positionOfAffectedBody,
        const Eigen::Vector3d& positionOfCentralBody = Eigen::Vector3d::Zero( ) );

//! Class for calculating third-body (gravitational) accelerations.
/*!
 *  Class for calculating third-body (gravitational accelerations),
 *  i.e. the gravitational acceleration on a body A, expressed in a frame fixed
 *  on a body B, due to a gravitating body C. The acceleration is calculated by
 *  subtracting the acceleration of the central body (B) due to body C
 *  from the direct acceleration of body C on body A.
 *  \tparam The gravitational acceleration class for which the third body acceleration is
 *  calculated (CentralGravitationalAccelerationModel,
 *  SphericalHarmonicsGravitationalAccelerationModel, ...)
 */
template< typename DirectAccelerationModelType,
          typename std::enable_if< is_direct_gravity_acceleration< DirectAccelerationModelType >::value, int >::type = 0 >
class ThirdBodyAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor for third body acceleration
    /*!
     *  Constructor, sets the two acceleration models (one direct, one on central body)
     *  \param accelerationModelForBodyUndergoingAcceleration Direct acceleration model on
     *  body undergoing acceleration (i.e. as expressed in an inertial frame)
     *  \param accelerationModelForCentralBody Acceleration model on central body
     *  (i.e. the body in a frame centered on which the third body acceleration is expressed)
     *  \param centralBodyName Name of the central body w.r.t. which the acceleration is computed.
     */
    ThirdBodyAcceleration(
            const std::shared_ptr< DirectAccelerationModelType >
            accelerationModelForBodyUndergoingAcceleration,
            const std::shared_ptr< DirectAccelerationModelType >
            accelerationModelForCentralBody,
            const std::string centralBodyName ):
        accelerationModelForBodyUndergoingAcceleration_(
            accelerationModelForBodyUndergoingAcceleration ),
        accelerationModelForCentralBody_( accelerationModelForCentralBody ),
        centralBodyName_( centralBodyName ){ }

    //! Update member variables to current state.
    /*!
     *  Update member variables to current state.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Update two constituent acceleration models.
            accelerationModelForBodyUndergoingAcceleration_->updateMembers( currentTime );
            accelerationModelForCentralBody_->updateMembers( currentTime );
            currentAcceleration_ = accelerationModelForBodyUndergoingAcceleration_->getAcceleration( ) -
                    accelerationModelForCentralBody_->getAcceleration( );
        }
    }

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    void resetCurrentTime( )
    {
        accelerationModelForBodyUndergoingAcceleration_->resetCurrentTime( );
        accelerationModelForCentralBody_->resetCurrentTime( );
    }

    //! Function to return the direct acceleration model on body undergoing acceleration.
    /*!
     *  Function to return the direct acceleration model on body undergoing acceleration.
     *  \return Direct acceleration model on body undergoing acceleration.
     */
    std::shared_ptr< DirectAccelerationModelType >
        getAccelerationModelForBodyUndergoingAcceleration( )
    {
        return accelerationModelForBodyUndergoingAcceleration_;
    }

    //! Function to return the acceleration model on central body
    /*!
     *  Function to return the acceleration model on central body
     *  \return Acceleration model on central body
     */
    std::shared_ptr< DirectAccelerationModelType >
        getAccelerationModelForCentralBody( )
    {
        return accelerationModelForCentralBody_;
    }

    //! Function to return the name of the central body w.r.t. which the acceleration is computed.
    /*!
     *  Function to return the name of the central body w.r.t. which the acceleration is computed.
     *  \return Name of the central body w.r.t. which the acceleration is computed.
     */
     std::string getCentralBodyName( )
     {
         return centralBodyName_;
     }

private:

    //! Direct acceleration model on body undergoing acceleration.
    /*!
     *  Direct acceleration model on body undergoing acceleration
     *  (i.e. as expressed in an inertial frame)
     */
    std::shared_ptr< DirectAccelerationModelType >
        accelerationModelForBodyUndergoingAcceleration_;

    //! Acceleration model on central body
    /*!
     *  Acceleration model on central body (i.e. the body in a frame centered on which the third
     *  body acceleration is expressed)
     */
    std::shared_ptr< DirectAccelerationModelType > accelerationModelForCentralBody_;

    //! Name of the central body w.r.t. which the acceleration is computed.
     std::string centralBodyName_;
};

//! Typedef for third body central gravity acceleration.
typedef ThirdBodyAcceleration< CentralGravitationalAccelerationModel3d >
ThirdBodyCentralGravityAcceleration;

//! Typedef for third body spherical harmonic gravity acceleration.
typedef ThirdBodyAcceleration< SphericalHarmonicsGravitationalAccelerationModel >
ThirdBodySphericalHarmonicsGravitationalAccelerationModel;

//! Typedef for third body mutual spherical harmonic gravity acceleration.
typedef ThirdBodyAcceleration< MutualSphericalHarmonicsGravitationalAccelerationModel >
ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel;


} // namespace gravitation

} // namespace tudat

#endif // TUDAT_THIRD_BODY_PERTURBATION_H
