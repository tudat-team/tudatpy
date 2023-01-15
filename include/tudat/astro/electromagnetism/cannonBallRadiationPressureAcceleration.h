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
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 */
//RP-OLD

#ifndef TUDAT_CANNON_BALL_RADIATION_PRESSURE_ACCELERATION_H
#define TUDAT_CANNON_BALL_RADIATION_PRESSURE_ACCELERATION_H 

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/electromagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute radiation pressure acceleration using a cannon-ball model.
/*!
 * Computes radiation pressure acceleration using a cannon-ball model, i.e. assuming force to be in
 * opposite direction of the vector to the source. This function is essentially a wrapper for the
 * function that computes the force.
 * opposite direction of the vector to the source.
 * \param radiationPressure Radiation pressure at target. N.B: the usual way of computing the
 *          radiation pressure at the target, in case the source is the Sun, is to take the
 *          radiation presssure at 1 AU and scale it based on the distance from the Sun     [N/m^2]
 * \param vectorToSource Vector pointing from target to source. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.                                                                             [-]
 * \param area Area on which radiation pressure is assumed to act.                            [m^2]
 * \param radiationPressureCoefficient Coefficient to scale effective force. Equal to 1 +
 *          emissivitty, assuming no diffuse reflection.                                        [-]
 * \param mass Mass of accelerated body.                                                       [kg]
 * \return Acceleration due to radiation pressure.                                          [m/s^2]
 * \sa computeCannonBallRadiationPressureForce().
 */
Eigen::Vector3d computeCannonBallRadiationPressureAcceleration(
        const double radiationPressure,
        const Eigen::Vector3d& vectorToSource,
        const double area,
        const double radiationPressureCoefficient,
        const double mass );

//! Cannon-ball Radiation pressure acceleration model class.
/*!
 * Class that can be used to compute the radiation pressure using a cannon-ball model, i.e.,
 * assuming force to be in opposite direction of the vector to the source.
 */
class CannonBallRadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef std::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef std::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          radiation pressure acceleration.
     * \param radiationPressureFunction Function returning current radiation pressure.     
     * \param radiationPressureCoefficientFunction Function returning current radiation pressure
     *          coefficient.
     * \param areaFunction Function returning current area assumed to undergo radiation pressure.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     */
    CannonBallRadiationPressureAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            DoubleReturningFunction radiationPressureFunction,
            DoubleReturningFunction radiationPressureCoefficientFunction,
            DoubleReturningFunction areaFunction,
            DoubleReturningFunction massFunction )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          radiationPressureCoefficientFunction_( radiationPressureCoefficientFunction ),
          areaFunction_( areaFunction ),
          massFunction_( massFunction )
    {
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position vectors and radiation pressure and
     * constant values for other parameters.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *              radiation pressure acceleration.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param radiationPressureCoefficient Constant radiation pressure coefficient.
     * \param area Constant area assumed to undergo radiation pressure.
     * \param mass Constant mass of body undergoing acceleration.
     */
    CannonBallRadiationPressureAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            DoubleReturningFunction radiationPressureFunction,
            const double radiationPressureCoefficient,
            const double area,
            const double mass )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          radiationPressureCoefficientFunction_(
              [ = ]( ){ return radiationPressureCoefficient; } ),
          areaFunction_( [ = ]( ){ return area; } ),
          massFunction_( [ = ]( ){ return mass; } )
    {
    }

    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration( ) function.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentVectorToSource_ = ( sourcePositionFunction_( )
                                       - acceleratedBodyPositionFunction_( ) ).normalized( );
            currentRadiationPressure_ = radiationPressureFunction_( );
            currentRadiationPressureCoefficient_ = radiationPressureCoefficientFunction_( );
            currentArea_ = areaFunction_( );
            currentMass_ = massFunction_( );
            currentAcceleration_ = computeCannonBallRadiationPressureAcceleration(
                        currentRadiationPressure_, currentVectorToSource_, currentArea_,
                        currentRadiationPressureCoefficient_, currentMass_ );
        }
    }

    //! Function to retrieve the function pointer returning mass of accelerated body.
    /*!
     * Function to retrieve the function pointer returning mass of accelerated body.
     * \return Function pointer returning mass of accelerated body.
     */
    std::function< double( ) > getMassFunction( )
    {
        return massFunction_;
    }

private:

    //! Function pointer returning position of source.
    /*!
     * Function pointer returning position of source (3D vector).
     */
    const Vector3dReturningFunction sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning radiation pressure.
    /*!
     * Function pointer returning radiation pressure [N/m^{2}].
     */
    const DoubleReturningFunction radiationPressureFunction_;

    //! Function pointer returning radiation pressure coefficient.
    /*!
     * Function pointer returning radiation pressure coefficient [-].
     */
    const DoubleReturningFunction radiationPressureCoefficientFunction_;

    //! Function pointer returning area on which radiation pressure is acting.
    /*!
     * Function pointer returning area on which radiation pressure is acting [m^{2}].
     */
    const DoubleReturningFunction areaFunction_;

    //! Function pointer returning mass of accelerated body.
    /*!
     * Function pointer returning mass of accelerated body [kg].
     */
    const DoubleReturningFunction massFunction_;

    //! Current vector from accelerated body to source.
    /*!
     * Current vector from accelerated body to source (3D vector).
     */
    Eigen::Vector3d currentVectorToSource_;

    //! Current radiation pressure.
    /*!
     * Current radiation pressure [N/m^{2}].
     */
    double currentRadiationPressure_;

    //! Current radiation pressure coefficient.
    /*!
     *  Current radiation pressure coefficient [-].
     */
    double currentRadiationPressureCoefficient_;

    //! Current area on which radiation pressure is acting.
    /*!
     * Current area on which radiation pressure is acting [m^{2}].
     */
    double currentArea_;

    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;
};

//! Typedef for shared-pointer to CannonBallRadiationPressureAcceleration.
typedef std::shared_ptr< CannonBallRadiationPressureAcceleration > CannonBallRadiationPressurePointer;

} // namespace electromagnetism
} // namespace tudat

#endif // TUDAT_CANNON_BALL_RADIATION_PRESSURE_ACCELERATION_H 
