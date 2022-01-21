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

#ifndef TUDAT_SOLAR_SAIL_ACCELERATION_H
#define TUDAT_SOLAR_SAIL_ACCELERATION_H

#include <functional>
#include <memory>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/electromagnetism/solarSailForce.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute solar sail acceleration with a non-ideal reflective model.
/*!
 * Compute solar sail acceleration with a non-ideal reflective model. The acceleration is computed from the force derived from
 * the solar sailing radiation pressure model, which is itself a function of the source power and geometry, as well as of the
 * solar sailing model.
 * \param frontEmissivityCoefficient Emissivity coefficient of the front of the sail.                    [-]
 * \param backEmissivityCoefficient Emissivity coefficient of the back of the sail.                      [-]
 * \param frontLambertianCoefficient Lambertian coefficient of the front of the sail.                    [-]
 * \param backLambertianCoefficient Lambertian coefficient of the back of the sail.                      [-]
 * \param reflectivityCoefficient Reflectivity coefficient of the sail.                                  [-]
 * \param specularReflectionCoefficient Specular reflection coefficient of the sail.                     [-]
 * \param normalisedVectorToSource Normalised vector pointing from target to source.
 *          N.B: this must be a unit vector!
 *          To compute the unit vector based on a given position vector, you can use the
 *          .normalize() or .normalized() member functions of an Eigen::Vector3d object.                 [-]
 * \param normalisedVelocityVector Normalised velocity vector of the spacecraft w.r.t. central body.     [-]
 * \param radiationPressure Radiation pressure at target.                                            [N/m^2]
 * \param area Area on which radiation pressure is assumed to act.                                     [m^2]
 * \param coneAngle Sail cone angle.                                                                   [rad]
 * \param clockAngle Sail clock angle.                                                                 [rad]
 * \param mass Mass of accelerated body.                                                                [kg]
 * \return Acceleration due to radiation pressure.                                                   [m/s^2]
 */
Eigen::Vector3d computeSolarSailAcceleration(
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const Eigen::Vector3d& normalisedVectorToSource,
        const Eigen::Vector3d& normalisedVelocityVector,
        const double radiationPressure,
        const double area,
        const double coneAngle,
        const double clockAngle,
        const double mass );

//! Solar sail acceleration model class.
/*!
 * Class for calculating the solar sail acceleration using a non-ideally reflective model.
 */
class SolarSailAcceleration: public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef std::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef std::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

    //! Typedef for Eigen::Vector6d-returning function.
    typedef std::function< Eigen::VectorXd( ) > VectorXdReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing acceleration.
     * \param acceleratedBodyVelocityFunction Function returning velocity of body undergoing acceleration.
     * \param centralBodyVelocityFunction Function returning velocity of central body.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param coneAngleFunction Function returning current cone angle.
     * \param clockAngleFunction Function returning current clock angle.
     * \param frontEmissivityCoefficientFunction Function returning current front emissivity coefficient.
     * \param backEmissivityCoefficientFunction Function returning current back emissivity coefficient.
     * \param frontLambertianCoefficientFunction Function returning current front Lambertian coefficient.
     * \param backLambertianCoefficientFunction Function returning current back Lambertian coefficient.
     * \param reflectivityCoefficientFunction Function returning current reflectivity coefficient.
     * \param specularReflectionCoefficientFunction Function returning current specular reflection coefficient.
     * \param areaFunction Function returning current area assumed to undergo radiation pressure.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     */
    SolarSailAcceleration(
            const Vector3dReturningFunction sourcePositionFunction,
            const Vector3dReturningFunction acceleratedBodyPositionFunction,
            const Vector3dReturningFunction acceleratedBodyVelocityFunction,
            const Vector3dReturningFunction centralBodyVelocityFunction,
            const DoubleReturningFunction radiationPressureFunction,
            const DoubleReturningFunction coneAngleFunction,
            const DoubleReturningFunction clockAngleFunction,
            const DoubleReturningFunction frontEmissivityCoefficientFunction,
            const DoubleReturningFunction backEmissivityCoefficientFunction,
            const DoubleReturningFunction frontLambertianCoefficientFunction,
            const DoubleReturningFunction backLambertianCoefficientFunction,
            const DoubleReturningFunction reflectivityCoefficientFunction,
            const DoubleReturningFunction specularReflectionCoefficientFunction,
            const DoubleReturningFunction areaFunction,
            const DoubleReturningFunction massFunction )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          centralBodyVelocityFunction_( centralBodyVelocityFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          coneAngleFunction_( coneAngleFunction ),
          clockAngleFunction_( clockAngleFunction ),
          frontEmissivityCoefficientFunction_( frontEmissivityCoefficientFunction ),
          backEmissivityCoefficientFunction_( backEmissivityCoefficientFunction ),
          frontLambertianCoefficientFunction_( frontLambertianCoefficientFunction ),
          backLambertianCoefficientFunction_( backLambertianCoefficientFunction ),
          reflectivityCoefficientFunction_( reflectivityCoefficientFunction ),
          specularReflectionCoefficientFunction_( specularReflectionCoefficientFunction ),
          areaFunction_( areaFunction ),
          massFunction_( massFunction )
    {
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position vectors and radiation pressure and
     * constant values for other parameters.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing acceleration.
     * \param acceleratedBodyVelocityFunction Function returning velocity vector of body undergoing acceleration.
     * \param centralBodyVelocityFunction Function returning velocity vector of central body.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param coneAngleFunction Function returning current sail cone angle.
     * \param clockAngleFunction Function returning current sail clock angle.
     * \param frontEmissivityCoefficient Constant front emissivity coefficient.
     * \param backEmissivityCoefficient Constant back emissivity coefficient.
     * \param frontLambertianCoefficient Constant front Lambertian coefficient.
     * \param backLambertianCoefficient Constant back Lambertian coefficient.
     * \param reflectivityCoefficient Constant reflectivity coefficient.
     * \param specularReflectionCoefficient Constant specular reflection coefficient
     * \param area Constant area assumed to undergo radiation pressure.
     * \param mass Constant mass of body undergoing acceleration.
     */
    SolarSailAcceleration(
            const Vector3dReturningFunction sourcePositionFunction,
            const Vector3dReturningFunction acceleratedBodyPositionFunction,
            const Vector3dReturningFunction acceleratedBodyVelocityFunction,
            const Vector3dReturningFunction centralBodyVelocityFunction,
            const DoubleReturningFunction radiationPressureFunction,
            const DoubleReturningFunction coneAngleFunction,
            const DoubleReturningFunction clockAngleFunction,
            const double frontEmissivityCoefficient,
            const double backEmissivityCoefficient,
            const double frontLambertianCoefficient,
            const double backLambertianCoefficient,
            const double reflectivityCoefficient,
            const double specularReflectionCoefficient,
            const double area,
            const double mass
            )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          centralBodyVelocityFunction_( centralBodyVelocityFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          coneAngleFunction_( coneAngleFunction ),
          clockAngleFunction_( clockAngleFunction ),
          frontEmissivityCoefficientFunction_( [ = ]( ){ return frontEmissivityCoefficient;} ),
          backEmissivityCoefficientFunction_( [ = ]( ){ return backEmissivityCoefficient;} ),
          frontLambertianCoefficientFunction_( [ = ]( ){ return frontLambertianCoefficient;} ),
          backLambertianCoefficientFunction_( [ = ]( ){ return backLambertianCoefficient;} ),
          reflectivityCoefficientFunction_( [ = ]( ){ return reflectivityCoefficient;} ),
          specularReflectionCoefficientFunction_( [ = ]( ){ return specularReflectionCoefficient;} ),
          areaFunction_( [ = ]( ){ return area;} ),
          massFunction_( [ = ]( ){ return mass;} )
    {
    }


    //! Constructor for setting up the acceleration model, with input SolarSailingRadiationPressureInterface (and massFunction)
    /*!
     *  Constructor for setting up the acceleration model, with input SolarSailingRadiationPressureInterface (and massFunction).
     *  The massFunction is only used to convert force to acceleration. All solar sailing radiation pressure properties are taken from the
     *  SolarSailingRadiationPressureInterface object.
     *  \param radiationPressureInterface Object in which solar sailing radiation pressure properties of accelerated body due to radiation
     *  from body causing acceleration is stored.
     *  \param massFunction Function returning the current mass of the body being accelerated.
     */
    SolarSailAcceleration( const std::shared_ptr< SolarSailingRadiationPressureInterface > radiationPressureInterface,
                           const std::function< double( ) > massFunction ):
        sourcePositionFunction_( std::bind( &RadiationPressureInterface::getCurrentSolarVector, radiationPressureInterface ) ),
        acceleratedBodyPositionFunction_( radiationPressureInterface->getTargetPositionFunction() ),
        acceleratedBodyVelocityFunction_( radiationPressureInterface->getTargetVelocityFunction() ),
        centralBodyVelocityFunction_( radiationPressureInterface->getCentralBodyVelocity() ),
        radiationPressureFunction_( std::bind( &RadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface ) ),
        coneAngleFunction_( std::bind( &SolarSailingRadiationPressureInterface::getCurrentConeAngle, radiationPressureInterface ) ),
        clockAngleFunction_( std::bind( &SolarSailingRadiationPressureInterface::getCurrentClockAngle, radiationPressureInterface ) ),
        frontEmissivityCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getFrontEmissivityCoefficient, radiationPressureInterface ) ),
        backEmissivityCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getBackEmissivityCoefficient, radiationPressureInterface ) ),
        frontLambertianCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getFrontLambertianCoefficient, radiationPressureInterface ) ),
        backLambertianCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getBackLambertianCoefficient, radiationPressureInterface ) ),
        reflectivityCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getReflectivityCoefficient, radiationPressureInterface ) ),
        specularReflectionCoefficientFunction_( std::bind( &SolarSailingRadiationPressureInterface::getSpecularReflectionCoefficient, radiationPressureInterface ) ),
        areaFunction_( std::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ) ),
        massFunction_( massFunction )
    {
    }


    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration( ) function.
     *
     * \param currentTime Time at which acceleration model is to be updated [s].
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentFrontEmissivityCoefficient_ = frontEmissivityCoefficientFunction_( );
            currentBackEmissivityCoefficient_ = backEmissivityCoefficientFunction_( );
            currentFrontLambertianCoefficient_ = frontLambertianCoefficientFunction_( );
            currentBackLambertianCoefficient_ = backLambertianCoefficientFunction_( );
            currentReflectivityCoefficient_ = reflectivityCoefficientFunction_( );
            currentSpecularReflectionCoefficient_ = specularReflectionCoefficientFunction_( );
            currentNormalizedVectorToSource_ = ( sourcePositionFunction_( ) - acceleratedBodyPositionFunction_( ) ).normalized( );
            currentNormalizedVelocityVector_= ( acceleratedBodyVelocityFunction_( ) - centralBodyVelocityFunction_( ) ).normalized( );
            currentRadiationPressure_ = radiationPressureFunction_( );
            currentArea_ = areaFunction_( );
            currentConeAngle_ = coneAngleFunction_( );
            currentClockAngle_ = clockAngleFunction_( );
            currentMass_ = massFunction_( );
            currentDistanceToSource_ = ( sourcePositionFunction_( ) - acceleratedBodyPositionFunction_( ) ).norm();
            currentVelocityWrtSource_ = ( acceleratedBodyVelocityFunction_( ) - centralBodyVelocityFunction_( ) ).norm( );
            currentAcceleration_ = computeSolarSailAcceleration(
                        currentFrontEmissivityCoefficient_, currentBackEmissivityCoefficient_,
                        currentFrontLambertianCoefficient_, currentBackLambertianCoefficient_,
                        currentReflectivityCoefficient_, currentSpecularReflectionCoefficient_,
                        currentNormalizedVectorToSource_, currentNormalizedVelocityVector_,
                        currentRadiationPressure_, currentArea_,
                        currentConeAngle_, currentClockAngle_, currentMass_ );
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

    //! Returns the current normalized vector from the accelerated body to the source body
    /*!
     *  Returns the current normalized vector from the accelerated body to the source body, as set by the last call to the updateMembers function
     *  \return Crrent normalized vector from the accelerated body to the source body [-].
     */
    Eigen::Vector3d getCurrentVectorToSource( )
    {
        return currentNormalizedVectorToSource_;
    }

    //! Returns the current normalized velocity vector vector of the accelerated body w.r.t. the central body.
    /*!
     *  Returns the current normalized velocity vector vector of the accelerated body w.r.t. the central body.
     *  \return Current normalized velocity vector of the accelerated body w.r.t. the central body [-].
     */
    Eigen::Vector3d getCurrentVelocityVector( )
    {
        return currentNormalizedVelocityVector_;
    }

    //! Returns the current distance from the accelerated body to the source body
    /*!
     *  Returns the current distance from the accelerated body to the source body, as set by the last call to the updateMembers function
     *  (i.e. norm of the current vector from accelerated body to the source)
     *  \return Current distance from the accelerated body to the source body [m].
     */
    double getCurrentDistanceToSource( )
    {
        return currentDistanceToSource_;
    }

    //! Returns the current velocity of the accelerated body w.r.t. the central body
    /*!
     *  Returns the current distance of the accelerated body w.r.t. the central body, as set by the last call to the updateMembers function
     *  (i.e. norm of the current velocity vector of the accelerated body w.r.t. central body)
     *  \return Current velocity of the accelerated body w.r.t. the central body [m/s].
     */
    double getCurrentVelocityWrtSource( )
    {
        return currentVelocityWrtSource_;
    }

    //! Returns the current radiation pressure at the accelerated body
    /*!
     *  Returns the current radiation pressure at the accelerated body, as set by the last call to the updateMembers function
     *  (i.e. incident flux, in W/m^2, divided by speed of light)
     *  \return Current radiation pressure at the accelerated body [N/m^2].
     */
    double getCurrentRadiationPressure( )
    {
        return currentRadiationPressure_;
    }

    //! Returns the current mass of the accelerated body
    /*!
     *  Returns the current mass of the accelerated body, as set by the last call to the updateMembers function
     *  \return Current mass of the accelerated body [kg].
     */
    double getCurrentMass( )
    {
        return currentMass_;
    }

    //! Returns the current area of the accelerated body
    /*!
     *  Returns the current area of the accelerated body, as set by the last call to the updateMembers function
     *  \return Current area of the accelerated body [m^2].
     */
    double getCurrentArea( )
    {
        return currentArea_;
    }



private:

    //! Function pointer returning position of source.
    /*!
     * Function pointer returning position of source (3D vector).
     */
    Vector3dReturningFunction sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning velocity of accelerated body.
    /*!
     * Function pointer returning velocity of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyVelocityFunction_;

    //! Function pointer returning velocity of central body.
    /*!
     * Function pointer returning velocity of central body (3D vector).
     */
    const Vector3dReturningFunction centralBodyVelocityFunction_;

    //! Function pointer returning radiation pressure.
    /*!
     * Function pointer returning radiation pressure [N/m^{2}].
     */
    const DoubleReturningFunction radiationPressureFunction_;

    //! Function pointer returning cone angle.
    /*!
     * Function pointer returning cone angle [rad].
     */
    const DoubleReturningFunction coneAngleFunction_;


    //! Function pointer returning clock angle.
    /*!
     * Function pointer returning clock angle [rad].
     */
    const DoubleReturningFunction clockAngleFunction_;


    //! Function pointer returning front emissivity coefficient.
    /*!
     * Function pointer returning front emissivity coefficient [-].
     */
    const DoubleReturningFunction frontEmissivityCoefficientFunction_;


    //! Function pointer returning back emissivity coefficient.
    /*!
     * Function pointer returning back emissivity coefficient [-].
     */
    const DoubleReturningFunction backEmissivityCoefficientFunction_;


    //! Function pointer returning front Lambertian coefficient.
    /*!
     * Function pointer returning front Lambertian coefficient [-].
     */
    const DoubleReturningFunction frontLambertianCoefficientFunction_;


    //! Function pointer returning back Lambertian coefficient.
    /*!
     * Function pointer returning back Lambertian coefficient [-].
     */
    const DoubleReturningFunction backLambertianCoefficientFunction_;


    //! Function pointer returning reflectivity coefficient.
    /*!
     * Function pointer returning reflectivity coefficient [-].
     */
    const DoubleReturningFunction reflectivityCoefficientFunction_;


    //! Function pointer returning specular reflection coefficient.
    /*!
     * Function pointer returning specular reflection coefficient [-].
     */
    const DoubleReturningFunction specularReflectionCoefficientFunction_;


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


    //! Current front emissivity coefficient.
    /*!
     * Current front emissivity coefficient [-].
     */
    double currentFrontEmissivityCoefficient_;

    //! Current back emissivity coefficient.
    /*!
     * Current back emissivity coefficient [-].
     */
    double currentBackEmissivityCoefficient_;

    //! Current front Lambertian coefficient.
    /*!
     * Current front Lambertian coefficient [-].
     */
    double currentFrontLambertianCoefficient_;

    //! Current back Lambertian coefficient.
    /*!
     * Current back Lambertian coefficient [-].
     */
    double currentBackLambertianCoefficient_;

    //! Current reflectivity coefficient.
    /*!
     * Current reflectivity coefficient [-].
     */
    double currentReflectivityCoefficient_;

    //! Current specular reflection coefficient.
    /*!
     * Current specular reflection coefficient [-].
     */
    double currentSpecularReflectionCoefficient_;

    //! Current normalised vector from accelerated body to source.
    /*!
     * Current normalised vector from accelerated body to source (3D vector) [-].
     */
    Eigen::Vector3d currentNormalizedVectorToSource_;

    //! Current normalised velocity unit vector of the propagated body w.r.t. central body.
    /*!
     * Current normalised velocity unit vector of the propagated body w.r.t. central body (3D vector) [-].
     */
    Eigen::Vector3d currentNormalizedVelocityVector_;

    //! Current radiation pressure.
    /*!
     * Current radiation pressure [N/m^{2}].
     */
    double currentRadiationPressure_;

    //! Current area on which radiation pressure is acting.
    /*!
     * Current area on which radiation pressure is acting [m^{2}].
     */
    double currentArea_;

    //! Current cone angle.
    /*!
     * Current cone angle [rad].
     */
    double currentConeAngle_;

    //! Current clock angle
    /*!
     * Current clock angle [rad].
     */
    double currentClockAngle_;

    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;

    //! Current distance between accelerated body and source body.
    /*!
     * Current distance between accelerated body and source body [m].
     */
    double currentDistanceToSource_;

    //! Current velocity of the accelerated body w.r.t. the central body.
    /*!
     * Current velocity of the accelerated body w.r.t. the central body [m/s].
     */
    double currentVelocityWrtSource_;

};

//! Typedef for shared-pointer to SolarSailAcceleration.
typedef std::shared_ptr< SolarSailAcceleration > SolarSailAccelerationPointer;

} // namespace electromagnetism
} // namespace tudat

#endif // TUDAT_SOLAR_SAIL_ACCELERATION_H
