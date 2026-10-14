/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#ifndef TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_MODEL_H
#define TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_MODEL_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityField.h"

namespace tudat
{

namespace gravitation
{

/*! Acceleration caused by a gravity field described by ellipsoidal harmonics.
 *
 * The acceleration is the gradient of the potential, computed in the body-fixed frame of the
 * body exerting it and then rotated into the frame the propagation runs in.
 */
class EllipsoidalHarmonicsGravitationalAccelerationModel : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
protected:
    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector3d& ) > StateFunction;

public:
    /*! Constructor.
     *
     * @param positionOfBodySubjectToAccelerationFunction Function returning the position of the
     *      body being accelerated.
     * @param gravityField Gravity field of the body exerting the acceleration.
     * @param gravitationalParameterFunction Function returning the gravitational parameter to
     *      use. This is the parameter of the body exerting the acceleration, or the sum of both
     *      bodies' parameters when the frame is fixed to the central body.
     * @param positionOfBodyExertingAccelerationFunction Function returning the position of the
     *      body exerting the acceleration; the origin by default.
     * @param rotationFromBodyFixedToIntegrationFrameFunction Function returning the rotation
     *      from the body-fixed frame of the body exerting the acceleration to the frame the
     *      propagation runs in.
     * @param isMutualAttractionUsed Whether the attraction of the accelerated body on the body
     *      exerting the acceleration is included, which is what makes the gravitational
     *      parameter above a sum.
     * @param updateGravitationalPotential Whether to also evaluate the potential on each
     *      update, which is only needed when it is written out as a dependent variable.
     */
    EllipsoidalHarmonicsGravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const std::shared_ptr< EllipsoidalHarmonicsGravityField > gravityField,
            const std::function< double( ) > gravitationalParameterFunction,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                    []( Eigen::Vector3d& input ) { input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    []( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            const bool updateGravitationalPotential = false ):
        subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ), gravityField_( gravityField ),
        gravitationalParameterFunction_( gravitationalParameterFunction ),
        sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
        rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
        isMutualAttractionUsed_( isMutualAttractionUsed ), currentPotential_( TUDAT_NAN ), updatePotential_( updateGravitationalPotential )
    {}

    //! Destructor.
    ~EllipsoidalHarmonicsGravitationalAccelerationModel( ) {}

    /*! Recomputes the acceleration for the given time.
     *
     * @param currentTime Time to update to. Nothing is recomputed if it is the time already set.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );

            subjectPositionFunction_( positionOfBodySubjectToAcceleration_ );
            sourcePositionFunction_( positionOfBodyExertingAcceleration_ );
            currentInertialRelativePosition_ = positionOfBodySubjectToAcceleration_ - positionOfBodyExertingAcceleration_;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;

            // The gravitational parameter is part of the field's own coefficients rather than a
            // separate factor, so wanting a different one, such as the sum of both bodies'
            // parameters, means scaling the whole field.
            const double scaling = gravitationalParameterFunction_( ) / gravityField_->getGravitationalParameter( );

            currentAccelerationInBodyFixedFrame_ = scaling * gravityField_->getGradientOfPotential( currentRelativePosition_ );
            currentAcceleration_ = rotationToIntegrationFrame_ * currentAccelerationInBodyFixedFrame_;

            if( updatePotential_ )
            {
                currentPotential_ = scaling * gravityField_->getGravitationalPotential( currentRelativePosition_ );
            }
        }
    }

    //! Position of the accelerated body relative to the accelerating one, in the body-fixed frame.
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Position of the accelerated body relative to the accelerating one, in the inertial frame.
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
    }

    //! Current rotation from the body-fixed frame to the propagation frame, as a quaternion.
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Current rotation from the body-fixed frame to the propagation frame, as a matrix.
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

    //! Function returning the gravitational parameter in use.
    std::function< double( ) > getGravitationalParameterFunction( )
    {
        return gravitationalParameterFunction_;
    }

    //! Gravity field of the body exerting the acceleration.
    std::shared_ptr< EllipsoidalHarmonicsGravityField > getGravityField( )
    {
        return gravityField_;
    }

    //! Current position of the accelerated body, in the inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodySubjectToAcceleration( )
    {
        return positionOfBodySubjectToAcceleration_;
    }

    //! Current position of the accelerating body, in the inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodyExertingAcceleration( )
    {
        return positionOfBodyExertingAcceleration_;
    }

    //! Whether the parameter in use is the sum of both bodies' gravitational parameters.
    bool getIsMutualAttractionUsed( )
    {
        return isMutualAttractionUsed_;
    }

    //! Current acceleration, in the body-fixed frame of the accelerating body.
    Eigen::Vector3d getAccelerationInBodyFixedFrame( )
    {
        return currentAccelerationInBodyFixedFrame_;
    }

    //! Potential at the last update, which is only kept if that was asked for.
    double getCurrentPotential( )
    {
        return currentPotential_;
    }

    //! Sets whether the potential is evaluated on each update.
    void resetUpdatePotential( const bool updatePotential )
    {
        updatePotential_ = updatePotential;
    }

private:
    //! Function returning the position of the accelerated body.
    const StateFunction subjectPositionFunction_;

    //! Gravity field of the body exerting the acceleration.
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > gravityField_;

    //! Function returning the gravitational parameter to use.
    const std::function< double( ) > gravitationalParameterFunction_;

    //! Function returning the position of the accelerating body.
    const StateFunction sourcePositionFunction_;

    //! Function returning the rotation from the body-fixed to the propagation frame.
    const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Whether the gravitational parameter in use is the sum of both bodies' parameters.
    const bool isMutualAttractionUsed_;

    Eigen::Vector3d positionOfBodySubjectToAcceleration_;

    Eigen::Vector3d positionOfBodyExertingAcceleration_;

    Eigen::Vector3d currentInertialRelativePosition_;

    Eigen::Vector3d currentRelativePosition_;

    Eigen::Quaterniond rotationToIntegrationFrame_;

    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    double currentPotential_;

    bool updatePotential_;
};

//! Pointer to an ellipsoidal harmonics gravitational acceleration model.
typedef std::shared_ptr< EllipsoidalHarmonicsGravitationalAccelerationModel > EllipsoidalHarmonicsGravitationalAccelerationModelPointer;

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_MODEL_H
