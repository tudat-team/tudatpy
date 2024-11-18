/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Precise computation of acceleration due to uniform ring or disk, Toshio Fukushima (2010), Celestial Mechanics
 *          and Dynamical Astronomy, 108:339–356.
 */

#ifndef TUDAT_RINGGRAVITYMODEL_H
#define TUDAT_RINGGRAVITYMODEL_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/ringGravityField.h"

namespace tudat
{
namespace gravitation
{

class RingGravitationalAccelerationModel: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{

protected:
    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector3d& ) > StateFunction;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of ring parameters.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and volume of the
     * body exerting the acceleration, polyhedron parameters, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter and volume. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param aRingRadius A (constant) ring radius [m].
     * \param ellipticIntegralSFromDAndB Flag indicating whether to compute S(m) from D(m) and B(m) (if true),
     *      or from K(m) and E(m) (if false). The former has a lower loss of accuracy due to numerical cancellation
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param ringCache Cache object for computing/retrieving repeated terms in ring potential gradient calculation.
     * \param updatePotential Flag indicating whether to update the gravitational potential when calling
     * the updateMembers function.
     */
    RingGravitationalAccelerationModel (
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double aRingRadius,
            const bool ellipticIntegralSFromDAndB,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                    [ ]( Eigen::Vector3d& input) { input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond ( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    [ ] ( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            const bool updateGravitationalPotential = false )
        : subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction_( [ = ]( ){ return aGravitationalParameter; } ),
          ringRadiusFunction_( [ = ]( ){ return aRingRadius; } ),
          ellipticIntegralSFromDAndB_( ellipticIntegralSFromDAndB ),
          sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          ringCache_( std::make_shared< RingGravityCache >( aRingRadius, ellipticIntegralSFromDAndB_ ) ),
          currentPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential )
    { }

    //! Constructor taking position-functions for bodies, and constant parameters of ring parameters.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and volume of the
     * body exerting the acceleration, polyhedron parameters, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter and volume. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param gravitationalParameterFunction Pointer to function returning the gravitational parameter.
     * \param ringRadiusFunction Pointer to function returning the ring radius [m].
     * \param ellipticIntegralSFromDAndB Flag indicating whether to compute S(m) from D(m) and B(m) (if true),
     *      or from K(m) and E(m) (if false). The former has a lower loss of accuracy due to numerical cancellation
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param ringCache Cache object for computing/retrieving repeated terms in ring potential gradient calculation.
     * \param updatePotential Flag indicating whether to update the gravitational potential when calling
     * the updateMembers function.
     */
    RingGravitationalAccelerationModel (
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const std::function< double() > gravitationalParameterFunction,
            const std::function< double() > ringRadiusFunction,
            const bool ellipticIntegralSFromDAndB,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                    [ ]( Eigen::Vector3d& input) { input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond ( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    [ ] ( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            const bool updateGravitationalPotential = false )
        : subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction_( gravitationalParameterFunction ),
          ringRadiusFunction_( ringRadiusFunction ),
          ellipticIntegralSFromDAndB_( ellipticIntegralSFromDAndB ),
          sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          ringCache_( std::make_shared< RingGravityCache >( ringRadiusFunction(), ellipticIntegralSFromDAndB_ ) ),
          currentPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential )
    { }

    ~RingGravitationalAccelerationModel( ){ }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class members of this class.
     * The potential and laplacian of potential are only updated if the associated flags indicate so.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {

            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );

            subjectPositionFunction_( positionOfBodySubjectToAcceleration_ );
            sourcePositionFunction_( positionOfBodyExertingAcceleration_ );
            currentInertialRelativePosition_ =
                positionOfBodySubjectToAcceleration_ - positionOfBodyExertingAcceleration_ ;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;

            ringCache_->setRingRadius( ringRadiusFunction_( ) );
            ringCache_->update( currentRelativePosition_ );

            // Compute the current acceleration
            currentAccelerationInBodyFixedFrame_ = computeRingGravitationalAcceleration(
                currentRelativePosition_,
                ringRadiusFunction_( ),
                gravitationalParameterFunction_(),
                ringCache_->getEllipticIntegralB( ),
                ringCache_->getEllipticIntegralE( ),
                ringCache_->getEllipticIntegralS( ) );

            currentAcceleration_ = rotationToIntegrationFrame_ * currentAccelerationInBodyFixedFrame_;

            // Compute the current gravitational potential
            if ( updatePotential_ )
            {
                currentPotential_ = computeRingGravitationalPotential(
                    currentRelativePosition_,
                    ringRadiusFunction_( ),
                    gravitationalParameterFunction_(),
                    ringCache_->getEllipticIntegralK( ) );
            }

        }
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
    //! fixed to body undergoing acceleration
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
    //! frame
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

    //! Function to return the function returning the relevant gravitational parameter.
    std::function< double( ) > getGravitationalParameterFunction( )
    {
        return gravitationalParameterFunction_;
    }

    //! Function to return the radius of the ring.
    std::function< double( ) > getRingRadiusFunction( )
    {
        return ringRadiusFunction_;
    }

    //! Function to return current position vector of body exerting gravitational acceleration in inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodySubjectToAcceleration( )
    {
        return positionOfBodySubjectToAcceleration_;
    }

    //! Function to return current position vector of body undergoing gravitational acceleration in inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodyExertingAcceleration( )
    {
        return positionOfBodyExertingAcceleration_;
    }

    //! Function to return the function returning position of body exerting acceleration.
    /*!
     * Function to return the function returning position of body exerting acceleration.
     * \return Function returning position of body exerting acceleration.
     */
    StateFunction getStateFunctionOfBodyExertingAcceleration( )
    { return sourcePositionFunction_; }

    //! Function to return the function returning position of body subject to acceleration.
    /*!
     * Function to return the function returning position of body subject to acceleration.
     * \return Function returning position of body subject to acceleration.
     */
    StateFunction getStateFunctionOfBodyUndergoingAcceleration( )
    { return subjectPositionFunction_; }

    //! Function to retrieve the ring cache for this acceleration.
    std::shared_ptr< RingGravityCache > getRingCache( )
    {
        return ringCache_;
    }

    //! Function to return the value of the current gravitational potential.
    double getCurrentPotential ( )
    { return currentPotential_; }

    //! Function to return the update potential flag.
    bool getUpdatePotential ( )
    { return updatePotential_; }

    //! Function to reset the update potential flag.
    void resetUpdatePotential ( bool updatePotential )
    { updatePotential_ = updatePotential; }

private:

    //! Pointer to function returning position of body subject to acceleration.
    const StateFunction subjectPositionFunction_;

    //! Function returning a gravitational parameter [m^3 s^-2].
    const std::function< double( ) > gravitationalParameterFunction_;

    //! Function returning the ring radius [m].
    const std::function< double() > ringRadiusFunction_;

    //! Flag indicating whether to compute S(m) from D(m) and B(m) (if true), or from K(m) and E(m) (if false)
    const bool ellipticIntegralSFromDAndB_;

    //! Pointer to function returning position of body exerting acceleration.
    const StateFunction sourcePositionFunction_;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Variable denoting whether mutual acceleration between bodies is included.
    bool isMutualAttractionUsed_;

    //!  Ring cache for this acceleration
    std::shared_ptr< RingGravityCache > ringCache_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    //! Current acceleration in frame fixed to body undergoing acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    //! Position of body subject to acceleration.
    Eigen::Vector3d positionOfBodySubjectToAcceleration_;

    //! Position of body exerting acceleration.
    Eigen::Vector3d positionOfBodyExertingAcceleration_;

    //! Current gravitational potential acting on the body undergoing acceleration, as computed by last call to
    //! updateMembers function
    double currentPotential_;

    //! Flag indicating whether to update the gravitational potential when calling the updateMembers function.
    bool updatePotential_;
};


} // namespace gravitation

} // namespace tudat

#endif //TUDAT_RINGGRAVITYMODEL_H
