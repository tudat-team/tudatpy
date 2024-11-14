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

#ifndef TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H
#define TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H

#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the matrix by which a quaternion vector is to be pre-multiplied to obtain this
//! quaternion's time-derivative.
/*!
 * Function to obtain the matrix by which a quaternion vector (representing body-fixed to inertial frame rotation)
 * is to be pre-multiplied to obtain this quaternion's time-derivative.
 * \param angularVelocityVectorInBodyFixedFrame  Current angular velocity vector of body, expressed in its
 * body-fixed frame.
 * \return Matrix by which a quaternion vector (representing body-fixed to inertial frame rotation) is to be
 * pre-multiplied to obtain this quaternion's time-derivative.
 */
Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Function to obtain the matrix by which an angular velocity vector is to be pre-multiplied to obtain the
//! quaternion's time-derivative.
/*!
 * Function to obtain the matrix by which an angular velocity vector of body (expressed in its body-fixed frame)
 * is to be pre-multiplied to obtain the quaternion's time-derivative.
 * \param quaternionVector Current quaternion vector, representing body-fixed to inertial frame rotation.
 * \return Matrix by which an angular velocity vector of body (expressed in its body-fixed frame) is to be
 * pre-multiplied to obtain the quaternion's time-derivative.
 */
Eigen::Matrix< double, 4, 3 > getAngularVelocityToQuaternionRateMatrix( const Eigen::Vector4d& quaternionVector );

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
/*!
 * Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 * \param currentQuaternionsToBaseFrame Quaternions (in vector representation) that define the rotation from body-fixed to inertial
 * frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateQuaternionDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
                                               const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Class for computing the state derivative for rotational dynamics of N bodies.
/*!
 *  Class for computing the state derivative for rotational dynamics of N bodies, using quaternion from body-fixed to inertial
 *  frame (in quaternion format) and angular velocity-vector of body expressed in body-fixed frame as the rotational state of a
 *  single body
 */
template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionQuaternionsStateDerivative: public RotationalMotionStateDerivative< StateScalarType, TimeType >
{
public:

    using SingleStateTypeDerivative< StateScalarType, TimeType >::postProcessState;

    //! Constructor.
    /*!
     * Constructor.
     * \param torqueModelsPerBody List of torque models (first map key body undergoing acceleration, second map key body exerting
     * acceleration).
     * \param bodiesToPropagate List of names of bodies for which rotational state is to be propagated.
     * \param bodyInertiaTensorFunctions List of functions returning inertia tensors of bodiesToPropagate (in same order).
     * \param bodyInertiaTensorTimeDerivativeFunctions List of functions returning time derivatives of inertia tensors of
     *  bodiesToPropagate (in same order). Default empty, denoting time-invariant inertia tensors.
     */
    RotationalMotionQuaternionsStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const std::vector< std::string >& bodiesToPropagate,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions =
            std::vector< std::function< Eigen::Matrix3d( ) > >( ) ):
        RotationalMotionStateDerivative< StateScalarType, TimeType >(
            torqueModelsPerBody, quaternions, bodiesToPropagate, bodyInertiaTensorFunctions,
            bodyInertiaTensorTimeDerivativeFunctions )
    { }

    //! Destructor
    ~RotationalMotionQuaternionsStateDerivative( ){ }

    //! Calculates the state derivative of the rotational motion of the system.
    /*!
     *  Calculates the state derivative of the rotational motion of the system at the given time and rotational state.
     *  \param time Time (seconds since reference epoch) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 7 * bodiesToPropagate_.size( ), containing rotation quaternion and
     *  angular velocity of the bodies being propagated. The order of the values is defined by the order of bodies in
     *  bodiesToPropagate_.
     *  \param stateDerivative Current state derivative (quaternion rate + angular acceleration) of system of bodies
     *  integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( stateOfSystemToBeIntegrated.rows( ), 1 );
        std::vector< Eigen::Vector3d > torquesActingOnBodies = this->sumTorquesPerBody( );

        for( unsigned int i = 0; i < torquesActingOnBodies.size( ); i++ )
        {
            Eigen::Matrix< StateScalarType, 4, 1 > currentQuaternions = stateOfSystemToBeIntegrated.block( i * 7, 0, 4, 1 );
            Eigen::Matrix< StateScalarType, 3, 1 > currentBodyFixedRotationRate = stateOfSystemToBeIntegrated.block( i * 7 + 4, 0, 3, 1 );

            stateDerivative.block( i * 7, 0, 4, 1 ) = calculateQuaternionDerivative(
                        currentQuaternions.template cast< double >( ), currentBodyFixedRotationRate.template cast< double >( ) ).
                    template cast< StateScalarType >( );
            stateDerivative.block( i * 7 + 4, 0, 3, 1 ) = evaluateRotationalEquationsOfMotion(
                        this->bodyInertiaTensorFunctions_.at( i )( ), torquesActingOnBodies.at( i ),
                        currentBodyFixedRotationRate.template cast< double >( ),
                        this->bodyInertiaTensorTimeDerivativeFunctions_.at( i )( ) ).template cast< StateScalarType >( );
        }
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param outputSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param internalSolution State in propagator-specific form (which is equal to outputSolution to conventional form for
     * this propagator)
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentLocalSolution State (internalSolution), converted to the 'conventional form',
     * which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSolution )
    {
        currentLocalSolution = internalSolution;
    }

    //! Function to process the state during propagation.
    /*!
     * Function to process the state during propagation. For quaternions, this function normalizes the quaternion vector
     * in case its magnitude differs from 1.0 by a value larger than the tolerance.
     * \param unprocessedState State computed after propagation.
     * \return Processed state (returned by reference).
     */
    void postProcessState( Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > unprocessedState )
    {
        // Loop over each body
        const double tolerance = 20.0 * std::numeric_limits< double >::epsilon( );
        Eigen::Matrix< StateScalarType, 4, 1 > quaternionsVector;
        StateScalarType quaternionsMagnitude;
        for( unsigned int i = 0; i < this->bodiesToPropagate_.size( ); i++ )
        {
            // Normalize quaternions
            quaternionsVector = unprocessedState.block( i * 7, 0, 4, 1 );
            quaternionsMagnitude = quaternionsVector.norm( );
            if ( std::fabs( quaternionsMagnitude - 1.0 ) >= tolerance )
            {
                // Normalize
                quaternionsVector /= quaternionsMagnitude;

                // Replace old quaternions with normalized quaternions
                unprocessedState.block( i * 7, 0, 4, 1 ) = quaternionsVector;
            }
        }
    }

    //! Function to return whether the state needs to be post-processed.
    /*!
     * Function to return whether the state needs to be post-processed. For quaternions this is true.
     * \return Boolean confirming that the state needs to be post-processed.
     */
    bool isStateToBePostProcessed( )
    {
        return true;
    }

private:

};


//extern template class RotationalMotionQuaternionsStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat

#endif // TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H
