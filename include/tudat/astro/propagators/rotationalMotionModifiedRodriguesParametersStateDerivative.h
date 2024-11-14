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

#ifndef TUDAT_ROTATIONAL_MOTION_MODIFIED_RODRIGUES_PARAMETERS_STATE_DERIVATIVE_H
#define TUDAT_ROTATIONAL_MOTION_MODIFIED_RODRIGUES_PARAMETERS_STATE_DERIVATIVE_H

#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/astro/basic_astro/attitudeElementConversions.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the time derivative of modified Rodrigues parameters (in vector representation) of body-fixed to inertial frame.
/*!
 * Function to obtain the time derivative of modified Rodrigues parameters (in vector representation) of body-fixed to inertial frame.
 * \param currentModiefiedRodriguesParametersToBaseFrame Modified Rodrigues parameters (in vector representation) that define the rotation
 * from body-fixed to inertial frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of modified Rodrigues parameters (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateModifiedRodriguesParametersDerivative(
        const Eigen::Vector4d& currentModiefiedRodriguesParametersToBaseFrame,
        const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Class for computing the state derivative for rotational dynamics of N bodies.
/*!
 *  Class for computing the state derivative for rotational dynamics of N bodies, using modified Rodrigues parameters from body-fixed to inertial
 *  frame (in modified Rodrigues parameters format) and angular velocity-vector of body expressed in body-fixed frame as the rotational state of a
 *  single body
 */
template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionModifiedRodriguesParametersStateDerivative: public RotationalMotionStateDerivative< StateScalarType, TimeType >
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
    RotationalMotionModifiedRodriguesParametersStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const std::vector< std::string >& bodiesToPropagate,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions =
            std::vector< std::function< Eigen::Matrix3d( ) > >( ) ):
        RotationalMotionStateDerivative< StateScalarType, TimeType >(
            torqueModelsPerBody, modified_rodrigues_parameters, bodiesToPropagate, bodyInertiaTensorFunctions,
            bodyInertiaTensorTimeDerivativeFunctions )
    { }

    //! Destructor
    ~RotationalMotionModifiedRodriguesParametersStateDerivative( ){ }

    //! Calculates the state derivative of the rotational motion of the system.
    /*!
     *  Calculates the state derivative of the rotational motion of the system at the given time and rotational state.
     *  \param time Time (seconds since reference epoch) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 7 * bodiesToPropagate_.size( ), containing rotation modified Rodrigues
     *  parameters and angular velocity of the bodies being propagated. The order of the values is defined by the order of
     *  bodies in bodiesToPropagate_.
     *  \param stateDerivative Current state derivative (modified Rodrigues parameters rate + angular acceleration) of
     *  system of bodies integrated numerically (returned by reference).
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
            Eigen::Matrix< StateScalarType, 4, 1 > currentModifiedRodriguesParameters = stateOfSystemToBeIntegrated.block( i * 7, 0, 4, 1 );
            Eigen::Matrix< StateScalarType, 3, 1 > currentBodyFixedRotationRate = stateOfSystemToBeIntegrated.block( i * 7 + 4, 0, 3, 1 );

            stateDerivative.block( i * 7, 0, 4, 1 ) = calculateModifiedRodriguesParametersDerivative(
                        currentModifiedRodriguesParameters.template cast< double >( ),
                        currentBodyFixedRotationRate.template cast< double >( ) ).
                    template cast< StateScalarType >( );
            stateDerivative.block( i * 7 + 4, 0, 3, 1 ) = evaluateRotationalEquationsOfMotion(
                        this->bodyInertiaTensorFunctions_.at( i )( ), torquesActingOnBodies.at( i ),
                        currentBodyFixedRotationRate.template cast< double >( ),
                        this->bodyInertiaTensorTimeDerivativeFunctions_.at( i )( ) ).template cast< StateScalarType >( );
        }
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form.
     * \param outputSolution State in 'conventional form'.
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form'.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( this->getPropagatedStateSize( ) );

        // Convert state to modified Rodrigues parameters for each body
        for( unsigned int i = 0; i < this->bodiesToPropagate_.size( ); i++ )
        {
            currentState.segment( i * 7, 4 ) =
                    orbital_element_conversions::convertQuaternionsToModifiedRodriguesParameterElements(
                        outputSolution.block( i * 7, 0, 4, 1 ).template cast< double >( ) ).template cast< StateScalarType >( );
            currentState.segment( i * 7 + 4, 3 ) = outputSolution.block( i * 7 + 4, 0, 3, 1 ); // rotational velocity is the same
        }

        return currentState;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form.
     * \param internalSolution State in propagator-specific form.
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentLocalSolution State (internalSolution), converted to the 'conventional form' (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSolution )
    {
        // Convert state to quaternions for each body
        for( unsigned int i = 0; i < this->bodiesToPropagate_.size( ); i++ )
        {
            currentLocalSolution.block( i * 7, 0, 4, 1 ) =
                    orbital_element_conversions::convertModifiedRodriguesParametersToQuaternionElements(
                        internalSolution.block( i * 7, 0, 4, 1 ).template cast< double >( ) ).template cast< StateScalarType >( );
            currentLocalSolution.block( i * 7 + 4, 0, 3, 1 ) = internalSolution.block( i * 7 + 4, 0, 3, 1 ); // rotational velocity is the same
        }
    }

    //! Function to process the state during propagation.
    /*!
     * Function to process the state during propagation. For modified Rodrigues parameters (MRP), this function converts to/from
     * shadow modified Rodrigues parameters (SMRP), in case the magnitude of the (S)MRP vector is larger than 1.0.
     * \param unprocessedState State computed after propagation.
     * \return Processed state (returned by reference).
     */
    void postProcessState( Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > unprocessedState )
    {
        // Loop over each body
        Eigen::Matrix< StateScalarType, 3, 1 > modifiedRodriguesParametersVector;
        StateScalarType modifiedRodriguesParametersMagnitude;
        for( unsigned int i = 0; i < this->bodiesToPropagate_.size( ); i++ )
        {
            // Convert to/from shadow modifed Rodrigues parameters (SMRP) (transformation is the same either way)
            modifiedRodriguesParametersVector = unprocessedState.block( i * 7, 0, 3, 1 );
            modifiedRodriguesParametersMagnitude = modifiedRodriguesParametersVector.norm( );
            if ( modifiedRodriguesParametersMagnitude >= 1.0 )
            {
                // Invert flag
                unprocessedState.block( i * 7 + 3, 0, 1, 1 ) = ( unprocessedState.block( i * 7 + 3, 0, 1, 1 ) -
                                                             Eigen::Matrix< StateScalarType, 1, 1 >::Ones( ) ).cwiseAbs( );

                // Convert to MRP/SMRP
                modifiedRodriguesParametersVector /= - modifiedRodriguesParametersMagnitude *
                        modifiedRodriguesParametersMagnitude;

                // Replace MRP with SMPR, or vice-versa
                unprocessedState.block( i * 7, 0, 3, 1 ) = modifiedRodriguesParametersVector;
            }
        }
    }

    //! Function to return whether the state needs to be post-processed.
    /*!
     * Function to return whether the state needs to be post-processed. For (shadow) modified Rodrigues parameters this is true.
     * \return Boolean confirming that the state needs to be post-processed.
     */
    bool isStateToBePostProcessed( )
    {
        return true;
    }

private:

};

//extern template class RotationalMotionModifiedRodriguesParametersStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat

#endif // TUDAT_ROTATIONAL_MOTION_MODIFIED_RODRIGUES_PARAMETERS_STATE_DERIVATIVE_H
