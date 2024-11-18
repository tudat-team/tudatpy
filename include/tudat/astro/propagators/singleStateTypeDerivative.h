/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_STATEDERIVATIVE_H
#define TUDAT_STATEDERIVATIVE_H

#include <map>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"
#include <tudat/basics/utilityMacros.h>

namespace tudat
{

namespace propagators
{

// Enum listing types of dynamics that can be numerically integrated
// ! @get_docstring(IntegratedStateType.__docstring__)

enum IntegratedStateType
{
    hybrid = 0,
    translational_state = 1,
    rotational_state = 2,
    body_mass_state = 3,
    custom_state = 4
};

std::string getIntegratedStateTypString( const IntegratedStateType stateType );

// Get size of state for single propagated state of given type.
/*
 * Get size of state for single propagated state of given type (i.e. 6 for translational state).
 * \param stateType Type of state
 * \return Size of single state.
 */
int getSingleIntegrationSize( const IntegratedStateType stateType );

// Get order of differential equation for governing equations of dynamics of given type.
/*
 * Get order of differential equation for governing equations of dynamics of given type (i.e. 2 for translational state).
 * \param stateType Type of state
 * \return Order of differential equations.
 */
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType );

// Function to get the size of the generalized acceleration for a given state type
/*
 * Function to get the size of the generalized acceleration (e.g. acceleration for translational dynamics, torque for rotational
 * dynamics, mass rate for mass) for a given state type
 * \param stateType State type for which generalized acceleration size is to be determined
 * \return Generalized acceleration size
 */
int getGeneralizedAccelerationSize( const IntegratedStateType stateType );

// Base class for calculating the state derivative model for a single type of dynamics.
/*
 *  Base class for calculating the state derivative model for a single
 *  type of dynamics (i.e. translational, rotational, etc.). Each type
 *  of dynamics requires its own derived class. Moreover, a specific
 *  type of propagator (Cowell, Encke, etc.  for translational
 *  dynamics) may require their own further derived class, depending
 *  on the exact requirements of such a propagator.
 */
template< typename StateScalarType = double, typename TimeType = double >
class SingleStateTypeDerivative
{
public:

    // Constructor.
    /*
     * Constructor.
     * \param integratedStateType Type of dynamics for whichh the state derivative is calculated.
     */
    SingleStateTypeDerivative( const IntegratedStateType integratedStateType ):
        integratedStateType_( integratedStateType )
    {
        if( isStateToBePostProcessed( ) )
        {
            unprocessedState_.setZero( getPropagatedStateSize( ) );
        }
    }

    // Virtual destructor.
    virtual ~SingleStateTypeDerivative( ){ }

    // Calculates the state derivative of the system of equations for the given type of dynamics
    /*
     * Calculates the state derivative of the system of equations for the given type of
     * dynamics. The environment and acceleration models (updateStateDerivativeModel) must be
     * updated before calling this function. It returns the state derivative in the form required
     * for the specific type of propagator used (defined by derived class).
     * \param time Time at which the state derivative is to be calculated.
     * \param stateOfSystemToBeIntegrated Current state of the system, in the form that the equations are propagated (i.e.
     * directly from numerical integrator)
     * \param stateDerivative Derivative of the state of the system, in the form that the equations are propagated
     * (i.e. to be piped directly to numerical integrator), returned by reference.
     */
    virtual void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative ) = 0;

    // Function to clear reference/cached values of state derivative model
    /*
     * Function to clear reference/cached values of state derivative model, such as the current time and/or state.
     * This function is to be implemented in each derived class
     */
    virtual void clearStateDerivativeModel( ) = 0;

    // Function to update the state derivative model to the current time.
    /*
     * Function to update the state derivative model (i.e. acceleration, torque, etc. models) to the
     * current time. Note that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function
     * \param currentTime Time to which the state derivative is to be updated.
     */
    virtual void updateStateDerivativeModel( const TimeType currentTime ) = 0;

    // Function to convert the propagator-specific form of the state to the conventional form in
    // the global frame.
    /*
     * Function to convert the propagator-specific form of the state to the conventional form in the
     * global frame.  The conventional form is one that is typically used to represent the current
     * state in the environment (e.g. Body class). For translational dynamics this is the Cartesian
     * position and velocity).  The inertial frame is typically the barycenter with J2000/ECLIPJ2000
     * orientation, but may differ depending on simulation settings
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration).
     * \param time Current time at which the state is valid.
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the 'conventional form' in inertial
     * coordinates, that can for instance be set directly  in the body object (returned by reference).
     */
    virtual void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton ) = 0;

    // Function to convert the state in the conventional form to the propagator-specific form.
    /*
     * Function to convert the state in the conventional form to the propagator-specific form.  The
     * conventional form is one that is typically used to represent the current state in the
     * environment (e.g. Body class). For translational dynamics this is the Cartesian position and
     * velocity).
     * \param outputSolution State in 'conventional form'
     * \param time Current time at which the state is valid.
     * \return State (outputSolution), converted to the 'propagator-specific form'
     */
    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time ) = 0;

    // Function to convert the propagator-specific form of the state to the conventional form.
    /*
     * Function to convert the propagator-specific form of the state to the conventional form. The
     * conventional form is one that is typically used to represent the current state in the
     * environment (e.g. Body class). For translational dynamics this is the Cartesian position and
     * velocity).  In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration).
     * \param time Current time at which the state is valid.
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the 'conventional form' (returned by
     * reference).
     */
    virtual void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton ) = 0;

    // Function to return the size of the conventional state handled by the object.
    /*
     * Function to return the size of the conventional state handled by the object. This is the size of the conventional
     * propagation state, e.g., size of Cartesian state for translational propagation.
     * \return Size of the state under consideration.
     */
    virtual int getConventionalStateSize( ) = 0;

    // Function to return the size of the propagated state handled by the object.
    /*
     * Function to return the size of the propagated state handled by the object. This is the size of the actual propagation
     * state, e.g., size of USM7 state for translational propagation.
     * \return Size of the propagated state under consideration.
     */
    virtual int getPropagatedStateSize( )
    {
        return getConventionalStateSize( );
    }

    // Function to return the type of dynamics for which the state derivative is calculated.
    /*
     * Function to return the type of dynamics for which the state derivative is calculated
     * \return Type of dynamics for which the state derivative is calculated.
     */
    IntegratedStateType getIntegratedStateType( )
    {
        return integratedStateType_;
    }

    // Function to process the state vector during propagation.
    /*
     * Function to process the state during propagation. Is especially useful for attitude states (e.g.,
     * normalization of quaternions and transformation to/from shadow attitude parameters).
     * \param unprocessedState State computed after propagation.
     */
    virtual void postProcessState( Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > unprocessedState )
    {
        TUDAT_UNUSED_PARAMETER( unprocessedState );
    }

    virtual void postProcessState( Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& unprocessedState )
    {
        unprocessedState_ = unprocessedState.block( 0, 0, getPropagatedStateSize( ), 1 );
        postProcessState( unprocessedState_.block( 0, 0, getPropagatedStateSize( ), 1 ) );
        unprocessedState.block( 0, 0, getPropagatedStateSize( ), 1 ) = unprocessedState_;
    }

    // Function to return whether the state needs to be post-processed.
    /*
     * Function to return whether the state needs to be post-processed. Default value is false.
     * \return Boolean informing whether the state needs to be post-processed.
     */
    virtual bool isStateToBePostProcessed( )
    {
        return false;
    }

protected:

    // Type of dynamics for which the state derivative is calculated.
    IntegratedStateType integratedStateType_;

    // Vector used during post-processing of state.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > unprocessedState_;
};

//extern template class SingleStateTypeDerivative< double, double >;


} // namespace propagators

} // namespace tudat


namespace std
{

// Hash for IntegratedStateType enum.
template< >
struct hash< tudat::propagators::IntegratedStateType >
{
    typedef tudat::propagators::IntegratedStateType argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const
    {
        using type = typename std::underlying_type<argument_type>::type;
        return std::hash< type >( )( static_cast< type >( x ) );
    }
};

} // namespace std


#endif // TUDAT_STATEDERIVATIVE_H
