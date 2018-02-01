/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <Eigen/Core>

namespace tudat
{

namespace propagators
{



//! Enum listing types of dynamics that can be numerically integrated
enum IntegratedStateType
{
    hybrid = 0,
    translational_state = 1,
    transational_state = translational_state,  // deprecated (typo)
    rotational_state = 2,
    body_mass_state = 3,
    custom_state = 4
};


//! Get size of state for single propagated state of given type.
/*!
 * Get size of state for single propagated state of given type (i.e. 6 for translational state).
 * \param stateType Type of state
 * \return Size of single state.
 */
int getSingleIntegrationSize( const IntegratedStateType stateType );

//! Get order of differential equation for governing equations of dynamics of given type.
/*!
 * Get order of differential equation for governing equations of dynamics of given type (i.e. 2 for translational state).
 * \param stateType Type of state
 * \return Order of differential equations.
 */
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType );


//! Base class for calculating the state derivative model for a single type of dynamics.
/*!
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

    //! Constructor.
    /*!
     * Constructor.
     * \param integratedStateType Type of dynamics for whichh the state derivative is calculated.
     */
    SingleStateTypeDerivative( const IntegratedStateType integratedStateType ):
        integratedStateType_( integratedStateType ){ }

    //! Virtual destructor.
    virtual ~SingleStateTypeDerivative( ){ }

    //! Calculates the state derivative of the system of equations for the given type of dynamics
    /*!
     * Calculates the state derivative of the system of equations for the given type of
     * dynamics. The environment and acceleration models (updateStateDerivativeModel) must be
     * updated before calling this function. It returns the state derivative in teh form required
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

    //! Function to clear reference/cached values of state derivative model
    /*!
     * Function to clear reference/cached values of state derivative model, such as the current time and/or state.
     * This function is to be implemented in each derived class
     */
    virtual void clearStateDerivativeModel( ) = 0;

    //! Function to update the state derivative model to the current time.
    /*!
     * Function to update the state derivative model (i.e. acceleration, torque, etc. models) to the
     * current time. Note that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function
     * \param currentTime Time to which the state derivative is to be updated.
     */
    virtual void updateStateDerivativeModel( const TimeType currentTime ) = 0;

    //! Function to convert the propagator-specific form of the state to the conventional form in
    //! the global frame.
    /*!
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


    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
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

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
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

    //! Function to return the size of the state handled by the object
    /*!
     * Function to return the size of the state handled by the object
     * \return Size of the state under consideration.
     */
    virtual int getStateSize( ) = 0;

    //! Function to return the type of dynamics for which the state derivative is calculated.
    /*!
     * Function to return the type of dynamics for which the state derivative is calculated
     * \return Type of dynamics for which the state derivative is calculated.
     */
    IntegratedStateType getIntegratedStateType( )
    {
        return integratedStateType_;
    }

protected:

    //! Type of dynamics for whichh the state derivative is calculated.
    IntegratedStateType integratedStateType_;

};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_STATEDERIVATIVE_H
