/*! \file integrator.h
 *    Header file that defines the base class for all integration methods
 *    included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 12
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 28 July, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *      Currently, individual boolean flags are used to check whether necessary
 *      set functions have been called to execute integration. In the future,
 *      if it turns out to be a performance bottleneck, it might be wise to
 *      consider using the STL bitset object to group these flags.
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100907    K. Kumar          File header and footer added.
 *      100908    K. Kumar          Edited to conform to protocols.
 *                                  Changed variable/function names to more
 *                                  suitable choices.
 *      100926    K. Kumar          Added solution to pointer-to-member-
 *                                  function problem and comments under
 *                                  Notes, removed fixed output interval
 *                                  option.
 *      100928    K. Kumar          Added missing comments.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class;
 *                                  implemented adaptor class to replace
 *                                  pointer-to-member functions.
 *      110203    J. Melman         Some minor comment corrections.
 *      110206    J. Melman         isSaveIntegrationHistory changed into
 *                                  isIntegrationHistoryToBeSaved.
 *      110207    K. Kumar          Path changed; added integrate() virtual
 *                                  function.
 *      110516    K. Kumar          Changed architecture so adaptor is not
 *                                  used. Global functions can no longer be
 *                                  used. StateDerivativeBase used to define
 *                                  classes with state derivative function.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

// Include statements.
#include <vector>
#include "basicMathematicsFunctions.h"
#include "state.h"
#include "stateDerivativeBase.h"

//! Integrator class.
/*!
 * Base class for all integration methods included in Tudat.
 */
class Integrator
{
public:

    //! Default constructor.
    /*!
     * Integrator class constructor.
     */
    Integrator( );

    //! Default destructor.
    /*!
     * Integrator class destructor.
     */
    virtual ~Integrator( );

    //! Set object containing state derivative.
    /*!
     * Sets class object containing the function to compute the state
     * derivative. The state derivative function must have the following form:
     * void computeStateDerivative( const double&, State*, State* ). The name
     * of the function must be 'computeStateDerivative'.
     * \param pointerToStateDerivative Pointer to class containing state
     *          derivative function.
     */
    void setObjectContainingStateDerivative( StateDerivativeBase*
                                             pointerToStateDerivative );

    //! Set initial state.
    /*!
     * Sets the initial state to be used during integration as a pointer to a
     * State object. This function must be called at least once for
     * integration to proceed.
     * \param pointerToInitialState Initial state given as pointer to State
     *          object.
     */
    void setInitialState( State* pointerToInitialState );

    //! Set initial stepsize.
    /*!
     * Sets the initial stepsize for all integration methods. This function
     * must be called at least once for integration to proceed.
     * \param initialStepsize Initial stepsize for integration method.
     */
    void setInitialStepsize( const double& initialStepsize );

    //! Set start of integration interval.
    /*!
     * Sets the start of the integration interval. This function must be called
     * at least once for integration to proceed.
     * \param integrationIntervalStart Start of integration interval.
     */
    void setIntegrationIntervalStart( const double& integrationIntervalStart );

    //! Set end of integration interval.
    /*!
     * Sets the end of the integration interval. This function must be called
     * at least once for integration to proceed.
     * \param integrationIntervalEnd End of integration interval.
     */
    void setIntegrationIntervalEnd( const double& integrationIntervalEnd );

    //! Get stepsize.
    /*!
     * Returns the stepsize.
     * \return Stepsize.
     */
    double& getStepsize( );

    //! Get number of integration steps.
    /*!
     * Returns the number of integration steps.
     * \return Number of integration steps.
     */
    unsigned int& getNumberOfIntegrationSteps( );

    //! Get start of integration interval.
    /*!
     * Returns the start of the integration interval.
     * \return Start of integration interval.
     */
    double& getIntegrationIntervalStart( );

    //! Get end of integration interval.
    /*!
     * Returns the end of the integration interval.
     * \return End of integration interval.
     */
    double& getIntegrationIntervalEnd( );

    //! Get initial state.
    /*!
     * Returns the initial state as a pointer to a State object.
     * \return Initial state given as pointer to State object.
     */
    State* getInitialState( );

    //! Get final state.
    /*!
     * Returns the final state as a pointer to a State object.
     * \return Final state given as pointer to State object.
     */
    State* getFinalState( );

    //! Integrate.
    /*!
     * This function executes an integration.
     */
    virtual void integrate( ) = 0;

protected:

    //! Dimension of current state.
    /*!
     * Dimension of current state.
     */
    unsigned int dimensionOfState_;

    //! Number of integration steps.
    /*!
     * Number of integration steps necessary to span integration interval.
     */
    unsigned int numberOfIntegrationSteps_;

    //! Stepsize.
    /*!
     * Stepsize for integration.
     */
    double stepsize_;

    //! Intial stepsize.
    /*!
     * Initial stepsize.
     */
    double initialStepsize_;

    //! Start of integration interval.
    /*!
     * Start of integration interval.
     */
    double integrationIntervalStart_;

    //! End of integration interval.
    /*!
     * End of integration interval.
     */
    double integrationIntervalEnd_;

    //! Size of integration interval.
    /*!
     * Size of integration interval.
     */
    double integrationInterval_;

    //! Current point in integration interval.
    /*!
     * Current point in integration interval. This is "time" in many
     * simulations.
     */
    double integrationIntervalCurrentPoint_;

    //! Stepsize of last step.
    /*!
     * Stepsize for last step of integration interval.
     */
    double lastStepStepsize_;

    //! Pointer to initial state.
    /*!
     * Pointer to initial state.
     */
    State* pointerToInitialState_;

    //! State derivative.
    /*!
     * A State object containing the state derivative.
     */
    State stateDerivative_;

    //! Final state.
    /*!
     * A State object containing the final state.
     */
    State finalState_;

    //! Vector containing current states.
    /*!
     * Vector containing the current states given as State objects.
     * It is a vector of states, since it has to be possible to contain more
     * than one state for the multi-step methods.  In the case of single-step
     * methods only the first entry of the vector is used.
     */
    std::vector< State > vectorOfCurrentStates_;

    //! Pointer to state derivative.
    /*!
     * Pointer to class object containing state derivative function.
     */
    StateDerivativeBase* pointerToStateDerivative_;

    //! Compute state derivative.
    /*!
     * Computes the state derivative by calling the object set containing the
     * state derivative function. This is simply an internal wrapper function.
     * \param integrationIntervalCurrentPoint Current point in integration
     *          interval.
     * \param pointerToState State given as pointer to State object.
     * \param pointerToStateDerivative State derivative given as pointer to
     *          State object. This is where the computed state derivative is
     *          stored.
     */
    void computeStateDerivative_( double& integrationIntervalCurrentPoint,
                                  State* pointerToState,
                                  State* pointerToStateDerivative );

private:

    //! Flag to indicate if stepsize is set.
    /*!
     * Flag to indicate if stepsize is set.
     */
    bool isStepsizeSet_;

    //! Flag to indicate if initial state is set.
    /*!
     * Flag to indicate if initial state is set.
     */
    bool isInitialStateSet_;

    //! Flag to indicate if start of integration interval is set.
    /*!
     * Flag to indicate if start of integration interval is set.
     */
    bool isIntegrationIntervalStartSet_;

    //! Flag to indicate if end of integration interval is set.
    /*!
     * Flag to indicate if end of integration interval is set.
     */
    bool isIntegrationIntervalEndSet_;

    //! Flag to indicate if state derivative is set.
    /*!
     * Flag to indicate if state derivative is set.
     */
    bool isStateDerivativeSet_;

    //! Compute internal derived integration parameters.
    /*!
     * Computes internal derived integration parameters necessary to perform
     * integration.
     */
    void computeInternalDerivedIntegrationParameters_( );
};

#endif // INTEGRATOR_H

// End of file.
