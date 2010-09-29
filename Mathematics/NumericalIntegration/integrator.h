/*! \file integrator.h
 *    Header file that defines the base class for all integration methods
 *    included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dikx@student.tudelft.nl
 *
 *    Date created      : 28 July, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *      Currently, individual boolean flags are used to check whether necessary
 *      set functions have been called to execute integration. In the future,
 *      if it turns out to be a performance bottleneck, it might be wise to
 *      consider using the STL bitset object to group these flags.
 *
 *      Currently, a workaround for a well-known pointer-to-member-function
 *      issue has been implemented to be able to use the Integrator class with
 *      global functions and member functions of the Propagator class, although
 *      this is not optimized for performance. In the future, the use of
 *      functors and templates should be implemented since this is C++-style
 *      and optimized for performance.
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
 *      YYMMDD    author              comment
 *      100907    K. Kumar            File header and footer added.
 *      100908    K. Kumar            Edited to conform to protocols
 *                                    Changed variable/function names to more
 *                                    suitable choices.
 *      100926    K. Kumar            Added solution to pointer-to-member-
 *                                    function problem and comments under
 *                                    Notes, removed fixed output interval
 *                                    option.
 *      100928    K. Kumar            Added missing comments.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor comment modifications.
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

// Include statements.
#include <map>
#include <vector>
#include "numericalPropagator.h"
#include "linearALgebra.h"

// Forward declartions.
class NumericalPropagator;

//! Integrator class.
/*!
 * Base class for all integration methods included in Tudat.
 */
class Integrator
{
public:

    // Definition of typedefs.
    // State derviative functions can be passed as pointers to global
    // functions or pointers to member functions of the NumericalPropagator
    // class.
    // When the integrator is to use a function to calculate derivatives from a
    // function which is neither global nor a member of NumericalPropagator,
    // a new (member)-function pointer is necessary.
    // A typedef should be added here for such a pointer.
    typedef void ( *pointerToVectorVectorTakingFunction ) ( Vector&, Vector& );

    typedef void ( NumericalPropagator::
                   *pointerToStateDerivativeFunction ) ( Vector&, Vector& );

    //! Integrator class constructor.
    /*!
     * Integrator class constructor.
     */
    Integrator( );

    //! Integrator class destructor.
    /*!
     * Integrator class destructor.
     */
    virtual ~Integrator( );

    //! Set address of function containing state derivatives.
    /*!
     * This function sets the address of the function that returns the state
     * derivatives to be used for integration. This function must be called at
     * least once for integration to proceed.
     * \param derivativeFunction Address of function containing state derivatives.
     */
    void setStateDerivativeFunction( pointerToVectorVectorTakingFunction
                                     derivativeFunction );

    //! Set address of function containing state derivatives.
    /*!
     * This function sets the address of the function that returns the state
     * derivatives to be used for integration. This function must be called at
     * least once for integration to proceed.
     * \param derivativeFunction Address of function containing state derivatives.
     * \param pointerToNumericalPropagator Pointer to object of
     * NumericalPropagator class containing state derivative function.
     */
    void setStateDerivativeFunction( pointerToStateDerivativeFunction
                                     derivativeFunction,
                                     NumericalPropagator*
                                     pointerToNumericalPropagator );

    //! Set initial state vector.
    /*!
     * This function sets the initial state vector to be used during
     * integration. This function must be called at least once for integration
     * to proceed.
     * \param initialStateVector Initial state vector.
     */
    void setInitialStateVector( Vector& initialStateVector );

    //! Set value of initial stepsize.
    /*!
     * This function sets the value of the initial stepsize for all integration
     * methods.
     * This function must be called at least once for integration to proceed.
     * \param stepsize Stepsize for fixed-stepsize integration methods.
     */
    void setInitialStepsize( const double& stepsize );

    //! Set start of integration interval.
    /*!
     * This function sets the value of the start of the integration interval.
     * This function must be called at least once for integration to proceed.
     * \param integrationIntervalStart Start of integration interval.
     */
    void setIntegrationIntervalStart( const double& integrationIntervalStart );

    //! Set end of integration interval.
    /*!
     * This function sets the value of the end of the integration interval.
     * This function must be called at least once for integration to proceed.
     * \param integrationIntervalEnd End of integration interval.
     */
    void setIntegrationIntervalEnd( const double& integrationIntervalEnd );

    //! Save integration history.
    /*!
     * This function sets whether integration history should be saved. By
     * default, it is not saved.
     * \param isSaveIntegrationHistory Boolean flag to indicate whether to
     *          save integration history.
     */
    void saveIntegrationHistory( const bool& isSaveIntegrationHistory );

    //! Get stepsize.
    /*!
     * This function returns the stepsize.
     * \return Stepsize.
     */
    double getStepsize( );

    //! Get number of integration steps.
    /*!
     * This function returns the number of integration steps.
     * \return Number of integration steps.
     */
    unsigned int getNumberOfIntegrationSteps( );

    //! Get start of integration interval.
    /*!
     * This function returns the value of the start of the integration
     * interval.
     * \return Start of integration interval.
     */
    double getIntegrationIntervalStart( );

    //! Get end of integration interval.
    /*!
     * This function returns the value of the end of the integration
     * interval.
     * \return End of integration interval.
     */
    double getIntegrationIntervalEnd( );

    //! Get initial state vector.
    /*!
     * This function returns address of the initial state vector.
     * \return Initial state vector.
     */
    Vector& getInitialStateVector( );

    //! Get final state vector.
    /*!
     * This function returns address of the final state vector.
     * \return Final state vector.
     */
    Vector& getFinalStateVector( );

    //! Get integration history.
    /*!
     * This function returns a map of the integration history.
     * \return Integration history.
     */
    std::map < double, Vector >& getIntegrationHistory( );

    //! Integrate.
    /*!
     * This function is called to perform an integration.
     */
    virtual void integrate( ) = 0;

protected:

    //! Flag to indicate if integration history should be saved.
    /*!
     * Flag to indicate if integration history should be saved.
     */
    bool isSaveIntegrationHistory_;

    //! Dimension of state vector.
    /*!
     * Dimension of state vector.
     */
    unsigned int dimensionOfCurrentStateVector_;

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
     * Current point in integration interval.
     */
    double integrationIntervalCurrentPoint_;

    //! Stepsize of last step.
    /*!
     * Value of stepsize for last step of integration interval.
     */
    double lastStepStepsize_;

    //! Initial state vector.
    /*!
     * Vector containing initial state.
     */
    Vector initialStateVector_;

    //! State derivatives vector.
    /*!
     * Vector containing state derivatives.
     */
    Vector stateDerivativeVector_;

    //! Final state vector.
    /*!
     * Vector containing final state.
     */
    Vector finalStateVector_;

    //! Vector containing current state vectors.
    /*!
     * Vector containing current state vectors.
     */
     std::vector< Vector > vectorOfCurrentStateVectors_;

     //! A map of integration history.
     /*!
      * A map of integration history with current point in integration interval
      * taken as key.
      */
     std::map < double, Vector > integrationHistory_;

    //! Pointer to function of type: void functionName( Vector, Vector ).
    /*!
     * Pointer to function type for state derivatives.
     */
    pointerToVectorVectorTakingFunction pointerToVectorVectorTakingFunction_;

    //! Pointer to state derivative function within Propagator class.
    /*!
     * Pointer to function type for state derivatives.
     */
    pointerToStateDerivativeFunction pointerToStateDerivativeFunction_;

    //! Pointer to NumericalPropagator class.
    /*!
     * Pointer to NumericalPropagator class.
     */
    NumericalPropagator* pointerToNumericalPropagator_;

    //! Compute state derivatives.
    /*!
     * This function computes state derivatives.
     * \param stateVector State vector.
     * \param stateDerivativeVector State derivative vector computed in
     * function.
     */
    void computeStateDerivatives_( Vector& stateVector,
                                   Vector& stateDerivativeVector );

    //! Store intermediate integration result.
    /*!
     * This function stores intermediate integration results obtained during
     * integration in integrationHistory.
    */
    void storeIntermediateIntegrationResult_( );

private:

    //! Flag to indicate if stepsize is set.
    /*!
     * Flag to indicate if stepsize is set.
     */
    bool isStepsizeSet_;

    //! Flag to indicate if initial state vector is set.
    /*!
     * Flag to indicate if initial state vector is set.
     */
    bool isInitialStateVectorSet_;

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

    //! Flag to indicate if state derivative function is set.
    /*!
     * Flag to indicate if state derivative function is set.
     */
    bool isStateDerivativeFunctionSet_;

    //! Compute internal derived integration parameters.
    /*!
     * Compute internal derived integration parameters necessary to perform
     * integration.
     */
    void computeInternalDerivedIntegrationParameters_( );
};

#endif // INTEGRATOR_H

// End of file.
