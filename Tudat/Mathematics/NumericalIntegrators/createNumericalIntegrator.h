/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CREATENUMERICALINTEGRATOR_H
#define TUDAT_CREATENUMERICALINTEGRATOR_H

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/euler.h"

#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{

namespace numerical_integrators
{

//! Enum to define available integrators.
enum AvailableIntegrators
{
    rungeKutta4,
    euler,
    rungeKuttaVariableStepSize
};

//! Class to define settings of numerical integrator
/*!
 *  Class to define settings of numerical integrator, for instance for use in numerical integration of equations of motion/
 *  variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
 *  require more settings to define have their own derived class (see below).
 */
template< typename TimeType = double >
class IntegratorSettings
{
public:

    //! Constructor
    /*!
     *  Constructor for integrator settings.
     *  \param integratorType Type of numerical integrator
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration
     *  for variable step size integrators.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     *  \param assessPropagationTerminationConditionDuringIntegrationSubsteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *  each integration step (`false`).
     */
    IntegratorSettings( const AvailableIntegrators integratorType, const TimeType initialTime,
                        const TimeType initialTimeStep,
                        const int saveFrequency = 1,
                        const bool assessPropagationTerminationConditionDuringIntegrationSubsteps = false )
        : integratorType_( integratorType ),
        initialTime_( initialTime ), initialTimeStep_( initialTimeStep ),
        saveFrequency_( saveFrequency ), assessPropagationTerminationConditionDuringIntegrationSubsteps_(
                                             assessPropagationTerminationConditionDuringIntegrationSubsteps ){ }

    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~IntegratorSettings( ) { }

    //! Type of numerical integrator
    /*!
     *  Type of numerical integrator, from enum of available integrators.
     */
    AvailableIntegrators integratorType_;

    //! Start time of numerical integration.
    /*!
     *  Start time (independent variable) of numerical integration.
     */
    TimeType initialTime_;

    //! Initial time step used in numerical integration
    /*!
     *  Initial time (independent variable) step used in numerical integration. Adapted during integration
     *  for variable step size integrators.
     */
    TimeType initialTimeStep_;

    //! Frequency which with to save numerical integration result.
    /*!
     *  Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     */
    int saveFrequency_;

    //! Whether the propagation termination conditions should be evaluated during the intermediate sub-steps.
    /*!
     * Whether the propagation termination conditions should be evaluated during the intermediate updates
     * performed by the integrator to compute the quantities necessary to integrate the state to a new epoch.
     * The default value is false, so the propagation termination condition is only checked at the end of each
     * integration step.
     */
    bool assessPropagationTerminationConditionDuringIntegrationSubsteps_;


};

//! Class to define settings of variable step RK numerical integrator
/*!
 *  Class to define settings of variable step RK  numerical integrator, for instance for use in numerical integration of equations of motion/
 *  variational equations.
 */
template< typename TimeType = double >
class RungeKuttaVariableStepSizeSettings: public IntegratorSettings< TimeType >
{
public:

    //! Constructor
    /*!
     *  Constructor for variable step RK integrator settings.
     *  \param integratorType Type of numerical integrator (must be an RK variable step type)
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration.
     *  Adapted during integration
     *  \param coefficientSet Coefficient set (butcher tableau) to use in integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *  comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     *  \param assessPropagationTerminationConditionDuringIntegrationSubsteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *  each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Maximum decrease factor in time step in subsequent iterations.
     */
    RungeKuttaVariableStepSizeSettings(
            const AvailableIntegrators integratorType,
            const TimeType initialTime,
            const TimeType initialTimeStep,
            const numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet,
            const TimeType minimumStepSize, const TimeType maximumStepSize,
            const TimeType relativeErrorTolerance = 1.0E-12,
            const TimeType absoluteErrorTolerance = 1.0E-12,
            const int saveFrequency = 1,
            const bool assessPropagationTerminationConditionDuringIntegrationSubsteps = false,
            const TimeType safetyFactorForNextStepSize = 0.8,
            const TimeType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeType minimumFactorDecreaseForNextStepSize = 0.1 ):
        IntegratorSettings< TimeType >( integratorType, initialTime, initialTimeStep, saveFrequency,
                                        assessPropagationTerminationConditionDuringIntegrationSubsteps ),
        coefficientSet_( coefficientSet ), minimumStepSize_( minimumStepSize ), maximumStepSize_( maximumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ), absoluteErrorTolerance_( absoluteErrorTolerance ),
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        maximumFactorIncreaseForNextStepSize_( maximumFactorIncreaseForNextStepSize ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ){ }

    //! Destructor
    /*!
     *  Destructor
     */
    ~RungeKuttaVariableStepSizeSettings( ){ }

    //! Type of numerical integrator (must be an RK variable step type)
    numerical_integrators::RungeKuttaCoefficients::CoefficientSets coefficientSet_;

    //! Minimum step size for integration.
    /*!
     *  Minimum step size for integration. Integration stops (exception thrown) if time step comes below this value.
     */
    TimeType minimumStepSize_;

    //! Maximum step size for integration.
    TimeType maximumStepSize_;

    //! Relative error tolerance for step size control
    TimeType relativeErrorTolerance_;

    //! Absolute error tolerance for step size control
    TimeType absoluteErrorTolerance_;

    //! Safety factor for step size control
    TimeType safetyFactorForNextStepSize_;

    //! Maximum increase factor in time step in subsequent iterations.
    TimeType maximumFactorIncreaseForNextStepSize_;

    //! Maximum decrease factor in time step in subsequent iterations.
    TimeType minimumFactorDecreaseForNextStepSize_;
};

//! Function to create a numerical integrator.
/*!
 *  Function to create a numerical integrator from given integrator settings, state derivative function and initial state.
 *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
 *  \param initialState Initial state for numerical integration
 *  \param integratorSettings Settings for numerical integrator.
 *  \return Numerical integrator object
 */
template< typename IndependentVariableType, typename DependentVariableType, typename TimeStepType = IndependentVariableType >
boost::shared_ptr< numerical_integrators::NumericalIntegrator< IndependentVariableType, DependentVariableType,
DependentVariableType, TimeStepType > > createIntegrator(
        boost::function< DependentVariableType(
            const IndependentVariableType, const DependentVariableType& ) > stateDerivativeFunction,
        const DependentVariableType initialState,
        boost::shared_ptr< IntegratorSettings< IndependentVariableType > > integratorSettings )

{    
    boost::shared_ptr< NumericalIntegrator
            < IndependentVariableType, DependentVariableType, DependentVariableType, TimeStepType > > integrator;

    // Retrieve requested type of integrator
    switch( integratorSettings->integratorType_ )
    {
    case euler:
        integrator = boost::make_shared< EulerIntegrator
                < IndependentVariableType, DependentVariableType, DependentVariableType, TimeStepType > >
                ( stateDerivativeFunction, integratorSettings->initialTime_, initialState ) ;
        break;
    case rungeKutta4:

        integrator = boost::make_shared< RungeKutta4Integrator
                < IndependentVariableType, DependentVariableType, DependentVariableType, TimeStepType > >
                ( stateDerivativeFunction, integratorSettings->initialTime_, initialState ) ;
        break;


    case rungeKuttaVariableStepSize:
    {
        // Check input consistency
        boost::shared_ptr< RungeKuttaVariableStepSizeSettings< IndependentVariableType > > variableStepIntegratorSettings =
                boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< IndependentVariableType > >(
                    integratorSettings );
        if( variableStepIntegratorSettings == NULL )
        {
           std::runtime_error( "Error, type of integrator settings (rungeKuttaVariableStepSize) not compatible with selected integrator (derived class of IntegratorSettings must be RungeKuttaVariableStepSizeSettings for this type)" );
        }
        else
        {
            // Get requested RK coefficients and create integrator.
            RungeKuttaCoefficients coefficients =  RungeKuttaCoefficients::get(
                        variableStepIntegratorSettings->coefficientSet_ );
            integrator = boost::make_shared<
                    RungeKuttaVariableStepSizeIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, TimeStepType > >
                    ( coefficients,
                      stateDerivativeFunction, integratorSettings->initialTime_, initialState,
                      static_cast< TimeStepType >( variableStepIntegratorSettings->minimumStepSize_ ),
                      static_cast< TimeStepType >( variableStepIntegratorSettings->maximumStepSize_ ),
                      variableStepIntegratorSettings->relativeErrorTolerance_,
                      variableStepIntegratorSettings->absoluteErrorTolerance_ ,
                      static_cast< TimeStepType >( variableStepIntegratorSettings->safetyFactorForNextStepSize_ ),
                      static_cast< TimeStepType >( variableStepIntegratorSettings->maximumFactorIncreaseForNextStepSize_ ),
                      static_cast< TimeStepType >( variableStepIntegratorSettings->minimumFactorDecreaseForNextStepSize_ ) );
        }
        break;
    }
    default:
        std::runtime_error(
                    "Error, integrator " +  std::to_string( integratorSettings->integratorType_ ) +
                    "not found. " );    }
    return integrator;
}

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_CREATENUMERICALINTEGRATOR_H
