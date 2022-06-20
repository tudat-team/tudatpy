/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <memory>
#include <boost/lexical_cast.hpp>

#include "tudat/basics/timeType.h"
#include "tudat/math/integrators/bulirschStoerVariableStepsizeIntegrator.h"
#include "tudat/math/integrators/numericalIntegrator.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/math/integrators/rungeKuttaFixedStepSizeIntegrator.h"
#include "tudat/math/integrators/euler.h"
#include "tudat/math/integrators/adamsBashforthMoultonIntegrator.h"
#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"

namespace tudat
{

namespace numerical_integrators
{

// Enum to define available integrators.
//! @get_docstring(AvailableIntegrators.__docstring__)
enum AvailableIntegrators
{
    euler,
    rungeKutta4,
    rungeKuttaFixedStepSize,
    rungeKuttaVariableStepSize,
    bulirschStoer,
    adamsBashforthMoulton
};

// Class to define settings of numerical integrator
/*
 *  Class to define settings of numerical integrator, for instance for use in numerical integration of equations of motion/
 *  variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
 *  require more settings to define have their own derived class (see below).
 */
template< typename IndependentVariableType = double >
class IntegratorSettings
{
public:

    // Constructor
    /*
     *  Constructor for integrator settings.
     *  \param integratorType Type of numerical integrator
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration
     *  for variable step size integrators.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *  each integration step (`false`).
     */
    IntegratorSettings( const AvailableIntegrators integratorType, const IndependentVariableType initialTime,
                        const IndependentVariableType initialTimeStep,
                        const int saveFrequency = 1,
                        const bool assessTerminationOnMinorSteps = false ) :
        integratorType_( integratorType ), initialTime_( initialTime ),
        initialTimeStep_( initialTimeStep ), saveFrequency_( saveFrequency ),
        assessTerminationOnMinorSteps_( assessTerminationOnMinorSteps )
    { }
    
    // Virtual destructor.
    /*
     *  Virtual destructor.
     */
    virtual ~IntegratorSettings( ) { }

    // Type of numerical integrator
    /*
     *  Type of numerical integrator, from enum of available integrators.
     */
    AvailableIntegrators integratorType_;

    // Start time of numerical integration.
    /*
     *  Start time (independent variable) of numerical integration.
     */
    IndependentVariableType initialTime_;

    // Initial time step used in numerical integration
    /*
     *  Initial time (independent variable) step used in numerical integration. Adapted during integration
     *  for variable step size integrators.
     */
    IndependentVariableType initialTimeStep_;

    // Frequency which with to save numerical integration result.
    /*
     *  Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     */
    int saveFrequency_;

    // Whether the propagation termination conditions should be evaluated during the intermediate sub-steps.
    /*
     * Whether the propagation termination conditions should be evaluated during the intermediate updates
     * performed by the integrator to compute the quantities necessary to integrate the state to a new epoch.
     * The default value is false, so the propagation termination condition is only checked at the end of each
     * integration step.
     */
    bool assessTerminationOnMinorSteps_;

};

// Class to define settings of fixed step RK numerical integrator.
/*
 *  Class to define settings of fixed step RK numerical integrator, that can take distinct Butcher tableaus.
 */
template< typename IndependentVariableType = double >
class RungeKuttaFixedStepSizeSettings: public IntegratorSettings< IndependentVariableType >
{
public:

    // Default constructor.
    /*
     *  Constructor for fixed step RK integrator base settings.
     *  \param integratorType Type of numerical integrator
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration
     *  for variable step size integrators.
     *  \param orderToUse Order of Butcher tableau to use (only used if variable step intragration coefficients are specified).
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *  time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *  conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *  each integration step (`false`).
     *  \param butcherTableau Runge-Kutta tableau to be used for fixed step RK integrator.
     */
    RungeKuttaFixedStepSizeSettings(
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const numerical_integrators::CoefficientSets coefficientSet,
            const RungeKuttaCoefficients::OrderEstimateToIntegrate orderToUse = RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false ):
        IntegratorSettings< IndependentVariableType >(
            rungeKuttaFixedStepSize, initialTime, initialTimeStep, saveFrequency,
            assessTerminationOnMinorSteps ),
            coefficientSet_(coefficientSet),
            orderToUse_(orderToUse)
    { }

    // Virtual destructor.
    /*
     *  Virtual destructor.
     */
    virtual ~RungeKuttaFixedStepSizeSettings( ) { }

    // Type of numerical integrator (must be an RK fixed step type).
    numerical_integrators::CoefficientSets coefficientSet_;

    // Order of Runge-Kutta method to be used.
    RungeKuttaCoefficients::OrderEstimateToIntegrate orderToUse_;
};

// Base class to define settings of variable step RK numerical integrator.
/*
 *  Base class to define settings of variable step RK numerical integrator. From this class, two classes are derived, to
 *  define the relative and absolute error tolerances. These can be defined as a scalar or as a vector.
 */
template< typename IndependentVariableType = double >
class RungeKuttaVariableStepSizeBaseSettings: public IntegratorSettings< IndependentVariableType >
{
public:

    // Default constructor.
    /*
     *  Constructor for variable step RK integrator base settings.
     *  \param areTolerancesDefinedAsScalar Boolean denoting whether the relative and absolute error tolerances are
     *      defined as a scalar. Alternatively, they can be defined as a vector.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param coefficientSet Coefficient set (butcher tableau) to use in integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control.
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Minimum decrease factor in time step in subsequent iterations.
     */
    RungeKuttaVariableStepSizeBaseSettings(
            const bool areTolerancesDefinedAsScalar,
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const numerical_integrators::CoefficientSets coefficientSet,
            const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType safetyFactorForNextStepSize = 0.8,
            const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
            const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
            const bool exceptionIfMinimumStepExceeded = true  ) :
        IntegratorSettings< IndependentVariableType >(
            rungeKuttaVariableStepSize, initialTime, initialTimeStep, saveFrequency,
            assessTerminationOnMinorSteps ),
        areTolerancesDefinedAsScalar_( areTolerancesDefinedAsScalar ), coefficientSet_( coefficientSet ),
        minimumStepSize_( minimumStepSize ), maximumStepSize_( maximumStepSize ),
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        maximumFactorIncreaseForNextStepSize_( maximumFactorIncreaseForNextStepSize ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ),
        exceptionIfMinimumStepExceeded_( exceptionIfMinimumStepExceeded )
    { }

    // Virtual destructor.
    /*
     *  Virtual destructor.
     */
    virtual ~RungeKuttaVariableStepSizeBaseSettings( ) { }

    // Boolean denoting whether integration error tolerances are defined as a scalar (or vector).
    bool areTolerancesDefinedAsScalar_;

    // Type of numerical integrator (must be an RK variable step type).
    numerical_integrators::CoefficientSets coefficientSet_;

    // Minimum step size for integration.
    /*
     *  Minimum step size for integration. Integration stops (exception thrown) if time step comes below this value.
     */
    IndependentVariableType minimumStepSize_;

    // Maximum step size for integration.
    IndependentVariableType maximumStepSize_;

    // Safety factor for step size control
    IndependentVariableType safetyFactorForNextStepSize_;

    // Maximum increase factor in time step in subsequent iterations.
    IndependentVariableType maximumFactorIncreaseForNextStepSize_;

    // Minimum decrease factor in time step in subsequent iterations.
    IndependentVariableType minimumFactorDecreaseForNextStepSize_;

    bool exceptionIfMinimumStepExceeded_;

};

// Class to define settings of variable step RK numerical integrator with scalar tolerances.
/*
 *  Class to define settings of variable step RK numerical integrator with scalar tolerances.
 */
template< typename IndependentVariableType = double >
class RungeKuttaVariableStepSizeSettingsScalarTolerances: public RungeKuttaVariableStepSizeBaseSettings< IndependentVariableType >
{
public:

    // Default constructor.
    /*
     *  Constructor for variable step RK integrator settings with scalar tolerances.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param coefficientSet Coefficient set (butcher tableau) to use in integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control, expressed as a scalar.
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control, expressed as a scalar.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control.
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Minimum decrease factor in time step in subsequent iterations.
     */
    RungeKuttaVariableStepSizeSettingsScalarTolerances(
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const numerical_integrators::CoefficientSets coefficientSet,
            const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
            const IndependentVariableType relativeErrorTolerance = 1.0E-12,
            const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType safetyFactorForNextStepSize = 0.8,
            const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
            const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
            const bool exceptionIfMinimumStepExceeded = true ) :
        RungeKuttaVariableStepSizeBaseSettings< IndependentVariableType >(
            true, initialTime, initialTimeStep, coefficientSet, minimumStepSize, maximumStepSize, saveFrequency,
            assessTerminationOnMinorSteps, safetyFactorForNextStepSize,
            maximumFactorIncreaseForNextStepSize, minimumFactorDecreaseForNextStepSize, exceptionIfMinimumStepExceeded ),
        relativeErrorTolerance_( relativeErrorTolerance ), absoluteErrorTolerance_( absoluteErrorTolerance )
    { }

    // Constructor.
    /*
     *  Constructor for variable step RK integrator settings with scalar tolerances (also requires the input of the integratorType,
     *  which has to be rungeKyttaVariableStepSize by definition).
     *  \param integratorType Type of numerical integrator.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param coefficientSet Coefficient set (butcher tableau) to use in integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control, expressed as a scalar.
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control, expressed as a scalar.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control.
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Minimum decrease factor in time step in subsequent iterations.
     */
    RungeKuttaVariableStepSizeSettingsScalarTolerances(
            const AvailableIntegrators integratorType,
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const numerical_integrators::CoefficientSets coefficientSet,
            const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
            const IndependentVariableType relativeErrorTolerance = 1.0E-12,
            const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType safetyFactorForNextStepSize = 0.8,
            const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
            const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
            const bool exceptionIfMinimumStepExceeded = true ) :
        RungeKuttaVariableStepSizeSettingsScalarTolerances(
            initialTime, initialTimeStep, coefficientSet, minimumStepSize, maximumStepSize,
            relativeErrorTolerance, absoluteErrorTolerance, saveFrequency,
            assessTerminationOnMinorSteps, safetyFactorForNextStepSize,
            maximumFactorIncreaseForNextStepSize, minimumFactorDecreaseForNextStepSize,
            exceptionIfMinimumStepExceeded )
    {
        // Give error if integrator type is wrong
        if ( integratorType != rungeKuttaVariableStepSize )
        {
            throw std::runtime_error( "Error while creating numerical integrator. The integrator settings are of class "
                                      "RungeKuttaVariableStepSizeSettingsScalarTolerances, but the input integrator type is not of type "
                                      "rungeKuttaVariableStepSize." );
        }

        // Warn of different constructor
        std::cerr << "Warning in numerical integrator. The Runge-Kutta variable step size integrator settings no longer require the "
                     "input of the integrator type, i.e., rungeKuttaVariableStepSize. You can thus safely remove the first input, especially "
                     "if you do not want this warning to show up again. Note that this constructor will be removed in a "
                     "future release." << std::endl;
    }

    // Destructor.
    /*
     *  Destructor.
     */
    ~RungeKuttaVariableStepSizeSettingsScalarTolerances( ) { }

    // Relative error tolerance for step size control.
    IndependentVariableType relativeErrorTolerance_;

    // Absolute error tolerance for step size control.
    IndependentVariableType absoluteErrorTolerance_;

};

// Alias for variable step RK numerical integrator with scalar tolerances (added for compatibility with old code).
template< typename IndependentVariableType = double >
using RungeKuttaVariableStepSizeSettings = RungeKuttaVariableStepSizeSettingsScalarTolerances< IndependentVariableType >;

// Class to define settings of variable step RK numerical integrator with vector tolerances.
/*
 *  Class to define settings of variable step RK numerical integrator with vector tolerances.
 */
template< typename IndependentVariableType = double >
class RungeKuttaVariableStepSizeSettingsVectorTolerances: public RungeKuttaVariableStepSizeBaseSettings< IndependentVariableType >
{
public:

    typedef Eigen::Matrix< IndependentVariableType, Eigen::Dynamic, Eigen::Dynamic > DependentVariableType;
    // Default constructor.
    /*
     *  Constructor for variable step RK integrator settings with vector tolerances.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param coefficientSet Coefficient set (butcher tableau) to use in integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control, expressed as a vector.
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control, expressed as a vector.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control.
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Minimum decrease factor in time step in subsequent iterations.
     */
    RungeKuttaVariableStepSizeSettingsVectorTolerances(
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const numerical_integrators::CoefficientSets coefficientSet,
            const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
            const DependentVariableType relativeErrorTolerance,
            const DependentVariableType absoluteErrorTolerance,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType safetyFactorForNextStepSize = 0.8,
            const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
            const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
            const bool exceptionIfMinimumStepExceeded = true ) :
        RungeKuttaVariableStepSizeBaseSettings< IndependentVariableType >(
            false, initialTime, initialTimeStep, coefficientSet, minimumStepSize, maximumStepSize, saveFrequency,
            assessTerminationOnMinorSteps, safetyFactorForNextStepSize,
            maximumFactorIncreaseForNextStepSize, minimumFactorDecreaseForNextStepSize, exceptionIfMinimumStepExceeded ),
        relativeErrorTolerance_( relativeErrorTolerance ), absoluteErrorTolerance_( absoluteErrorTolerance )
    { }

    // Destructor.
    /*
     *  Destructor.
     */
    ~RungeKuttaVariableStepSizeSettingsVectorTolerances( ) { }

    // Relative error tolerance for step size control.
    DependentVariableType relativeErrorTolerance_;

    // Absolute error tolerance for step size control.
    DependentVariableType absoluteErrorTolerance_;

};

template< typename IndependentVariableType = double >
class BulirschStoerIntegratorSettings: public IntegratorSettings< IndependentVariableType >
{
public:

    // Constructor.
    /*
     *  Constructor for variable step RK integrator settings.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param extrapolationSequence Type of sequence that is to be used for Bulirsch-Stoer integrator.
     *  \param maximumNumberOfSteps Number of entries in the sequence, e.g. number of integrations used for a single
     *      extrapolation.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control.
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control.
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param safetyFactorForNextStepSize Safety factor for step size control.
     *  \param maximumFactorIncreaseForNextStepSize Maximum increase factor in time step in subsequent iterations.
     *  \param minimumFactorDecreaseForNextStepSize Minimum decrease factor in time step in subsequent iterations.
     */
    BulirschStoerIntegratorSettings(
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const ExtrapolationMethodStepSequences extrapolationSequence,
            const unsigned int maximumNumberOfSteps,
            const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
            const IndependentVariableType relativeErrorTolerance = 1.0E-12,
            const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType safetyFactorForNextStepSize = 0.7,
            const IndependentVariableType maximumFactorIncreaseForNextStepSize = 10.0,
            const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1 ):
        IntegratorSettings< IndependentVariableType >(
            bulirschStoer, initialTime, initialTimeStep, saveFrequency,
            assessTerminationOnMinorSteps ),
        extrapolationSequence_( extrapolationSequence ), maximumNumberOfSteps_( maximumNumberOfSteps ),
        minimumStepSize_( minimumStepSize ), maximumStepSize_( maximumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ), absoluteErrorTolerance_( absoluteErrorTolerance ),
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        maximumFactorIncreaseForNextStepSize_( maximumFactorIncreaseForNextStepSize ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ){ }

    // Destructor.
    /*
     *  Destructor.
     */
    ~BulirschStoerIntegratorSettings( ){ }

    // Type of sequence that is to be used for Bulirsch-Stoer integrator
    ExtrapolationMethodStepSequences extrapolationSequence_;

    // Number of entries in the sequence, e.g. number of integrations used for a single extrapolation.
    unsigned int maximumNumberOfSteps_;

    // Minimum step size for integration.
    /*
     *  Minimum step size for integration. Integration stops (exception thrown) if time step comes below this value.
     */
    const IndependentVariableType minimumStepSize_;

    // Maximum step size for integration.
    const IndependentVariableType maximumStepSize_;

    // Relative error tolerance for step size control
    const IndependentVariableType relativeErrorTolerance_;

    // Absolute error tolerance for step size control
    const IndependentVariableType absoluteErrorTolerance_;

    // Safety factor for step size control
    const IndependentVariableType safetyFactorForNextStepSize_;

    // Maximum increase factor in time step in subsequent iterations.
    const IndependentVariableType maximumFactorIncreaseForNextStepSize_;

    // Minimum decrease factor in time step in subsequent iterations.
    const IndependentVariableType minimumFactorDecreaseForNextStepSize_;

};

// Class to define settings of variable step ABAM numerical integrator
/*
 *  Class to define settings of variable step ABAM  numerical integrator,
 *  for instance for use in numerical integration of equations of motion/
 *  variational equations.
 */
template< typename IndependentVariableType = double >
class AdamsBashforthMoultonSettings: public IntegratorSettings< IndependentVariableType >
{
public:

    // Constructor
    /*
     *  Constructor for variable step RK integrator settings.
     *  \param initialTime Start time (independent variable) of numerical integration.
     *  \param initialTimeStep Initial time (independent variable) step used in numerical integration. Adapted during integration.
     *  \param minimumStepSize Minimum step size for integration. Integration stops (exception thrown) if time step
     *      comes below this value.
     *  \param maximumStepSize Maximum step size for integration.
     *  \param relativeErrorTolerance Relative error tolerance for step size control.
     *  \param absoluteErrorTolerance Absolute error tolerance for step size control.
     *  \param minimumOrder Minimum order of integrator (default 6).
     *  \param maximumOrder Maximum order of integrator (default 11).
     *  \param saveFrequency Frequency at which to save the numerical integrated states (in units of i.e. per n integration
     *      time steps, with n = saveFrequency).
     *  \param assessTerminationOnMinorSteps Whether the propagation termination
     *      conditions should be evaluated during the intermediate sub-steps of the integrator (`true`) or only at the end of
     *      each integration step (`false`).
     *  \param bandwidth Maximum error factor for doubling the stepsize (default: 200).
     */
    AdamsBashforthMoultonSettings(
            const IndependentVariableType initialTime,
            const IndependentVariableType initialTimeStep,
            const IndependentVariableType minimumStepSize,
            const IndependentVariableType maximumStepSize,
            const IndependentVariableType relativeErrorTolerance = 1.0E-12,
            const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
            const int minimumOrder = 6,
            const int maximumOrder = 11,
            const int saveFrequency = 1,
            const bool assessTerminationOnMinorSteps = false,
            const IndependentVariableType bandwidth = 200. ):
        IntegratorSettings< IndependentVariableType >(
            adamsBashforthMoulton, initialTime, initialTimeStep, saveFrequency,
            assessTerminationOnMinorSteps ),
        minimumStepSize_( minimumStepSize ), maximumStepSize_( maximumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ), absoluteErrorTolerance_( absoluteErrorTolerance ),
        minimumOrder_( minimumOrder ), maximumOrder_( maximumOrder ),
        bandwidth_( bandwidth ) { }

    // Destructor
    /*
     *  Destructor
     */
    ~AdamsBashforthMoultonSettings( ){ }

    // Minimum step size for integration.
    /*
     *  Minimum step size for integration. Integration stops (exception thrown) if time step comes below this value.
     */
    IndependentVariableType minimumStepSize_;

    // Maximum step size for integration.
    IndependentVariableType maximumStepSize_;

    // Relative error tolerance for step size control
    IndependentVariableType relativeErrorTolerance_;

    // Absolute error tolerance for step size control
    IndependentVariableType absoluteErrorTolerance_;

    // Minimum order of integrator
    const int minimumOrder_;

    // Maximum order of integrator
    const int maximumOrder_;

    // Safety factor for step size control
    IndependentVariableType bandwidth_;

};

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > eulerSettings(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false )
{
    return std::make_shared< IntegratorSettings< IndependentVariableType > >(
                euler, initialTime, initialTimeStep, saveFrequency, assessTerminationOnMinorSteps );
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > rungeKutta4Settings(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false )
{
    return std::make_shared< IntegratorSettings< IndependentVariableType > >(
                rungeKutta4, initialTime, initialTimeStep, saveFrequency, assessTerminationOnMinorSteps );
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > rungeKuttaFixedStepSettings(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const numerical_integrators::CoefficientSets coefficientSet,
        const RungeKuttaCoefficients::OrderEstimateToIntegrate orderToUse = RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false )
{
    return std::make_shared< RungeKuttaFixedStepSizeSettings< IndependentVariableType > >(
                initialTime, initialTimeStep, coefficientSet, orderToUse, saveFrequency, assessTerminationOnMinorSteps );
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > rungeKuttaVariableStepSettingsScalarTolerances(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const numerical_integrators::CoefficientSets coefficientSet,
        const IndependentVariableType minimumStepSize,
        const IndependentVariableType maximumStepSize,
        const IndependentVariableType& relativeErrorTolerance,
        const IndependentVariableType& absoluteErrorTolerance,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false,
        const IndependentVariableType safetyFactorForNextStepSize = 0.8,
        const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
        const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
        const bool exceptionIfMinimumStepExceeded = true )
{
    return std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances<
            IndependentVariableType > >(
                initialTime, initialTimeStep,
                coefficientSet, minimumStepSize, maximumStepSize,
                relativeErrorTolerance, absoluteErrorTolerance,
                saveFrequency, assessTerminationOnMinorSteps, safetyFactorForNextStepSize,
                maximumFactorIncreaseForNextStepSize, minimumFactorDecreaseForNextStepSize,
                exceptionIfMinimumStepExceeded );
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > rungeKuttaVariableStepSettingsVectorTolerances(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const numerical_integrators::CoefficientSets coefficientSet,
        const IndependentVariableType minimumStepSize,
        const IndependentVariableType maximumStepSize,
        const Eigen::Matrix< IndependentVariableType, Eigen::Dynamic, Eigen::Dynamic > relativeErrorTolerance,
        const Eigen::Matrix< IndependentVariableType, Eigen::Dynamic, Eigen::Dynamic > absoluteErrorTolerance,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false,
        const IndependentVariableType safetyFactorForNextStepSize = 0.8,
        const IndependentVariableType maximumFactorIncreaseForNextStepSize = 4.0,
        const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1,
        const bool exceptionIfMinimumStepExceeded = true )
{
    auto settings = std::make_shared< RungeKuttaVariableStepSizeSettingsVectorTolerances<
            IndependentVariableType > >(
                initialTime, initialTimeStep,
                coefficientSet, minimumStepSize, maximumStepSize,
                relativeErrorTolerance, absoluteErrorTolerance,
                saveFrequency, assessTerminationOnMinorSteps, safetyFactorForNextStepSize,
                maximumFactorIncreaseForNextStepSize, minimumFactorDecreaseForNextStepSize,
                exceptionIfMinimumStepExceeded );

    return settings;
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > bulirschStoerIntegratorSettings(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const ExtrapolationMethodStepSequences extrapolationSequence,
        const unsigned int maximumNumberOfSteps,
        const IndependentVariableType minimumStepSize, const IndependentVariableType maximumStepSize,
        const IndependentVariableType relativeErrorTolerance = 1.0E-12,
        const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false,
        const IndependentVariableType safetyFactorForNextStepSize = 0.7,
        const IndependentVariableType maximumFactorIncreaseForNextStepSize = 10.0,
        const IndependentVariableType minimumFactorDecreaseForNextStepSize = 0.1 )
{
    return std::make_shared< BulirschStoerIntegratorSettings< IndependentVariableType > >(
                initialTime, initialTimeStep,
                extrapolationSequence, maximumNumberOfSteps,
                minimumStepSize, maximumStepSize,
                relativeErrorTolerance, absoluteErrorTolerance,
                saveFrequency,  assessTerminationOnMinorSteps,
                safetyFactorForNextStepSize,
                maximumFactorIncreaseForNextStepSize,
                minimumFactorDecreaseForNextStepSize );
}

template< typename IndependentVariableType = double >
inline std::shared_ptr< IntegratorSettings< IndependentVariableType > > adamsBashforthMoultonSettings(
        const IndependentVariableType initialTime,
        const IndependentVariableType initialTimeStep,
        const IndependentVariableType minimumStepSize,
        const IndependentVariableType maximumStepSize,
        const IndependentVariableType relativeErrorTolerance = 1.0E-12,
        const IndependentVariableType absoluteErrorTolerance = 1.0E-12,
        const int minimumOrder = 6,
        const int maximumOrder = 11,
        const int saveFrequency = 1,
        const bool assessTerminationOnMinorSteps = false,
        const IndependentVariableType bandwidth = 200. )
{
    return std::make_shared< AdamsBashforthMoultonSettings< IndependentVariableType > >(
                initialTime, initialTimeStep,
                minimumStepSize, maximumStepSize,
                relativeErrorTolerance, absoluteErrorTolerance,
                minimumOrder, maximumOrder,
                saveFrequency,assessTerminationOnMinorSteps, bandwidth );
}

// Function to create a numerical integrator.
/*
 *  Function to create a numerical integrator from given integrator settings, state derivative function and initial state.
 *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
 *  \param initialState Initial state for numerical integration.
 *  \param integratorSettings Settings for numerical integrator.
 *  \return Numerical integrator object.
 */
template< typename IndependentVariableType, typename DependentVariableType,
          typename IndependentVariableStepType = IndependentVariableType >
std::shared_ptr< numerical_integrators::NumericalIntegrator< IndependentVariableType, DependentVariableType,
DependentVariableType, IndependentVariableStepType > > createIntegrator(
        std::function< DependentVariableType(
            const IndependentVariableType, const DependentVariableType& ) > stateDerivativeFunction,
        const DependentVariableType initialState,
        std::shared_ptr< IntegratorSettings< IndependentVariableType > > integratorSettings )
{
    // Declare eventual output
    std::shared_ptr< NumericalIntegrator
            < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > > integrator;

    // Retrieve requested type of integrator
    switch( integratorSettings->integratorType_ )
    {
    case euler:
    {
        // Create Euler integrator
        integrator = std::make_shared< EulerIntegrator
                < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                ( stateDerivativeFunction, integratorSettings->initialTime_, initialState ) ;
        break;
    }
    case rungeKutta4:
    {
        // Create Runge-Kutta 4 integrator
        integrator = std::make_shared< RungeKutta4Integrator
                < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                ( stateDerivativeFunction, integratorSettings->initialTime_, initialState ) ;
        break;
    }
    case rungeKuttaFixedStepSize:
    {
        // Cast integrator
        std::shared_ptr< RungeKuttaFixedStepSizeSettings< IndependentVariableType > >
                fixedStepIntegratorSettings = std::dynamic_pointer_cast< RungeKuttaFixedStepSizeSettings<
                IndependentVariableType > >( integratorSettings );

        // Create Runge-Kutta fixed step integrator
        integrator = std::make_shared< RungeKuttaFixedStepSizeIntegrator
                < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                ( stateDerivativeFunction, fixedStepIntegratorSettings->initialTime_, initialState, fixedStepIntegratorSettings->coefficientSet_, fixedStepIntegratorSettings->orderToUse_) ;
        break;
    }
    case rungeKuttaVariableStepSize:
    {
        // Cast integrator
        std::shared_ptr< RungeKuttaVariableStepSizeBaseSettings< IndependentVariableType > >
                variableStepIntegratorSettings = std::dynamic_pointer_cast< RungeKuttaVariableStepSizeBaseSettings<
                IndependentVariableType > >( integratorSettings );

        if( std::fabs( static_cast<double>(variableStepIntegratorSettings->initialTimeStep_) ) <
                std::fabs( static_cast<double>(variableStepIntegratorSettings->minimumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making RK variable step-size integrator: initial step size is smaller than minimum step" );
        }


        if( std::fabs( static_cast<double>(variableStepIntegratorSettings->initialTimeStep_) ) >
                std::fabs( static_cast<double>(variableStepIntegratorSettings->maximumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making RK variable step-size integrator: initial step size is larger than maximum step" );
        }

        // Check input consistency
        if ( variableStepIntegratorSettings == nullptr )
        {
            throw std::runtime_error( "Error, type of integrator settings (rungeKuttaVariableStepSize) not compatible with "
                                      "selected integrator (derived class of IntegratorSettings must be "
                                      "RungeKuttaVariableStepSizeBaseSettings for this type)." );
        }

        // Get requested RK coefficients
        RungeKuttaCoefficients coefficients = RungeKuttaCoefficients::get( variableStepIntegratorSettings->coefficientSet_ );

        // Check which constructor is being used
        if ( variableStepIntegratorSettings->areTolerancesDefinedAsScalar_ )
        {
            // Settings with scalar tolerances
            std::shared_ptr< RungeKuttaVariableStepSizeSettingsScalarTolerances< IndependentVariableType > >
                    scalarTolerancesIntegratorSettings = std::dynamic_pointer_cast<
                    RungeKuttaVariableStepSizeSettingsScalarTolerances< IndependentVariableType > >( variableStepIntegratorSettings );

            // Check input consistency
            if ( scalarTolerancesIntegratorSettings == nullptr )
            {
                throw std::runtime_error( "Error while creating Runge-Kutta variable step size integrator. Input class must be of "
                                          "RungeKuttaVariableStepSizeSettingsScalarTolerances type." );
            }

            // Create Runge-Kutta integrator with scalar tolerances
            integrator = std::make_shared< RungeKuttaVariableStepSizeIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                    ( coefficients, stateDerivativeFunction, integratorSettings->initialTime_, initialState,
                      static_cast< IndependentVariableStepType >( scalarTolerancesIntegratorSettings->minimumStepSize_ ),
                      static_cast< IndependentVariableStepType >( scalarTolerancesIntegratorSettings->maximumStepSize_ ),
                      static_cast< typename DependentVariableType::Scalar >( scalarTolerancesIntegratorSettings->relativeErrorTolerance_ ),
                      static_cast< typename DependentVariableType::Scalar >( scalarTolerancesIntegratorSettings->absoluteErrorTolerance_ ),
                      static_cast< IndependentVariableStepType >( scalarTolerancesIntegratorSettings->safetyFactorForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( scalarTolerancesIntegratorSettings->maximumFactorIncreaseForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( scalarTolerancesIntegratorSettings->minimumFactorDecreaseForNextStepSize_ ),
                      nullptr, scalarTolerancesIntegratorSettings->exceptionIfMinimumStepExceeded_ );
        }
        else
        {
            // Settings with vector tolerances
            std::shared_ptr< RungeKuttaVariableStepSizeSettingsVectorTolerances< IndependentVariableType > >
                    vectorTolerancesIntegratorSettings = std::dynamic_pointer_cast<
                    RungeKuttaVariableStepSizeSettingsVectorTolerances< IndependentVariableType > >(
                        variableStepIntegratorSettings );

            // Check input consistency
            if ( vectorTolerancesIntegratorSettings == nullptr )
            {
                throw std::runtime_error( "Error while creating Runge-Kutta variable step size integrator. Input class must be of "
                                          "RungeKuttaVariableStepSizeSettingsVectorTolerances type with suitable independent variable type." );
            }

            // Check that sizes of tolerances and initial state match
            auto relativeErrorTolerance = vectorTolerancesIntegratorSettings->relativeErrorTolerance_.template cast<
                    typename DependentVariableType::Scalar >( );
            auto absoluteErrorTolerance = vectorTolerancesIntegratorSettings->absoluteErrorTolerance_.template cast<
                    typename DependentVariableType::Scalar >( );
            if ( ( relativeErrorTolerance.rows( ) != initialState.rows( ) ) ||
                 ( relativeErrorTolerance.cols( ) != initialState.cols( ) ) ||
                 ( absoluteErrorTolerance.rows( ) != initialState.rows( ) ) ||
                 ( absoluteErrorTolerance.cols( ) != initialState.cols( ) ) )
            {
                throw std::runtime_error( "Error while creating Runge-Kutta variable step size integrator. The sizes of the "
                                          "relative and absolute tolerance vectors do not match the size of the initial state. "
                                          "This could be the case if you are propagating more than just one state, e.g., translational "
                                          "and/or rotational dynamics and mass, or more than one body." );
            }

            // Create Runge-Kutta integrator with vector tolerances
            integrator = std::make_shared< RungeKuttaVariableStepSizeIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                    ( coefficients, stateDerivativeFunction, integratorSettings->initialTime_, initialState,
                      static_cast< IndependentVariableStepType >( vectorTolerancesIntegratorSettings->minimumStepSize_ ),
                      static_cast< IndependentVariableStepType >( vectorTolerancesIntegratorSettings->maximumStepSize_ ),
                      relativeErrorTolerance, absoluteErrorTolerance,
                      static_cast< IndependentVariableStepType >( vectorTolerancesIntegratorSettings->safetyFactorForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( vectorTolerancesIntegratorSettings->maximumFactorIncreaseForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( vectorTolerancesIntegratorSettings->minimumFactorDecreaseForNextStepSize_ ),
                      nullptr, vectorTolerancesIntegratorSettings->exceptionIfMinimumStepExceeded_  );
        }
        break;
    }
    case bulirschStoer:
    {
        // Check input consistency
        std::shared_ptr< BulirschStoerIntegratorSettings< IndependentVariableType > > bulirschStoerIntegratorSettings =
                std::dynamic_pointer_cast< BulirschStoerIntegratorSettings< IndependentVariableType > >(
                    integratorSettings );


        if( std::fabs( static_cast<double>(bulirschStoerIntegratorSettings->initialTimeStep_) ) <
                std::fabs( static_cast<double>(bulirschStoerIntegratorSettings->minimumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making BS variable step-size integrator: initial step size is smaller than minimum step" );
        }


        if( std::fabs( static_cast<double>(bulirschStoerIntegratorSettings->initialTimeStep_)  ) >
                std::fabs( static_cast<double>(bulirschStoerIntegratorSettings->maximumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making BS variable step-size integrator: initial step size is larger than maximum step" );
        }

        // Check that integrator type has been cast properly
        if ( bulirschStoerIntegratorSettings == nullptr )
        {
            throw std::runtime_error( "Error, type of integrator settings (bulirschStoer) not compatible with "
                                      "selected integrator (derived class of IntegratorSettings must be BulirschStoerIntegratorSettings "
                                      "for this type)." );
        }
        else
        {
            integrator = std::make_shared< BulirschStoerVariableStepSizeIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                    ( getBulirschStoerStepSequence( bulirschStoerIntegratorSettings->extrapolationSequence_,
                                                    bulirschStoerIntegratorSettings->maximumNumberOfSteps_ ),
                      stateDerivativeFunction, integratorSettings->initialTime_, initialState,
                      static_cast< IndependentVariableStepType >( bulirschStoerIntegratorSettings->minimumStepSize_ ),
                      static_cast< IndependentVariableStepType >( bulirschStoerIntegratorSettings->maximumStepSize_ ),
                      bulirschStoerIntegratorSettings->relativeErrorTolerance_,
                      bulirschStoerIntegratorSettings->absoluteErrorTolerance_,
                      static_cast< IndependentVariableStepType >( bulirschStoerIntegratorSettings->safetyFactorForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( bulirschStoerIntegratorSettings->maximumFactorIncreaseForNextStepSize_ ),
                      static_cast< IndependentVariableStepType >( bulirschStoerIntegratorSettings->minimumFactorDecreaseForNextStepSize_ ) );
        }
        break;
    }
    case adamsBashforthMoulton:
    {
        // Check input consistency
        std::shared_ptr< AdamsBashforthMoultonSettings< IndependentVariableType > > adamsBashforthMoultonIntegratorSettings =
                std::dynamic_pointer_cast< AdamsBashforthMoultonSettings< IndependentVariableType > >(
                    integratorSettings );

        if( std::fabs( static_cast<double>(adamsBashforthMoultonIntegratorSettings->initialTimeStep_) ) <
                std::fabs( static_cast<double>(adamsBashforthMoultonIntegratorSettings->minimumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making ABM variable step-size integrator: initial step size is smaller than minimum step" );
        }


        if( std::fabs( static_cast<double>(adamsBashforthMoultonIntegratorSettings->initialTimeStep_) )  >
            std::fabs( static_cast<double>(adamsBashforthMoultonIntegratorSettings->maximumStepSize_) ) )
        {
            throw std::runtime_error( "Error when making ABM variable step-size integrator: initial step size is larger than maximum step" );
        }

        // Check that integrator type has been cast properly
        if ( adamsBashforthMoultonIntegratorSettings == nullptr )
        {
            throw std::runtime_error( "Error, type of integrator settings (AdamsBashforthMoultonSettings) not compatible with "
                                      "selected integrator (derived class of IntegratorSettings must be AdamsBashforthMoultonSettings "
                                      "for this type)." );
        }
        else
        {
            // Create integrator
            integrator = std::make_shared< AdamsBashforthMoultonIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >
                    ( stateDerivativeFunction, integratorSettings->initialTime_, initialState,
                      static_cast< IndependentVariableStepType >( adamsBashforthMoultonIntegratorSettings->minimumStepSize_ ),
                      static_cast< IndependentVariableStepType >( adamsBashforthMoultonIntegratorSettings->maximumStepSize_ ),
                      adamsBashforthMoultonIntegratorSettings->relativeErrorTolerance_,
                      adamsBashforthMoultonIntegratorSettings->absoluteErrorTolerance_ ,
                      static_cast< IndependentVariableStepType >( adamsBashforthMoultonIntegratorSettings->bandwidth_ ) );

            std::dynamic_pointer_cast< AdamsBashforthMoultonIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >(
                        integrator )->setMinimumOrder( adamsBashforthMoultonIntegratorSettings->minimumOrder_ );
            std::dynamic_pointer_cast< AdamsBashforthMoultonIntegrator
                    < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >(
                        integrator )->setMaximumOrder( adamsBashforthMoultonIntegratorSettings->maximumOrder_ );

            if ( adamsBashforthMoultonIntegratorSettings->minimumOrder_ == adamsBashforthMoultonIntegratorSettings->maximumOrder_ )
            {
                std::dynamic_pointer_cast< AdamsBashforthMoultonIntegrator
                        < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >(
                            integrator )->setFixedOrder( true );
            }

            if ( adamsBashforthMoultonIntegratorSettings->minimumStepSize_ ==
                 adamsBashforthMoultonIntegratorSettings->maximumStepSize_ )
            {
                std::dynamic_pointer_cast< AdamsBashforthMoultonIntegrator
                        < IndependentVariableType, DependentVariableType, DependentVariableType, IndependentVariableStepType > >(
                            integrator )->setFixedStepSize( true );
            }
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, integrator " +  std::to_string( integratorSettings->integratorType_ ) + " not found." );
    }

    // Check that assignment of integrator went well
    if ( integrator == nullptr )
    {
        throw std::runtime_error( "Error while creating integrator. The resulting integrator pointer is null." );
    }

    // Give back integrator
    return integrator;
}


extern template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::VectorXd,
                                                                             Eigen::VectorXd, double > > createIntegrator< double, Eigen::VectorXd, double >(
        std::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) > stateDerivativeFunction,
        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

extern template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
                                                                             Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double > > createIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
            const double, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

extern template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::VectorXd,
                                                                             Eigen::VectorXd, long double > > createIntegrator< Time, Eigen::VectorXd, long double >(
        std::function< Eigen::VectorXd( const Time, const Eigen::VectorXd& ) > stateDerivativeFunction,
        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );

extern template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
                                                                             Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double > > createIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
            const Time, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );
} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_CREATENUMERICALINTEGRATOR_H
