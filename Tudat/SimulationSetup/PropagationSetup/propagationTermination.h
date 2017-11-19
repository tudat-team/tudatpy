/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
#define TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/SimulationSetup/PropagationSetup/propagationOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace propagators
{

//! Possible events that can trigger the termination of a propagation
enum PropagationTerminationReason
{
    propagation_never_run,
    unknown_propagation_termination_reason,
    termination_condition_reached,
    runtime_error_caught_in_propagation
};

//! Base class for checking whether the numerical propagation is to be stopped at current time step or not
/*!
 *  Base class for checking whether the numerical propagation is to be stopped at current time step or not. Derived
 *  classes implement the various types of conditions (and associated threshold values) under which the propagation is to
 *  be stopped.
 */
class PropagationTerminationCondition
{
public:

    //! Constructor
    PropagationTerminationCondition(
            const PropagationTerminationTypes terminationType,
            const bool terminateExactlyOnFinalCondition = false,
            const double terminationTolerance = TUDAT_NAN ):
        terminationType_( terminationType ){ }

    //! Destructor
    virtual ~PropagationTerminationCondition( ){ }

    //! (Pure virtual) function to check whether the propagation should be stopped
    /*!
     * (Pure virtual) function to check whether the propagation should be stopped. Note that the accelerations and
     * environment must be updated (done automatically during numerical propagation) to check the stopping condition.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    virtual bool checkStopCondition( const double time, const double cpuTime ) = 0;

    virtual PropagationTerminationTypes getTerminationType( )
    {
        return terminationType_;
    }

protected:

    PropagationTerminationTypes terminationType_;

};

//! Class for stopping the propagation after a fixed amount of time (i.e. for certain independent variable value)
class FixedTimePropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stopTime Time at which the propagation is to stop.
     * \param propagationDirectionIsPositive Boolean denoting whether propagation is forward (if true) or backwards
     * (if false) in time.
     */
    FixedTimePropagationTerminationCondition(
            const double stopTime,
            const bool propagationDirectionIsPositive,
            const bool terminateExactlyOnFinalCondition = false,
            const double terminationTolerance = TUDAT_NAN  ):
        PropagationTerminationCondition( time_stopping_condition, terminateExactlyOnFinalCondition, terminationTolerance ),
        stopTime_( stopTime ),
        propagationDirectionIsPositive_( propagationDirectionIsPositive ){ }


    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the stopTime_ has been reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

private:

    //! Time at which the propagation is to stop.
    double stopTime_;

    //!  Boolean denoting whether propagation is forward (if true) or backwards (if false) in time.
    bool propagationDirectionIsPositive_;
};

//! Class for stopping the propagation after a fixed amount of CPU time
class FixedCPUTimePropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param cpuStopTime CPU time at which the propagation is to stop.
     */
    FixedCPUTimePropagationTerminationCondition( const double cpuStopTime ) :
        PropagationTerminationCondition( cpu_time_stopping_condition, false, TUDAT_NAN ),
        cpuStopTime_( cpuStopTime ) { }


    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the stopTime_ has been reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

private:

    //! Time at which the propagation is to stop.
    double cpuStopTime_;

};

//! Class for stopping the propagation when a dependent variable reaches a given value (either upper or lower bound)
class SingleVariableLimitPropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariableSettings Settings for dependent variable that is to be checked
     * \param variableRetrievalFuntion Function returning the dependent variable.
     * \param limitingValue Value at which the propagation is to be stopped
     * \param useAsLowerBound Boolean denoting whether the propagation should stop if the dependent variable goes below
     * (if true) or above (if false) limitingValue
     */
    SingleVariableLimitPropagationTerminationCondition(
            const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const boost::function< double( ) > variableRetrievalFuntion,
            const double limitingValue,
            const bool useAsLowerBound,
            const bool terminateExactlyOnFinalCondition = false,
            const double terminationTolerance = TUDAT_NAN,
            const boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = NULL ):
        PropagationTerminationCondition(
            dependent_variable_stopping_condition, terminateExactlyOnFinalCondition, terminationTolerance ),
        dependentVariableSettings_( dependentVariableSettings ), variableRetrievalFuntion_( variableRetrievalFuntion ),
        limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ),
    terminationRootFinderSettings_( terminationRootFinderSettings ){ }

    //! Destructor.
    ~SingleVariableLimitPropagationTerminationCondition( ){ }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the given dependent variable has been
     * reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    double getStopConditionError( )
    {
         return variableRetrievalFuntion_( ) - limitingValue_;
    }

    boost::shared_ptr< root_finders::RootFinderSettings > getTerminationRootFinderSettings( )
    {
        return terminationRootFinderSettings_;
    }

private:

    //! Settings for dependent variable that is to be checked
    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Function returning the dependent variable.
    boost::function< double( ) > variableRetrievalFuntion_;

    //! Value at which the propagation is to be stopped
    double limitingValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below
    //! (if true) or above (if false) limitingValue
    bool useAsLowerBound_;

    boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings_;
};

//! Class for stopping the propagation when one or all of a given set of stopping conditions is reached.
class HybridPropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param propagationTerminationCondition List of termination conditions that are checked when calling
     * checkStopCondition is called.
     * \param fulFillSingleCondition Boolean denoting whether a single (if true) or all (if false) of the entries in the
     * propagationTerminationCondition_ should return true from the checkStopCondition function to stop the propagation.
     */
    HybridPropagationTerminationCondition(
            const std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition,
            const bool fulFillSingleCondition = 0 ):
        PropagationTerminationCondition( hybrid_stopping_condition, false, TUDAT_NAN ),
        propagationTerminationCondition_( propagationTerminationCondition ),
        fulFillSingleCondition_( fulFillSingleCondition ){ }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. one or all (depending on value of
     * fulFillSingleCondition_) of the stopping conditions are fulfilled.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    std::vector< boost::shared_ptr< PropagationTerminationCondition > > getPropagationTerminationConditions( )
    {
        return propagationTerminationCondition_;
    }

    bool getFulFillSingleCondition( )
    {
        return fulFillSingleCondition_;
    }

private:

    //! List of termination conditions that are checked when calling checkStopCondition is called.
    std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition_;

    //!  Boolean denoting whether a single (if true) or all (if false) of the entries in the propagationTerminationCondition_
    //!  should return true from the checkStopCondition function to stop the propagation.
    bool fulFillSingleCondition_;
};

//! Function to create propagation termination conditions from associated settings
/*!
 * Function to create propagation termination conditions from associated settings
 * \param terminationSettings Settings for propagation termination conditions
 * \param bodyMap List of body objects that contains all environment models
 * \param initialTimeStep Time step at first call of numerical integration.
 * \return Object used to check whether propagation is to be stopped or not.
 */
boost::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep );

} // namespace propagators

} // namespace tudat


#endif // TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
