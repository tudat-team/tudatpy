/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <memory>

#include "tudat/simulation/propagation_setup/propagationOutput.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

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
    runtime_error_caught_in_propagation,
    nan_or_inf_detected_in_state
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
    /*!
     * Constructor
     * \param terminationType Type of termination condition
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTerminationCondition(
            const PropagationTerminationTypes terminationType,
            const bool checkTerminationToExactCondition = false ):
        terminationType_( terminationType ), checkTerminationToExactCondition_( checkTerminationToExactCondition ){ }

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

    //! Function to retrieve type of termination condition
    /*!
     *  Function to retrieve type of termination condition
     *  \return Type of termination condition
     */
    virtual PropagationTerminationTypes getTerminationType( )
    {
        return terminationType_;
    }

    //! Function to retrieve boolean to denote whether the propagation is to terminate exactly on the final condition
    /*!
     *  Function to retrieve boolean to denote whether the propagation is to terminate exactly on the final condition
     *  \return Boolean to denote whether the propagation is to terminate exactly on the final condition
     */
    bool getcheckTerminationToExactCondition( )
    {
        return checkTerminationToExactCondition_;
    }

protected:

    //! Type of termination condition
    PropagationTerminationTypes terminationType_;

    //! Boolean to denote whether the propagation is to terminate exactly on the final condition, or whether it is to terminate
    //! on the first step where it is violated.
    bool checkTerminationToExactCondition_;

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
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    FixedTimePropagationTerminationCondition(
            const double stopTime,
            const bool propagationDirectionIsPositive,
            const bool checkTerminationToExactCondition = false ):
        PropagationTerminationCondition( time_stopping_condition, checkTerminationToExactCondition ),
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

    //! Function to retrieve time at which the propagation is to stop.
    /*!
     *  Function to retrieve time at which the propagation is to stop.
     *  \return Type of termination condition
     */
    double getStopTime( )
    {
        return stopTime_;
    }

private:

    //! Time at which the propagation is to stop.
    double stopTime_;

private:

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
        PropagationTerminationCondition( cpu_time_stopping_condition, false ),
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
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     * \param terminationRootFinderSettings Settings to create root finder used to converge on exact final condition.
     */
    SingleVariableLimitPropagationTerminationCondition(
            const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const std::function< double( ) > variableRetrievalFuntion,
            const double limitingValue,
            const bool useAsLowerBound,
            const bool checkTerminationToExactCondition = false,
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = nullptr ):
        PropagationTerminationCondition(
            dependent_variable_stopping_condition, checkTerminationToExactCondition ),
        dependentVariableSettings_( dependentVariableSettings ), variableRetrievalFunction_( variableRetrievalFuntion ),
        limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ),
        terminationRootFinderSettings_( terminationRootFinderSettings )
    {
        if( ( checkTerminationToExactCondition == false ) && ( terminationRootFinderSettings != nullptr ) )
        {
            std::cerr << "Warning, root finder provided to SingleVariableLimitPropagationTerminationCondition, "
                         "but termination on final conditions set to false." << std::endl;
        }
        if( ( checkTerminationToExactCondition ) && doesRootFinderRequireDerivatives( terminationRootFinderSettings ) )
        {
            throw std::runtime_error( "Error when setting exact dependent variable termination, requested root finder "
                                      "requires derivatives; not available in state derivative model." );
        }
    }

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

    //! Function to return current difference between termination variable, and the value at which the propagation must terminate
    /*!
     * Function to return current difference between termination variable, and the value at which the propagation must terminate
     * \return Current difference between termination variable, and the value at which the propagation must terminate
     */
    double getStopConditionError( )
    {
        return variableRetrievalFunction_( ) - limitingValue_;
    }

    //! Function to retrieve settings to create root finder used to converge on exact final condition.
    /*!
     *  Function to retrieve settings to create root finder used to converge on exact final condition.
     *  \return Settings to create root finder used to converge on exact final condition.
     */
    std::shared_ptr< root_finders::RootFinderSettings > getTerminationRootFinderSettings( )
    {
        return terminationRootFinderSettings_;
    }

private:

    //! Settings for dependent variable that is to be checked
    std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Function returning the dependent variable.
    std::function< double( ) > variableRetrievalFunction_;

    //! Value at which the propagation is to be stopped
    double limitingValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below
    //! (if true) or above (if false) limitingValue
    bool useAsLowerBound_;

    //! Settings to create root finder used to converge on exact final condition.
    std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings_;
};

//! Class for stopping the propagation with custom stopping function.
class CustomTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param checkStopCondition Custom function to check for the attainment of the propagation stopping conditions.
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    CustomTerminationCondition(
            std::function< bool( const double ) >& checkStopCondition,
            const bool checkTerminationToExactCondition = false ):
        PropagationTerminationCondition( custom_stopping_condition, checkTerminationToExactCondition ),
        checkStopCondition_( checkStopCondition )
    { }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped via the user-provided function.
     * \param time Current time in propagation.
     * \param cpuTime Current CPU time in propagation.
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime )
    {
        TUDAT_UNUSED_PARAMETER( cpuTime );
        return checkStopCondition_( time );
    }

private:

    //! Custom temination function.
    std::function< bool( const double ) > checkStopCondition_;

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
     * \param fulfillSingleCondition Boolean denoting whether a single (if true) or all (if false) of the entries in the
     * propagationTerminationCondition_ should return true from the checkStopCondition function to stop the propagation
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    HybridPropagationTerminationCondition(
            const std::vector< std::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition,
            const bool fulfillSingleCondition = 0,
            const bool checkTerminationToExactCondition = 0 ):
        PropagationTerminationCondition( hybrid_stopping_condition, checkTerminationToExactCondition ),
        propagationTerminationCondition_( propagationTerminationCondition ),
        fulfillSingleCondition_( fulfillSingleCondition )
    {
        isConditionMetWhenStopping_.resize( propagationTerminationCondition.size( ) );
    }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. one or all (depending on value of
     * fulfillSingleCondition_) of the stopping conditions are fulfilled.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    //! Function to retrieve list of termination conditions that are checked when calling checkStopCondition is called.
    /*!
     *  Function to retrieve list of termination conditions that are checked when calling checkStopCondition is called.
     *  \return List of termination conditions that are checked when calling checkStopCondition is called.
     */
    std::vector< std::shared_ptr< PropagationTerminationCondition > > getPropagationTerminationConditions( )
    {
        return propagationTerminationCondition_;
    }

    //! Function to retrieve whether all or a single termination condition should be met
    /*!
     *  Function to retrieve whether all or a single termination condition should be met
     *  \return Boolean denoting whether a single (if true) or all (if false) of the entries in the
     *  propagationTerminationCondition_ should return true from the checkStopCondition function to stop the propagation.
     */
    bool getFulfillSingleCondition( )
    {
        return fulfillSingleCondition_;
    }

    std::vector< bool > getIsConditionMetWhenStopping( )
    {
        return isConditionMetWhenStopping_;
    }

private:

    //! List of termination conditions that are checked when calling checkStopCondition is called.
    std::vector< std::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition_;

    //!  Boolean denoting whether a single (if true) or all (if false) of the entries in the propagationTerminationCondition_
    //!  should return true from the checkStopCondition function to stop the propagation.
    bool fulfillSingleCondition_;

    std::vector< bool > isConditionMetWhenStopping_;

};

//! Function to create propagation termination conditions from associated settings
/*!
 * Function to create propagation termination conditions from associated settings
 * \param terminationSettings Settings for propagation termination conditions
 * \param bodies List of body objects that contains all environment models
 * \param initialTimeStep Time step at first call of numerical integration.
 * \return Object used to check whether propagation is to be stopped or not.
 */
//! Function to create propagation termination conditions from associated settings
template< typename StateScalarType = double , typename TimeType = double >
std::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const double initialTimeStep,
        const std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr
        < SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr
                < SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition;

    // Check termination type.
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
    {
        // Create stopping time termination condition.
        std::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
                std::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = std::make_shared< FixedTimePropagationTerminationCondition >(
                    timeTerminationSettings->terminationTime_, ( initialTimeStep > 0 ),
                    timeTerminationSettings->checkTerminationToExactCondition_ );
        break;
    }
    case cpu_time_stopping_condition:
    {
        // Create stopping time termination condition.
        std::shared_ptr< PropagationCPUTimeTerminationSettings > cpuTimeTerminationSettings =
                std::dynamic_pointer_cast< PropagationCPUTimeTerminationSettings >( terminationSettings );
        propagationTerminationCondition = std::make_shared< FixedCPUTimePropagationTerminationCondition >(
                    cpuTimeTerminationSettings->cpuTerminationTime_ );
        break;
    }
    case dependent_variable_stopping_condition:
    {
        std::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                std::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );

        // Get dependent variable function
        std::function< double( ) > dependentVariableFunction;
        if( getDependentVariableSaveSize( dependentVariableTerminationSettings->dependentVariableSettings_ ) == 1 )
        {
            dependentVariableFunction =
                    getDoubleDependentVariableFunction(
                        dependentVariableTerminationSettings->dependentVariableSettings_, bodies, stateDerivativeModels );
        }
        else
        {
            throw std::runtime_error( "Error, cannot make stopping condition from vector dependent variable" );
        }

        // Create dependent variable termination condition.
        propagationTerminationCondition = std::make_shared< SingleVariableLimitPropagationTerminationCondition >(
                    dependentVariableTerminationSettings->dependentVariableSettings_,
                    dependentVariableFunction, dependentVariableTerminationSettings->limitValue_,
                    dependentVariableTerminationSettings->useAsLowerLimit_,
                    dependentVariableTerminationSettings->checkTerminationToExactCondition_,
                    dependentVariableTerminationSettings->terminationRootFinderSettings_ );
        break;
    }
    case custom_stopping_condition:
    {
        std::shared_ptr< PropagationCustomTerminationSettings > customTerminationSettings =
                std::dynamic_pointer_cast< PropagationCustomTerminationSettings >( terminationSettings );

        // Create dependent variable termination condition.
        propagationTerminationCondition = std::make_shared< CustomTerminationCondition >(
                    customTerminationSettings->checkStopCondition_,
                    customTerminationSettings->checkTerminationToExactCondition_ );
        break;
    }
    case hybrid_stopping_condition:
    {
        std::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                std::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );

        // Recursively create termination condition list.
        std::vector< std::shared_ptr< PropagationTerminationCondition > > propagationTerminationConditionList;
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            propagationTerminationConditionList.push_back(
                        createPropagationTerminationConditions(
                            hybridTerminationSettings->terminationSettings_.at( i ),
                            bodies, initialTimeStep, stateDerivativeModels ) );
        }
        propagationTerminationCondition = std::make_shared< HybridPropagationTerminationCondition >(
                    propagationTerminationConditionList, hybridTerminationSettings->fulfillSingleCondition_,
                    hybridTerminationSettings->checkTerminationToExactCondition_ );
        break;
    }
    default:
        std::string errorMessage = "Error, stopping condition type " + std::to_string(
                    terminationSettings->terminationType_ ) + " not recognized when making stopping conditions object";
        throw std::runtime_error( errorMessage );
    }
    return propagationTerminationCondition;
}
//! Class for storing details on the propagation termination
class PropagationTerminationDetails
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param propagationTerminationReason Reason for termination
     * \param terminationOnExactCondition True if exact termination condition is used, false if not, -1 if neither is relevant
     */
    PropagationTerminationDetails( const PropagationTerminationReason propagationTerminationReason,
                                   const bool terminationOnExactCondition = 0 ):
        propagationTerminationReason_( propagationTerminationReason ),
        terminationOnExactCondition_( terminationOnExactCondition ){ }

    //! Destructor
    virtual ~PropagationTerminationDetails( ) { }

    //! Function to retrieve reason for termination
    /*!
     * Function to retrieve reason for termination
     * \return Reason for termination
     */
    PropagationTerminationReason getPropagationTerminationReason( )
    {
        return propagationTerminationReason_;
    }

    //! Function to retrieve boolean to denote whether exact termination conditions are used.
    /*!
     * Function to retrieve boolean to denote whether exact termination conditions are used.
     * \return Boolean to denote whether exact termination conditions are used.
     */
    bool getTerminationOnExactCondition( )
    {
        return terminationOnExactCondition_;
    }

protected:

    //! Reason for termination
    PropagationTerminationReason propagationTerminationReason_;

    //! Boolean to denote whether exact termination conditions are used.
    /*!
     *  Boolean to denote whether exact termination conditions are used. True if exact termination condition is used,
     *  false if not, -1 if neither is relevant.
     */
    bool terminationOnExactCondition_;

};

//! Class for storing details on the propagation termination when using hybrid termination conditions
class PropagationTerminationDetailsFromHybridCondition: public PropagationTerminationDetails
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param terminationOnExactCondition True if exact termination condition is used, false if not.
     * \param terminationCondition Hybrid termination conditions that were used
     */
    PropagationTerminationDetailsFromHybridCondition(
            const bool terminationOnExactCondition,
            const std::shared_ptr< HybridPropagationTerminationCondition > terminationCondition ):
        PropagationTerminationDetails( termination_condition_reached, terminationOnExactCondition ),
        isConditionMetWhenStopping_( terminationCondition->getIsConditionMetWhenStopping( ) ){ }

    //! Destructor
    ~PropagationTerminationDetailsFromHybridCondition( ) { }

    //! Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    /*!
     * Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     * \return List of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     */
    std::vector< bool > getWasConditionMetWhenStopping( )
    {
        if( terminationOnExactCondition_ )
        {
            std::cerr << "Warning when retrieving list of conditions that were met in hybrid propagation "
                         "termination details. propagation was terminated on exact conditions using root finder, "
                         "list of conditions may not be reliable" << std::endl;
        }
        return isConditionMetWhenStopping_;
    }

private:

    //! List of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    std::vector< bool > isConditionMetWhenStopping_;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
