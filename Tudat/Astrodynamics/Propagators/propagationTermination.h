/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include "Tudat/Astrodynamics/Propagators/propagationOutput.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

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
    PropagationTerminationCondition( ){ }

    //! Destructor
    virtual ~PropagationTerminationCondition( ){ }

    //! (Pure virtual) function to check whether the propagation should be stopped
    /*!
     * (Pure virtual) function to check whether the propagation should be stopped. Note that the accelerations and
     * environment must be updated (done automatically during numerical propagation) to check the stopping condition.
     * \param time Current time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    virtual bool checkStopCondition( const double time ) = 0;
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
            const bool propagationDirectionIsPositive ):
        stopTime_( stopTime ),
        propagationDirectionIsPositive_( propagationDirectionIsPositive ){ }


    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the stopTime_ has been reached or not.
     * \param time Current time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time );

private:

    //! Time at which the propagation is to stop.
    double stopTime_;

    //!  Boolean denoting whether propagation is forward (if true) or backwards (if false) in time.
    bool propagationDirectionIsPositive_;
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
            const bool useAsLowerBound ):
        dependentVariableSettings_( dependentVariableSettings ), variableRetrievalFuntion_( variableRetrievalFuntion ),
        limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ){ }

    //! Destructor.
    ~SingleVariableLimitPropagationTerminationCondition( ){ }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the given dependent variable has been
     * reached or not.
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time );

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
        propagationTerminationCondition_( propagationTerminationCondition ),
        fulFillSingleCondition_( fulFillSingleCondition ){ }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. one or all (depending on value of
     * fulFillSingleCondition_) of the stopping conditions are fulfilled.
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time );

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
