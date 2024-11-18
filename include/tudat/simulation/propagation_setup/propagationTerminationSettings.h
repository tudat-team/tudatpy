/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONTERMINATIONSETTINGS_H
#define TUDAT_PROPAGATIONTERMINATIONSETTINGS_H

#include <vector>

#include <memory>

#include "tudat/math/root_finders/createRootFinder.h"

namespace tudat
{

namespace propagators
{

class SingleDependentVariableSaveSettings;

//! Enum listing the available types of propagation termination settings.
enum PropagationTerminationTypes
{
    time_stopping_condition = 0,
    cpu_time_stopping_condition = 1,
    dependent_variable_stopping_condition = 2,
    hybrid_stopping_condition = 3,
    custom_stopping_condition = 4,
    non_sequential_stopping_condition = 5
};


//! Base class for defining propagation termination settings.
/*!
 *  Base class for defining propagation termination settings, i.e. conditions on which the porpagation is deemed to be
 *  'finished'. Each particular type of stopping condition requires a different derived class.
 */
class PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param terminationType Type of stopping condition that is to be used.
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTerminationSettings( const PropagationTerminationTypes terminationType,
                                    const bool checkTerminationToExactCondition = false ):
        terminationType_( terminationType ), checkTerminationToExactCondition_( checkTerminationToExactCondition ){ }

    //! Destructor
    virtual ~PropagationTerminationSettings( ){ }

    //! Type of stopping condition that is to be used.
    PropagationTerminationTypes terminationType_;

    //! Boolean to denote whether the propagation is to terminate exactly on the final condition, or whether it is to terminate
    //! on the first step where it is violated.
    bool checkTerminationToExactCondition_;

};

//! Class for propagation stopping conditions settings: stopping the propagation after a fixed amount of time
/*!
 *  Class for propagation stopping conditions settings: stopping the propagation after a fixed amount of time. Note that the
 *  propagator will finish a given time step, slightly surpassing the defined final time.
 */
class PropagationTimeTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param terminationTime Maximum time for the propagation, upon which the propagation is to be stopped
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTimeTerminationSettings( const double terminationTime,
                                        const bool checkTerminationToExactCondition = false ):
        PropagationTerminationSettings( time_stopping_condition, checkTerminationToExactCondition ),
        terminationTime_( terminationTime ){ }

    //! Destructor
    ~PropagationTimeTerminationSettings( ){ }

    //! Maximum time for the propagation, upon which the propagation is to be stopped
    double terminationTime_;

};

//! Class for propagation stopping conditions settings: stopping the propagation after a fixed amount of CPU time
/*!
 *  Class for propagation stopping conditions settings: stopping the propagation after a fixed amount of CPU time.
 *  Note that the propagator will finish a given time step, slightly surpassing the defined final CPU time.
 */
class PropagationCPUTimeTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param cpuTerminationTime Maximum cpu time for the propagation, upon which the propagation is to be stopped
     */
    PropagationCPUTimeTerminationSettings( const double cpuTerminationTime ):
        PropagationTerminationSettings( cpu_time_stopping_condition ),
        cpuTerminationTime_( cpuTerminationTime ){ }

    //! Destructor
    ~PropagationCPUTimeTerminationSettings( ){ }

    //! Maximum cpu time for the propagation, upon which the propagation is to be stopped
    double cpuTerminationTime_;

};

//! Class for propagation stopping conditions settings: stopping the propagation after a given dependent variable reaches a
//! certain value.
/*!
 *  Class for propagation stopping conditions settings: stopping the propagation after a given dependent variable reaches a
 *  certain value. The limit value may be set as both an upper or lower bound (i.e. the propagation continues while the
 *  value is below or above some given value).
 *  Note that the propagator will finish a given time step, slightly surpassing the defined limit value of the dependent
 *  variable
 */
class PropagationDependentVariableTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariableSettings Settings for dependent variable that is to be checked
     * \param limitValue Value at which the propagation is to be stopped
     * \param useAsLowerLimit Boolean denoting whether the propagation should stop if the dependent variable goes below
     * (if true) or above (if false) limitingValue
     * \param checkTerminationToExactCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     * \param terminationRootFinderSettings Settings to create root finder used to converge on exact final condition.
     */
    PropagationDependentVariableTerminationSettings(
            const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double limitValue,
            const bool useAsLowerLimit,
            const bool checkTerminationToExactCondition = false,
            const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = nullptr ):
        PropagationTerminationSettings(
            dependent_variable_stopping_condition, checkTerminationToExactCondition ),
        dependentVariableSettings_( dependentVariableSettings ),
        limitValue_( limitValue ), useAsLowerLimit_( useAsLowerLimit ),
        terminationRootFinderSettings_( terminationRootFinderSettings )
    {
        if( checkTerminationToExactCondition_ && ( terminationRootFinderSettings_ == nullptr ) )
        {
            throw std::runtime_error( "Error when defining exact dependent variable propagation termination settings. Root finder not defined." );
        }
    }

    //! Destructor
    ~PropagationDependentVariableTerminationSettings( ){ }

    //! Settings for dependent variable that is to be checked
    std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Value at which the propagation is to be stopped
    double limitValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below (if true) or above
    //! (if false) limitingValue
    bool useAsLowerLimit_;

    //! Settings to create root finder used to converge on exact final condition.
    std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings_;
};

//! Class for propagation stopping conditions settings: stopping the propagation based on custom requirements
/*!
 *  Class for propagation stopping conditions settings: stopping the propagation after a fixed amount of time. Note that the
 *  propagator will finish a given time step, slightly surpassing the defined final time.
 */
class PropagationCustomTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param checkStopCondition Function that takes the current time as input and outputs whether the propagation should be
     *      stopped.
     */
    PropagationCustomTerminationSettings( const std::function< bool( const double ) >& checkStopCondition ):
        PropagationTerminationSettings( custom_stopping_condition ),
        checkStopCondition_( [=]( const double time, const Eigen::MatrixXd&  ){ return checkStopCondition( time ); } ){ }

    PropagationCustomTerminationSettings( const std::function< bool( const double, const Eigen::MatrixXd&  ) >& checkStopCondition ):
        PropagationTerminationSettings( custom_stopping_condition ),
        checkStopCondition_( checkStopCondition ){ }

    //! Destructor
    ~PropagationCustomTerminationSettings( ){ }

    //! Custom temination function.
    std::function< bool( const double, const Eigen::MatrixXd& ) > checkStopCondition_;

};

//! Class for propagation stopping conditions settings: combination of other stopping conditions.
/*!
 *  Class for propagation stopping conditions settings: combination of other stopping conditions. This class can be used
 *  to define that all of any number of conditions must be met, or that a single of these settings must be met to
 *  stop the propagation.
 */
class PropagationHybridTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * \brief PropagationHybridTerminationSettings
     * \param terminationSettings List of termination settings for which stopping conditions are created.
     * \param fulfillSingleCondition Boolean denoting whether a single (if true) or all (if false) of the conditions
     * defined by the entries in the terminationSettings list should be met.
     */
    PropagationHybridTerminationSettings(
            const std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettings,
            const bool fulfillSingleCondition = false ):
        PropagationTerminationSettings( hybrid_stopping_condition ),
        terminationSettings_( terminationSettings ),
        fulfillSingleCondition_( fulfillSingleCondition )
    {
        checkTerminationToExactCondition_ = false;
        for( unsigned int i = 0; i < terminationSettings_.size( ); i++ )
        {
            if( terminationSettings_.at( i )->checkTerminationToExactCondition_ )
            {
                checkTerminationToExactCondition_ = true;
            }
        }
    }

    //! Destructor
    ~PropagationHybridTerminationSettings( ){ }

    //! List of termination settings for which stopping conditions are created.
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettings_;

    //! Boolean denoting whether a single (if true) or all (if false) of the conditions
    //! defined by the entries in the terminationSettings list should be met.
    bool fulfillSingleCondition_;

};

//! Class for propagation stopping conditions settings: combination of two stopping conditions.
/*!
 *  Class for propagation stopping conditions settings: combination of two stopping conditions. This class should be used
 *  to define a "forward" and a "backward" propagation legs.
 */
class NonSequentialPropagationTerminationSettings: public PropagationTerminationSettings
{
public:

    //! Constructor
    /*!
     * \brief NonSequentialPropagationTerminationSettings
     * \param forwardTerminationSettings Termination settings for forward propagation.
     * \param backwardTerminationSettings Termination settings for backward propagation.
     */
    NonSequentialPropagationTerminationSettings(
            const std::shared_ptr< PropagationTerminationSettings > forwardTerminationSettings,
            const std::shared_ptr< PropagationTerminationSettings > backwardTerminationSettings ):
            PropagationTerminationSettings( non_sequential_stopping_condition ),
            forwardTerminationSettings_( forwardTerminationSettings ),
            backwardTerminationSettings_( backwardTerminationSettings )
    {
        checkTerminationToExactCondition_ = forwardTerminationSettings_->checkTerminationToExactCondition_;
    }

    //! Destructor
    ~NonSequentialPropagationTerminationSettings( ){ }

    //! Termination settings for forward propagation leg.
    std::shared_ptr< PropagationTerminationSettings > forwardTerminationSettings_;

    //! Termination settings for backward propagation leg.
    std::shared_ptr< PropagationTerminationSettings > backwardTerminationSettings_;

};

inline std::shared_ptr< PropagationTerminationSettings > propagationDependentVariableTerminationSettings(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const double limitValue,
        const bool useAsLowerLimit,
        const bool checkTerminationToExactCondition = false,
        const std::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = nullptr )
{
    return std::make_shared< PropagationDependentVariableTerminationSettings >(
                dependentVariableSettings, limitValue, useAsLowerLimit, checkTerminationToExactCondition,
                terminationRootFinderSettings );
}

inline std::shared_ptr< PropagationTerminationSettings > propagationTimeTerminationSettings(
        const double terminationTime,
        const bool checkTerminationToExactCondition = false  )
{
    return std::make_shared< PropagationTimeTerminationSettings >(
                terminationTime, checkTerminationToExactCondition );
}

inline std::shared_ptr< PropagationTerminationSettings > propagationCPUTimeTerminationSettings(
        const double cpuTerminationTime )
{
    return std::make_shared< PropagationCPUTimeTerminationSettings >(
                cpuTerminationTime );
}

inline std::shared_ptr< PropagationTerminationSettings > propagationHybridTerminationSettings(
        const std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettings,
        const bool fulfillSingleCondition = false )
{
    return std::make_shared< PropagationHybridTerminationSettings >(
                terminationSettings, fulfillSingleCondition );
}

inline std::shared_ptr< PropagationTerminationSettings > popagationCustomTerminationSettings(
        const std::function< bool( const double ) > checkStopCondition )
{
    return std::make_shared< PropagationCustomTerminationSettings >(
                checkStopCondition );
}

inline std::shared_ptr< PropagationTerminationSettings > popagationCustomTerminationSettingsFromFullState(
    const std::function< bool( const double, const Eigen::MatrixXd& ) > checkStopCondition )
{
    return std::make_shared< PropagationCustomTerminationSettings >(
        checkStopCondition );
}

inline std::shared_ptr< PropagationTerminationSettings > nonSequentialPropagationTerminationSettings(
        const std::shared_ptr< PropagationTerminationSettings > forwardTerminationSettings,
        const std::shared_ptr< PropagationTerminationSettings > backwardTerminationSettings )
{
    return std::make_shared< NonSequentialPropagationTerminationSettings >(
            forwardTerminationSettings, backwardTerminationSettings );
}



} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONTERMINATIONSETTINGS_H
