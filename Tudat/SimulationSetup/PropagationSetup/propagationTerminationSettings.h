/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/createRootFinder.h"

namespace tudat
{

namespace propagators
{

class SingleDependentVariableSaveSettings;

//! Enum listing the available types of propagation termination settings.
enum PropagationTerminationTypes
{
    time_stopping_condition,
    cpu_time_stopping_condition,
    dependent_variable_stopping_condition,
    hybrid_stopping_condition
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
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTerminationSettings( const PropagationTerminationTypes terminationType,
                                    const bool terminateExactlyOnFinalCondition = false ):
        terminationType_( terminationType ), terminateExactlyOnFinalCondition_( terminateExactlyOnFinalCondition ){ }

    //! Destructor
    virtual ~PropagationTerminationSettings( ){ }

    //! Type of stopping condition that is to be used.
    PropagationTerminationTypes terminationType_;

    //! Boolean to denote whether the propagation is to terminate exactly on the final condition, or whether it is to terminate
    //! on the first step where it is violated.
    bool terminateExactlyOnFinalCondition_;
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
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTimeTerminationSettings( const double terminationTime,
                                        const bool terminateExactlyOnFinalCondition = false ):
        PropagationTerminationSettings( time_stopping_condition, terminateExactlyOnFinalCondition ),
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
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     * \param terminationRootFinderSettings Settings to create root finder used to converge on exact final condition.
     */
    PropagationDependentVariableTerminationSettings(
            const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double limitValue,
            const bool useAsLowerLimit,
            const bool terminateExactlyOnFinalCondition = false,
            const boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = NULL ):
        PropagationTerminationSettings(
            dependent_variable_stopping_condition, terminateExactlyOnFinalCondition ),
        dependentVariableSettings_( dependentVariableSettings ),
        limitValue_( limitValue ), useAsLowerLimit_( useAsLowerLimit ),
        terminationRootFinderSettings_( terminationRootFinderSettings )
    {
        if( terminateExactlyOnFinalCondition_ && ( terminationRootFinderSettings_ == NULL ) )
        {
            throw std::runtime_error( "Error when defining exavct dependent variable propagation termination settings. Root finder not defined" );
        }
    }

    //! Destructor
    ~PropagationDependentVariableTerminationSettings( ){ }

    //! Settings for dependent variable that is to be checked
    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Value at which the propagation is to be stopped
    double limitValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below (if true) or above
    //! (if false) limitingValue
    bool useAsLowerLimit_;

    //! Settings to create root finder used to converge on exact final condition.
    boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings_;
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
     * \param fulFillSingleCondition Boolean denoting whether a single (if true) or all (if false) of the conditions
     * defined by the entries in the terminationSettings list should be met.
     */
    PropagationHybridTerminationSettings(
            const std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings,
            const bool fulFillSingleCondition = 0 ):
        PropagationTerminationSettings( hybrid_stopping_condition ),
        terminationSettings_( terminationSettings ),
        fulFillSingleCondition_( fulFillSingleCondition )
    {
        for( unsigned int i = 0; i < terminationSettings_.size( ); i++ )
        {
            if( i == 0 )
            {
                terminateExactlyOnFinalCondition_ = terminationSettings_.at( 0 )->terminateExactlyOnFinalCondition_;
            }
            else if( terminationSettings_.at( i )->terminateExactlyOnFinalCondition_ != terminateExactlyOnFinalCondition_ )
            {
                throw std::runtime_error( "Error in hybrid termination settings, terminateExactlyOnFinalCondition_ is inconsistent" );
            }
        }
    }

    //! Destructor
    ~PropagationHybridTerminationSettings( ){ }

    //! List of termination settings for which stopping conditions are created.
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings_;

    //! Boolean denoting whether a single (if true) or all (if false) of the conditions
    //! defined by the entries in the terminationSettings list should be met.
    bool fulFillSingleCondition_;
};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONTERMINATIONSETTINGS_H
