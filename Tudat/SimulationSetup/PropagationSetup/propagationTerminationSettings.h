/*    Copyright (c) 2010-2016, Delft University of Technology
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

namespace tudat
{

namespace propagators
{

class SingleDependentVariableSaveSettings;

//! Enum listing the available types of propagation termination settings.
enum PropagationTerminationTypes
{
    time_stopping_condition,
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
     */
    PropagationTerminationSettings( const PropagationTerminationTypes terminationType ):
        terminationType_( terminationType ){ }

    //! Destructor
    virtual ~PropagationTerminationSettings( ){ }

    //! Type of stopping condition that is to be used.
    PropagationTerminationTypes terminationType_;
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
     */
    PropagationTimeTerminationSettings( const double terminationTime ):
        PropagationTerminationSettings( time_stopping_condition ),
        terminationTime_( terminationTime ){ }

    //! Destructor
    ~PropagationTimeTerminationSettings( ){ }

    //! Maximum time for the propagation, upon which the propagation is to be stopped
    double terminationTime_;
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
     */
    PropagationDependentVariableTerminationSettings(
            const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double limitValue,
            const bool useAsLowerLimit ):
        PropagationTerminationSettings( dependent_variable_stopping_condition ),
        dependentVariableSettings_( dependentVariableSettings ),
        limitValue_( limitValue ), useAsLowerLimit_( useAsLowerLimit ){ }

    //! Destructor
    ~PropagationDependentVariableTerminationSettings( ){ }

    //! Settings for dependent variable that is to be checked
    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Value at which the propagation is to be stopped
    double limitValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below (if true) or above
    //! (if false) limitingValue
    bool useAsLowerLimit_;
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
        fulFillSingleCondition_( fulFillSingleCondition ){ }

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
