#ifndef TUDAT_PROPAGATIONTERMINATIONSETTINGS_H
#define TUDAT_PROPAGATIONTERMINATIONSETTINGS_H

#include <vector>

#include <boost/shared_ptr.hpp>

namespace tudat
{

namespace propagators
{

//! Enum listing the available types of propagation termination settings.
enum PropagationTerminationTypes
{
    time_stopping_condition,
    dependent_variable_stopping_condition,
    hybrid_stopping_condition
};


//! Based class for defining propagation termination settings
class PropagationTerminationSettings
{
public:
    PropagationTerminationSettings( const PropagationTerminationTypes terminationType ):
        terminationType_( terminationType ){ }

    virtual ~PropagationTerminationSettings( ){ }

    PropagationTerminationTypes terminationType_;
};

class PropagationTimeTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationTimeTerminationSettings( const double terminationTime ):
        PropagationTerminationSettings( time_stopping_condition ),
        terminationTime_( terminationTime ){ }

    ~PropagationTimeTerminationSettings( ){ }

    double terminationTime_;
};

class PropagationDependentVariableTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationDependentVariableTerminationSettings(
            const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const double limitValue,
            const bool useAsLowerLimit ):
        PropagationTerminationSettings( dependent_variable_stopping_condition ),
        dependentVariableSettings_( dependentVariableSettings ),
        limitValue_( limitValue ), useAsLowerLimit_( useAsLowerLimit ){ }

    virtual ~PropagationDependentVariableTerminationSettings( ){ }

    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    double limitValue_;

    bool useAsLowerLimit_;
};

class PropagationHybridTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationHybridTerminationSettings(
            const std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings,
            const bool fulFillSingleCondition = 0 ):
        PropagationTerminationSettings( hybrid_stopping_condition ),
        terminationSettings_( terminationSettings ),
        fulFillSingleCondition_( fulFillSingleCondition ){ }

    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings_;

    bool fulFillSingleCondition_;
};

}

}

#endif // TUDAT_PROPAGATIONTERMINATIONSETTINGS_H
