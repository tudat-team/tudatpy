#ifndef TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
#define TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Propagators/propagationOutput.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

class PropagationTerminationCondition
{
public:
    PropagationTerminationCondition( ){ }

    virtual ~PropagationTerminationCondition( ){ }

    virtual bool checkStopCondition( const double time ) = 0;
};

class FixedTimePropagationTerminationCondition: public PropagationTerminationCondition
{
public:
    FixedTimePropagationTerminationCondition(
            const double stopTime,
            const bool propagationDirectionIsPositive ):
        stopTime_( stopTime ),
        propagationDirectionIsPositive_( propagationDirectionIsPositive ){ }


    bool checkStopCondition( const double time );

private:
    double stopTime_;

    bool propagationDirectionIsPositive_;
};

class SingleVariableLimitPropagationTerminationCondition: public PropagationTerminationCondition
{
public:
    SingleVariableLimitPropagationTerminationCondition(
            const std::pair< PropagationDependentVariables, std::string > variableType,
            const boost::function< double( ) > variableRetrievalFuntion,
            const double limitingValue,
            const bool useAsLowerBound ):
    variableType_( variableType ), variableRetrievalFuntion_( variableRetrievalFuntion ),
    limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ){ }

    virtual ~SingleVariableLimitPropagationTerminationCondition( ){ }

    bool checkStopCondition( const double time );

private:
    std::pair< PropagationDependentVariables, std::string > variableType_;

    boost::function< double( ) > variableRetrievalFuntion_;

    double limitingValue_;

    bool useAsLowerBound_;
};

class HybridPropagationTerminationCondition: public PropagationTerminationCondition
{
public:
    HybridPropagationTerminationCondition(
            const std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition,
            const bool fulFillSingleCondition = 0 ):
        propagationTerminationCondition_( propagationTerminationCondition ), fulFillSingleCondition_( fulFillSingleCondition ){ }

    bool checkStopCondition( const double time );

private:

    std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition_;

    bool fulFillSingleCondition_;
};


boost::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep );

}

}

#endif // TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
