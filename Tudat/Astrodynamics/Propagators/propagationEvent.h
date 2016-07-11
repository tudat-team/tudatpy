#ifndef PROPAGATIONEVENT_H
#define PROPAGATIONEVENT_H

#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Astrodynamics/Propagators/propagationTermination.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

enum PropagationEventType
{
    impulsive_thrust,
    staging,
    acceleration_model_change
};

class PropagationEvent
{
public:
    PropagationEvent(
            const PropagationEventType propagationEventType,
            const boost::shared_ptr< PropagationTerminationCondition > eventCondition ):
        propagationEventType_( propagationEventType ),
    eventCondition_( eventCondition ){ }

    PropagationEventType getPropagationEventType( )
    {
        return propagationEventType_;
    }

protected:
    PropagationEventType propagationEventType_;

    boost::shared_ptr< PropagationTerminationCondition > eventCondition_;

};


template< typename TimeType = double, typename StateScalarType = double >
void triggerEvent(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > stateDerivativeModel,
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState,
        const boost::shared_ptr< PropagationEvent > propagationEvent )
{

}

}

}

#endif // PROPAGATIONEVENT_H
