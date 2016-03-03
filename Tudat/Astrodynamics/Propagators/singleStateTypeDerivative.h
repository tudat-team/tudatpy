#ifndef STATEDERIVATIVE_H
#define STATEDERIVATIVE_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType = double, typename TimeType = double >
class SingleStateTypeDerivative
{
public:
    SingleStateTypeDerivative( const IntegratedStateType integratedStateType ):
        integratedStateType_( integratedStateType ){ }

    virtual ~SingleStateTypeDerivative( ){ }

    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated ) =  0;

    virtual void updateStateDerivativeModel( const TimeType currentTime ) = 0;

    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time ) = 0;

    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time ) = 0;

    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time ) = 0;

    virtual int getStateSize( ) = 0;

    IntegratedStateType getIntegratedStateType( )
    {
        return integratedStateType_;
    }

protected:
    IntegratedStateType integratedStateType_;

};

}

}

#endif // STATEDERIVATIVE_H
