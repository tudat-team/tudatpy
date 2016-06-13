#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

//! Get size of state for single propagated state of given type.
int getSingleIntegrationSize( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 6;
        break;
    default:
        std::cerr<<"Did not recognize state type when getting size"<<std::endl;
    }
    return singleStateSize;
}

//! Get order of differential equation for governing equations of dynamics of given type.
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 2;
        break;
    default:
        std::cerr<<"Did not recognize state type when getting order"<<std::endl;
    }
    return singleStateSize;
}


}

}
