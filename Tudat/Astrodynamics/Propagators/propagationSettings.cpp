#include "Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

int getSingleIntegrationSize( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case hybrid:
        std::cerr<<"Error, cannot return size of hybrid state type"<<std::endl;
    case transational_state:
        singleStateSize = 6;
        break;
    case rotational_state:
        singleStateSize = 7;
        break;
    case proper_time:
        singleStateSize = 1;
        break;
    }
    return singleStateSize;
}

int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case hybrid:
        std::cerr<<"Error, cannot return order of hybrid state type"<<std::endl;
    case transational_state:
        singleStateSize = 2;
        break;
    case rotational_state:
        std::cerr<<"Error, cannot yet return order of rotational state"<<std::endl;
        break;
    case proper_time:
        singleStateSize = 1;
        break;
    }
    return singleStateSize;
}


}

}
