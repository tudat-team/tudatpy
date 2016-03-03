#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

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
    default:
        std::cerr<<"Error, did not recognize state type "<<stateType<<" when getting integration size"<<std::endl;
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
    default:
        std::cerr<<"Error, did not recognize state type "<<stateType<<" when getting equation order"<<std::endl;
    }
    return singleStateSize;
}


}

}
