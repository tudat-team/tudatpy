#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case initial_body_state:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}


std::string getParameterType( EstimatebleParametersEnum parameter )
{
    std::string parameterType;
    switch( parameter )
    {
    case gravitational_parameter:
        parameterType = "double";
        break;
    case constant_drag_coefficient:
        parameterType = "double";
        break;
    case radiation_pressure_coefficient:
        parameterType = "double";
        break;
    default:
        std::cerr<<"Error, parameter type "<<parameter<<" not found when getting parameter type"<<std::endl;
    }
    return parameterType;
}


}

}

