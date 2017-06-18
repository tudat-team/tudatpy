#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"

namespace tudat
{

namespace basic_astrodynamics
{


AvailableTorque getTorqueModelType(
        boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel )
{
    AvailableTorque torqueType = underfined_torque;
    if( false )
    {

    }
    else
    {
        std::cerr<<"Error, could not identify torque type"<<std::endl;
    }
    return torqueType;
}

}

}

