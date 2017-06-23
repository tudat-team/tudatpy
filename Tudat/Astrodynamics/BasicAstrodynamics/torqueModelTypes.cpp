/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to identify the derived class type of a torque model.
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

