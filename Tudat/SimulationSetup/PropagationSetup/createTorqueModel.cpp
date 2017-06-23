/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create torque model object.
boost::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const boost::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel;

    switch( torqueSettings->torqueType_ )
    {
    default:
        throw std::runtime_error(
                    "Error, did not recognize type " + boost::lexical_cast< std::string >( torqueSettings->torqueType_ ) +
                    " when making torque model" );
    }

    return torqueModel;
}


}

}

