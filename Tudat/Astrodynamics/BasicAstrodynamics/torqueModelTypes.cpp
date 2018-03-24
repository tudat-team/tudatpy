/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicTorque.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to identify the derived class type of a torque model.
AvailableTorque getTorqueModelType(
        boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel )
{
    AvailableTorque torqueType = underfined_torque;
    if( boost::dynamic_pointer_cast< gravitation::SecondDegreeGravitationalTorqueModel >( torqueModel ) != NULL )
    {
        torqueType = second_order_gravitational_torque;
    }
    else if( boost::dynamic_pointer_cast< aerodynamics::AerodynamicTorque >( torqueModel ) != NULL )
    {
        torqueType = aerodynamic_torque;
    }
    else
    {
        std::cerr << "Error, could not identify torque type" << std::endl;
    }
    return torqueType;
}

//! Function to get a string representing a 'named identification' of an torque type
std::string getTorqueModelName( const AvailableTorque torqueType )
{

    std::string torqueName;
    switch( torqueType )
    {
    case second_order_gravitational_torque:
        torqueName = "second-order gravitational torque ";
        break;
    case aerodynamic_torque:
        torqueName = "aerodynamic torque ";
        break;
    default:
        std::string errorMessage = "Error, torque type " +
                std::to_string( torqueType ) +
                "not found when retrieving torque name ";
        throw std::runtime_error( errorMessage );
    }
    return torqueName;
}

//! Function to get all torque models of a given type from a list of models
std::vector< boost::shared_ptr< TorqueModel > > getTorqueModelsOfType(
        const std::vector< boost::shared_ptr< TorqueModel > >& fullList,
        const AvailableTorque modelType )
{
    std::vector< boost::shared_ptr< TorqueModel > > torqueList;
    for( unsigned int i = 0; i < fullList.size( ); i++ )
    {
        if( getTorqueModelType( fullList.at( i ) ) == modelType )
        {
            torqueList.push_back( fullList.at( i  ) );
        }
    }
    return torqueList;
}

}

}

