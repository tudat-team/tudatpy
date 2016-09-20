/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#endif
#include "Tudat/SimulationSetup/createRotationModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a rotation model.
boost::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const boost::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body )
{
    using namespace tudat::ephemerides;

    // Declare return object.
    boost::shared_ptr< RotationalEphemeris > rotationalEphemeris;

    // Check which type of rotation model is to be created.
    switch( rotationModelSettings->getRotationType( ) )
    {
    case simple_rotation_model:
    {
        // Check whether settings for simple rotation model are consistent with its type.
        boost::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
                boost::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
        if( simpleRotationSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected simple rotation model settings for " + body );
        }
        else
        {
            // Create and initialize simple rotation model.
            rotationalEphemeris = boost::make_shared< SimpleRotationalEphemeris >(
                        simpleRotationSettings->getInitialOrientation( ),
                        simpleRotationSettings->getRotationRate( ),
                        simpleRotationSettings->getInitialTime( ),
                        basic_astrodynamics::JULIAN_DAY_ON_J2000,
                        simpleRotationSettings->getOriginalFrame( ),
                        simpleRotationSettings->getTargetFrame( ) );
        }
        break;
    }
    #if USE_CSPICE
    case spice_rotation_model:
    {
        // Create rotational ephemeris directly from Spice.
        rotationalEphemeris = boost::make_shared< SpiceRotationalEphemeris >(
                    rotationModelSettings->getOriginalFrame( ),
                    rotationModelSettings->getTargetFrame( ) );
        break;
    }
    #endif
    default:
        throw std::runtime_error(
                 "Error, did not recognize rotation model settings type " +
                  boost::lexical_cast< std::string >( rotationModelSettings->getRotationType( ) ) );
    }

    return rotationalEphemeris;
}

} // namespace simulation_setup

} // namespace tudat
