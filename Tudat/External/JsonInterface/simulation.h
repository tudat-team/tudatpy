/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SIMULATION_H
#define TUDAT_JSONINTERFACE_SIMULATION_H

#include <boost/filesystem.hpp>

#include "json/src/json.hpp"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "utilities.h"
#include "integratorSettings.h"


namespace tudat
{

namespace json_interface
{

typedef boost::filesystem::path path;


//! -DOC
template< typename TimeType = double >
class Simulation
{
private:
    //! Absolute path to the input file.
    path inputPath;

    //! JSON object with all the settings from the input file.
    json settings;


public:
    //! Set with filenames of the Spice kernels located in "Tudat/External/SpiceInterface/Kernels/"
    //! to be used for the propagation.
    std::set< std::string > spiceKernels;

    //! Whether to preload the ephemeris of the celestial bodies for the simulation period,
    //! or to retrieve this directly from Spice during the propagation at each integration step.
    //! Preloading Spice data generally results in faster propagations, unless:
    //! * The simulation ends much earlier than the specified maximum simulation end epoch.
    //! * The integrator step-size is very large (in the order of several hours or days).
    bool preloadSpiceData;

    //! Initial epoch for the simulation.
    TimeType startEpoch;

    //! Maximum end epoch for the simulation.
    TimeType endEpoch;

    //! Global frame origin.
    std::string globalFrameOrigin;

    //! Global frame orientation.
    std::string globalFrameOrientation;

    //! Set with the names of all the bodies.
    std::set< std::string > bodies;

    //! Set with the names of the celestial bodies (handled by Spice).
    std::set< std::string > celestialBodies;

    //! Set with the names of the central bodies.
    std::set< std::string > centralBodies;

    //! Set with the names of the bodies to be propagated.
    std::set< std::string > bodiesToPropagate;

    //! Body map.
    simulation_setup::NamedBodyMap bodyMap;

    //! Integrator settings.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings;


    //! -DOC
    Simulation( std::string inputFile )
    {
        setInputFile( inputFile );
        reset( );
    }


    //! -DOC
    void setInputFile( std::string inputFile )
    {
        // Add .json extension if necessary
        std::string inputFileWithExtension = inputFile;
        if ( inputFileWithExtension.find(".json") == std::string::npos )
        {
            inputFileWithExtension += ".json";
        }

        // Determine the absolute path of the directory where the input file is located
        path inputDirectory = boost::filesystem::canonical( inputFileWithExtension ).parent_path();

        // Determine the filename (without file extension) of the input file
        // path inputFilename = path( inputFile ).stem();

        // Build absolute path
        inputPath = inputDirectory / path( inputFileWithExtension ).filename();

        // Check that input file exists
        if ( ! boost::filesystem::exists( inputPath ) )
        {
            std::cerr << "The input file \"" << inputPath << "\" does not exist." << std::endl;
            std::terminate();
        }

        // Read and parse input file
        settings = json::parse( std::ifstream( inputPath.string( ) ) );
    }


    //! Update/create all the objects from the JSON data before the simulation can be run.
    void reset( )
    {
        resetGeneral( );
        resetSpice( );
        resetBodies( );
        resetAccelerations( );
        resetPropagators( );
        resetIntegrator();
        // resetOutput( );
    }


    //! -DOC
    void resetGeneral( )
    {
        // Start and end epochs
        startEpoch = getValue< TimeType >( settings, { "simulation", "startEpoch" } );
        endEpoch = getValue< TimeType >( settings, { "simulation", "endEpoch" } );

        // Global fram origin and orientation
        globalFrameOrigin = getValue< std::string >( settings, { "simulation", "globalFrameOrigin" } );
        globalFrameOrientation = getValue< std::string >( settings, { "simulation", "globalFrameOrientation" } );
    }


    //! -DOC
    void resetSpice( )
    {
        // Spice
        spiceKernels = getValue< std::set< std::string > >( settings, { "simulation", "spiceKernels" }, { } );
        preloadSpiceData = getValue< bool >( settings, { "simulation", "preloadSpiceData" }, true );

        // Load requested Spice kernels.
        for ( std::string kernelFile : spiceKernels )
        {
            spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + kernelFile );
        }
    }


    //! -DOC
    void resetBodies( )
    {
        using namespace simulation_setup;

        // All bodies
        bodies = { };
        for ( auto ent : getValue< std::map< std::string, json > >( settings, "bodies" ) )
        {
            bodies.insert( ent.first );
        }

        // Bodies to propagate
        std::vector< json > propagators = getValue< std::vector< json > >( settings, "propagators" );
        bodiesToPropagate = { };
        for ( json propagator : propagators )
        {
            std::set< std::string > bodies = getValue< std::set< std::string > >( propagator, "bodies" );
            for ( std::string body : bodies )
            {
                if ( bodies.count( body ) )
                {
                    bodiesToPropagate.insert( body );
                }
                else
                {
                    throw std::runtime_error( "Requested to propagate body named \"" + body +
                                              "\", but no body with this name was found.");
                }
            }
        }

        // Celestial (and central) bodies
        celestialBodies = { };
        centralBodies = { };
        for ( std::string bodyName : bodies )
        {
            if ( bodiesToPropagate.count( bodyName ) == 0 )
            {
                celestialBodies.insert( bodyName );
                if ( getValue( settings, { "bodies", bodyName, "isCentralBody" }, false ) )
                {
                    centralBodies.insert( bodyName );
                }
            }
        }

        // Get default settings for celestial bodies.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
        if ( preloadSpiceData )
        {
            bodySettings = getDefaultBodySettings( set2vector( celestialBodies ),
                                                   startEpoch - 300.0, endEpoch + 300.0 );
        }
        else
        {
            bodySettings = getDefaultBodySettings( set2vector( celestialBodies ) );
        }

        // Overwrite values that were not set to be default / were not omitted.
        // ... such as µ

        // Set celestial bodies ephemeris and rotational models.
        for ( std::string celestialBody : celestialBodies )
        {
            bodySettings[ celestialBody ]->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
            bodySettings[ celestialBody ]->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
        }

        // Create celestial bodies.
        bodyMap = simulation_setup::createBodies( bodySettings );

        // Create bodies to be propagated.
        for ( std::string bodyName : bodiesToPropagate )
        {
            boost::shared_ptr< Body > body = boost::make_shared< Body >( );
            body->setConstantBodyMass( getValue< double >( settings, { "bodies", bodyName, "mass" } ) );
            bodyMap[ bodyName ] = body;
            // ...
        }

    }


    //! -DOC
    void resetAccelerations( )
    {

    }


    //! -DOC
    void resetPropagators( )
    {

    }


    //! -DOC
    void resetIntegrator( )
    {
        // Integrator settings
        json jsonIntegratorSettings = getValue< json >( settings, { "integrator" } );
        if ( getValuePointer< TimeType >( jsonIntegratorSettings, { "initialTime" } ) == NULL )
        {
            jsonIntegratorSettings[ "initialTime" ] = startEpoch;
        }
        integratorSettings = createIntegratorSettings< TimeType >( jsonIntegratorSettings );
    }



    //! DOC
    void run( )
    {

    }


    //! DOC
    void exportResults( )
    {

    }


    //! DOC
    json getOriginalSettings( )
    {
        return settings;
    }


};


//! Function to create a `json` object from a `Simulation` object.
//! Called automatically by `nlohmann::json` when writing `json( simulation )`.
template< typename TimeType = double >
void to_json( json& j, const Simulation< TimeType >& simulation ) {
    // Initialise
    j = json();

    // Simulation
    json jSimulation;
    jSimulation[ "startEpoch" ] = simulation.startEpoch;
    jSimulation[ "endEpoch" ] = simulation.endEpoch;
    j[ "simulation" ] = jSimulation;

    // Integrator
    j[ "integrator" ] = json( simulation.integratorSettings );
}


} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SIMULATION_H
