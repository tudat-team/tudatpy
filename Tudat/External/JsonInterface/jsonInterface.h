/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_H
#define TUDAT_JSONINTERFACE_H

#include <boost/filesystem.hpp>

#include "json/src/json.hpp"

#include "integratorSettings.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


namespace tudat
{

namespace json_interface
{

using json = nlohmann::json;
typedef boost::filesystem::path path;

class Simulation
{
private:
    //! Absolute path to the input file.
    path inputPath;

    //! JSON object with all the settings from the input file.
    json settings;

    void parseInputFile( )
    {
        /// Read and parse
        settings = json::parse( std::ifstream( inputPath.string( ) ) );


        /// Validate values and update members

        // Spice kernels
        spiceKernels = settings[ "simulation" ][ "spiceKernels" ];

        // Start and end epochs
        startEpoch = settings[ "simulation" ][ "startEpoch" ];
        endEpoch = settings[ "simulation" ][ "endEpoch" ];
        BOOST_ASSERT_MSG( endEpoch > startEpoch, "The end epoch must be after the start epoch." );

        // Global fram origin and orientation
        globalFrameOrigin = settings[ "simulation" ][ "globalFrameOrigin" ];
        globalFrameOrientation = settings[ "simulation" ][ "globalFrameOrientation" ];

        // All bodies
        bodies = { };
        std::map< std::string, json > bodySettings = settings[ "bodies" ];
        for ( auto ent : bodySettings )
        {
            bodies.insert( ent->first );
        }

        // Bodies to propagate
        std::vector< json > propagators = settings[ "propagators" ];
        bodiesToPropagate = { };
        for ( json propagator : propagators )
        {
            std::set< std::string > bodies = propagator[ "bodies" ];
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
        for ( std::string body : bodies )
        {
            if ( bodiesToPropagate.count( body ) == 0 )
            {
                celestialBodies.insert( body );
                if ( bodySettings[ body ][ "isCentralBody" ] )
                {
                    centralBodies.insert( body );
                }
            }
        }

    }


public:
    //! Set with filenames of the Spice kernels located in "Tudat/External/SpiceInterface/Kernels/"
    //! to be used for the propagation.
    std::set< std::string > spiceKernels;

    //! Initial epoch for the simulation.
    double startEpoch;

    //! Maximum end epoch for the simulation.
    double endEpoch;

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


    Simulation( std::string inputFile )
    {
        setInputFile( inputFile );
        setUp( );
    }


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

        // Read and parse input file and update members
        parseInputFile( );
    }


    //! Update/create all the objects needed before the simulation can be run.
    void setUp( )
    {
        loadSpiceKernels( );
        setUpEnvironment( );
        createBodies( );
        setUpPropagators( );
        setUpIntegrator( );
    }


    //! Load requested Spice kernels.
    void loadSpiceKernels( )
    {
        for ( std::string kernelFile : spiceKernels )
        {
            spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + kernelFile );
        }
    }


    void setUpEnvironment( )
    {

    }


    void createBodies( )
    {
        using namespace simulation_setup;

        // Get default settings for celestial bodies.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( celestialBodies, startEpoch - 300.0, endEpoch + 300.0 );

        // Overwrite values that were not set to be default.
        // ...

        // Set celestial bodies ephemeris and rotational models.
        for ( std::string celestialBody : celestialBodies )
        {
            bodySettings[ celestialBody ]->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
            bodySettings[ celestialBody ]->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
        }

        // Create celestial bodies.
        bodyMap = createBodies( bodySettings );

        // Create bodies to be propagated.
        for ( std::string bodyName : bodiesToPropagate )
        {
            bodyMap[ bodyName ] = boost::make_shared< simulation_setup::Body >( );
            bodyMap[ bodyName ]->setConstantBodyMass( settings[ "bodies" ][ bodyName ][ "mass" ] );
        }


    }


    void setUpPropagators( )
    {

    }


    void setUpIntegrator( )
    {

    }


    void run( )
    {

    }


    void exportResults( )
    {

    }


};


} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_H
