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

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "jsonInterface.h"

#include "Body/body.h"
#include "integrator.h"

namespace tudat
{

namespace json_interface
{

typedef boost::filesystem::path path;


//! -DOC
template< typename TimeType = double, typename StateScalarType = double >
class Simulation
{
public:
    //! Vector with filenames of the Spice kernels located in "Tudat/External/SpiceInterface/Kernels/"
    //! to be used for the propagation.
    std::vector< std::string > spiceKernels;

    //! Whether to preload the ephemeris of the celestial bodies for the simulation period,
    //! or to retrieve this directly from Spice during the propagation at each integration step.
    //! Preloading Spice data generally results in faster propagations, unless:
    //! * The simulation ends much earlier than the specified maximum simulation end epoch.
    //! * The integrator step-size is very large (in the order of several hours or days).
    bool preloadSpiceData;

    //! Offset for the
    std::pair< TimeType, TimeType > spiceIntervalOffsets;

    //! Initial epoch for the simulation.
    TimeType startEpoch;

    //! Maximum end epoch for the simulation.
    TimeType endEpoch;

    //! Global frame origin.
    std::string globalFrameOrigin;

    //! Global frame orientation.
    std::string globalFrameOrientation;

    //! Vector with the names of all the bodies.
    std::vector< std::string > bodies;

    //! Vector with the names of the celestial bodies (handled by Spice).
    std::vector< std::string > celestialBodies;

    //! Vector with the names of the bodies to be propagated.
    std::vector< std::string > bodiesToPropagate;

    //! Body map.
    simulation_setup::NamedBodyMap bodyMap;

    //! Body settings.
    std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap;

    //! Integrator settings.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings;

    //! Integrator settings.
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< TimeType > > propagatorSettings;

    //! Dynamics simulator.
    boost::shared_ptr< propagators::DynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator;

    //! -DOC
    Simulation( const std::string& inputFile )
    {
        setInputFile( inputFile );
        reset( );
    }


    //! -DOC
    void setInputFile( std::string inputFile )
    {
        // Add .json extension if necessary
        if ( inputFile.find(".json") == std::string::npos )
        {
            inputFile += ".json";
        }

        // Determine the absolute path of the directory where the input file is located
        path inputDirectory = boost::filesystem::canonical( inputFile ).parent_path( );

        // Determine the filename (without file extension) of the input file
        // path inputFilename = path( inputFile ).stem();

        // Build absolute path
        inputPath = inputDirectory / path( inputFile ).filename();

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
        resetOutput( );
    }


    //! -DOC
    void resetGeneral( )
    {
        // Start and end epochs
        startEpoch = getNumber< TimeType >( settings, KeyTrees::Simulation::startEpoch );
        endEpoch = getNumber< TimeType >( settings, KeyTrees::Simulation::endEpoch );

        // Global frame origin and orientation
        globalFrameOrigin = getValue< std::string >( settings, KeyTrees::Simulation::globalFrameOrigin );
        globalFrameOrientation = getValue< std::string >( settings, KeyTrees::Simulation::globalFrameOrientation );
    }


    //! -DOC
    void resetSpice( )
    {
        spiceKernels = getValue< std::vector< std::string > >( settings, KeyTrees::Simulation::spiceKernels, { } );
        preloadSpiceData = getValue< bool >( settings, KeyTrees::Simulation::preloadSpiceData, true );
        spiceIntervalOffsets = preloadSpiceData ?
                    std::make_pair( -300.0, 300.0 ) : std::make_pair( TUDAT_NAN, TUDAT_NAN );

        // Clear all loaded kernels
        spice_interface::clearSpiceKernels( );

        // Load requested Spice kernels
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
        bodies.clear( );
        for ( auto ent : getValue< std::map< std::string, json > >( settings, Keys::bodies ) )
        {
            bodies.push_back( ent.first );
        }

        // Bodies to propagate
        std::vector< json > propagators = getValue< std::vector< json > >( settings, Keys::propagators );
        bodiesToPropagate.clear( );
        for ( json propagator : propagators )
        {
            std::vector< std::string > requestedBodies =
                    getValue< std::vector< std::string > >( propagator, Keys::Propagator::bodiesToPropagate );
            for ( std::string bodyName : requestedBodies )
            {
                if ( contains( bodies, bodyName ) )
                {
                    bodiesToPropagate.push_back( bodyName );
                }
                else
                {
                    throw std::runtime_error( "Requested to propagate body named \"" + bodyName +
                                              "\", but no body with this name was found." );
                }
            }
        }

        // Celestial bodies
        celestialBodies.clear( );
        for ( std::string bodyName : bodies )
        {
            if ( ! contains( bodiesToPropagate, bodyName ) )
            {
                celestialBodies.push_back( bodyName );
            }
        }

        // Get default settings for celestial bodies.
        bodySettingsMap = getDefaultBodySettings(
                    celestialBodies, startEpoch - spiceIntervalOffsets.first, endEpoch + spiceIntervalOffsets.second );

        // Set celestial bodies ephemeris and rotational models.
        for ( std::string celestialBody : celestialBodies )
        {
            bodySettingsMap[ celestialBody ]->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
            bodySettingsMap[ celestialBody ]->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
        }

        // Update celestial body settings from the JSON settings.
        for ( std::string bodyName : celestialBodies )
        {
            updateBodySettings( bodySettingsMap, bodyName, getValue< json >( settings, { Keys::bodies, bodyName } ) );
        }

        // Create celestial bodies.
        bodyMap = createBodies( bodySettingsMap );

        // Create rest of body settings (and corresponding bodies) once that celestial bodies have been created.
        for ( std::string bodyName : bodies )
        {
            if ( ! contains( celestialBodies, bodyName) )
            {
                updateBodySettings(
                            bodySettingsMap, bodyName, getValue< json >( settings, { Keys::bodies, bodyName } ) );
                addBody( bodyMap, bodyName, bodySettingsMap[ bodyName ] );
            }
        }

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, globalFrameOrigin, globalFrameOrientation );
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
        json jsonIntegratorSettings = getValue< json >( settings, Keys::integrator );
        if ( ! getNumberPointer< TimeType >( jsonIntegratorSettings, Keys::Integrator::initialTime ) )
        {
            jsonIntegratorSettings[ Keys::Integrator::initialTime ] = startEpoch;
        }
        integratorSettings = createIntegratorSettings< TimeType >( jsonIntegratorSettings );
    }


    //! -DOC
    void resetOutput( )
    {

    }


    //! -DOC
    void run( )
    {
        using namespace propagators;

        // Create simulation object
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, propagatorSettings, false );

        // FIXME: MultiArc

        boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > singleArcDynamicsSimulator
                = boost::dynamic_pointer_cast< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    dynamicsSimulator );
        if ( singleArcDynamicsSimulator )
        {
            singleArcDynamicsSimulator->integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
            return;
        }

        boost::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > multiArcDynamicsSimulator
                = boost::dynamic_pointer_cast< MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                    dynamicsSimulator );
        if ( multiArcDynamicsSimulator )
        {
            throw std::runtime_error( "MultiArcDynamicsSimulator not supported by JSON interface." );
        }
    }


    //! -DOC
    void exportResults( )
    {

    }


    //! -DOC
    json getOriginalSettings( )
    {
        return settings;
    }


protected:



private:
    //! Absolute path to the input file.
    path inputPath;

    //! JSON object with all the settings from the input file.
    json settings;


};


//! Function to create a `json` object from a `Simulation` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( simulation )`.
template< typename TimeType = double, typename StateScalarType = double >
void to_json( json& jsonObject, const Simulation< TimeType, StateScalarType >& simulation )
{
    // Initialise
    jsonObject = json( );

    // Simulation (general settings)
    json jsonSimulation;
    jsonSimulation[ Keys::Simulation::startEpoch ] = simulation.startEpoch;
    jsonSimulation[ Keys::Simulation::endEpoch ] = simulation.endEpoch;
    jsonSimulation[ Keys::Simulation::globalFrameOrigin ] = simulation.globalFrameOrigin;
    jsonSimulation[ Keys::Simulation::globalFrameOrientation ] = simulation.globalFrameOrientation;
    jsonSimulation[ Keys::Simulation::spiceKernels ] = simulation.spiceKernels;
    jsonSimulation[ Keys::Simulation::preloadSpiceData ] = simulation.preloadSpiceData;
    jsonObject[ Keys::simulation ] = jsonSimulation;

    // Bodies
    json jsonBodies;
    for ( auto entry : simulation.bodySettingsMap )
    {
        jsonBodies[ entry.first ] = entry.second;
    }
    jsonObject[ Keys::bodies ] = jsonBodies;

    // Integrator
    jsonObject[ Keys::integrator ] = simulation.integratorSettings;
}


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SIMULATION_H
