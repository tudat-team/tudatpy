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

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "Support/modular.h"
#include "Support/valueAccess.h"
#include "Support/valueConversions.h"

#include "Environment/spice.h"
#include "Environment/body.h"
#include "Propagation/variable.h"
#include "Propagation/propagator.h"
#include "Mathematics/integrator.h"
#include "Propagation/export.h"
#include "Support/validation.h"

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

    //! -DOC
    Simulation( const std::string& inputFile )
    {
        setInputFile( inputFile );
    }

    //! -DOC
    void setInputFile( const std::string& inputFile )
    {
        inputFilePath = getPathForJSONFile( inputFile );
        boost::filesystem::current_path( inputFilePath.parent_path( ) );

        jsonObject = getParsedModularJSON( inputFilePath );
        originalJsonObject = jsonObject;

        // std::cout << originalJsonObject.dump( 2 ) << std::endl;
        // throw;

        updateSettingsFromJSONObject( );
    }

    //! Synchronize JSON object and members.
    /*!
     * Synchronize JSON object and members.
     *
     * To be called by the user after he modifies any of the public members manually.
     * Only necessary if then he needs to access some specific settings that may depend on the modified settings.
     * If the user modifies the settings, it is not necessary to call this method:
     * - Before running the simulation (the `run` mehtod always calls `sync` before starting the propagation).
     * - Before converting the simulation to `json` (`sync` will be called before accessing the object's settings).
     *
     * This is how this method works:
     * 1. Update `jsonObject` from the data contained in the public settings objects.
     * 2. Update the settings objects from the data contained in `jsonObject`.
     *
     * FIXME: If startEpoch is read by IntegratorSettings, then startEpoch changed, after the sync
     * IntegratorSettings will have the old start epoch... Because it is not undefined inside IntegratorSettings
     * anymore, so it will not try to access the key startEpoch of the root JSON object...
     */
    void sync( )
    {
        updateJSONObjectFromSettings( );
        updateSettingsFromJSONObject( );
    }

    //! -DOC
    virtual void run( )
    {
        // FIXME: MultiArc
        using namespace propagators;

        // FIXME: ? sync( );

        // Create simulation object
        dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodyMap, integratorSettings, propagationSettings, false );

        dynamicsSimulator->integrateEquationsOfMotion( propagationSettings->getInitialStates( ) );
    }

    //! -DOC
    virtual void exportResults( )
    {
        // FIXME: MultiArc
        using namespace propagators;
        using namespace input_output;
        using namespace simulation_setup;

        boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > singleArcDynamicsSimulator
                = boost::dynamic_pointer_cast< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    dynamicsSimulator );
        if ( singleArcDynamicsSimulator )
        {
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > statesHistory =
                    singleArcDynamicsSimulator->getEquationsOfMotionNumericalSolution( );
            std::map< TimeType, Eigen::VectorXd > dependentVariables =
                    singleArcDynamicsSimulator->getDependentVariableHistory( );

            // Cast TimeType and StateVariableType to double (only for exporting as table columns)
            std::map< TimeType, Eigen::VectorXd > epochs, states;
            for ( auto entry : statesHistory )
            {
                epochs[ entry.first ] = ( Eigen::VectorXd( 1 ) << static_cast< double >( entry.first ) ).finished( );
                states[ entry.first ] = entry.second.template cast< double >( );
            }

            for ( boost::shared_ptr< ExportSettings > exportSettings : exportSettingsVector )
            {
                if ( exportSettings->onlyFinalStep )
                {
                    reduceToLast( epochs );
                    reduceToLast( states );
                    reduceToLast( dependentVariables );
                }

                // Determine number of columns (not including first column = epoch).
                unsigned int cols = 0;
                for ( boost::shared_ptr< VariableSettings > variable : exportSettings->variables )
                {
                    switch ( variable->variableType_ )
                    {
                    case epochVariable:
                    {
                        cols += 1;
                        break;
                    }
                    case stateVariable:
                    {
                        cols += 6;  // FIXME: will depend on which states have been requested (satellite's? also mass?)
                        break;
                    }
                    case dependentVariable:
                    {
                        const boost::shared_ptr< SingleDependentVariableSaveSettings > depVariable =
                                boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
                        enforceNonNullPointer( depVariable );
                        cols += getDependentVariableSize( depVariable->dependentVariableType_ );
                        break;
                    }
                    default:
                    {
                        throw std::runtime_error( "Could not export variable of unsupported type." );
                    }
                    }
                }

                // Concatenate requested results
                std::map< TimeType, Eigen::VectorXd > results;
                for ( auto entry : epochs )
                {
                    unsigned int currentIndex = 0;

                    const TimeType epoch = entry.first;
                    Eigen::VectorXd result( cols );
                    for ( boost::shared_ptr< VariableSettings > variable : exportSettings->variables )
                    {
                        unsigned int size;
                        switch ( variable->variableType_ )
                        {
                        case epochVariable:
                        {
                            size = 1;
                            result.segment( currentIndex, size ) = entry.second;
                            break;
                        }
                        case stateVariable:
                        {
                            size = 6;  // FIXME
                            result.segment( currentIndex, size ) = states.at( epoch );
                            break;
                        }
                        case dependentVariable:
                        {
                            const boost::shared_ptr< SingleDependentVariableSaveSettings > depVariable =
                                    boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
                            enforceNonNullPointer( depVariable );
                            size = getDependentVariableSize( depVariable->dependentVariableType_ );
                            unsigned int index =
                                    getKeyWithValue( singleArcDynamicsSimulator->getDependentVariableIds( ),
                                                     getDependentVariableId( depVariable ) );
                            result.segment( currentIndex, size ) =
                                    dependentVariables.at( epoch ).segment( index, size );
                            break;
                        }
                        default:
                            throw std::runtime_error( "Could not export variable of unsupported type." );
                        }
                        currentIndex += size;
                    }
                    results[ epoch ] = result;
                }

                // Write propagation history to file.
                writeDataMapToTextFile( results,
                                        exportSettings->outputFile.filename( ).string( ),
                                        exportSettings->outputFile.parent_path( ).string( ),
                                        "",
                                        exportSettings->numericalPrecision,
                                        exportSettings->numericalPrecision,
                                        "," );
            }
        }
    }

    //! -DOC
    json getAsJSON( )
    {
        // sync( );
        updateJSONObjectFromSettings( );
        return jsonObject;
    }

    //! -DOC
    void exportAsJSON( const path& exportPath, const unsigned int tabSize = 2 )
    {
        std::ofstream outputFile( exportPath.string( ) );
        outputFile << getAsJSON( ).dump( tabSize );
        outputFile.close( );
    }

    //! -DOC
    void exportAsJSON( const unsigned int tabSize = 2 )
    {
        std::string filename = inputFilePath.stem( ).string( ) + ".populated" + inputFilePath.extension( ).string( );
        exportAsJSON( inputFilePath.parent_path( ) / filename, tabSize );
    }

    //! -DOC
    json getOriginalJSONObject( )
    {
        return originalJsonObject;
    }


    /// Data contained in the JSON file:

    //! Initial epoch for the simulation.
    TimeType startEpoch;

    //! Maximum end epoch for the simulation.
    TimeType endEpoch;

    //! Global frame origin.
    std::string globalFrameOrigin;

    //! Global frame orientation.
    std::string globalFrameOrientation;

    //! Spice settings ( NULL if Spice is not used ).
    boost::shared_ptr< simulation_setup::SpiceSettings > spiceSettings;

    //! Map of body settings.
    std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap;

    //! Map of variable settings.
    std::unordered_map< std::string, boost::shared_ptr< propagators::VariableSettings > > variableSettingsMap;

    //! Propagation settings.
    boost::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagationSettings;

    //! Integrator settings.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings;

    //! Vector of export settings.
    std::vector< boost::shared_ptr< simulation_setup::ExportSettings > > exportSettingsVector;

    //! Validation settings.
    boost::shared_ptr< ValidationSettings > validationSettings;


protected:

    //! -DOC
    virtual void resetGeneral( )
    {
        // Start and end epochs
        startEpoch = getEpoch< TimeType >( jsonObject, Keys::startEpoch );
        endEpoch = getEpoch< TimeType >( jsonObject, Keys::endEpoch );

        // Global frame origin and orientation
        globalFrameOrigin = getValue< std::string >( jsonObject, Keys::globalFrameOrigin );
        globalFrameOrientation = getValue< std::string >( jsonObject, Keys::globalFrameOrientation );
    }

    //! -DOC
    virtual void resetSpice( )
    {
        spiceSettings = NULL;
        updateFromJSONIfDefined( spiceSettings, jsonObject, Keys::spice );
        if ( spiceSettings )
        {
            spice_interface::clearSpiceKernels( );
            for ( const path kernel : spiceSettings->kernels_ )
            {
                spice_interface::loadSpiceKernelInTudat( kernel.string( ) );
            }
        }
    }

    //! -DOC
    virtual void resetBodies( )
    {
        using namespace simulation_setup;

        bodySettingsMap.clear( );

        std::map< std::string, json > jsonBodySettingsMap =
                getValue< std::map< std::string, json > >( jsonObject, Keys::bodies );

        std::vector< std::string > defaultBodyNames;
        for ( auto entry : jsonBodySettingsMap )
        {
            const std::string bodyName = entry.first;
            if ( getValue( jsonObject, Keys::bodies / bodyName / Keys::Body::useDefaultSettings, false ) )
            {
                defaultBodyNames.push_back( bodyName );
            }
        }

        // Create map with default body settings.
        if ( ! defaultBodyNames.empty( ) )
        {
            if ( spiceSettings )
            {
                if ( spiceSettings->preloadKernels_ )
                {
                    bodySettingsMap = getDefaultBodySettings( defaultBodyNames,
                                                              startEpoch + spiceSettings->preloadOffsets_.first,
                                                              endEpoch + spiceSettings->preloadOffsets_.second );
                }
                else
                {
                    bodySettingsMap = getDefaultBodySettings( defaultBodyNames );
                }
            }
            else
            {
                throw std::runtime_error(
                            "Could not get default bodies settings because no Spice settings were found." );
            }
        }

        // Get body settings from JSON.
        for ( auto entry : jsonBodySettingsMap )
        {
            const std::string bodyName = entry.first;
            const json jsonBodySettings = jsonBodySettingsMap[ bodyName ];
            if ( bodySettingsMap.count( bodyName ) )
            {
                // Reset ephemeris and rotational models frames.
                boost::shared_ptr< BodySettings >& bodySettings = bodySettingsMap[ bodyName ];
                bodySettings->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
                bodySettings->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
                // Update body settings from JSON.
                updateBodySettings( bodySettings, jsonBodySettings );
            }
            else
            {
                // Create body settings from JSON.
                bodySettingsMap[ bodyName ] = createBodySettings( jsonBodySettings );
            }
        }

        // Create bodies.
        bodyMap = createBodies( bodySettingsMap );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, globalFrameOrigin, globalFrameOrientation );
    }

    //! -DOC
    virtual void resetVariables( )
    {
        variableSettingsMap.clear( );
        updateFromJSONIfDefined( variableSettingsMap, jsonObject, Keys::variables );
        // FIXME: modify jsonObject to define "save" for the corresponding propagators
    }

    //! -DOC
    virtual void resetPropagation( )
    {
        updateFromJSON( propagationSettings, jsonObject, Keys::propagation );
        propagationSettings->createIntegratedStateModels( bodyMap );
    }

    //! -DOC
    virtual void resetIntegrator( )
    {
        updateFromJSON( integratorSettings, jsonObject, Keys::integrator );
    }

    //! -DOC
    virtual void resetExport( )
    {
        exportSettingsVector.clear( );
        updateFromJSONIfDefined( exportSettingsVector, jsonObject, Keys::xport );
    }

    //! -DOC
    virtual void resetValidation( )
    {
        validationSettings = boost::make_shared< ValidationSettings >( );
        updateFromJSONIfDefined( validationSettings, jsonObject, Keys::validation );
    }


private:

    //! Update the JSON object with all the data from the current settings (objests).
    void updateJSONObjectFromSettings( )
    {
        jsonObject = *this;
    }

    //! Update all the settings (objects) from the JSON object.
    void updateSettingsFromJSONObject( )
    {
        clearAccessHistory( );

        resetGeneral( );
        resetSpice( );
        resetBodies( );
        resetVariables( );
        resetPropagation( );
        resetIntegrator( );
        resetExport( );
        resetValidation( );

        checkUnusedKeys( jsonObject, validationSettings->unusedKey );
    }

    //! Absolute path to the input file.
    path inputFilePath;

    //! Original JSON object with all the settings read directly from the input file.
    json originalJsonObject;

    //! JSON object with the current settings.
    json jsonObject;

    //! Body map.
    simulation_setup::NamedBodyMap bodyMap;

    //! Dynamics simulator.
    boost::shared_ptr< propagators::DynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator;

};


//! Function to create a `json` object from a `Simulation` object.
template< typename TimeType = double, typename StateScalarType = double >
void to_json( json& jsonObject, const Simulation< TimeType, StateScalarType >& simulation )
{
    jsonObject.clear( );

    // assignIfNot( jsonObject, Keys::simulationType, simulation.type, customSimulation );
    jsonObject[ Keys::startEpoch ] = simulation.startEpoch;
    jsonObject[ Keys::endEpoch ] = simulation.endEpoch;
    jsonObject[ Keys::globalFrameOrigin ] = simulation.globalFrameOrigin;
    jsonObject[ Keys::globalFrameOrientation ] = simulation.globalFrameOrientation;
    assignIfNotNull( jsonObject, Keys::spice, simulation.spiceSettings );
    jsonObject[ Keys::bodies ] = simulation.bodySettingsMap;
    // assignIfNotEmpty( jsonObject, Keys::variables, simulation.variableSettingsMap );
    jsonObject[ Keys::propagation ] = simulation.propagationSettings;
    jsonObject[ Keys::integrator ] = simulation.integratorSettings;
    assignIfNotEmpty( jsonObject, Keys::xport, simulation.exportSettingsVector );
    jsonObject[ Keys::validation ] = simulation.validationSettings;
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SIMULATION_H
