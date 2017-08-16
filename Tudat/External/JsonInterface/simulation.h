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
#include "Propagation/propagator.h"
#include "Mathematics/integrator.h"
#include "Propagation/export.h"
#include "Support/options.h"

namespace tudat
{

namespace json_interface
{

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
        inputFilePath_ = getPathForJSONFile( inputFile );
        boost::filesystem::current_path( inputFilePath_.parent_path( ) );

        jsonObject_ = getParsedModularJSON( inputFilePath_ );
        originalJsonObject_ = jsonObject_;

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
        using namespace propagators;

        // FIXME: ? sync( );

        if ( ! applicationOptions_->populatedFile_.empty( ) )
        {
            exportAsJSON( applicationOptions_->populatedFile_ );
        }

        // Create dynamics simulator object
        dynamicsSimulator_ = boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodyMap_, integratorSettings_, propagationSettings_, false );

        if ( applicationOptions_->notifyOnPropagationStart_ )
        {
            std::cout << "Propagation of file " << inputFilePath_ << " started." << std::endl;
        }

        // Run simulation
        dynamicsSimulator_->integrateEquationsOfMotion( propagationSettings_->getInitialStates( ) );

        if ( applicationOptions_->notifyOnPropagationTermination_ )
        {
            std::cout << "Propagation of file " << inputFilePath_ << " terminated with no errors.\n" << std::endl;
        }

        // FIXME: MultiArc
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
                    dynamicsSimulator_ );
        if ( singleArcDynamicsSimulator )
        {
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > statesHistory =
                    singleArcDynamicsSimulator->getEquationsOfMotionNumericalSolution( );
            std::map< TimeType, Eigen::VectorXd > dependentVariables =
                    singleArcDynamicsSimulator->getDependentVariableHistory( );

            for ( boost::shared_ptr< ExportSettings > exportSettings : exportSettingsVector_ )
            {
                // Determine number of columns (not including first column = epoch).
                unsigned int cols = 0;
                for ( boost::shared_ptr< VariableSettings > variable : exportSettings->variables )
                {
                    switch ( variable->variableType_ )
                    {
                    case independentVariable:
                    {
                        cols += 1;
                        break;
                    }
                    case stateVariable:
                    {
                        cols += statesHistory.begin( )->second.rows( );
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
                for ( auto it = statesHistory.begin( ); it != statesHistory.end( ); ++it )
                {
                    if ( ( it == statesHistory.begin( ) && ! exportSettings->onlyInitialStep &&
                           exportSettings->onlyFinalStep ) ||
                         ( it == --statesHistory.end( ) && ! exportSettings->onlyFinalStep &&
                           exportSettings->onlyInitialStep )
                         || ( ( it != statesHistory.begin( ) && it != --statesHistory.end( ) ) &&
                              ( exportSettings->onlyInitialStep || exportSettings->onlyFinalStep ) ) )
                    {
                        continue;
                    }
                    unsigned int currentIndex = 0;

                    const TimeType epoch = it->first;
                    Eigen::VectorXd result( cols );
                    for ( boost::shared_ptr< VariableSettings > variable : exportSettings->variables )
                    {
                        unsigned int size;
                        switch ( variable->variableType_ )
                        {
                        case independentVariable:
                        {
                            size = 1;
                            result.segment( currentIndex, size ) =
                                    ( Eigen::VectorXd( 1 ) << static_cast< double >( epoch ) ).finished( );
                            break;
                        }
                        case stateVariable:
                        {
                            size = it->second.rows( );
                            result.segment( currentIndex, size ) = it->second.template cast< double >( );
                            break;
                        }
                        case dependentVariable:
                        {
                            const boost::shared_ptr< SingleDependentVariableSaveSettings > depVariable =
                                    boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
                            enforceNonNullPointer( depVariable );
                            size = getDependentVariableSize( depVariable->dependentVariableType_ );
                            const std::string varID = getDependentVariableId( depVariable );
                            try
                            {
                                unsigned int index = getKeyWithValue(
                                            singleArcDynamicsSimulator->getDependentVariableIds( ), varID );
                                result.segment( currentIndex, size ) =
                                        dependentVariables.at( epoch ).segment( index, size );
                            }
                            catch ( ... )
                            {
                                std::cerr << "Could not export the results for variable \"" << varID << "\" "
                                          << "because none of the propagators was configured to compute this variable."
                                          << std::endl;
                                throw;
                            }
                            break;
                        }
                        default:
                            throw std::runtime_error( "Could not export variable of unsupported type." );
                        }
                        currentIndex += size;
                    }
                    results[ epoch ] = result;
                }

                if ( exportSettings->epochsInFirstColumn )
                {
                    // Write results map to file.
                    writeDataMapToTextFile( results,
                                            exportSettings->outputFile,
                                            "",
                                            exportSettings->numericalPrecision );
                }
                else
                {
                    // Write results matrix to file.
                    Eigen::MatrixXd resultsMatrix( results.size( ), cols );
                    int currentRow = 0;
                    for ( auto entry : results )
                    {
                        resultsMatrix.row( currentRow++ ) = entry.second.transpose( );
                    }
                    writeMatrixToFile( resultsMatrix,
                                       exportSettings->outputFile.filename( ).string( ),
                                       exportSettings->numericalPrecision,
                                       exportSettings->outputFile.parent_path( ) );
                }
            }
        }
    }

    //! -DOC
    json getAsJSON( )
    {
        // sync( );
        updateJSONObjectFromSettings( );
        return jsonObject_;
    }

    //! -DOC
    void exportAsJSON( const path& exportPath, const unsigned int tabSize = 2 )
    {
        std::ofstream outputFile( exportPath.string( ) );
        outputFile << getAsJSON( ).dump( tabSize );
        outputFile.close( );
    }

    //! -DOC
    json getOriginalJSONObject( )
    {
        return originalJsonObject_;
    }


    // Data contained in the JSON file:

    //! Integrator settings.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Spice settings ( NULL if Spice is not used ).
    boost::shared_ptr< simulation_setup::SpiceSettings > spiceSettings_;

    //! Map of body settings.
    std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap_;

    //! Propagation settings.
    boost::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagationSettings_;

    //! Vector of export settings.
    std::vector< boost::shared_ptr< simulation_setup::ExportSettings > > exportSettingsVector_;

    //! Application options.
    boost::shared_ptr< ApplicationOptions > applicationOptions_;


protected:

    //! -DOC
    virtual void resetIntegrator( )
    {
        updateFromJSON( integratorSettings_, jsonObject_, Keys::integrator );
    }

    //! -DOC
    virtual void resetSpice( )
    {
        spiceSettings_ = NULL;
        updateFromJSONIfDefined( spiceSettings_, jsonObject_, Keys::spice );
        if ( spiceSettings_ )
        {
            spice_interface::clearSpiceKernels( );
            for ( const path kernel : spiceSettings_->kernels_ )
            {
                spice_interface::loadSpiceKernelInTudat( kernel.string( ) );
            }
        }
    }

    //! -DOC
    virtual void resetBodies( )
    {
        using namespace simulation_setup;

        bodySettingsMap_.clear( );

        std::map< std::string, json > jsonBodySettingsMap =
                getValue< std::map< std::string, json > >( jsonObject_, Keys::bodies );

        std::vector< std::string > defaultBodyNames;
        for ( auto entry : jsonBodySettingsMap )
        {
            const std::string bodyName = entry.first;
            if ( getValue( jsonObject_, Keys::bodies / bodyName / Keys::Body::useDefaultSettings, false ) )
            {
                defaultBodyNames.push_back( bodyName );
            }
        }

        // Create map with default body settings.
        if ( ! defaultBodyNames.empty( ) )
        {
            if ( spiceSettings_ )
            {
                if ( spiceSettings_->preloadKernels_ )
                {
                    bodySettingsMap_ = getDefaultBodySettings( defaultBodyNames,
                                                               integratorSettings_->initialTime_
                                                               + spiceSettings_->preloadOffsets_.first,
                                                               getEpoch< TimeType >( jsonObject_, Keys::endEpoch )
                                                               + spiceSettings_->preloadOffsets_.second );
                }
                else
                {
                    bodySettingsMap_ = getDefaultBodySettings( defaultBodyNames );
                }
            }
            else
            {
                throw std::runtime_error(
                            "Could not get default bodies settings because no Spice settings were found." );
            }
        }

        // Global frame origin and orientation
        const std::string globalFrameOrigin = getValue< std::string >( jsonObject_, Keys::globalFrameOrigin, "SSB" );
        const std::string globalFrameOrientation =
                getValue< std::string >( jsonObject_, Keys::globalFrameOrientation, "ECLIPJ2000" );

        // Get body settings from JSON.
        for ( auto entry : jsonBodySettingsMap )
        {
            const std::string bodyName = entry.first;
            const json jsonBodySettings = jsonBodySettingsMap[ bodyName ];
            if ( bodySettingsMap_.count( bodyName ) )
            {
                boost::shared_ptr< BodySettings >& bodySettings = bodySettingsMap_[ bodyName ];
                // Reset ephemeris and rotational models frames.
                if ( bodySettings->ephemerisSettings )
                {
                    bodySettings->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
                }
                if ( bodySettings->rotationModelSettings )
                {
                    bodySettings->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
                }
                // Update body settings from JSON.
                updateBodySettings( bodySettings, jsonBodySettings );
            }
            else
            {
                // Create body settings from JSON.
                bodySettingsMap_[ bodyName ] = createBodySettings( jsonBodySettings );
            }
        }

        // Create bodies.
        bodyMap_ = createBodies( bodySettingsMap_ );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap_, globalFrameOrigin, globalFrameOrientation );
    }

    //! -DOC
    virtual void resetPropagation( )
    {
        using namespace propagators;

        updateFromJSON( propagationSettings_, jsonObject_, Keys::propagation );

        // Try to get initial states from bodies' epehemeris if not provided
        if ( propagationSettings_->getInitialStates( ).rows( ) == 0 )
        {
            boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagator =
                    boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                        propagationSettings_ );
            if ( translationalPropagator )
            {
                try
                {
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
                            getInitialStatesOfBodies( translationalPropagator->bodiesToIntegrate_,
                                                      translationalPropagator->centralBodies_,
                                                      bodyMap_, integratorSettings_->initialTime_ );
                    translationalPropagator->resetInitialStates( systemInitialState );
                }
                catch ( ... ) { }
            }
        }
        if ( propagationSettings_->getInitialStates( ).rows( ) == 0 )
        {
            std::cerr << "Could not determine initial integration state. "
                      << "Please provide it to the propagator settings using the key \""
                      << Keys::Propagator::initialStates << "\"." << std::endl;
            throw;
        }

        // Create integrated state models (acceleration, mass-rate, rotational models)
        propagationSettings_->createIntegratedStateModels( bodyMap_ );
    }

    //! -DOC
    virtual void resetExport( )
    {
        exportSettingsVector_.clear( );
        updateFromJSONIfDefined( exportSettingsVector_, jsonObject_, Keys::xport );
    }

    //! -DOC
    virtual void resetApplicationOptions( )
    {
        applicationOptions_ = boost::make_shared< ApplicationOptions >( );
        updateFromJSONIfDefined( applicationOptions_, jsonObject_, Keys::options );
    }


private:

    //! Update the JSON object with all the data from the current settings (objests).
    void updateJSONObjectFromSettings( )
    {
        jsonObject_ = *this;
    }

    //! Update all the settings (objects) from the JSON object.
    void updateSettingsFromJSONObject( )
    {
        clearAccessHistory( );

        resetIntegrator( );  // Integrator first, because then integrator->initialTime is used by other methods
        resetSpice( );
        resetBodies( );
        resetPropagation( );
        resetExport( );
        resetApplicationOptions( );

        checkUnusedKeys( jsonObject_, applicationOptions_->unusedKey_ );
    }

    //! Absolute path to the input file.
    path inputFilePath_;

    //! Original JSON object with all the settings read directly from the input file.
    json originalJsonObject_;

    //! JSON object with the current settings.
    json jsonObject_;

    //! Body map.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Dynamics simulator.
    boost::shared_ptr< propagators::DynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

};


//! Function to create a `json` object from a Simulation object.
template< typename TimeType = double, typename StateScalarType = double >
void to_json( json& jsonObject, const Simulation< TimeType, StateScalarType >& simulation )
{
    jsonObject.clear( );

    // assignIfNot( jsonObject, Keys::simulationType, simulation.type, customSimulation );
    jsonObject[ Keys::integrator ] = simulation.integratorSettings_;
    assignIfNotNull( jsonObject, Keys::spice, simulation.spiceSettings_ );
    jsonObject[ Keys::bodies ] = simulation.bodySettingsMap_;
    jsonObject[ Keys::propagation ] = simulation.propagationSettings_;
    assignIfNotEmpty( jsonObject, Keys::xport, simulation.exportSettingsVector_ );
    jsonObject[ Keys::options ] = simulation.applicationOptions_;
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SIMULATION_H
