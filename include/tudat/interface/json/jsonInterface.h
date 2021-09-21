/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/simulation.h"

#include "support/deserialization.h"
#include "support/valueAccess.h"
#include "support/valueConversions.h"

#include "tudat/interface/json/environment/spice.h"
#include "tudat/interface/json/environment/body.h"
#include "tudat/interface/json/propagation/propagator.h"
#include "tudat/interface/json/math/integrator.h"
#include "tudat/interface/json/propagation/export.h"
#include "support/options.h"

namespace tudat
{

namespace json_interface
{

enum JsonSimulationTypes
{
    equations_of_motion_propagation,
    variational_equations_propagation,
    state_estimation
};


//! Map of `ObservableType` string representations.
static std::map< JsonSimulationTypes, std::string > simulationTypes =
{
    { equations_of_motion_propagation, "EoM" },
    { variational_equations_propagation, "Variational" },
    { state_estimation, "Estimation" }
};

//! Map of `ObservableType` string representations.
static std::map< std::string, JsonSimulationTypes > simulationTypesInverse =
{
    { "EoM", equations_of_motion_propagation },
    { "Variational", variational_equations_propagation },
    { "Estimation", state_estimation }

};

//! Convert `ObservableType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const JsonSimulationTypes& simulationType )
{
    jsonObject = json_interface::stringFromEnum( simulationType, simulationTypes );
}

//! Convert `json` to `ObservableType`.
inline void from_json( const nlohmann::json& jsonObject, JsonSimulationTypes& simulationType )
{
    simulationType = json_interface::enumFromString( jsonObject, simulationTypes );
}

//! Class for managing JSON-based simulations.
template< typename TimeType = double, typename StateScalarType = double >
class JsonSimulationManager
{

public:
    //! Constructor from JSON file.
    /*!
     * Constructor.
     * \param inputFilePath Path to the root JSON input file. Can be absolute or relative (to the working directory).
     * \param initialClockTime Initial clock time from which the cumulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    JsonSimulationManager(
            const std::string& inputFilePath,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
        : initialClockTime_( initialClockTime )
    {
        resetInputFile( inputFilePath );
    }

    //! Constructor from JSON object.
    /*!
     * Constructor.
     * \param jsonObject The root JSON object.
     * \param initialClockTime Initial clock time from which the cumulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    JsonSimulationManager(
            const nlohmann::json& jsonObject,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
        : initialClockTime_( initialClockTime )
    {
        resetJsonObject( jsonObject );
    }

    virtual ~JsonSimulationManager( ){ }

    //! Reset the root JSON input file.
    /*!
     * Reset the root JSON input file.
     * \param inputFilePath Path to the root JSON input file. Can be absolute or relative (to the working directory).
     */
    void resetInputFile( const std::string& inputFilePath )
    {
        inputFilePath_ = getPathForJSONFile( inputFilePath );
        boost::filesystem::current_path( inputFilePath_.parent_path( ) );
        resetJsonObject( getDeserializedJSON( inputFilePath_ ) );
    }

    //! Reset the `json` object.
    /*!
     * Reset the `json` object.
     * \param jsonObject The new `json` object to be used for creating the simulation settings.
     */
    void resetJsonObject( const nlohmann::json& jsonObject )
    {
        jsonObject_ = jsonObject;
        originalJsonObject_ = jsonObject_;
    }

    virtual void createSimulationObjects( )
    {
        resetDynamicsSimulator( );
    }

    virtual void updateSettings( )
    {
        // Clear global variable keeping track of the keys that have been accessed
        clearAccessHistory( );

        if ( profiling )
        {
            std::cout << "parse: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }

        simulationType_ =
                simulationTypesInverse.at( getValue< std::string >( jsonObject_, Keys::simulationType, "EoM" ) );

        resetIntegratorSettings( );
        resetSpice( );
        resetBodies( );              // must be called after resetIntegratorSettings and resetSpice
        resetExportSettings( );
        resetPropagatorSettings( );  // must be called after resetBodies and resetExportSettings
        resetApplicationOptions( );
        if( simulationType_ == equations_of_motion_propagation )
        {
            createSimulationObjects( );
        }
    }

    //! Run the propagation.
    /*!
     * @copybrief runPropagation
     * <br/>
     * Before running the simulation, the JSON representation of `this` will be exported if requested in
     * applicationOptions_, and a message will be printed if requested in applicationOptions_.
     * <br/>
     * If some of the keys in jsonObject_ haven't been used, a message may be printed or an error may be thrown
     * depending on applicationOptions_.
     * <br/>
     * After running the simulation, a message will be printed if requested in applicationOptions_.
     */
    virtual void runPropagation( )
    {
        // Check if any keys in jsonObject_ haven't been used
        checkUnusedKeys( jsonObject_, applicationOptions_->unusedKey_ );

        if ( profiling )
        {
            std::cout << "checkUnusedKeys: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }

        // Export full settings JSON file if requested
        if ( ! applicationOptions_->fullSettingsFile_.empty( ) )
        {
            exportAsJson( applicationOptions_->fullSettingsFile_ );
        }

        // Print message on propagation start if requested
        if ( applicationOptions_->notifyOnPropagationStart_ )
        {
            std::cout << "Propagation of file " << inputFilePath_ << " started." << std::endl;
        }

        runJsonSimulation( );

        // Print message on propagation termination if requested
        if ( applicationOptions_->notifyOnPropagationTermination_ )
        {
            if ( dynamicsSimulator_->integrationCompletedSuccessfully( ) )
            {
                std::cout << "SUCCESS: propagation of file " << inputFilePath_ << " terminated with no errors."
                          << std::endl;
            }
            else
            {
                std::cout << "FAILURE: propagation of file " << inputFilePath_ << " terminated with errors."
                          << std::endl;
            }
        }

        if ( profiling )
        {
            std::cout << "run: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Export the results of the dynamics simulation according to the export settings.
    /*!
     * @copybrief exportResults
     */
    virtual void exportResults( )
    {
        if ( applicationOptions_->tagOutputFilesIfPropagationFails_ &&
             !dynamicsSimulator_->integrationCompletedSuccessfully( ) )
        {
            // Add header "FAILURE" to output files
            for ( std::shared_ptr< ExportSettings >& exportSettings : exportSettingsVector_ )
            {
                exportSettings->header_ = "FAILURE\n" + exportSettings->header_;
            }
        }

        exportResultsOfDynamicsSimulator( dynamicsSimulator_, exportSettingsVector_,
                                          !( simulationType_ == equations_of_motion_propagation ) );

        if ( profiling )
        {
            std::cout << "exportResults: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    virtual void runJsonSimulation( )
    {
        dynamicsSimulator_->integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
    }


    //! Export `this` as a `json` object.
    /*!
     * @copybrief getAsJson
     * \return JSON representation of `this`.
     */
    nlohmann::json getAsJson( )
    {
        updateJsonObjectFromSettings( );
        return jsonObject_;
    }

    //! Export `this` as a `json` object to the file \p exportPath.
    /*!
     * @copybrief exportAsJson
     * \param exportPath Path to which the contents of the JSON object are to be exported.
     * \param tabSize Size of tabulations in the exported file (default = 2, i.e. 2 spaces). If set to 0, no
     * indentantion or line breaks will be used, so the exported file will contain just one line.
     */
    void exportAsJson( const boost::filesystem::path& exportPath, const unsigned int tabSize = 2 )
    {
        if ( ! boost::filesystem::exists( exportPath.parent_path( ) ) )
        {
            boost::filesystem::create_directories( exportPath.parent_path( ) );
        }
        std::ofstream outputFile( exportPath.string( ) );
        outputFile << getAsJson( ).dump( tabSize );
        outputFile.close( );
    }

    //! Get original JSON object (defined at construction or last time setInputFile was called).
    /*!
     * @copybrief getOriginalJsonObject
     * \remark The returned JSON object is the result of combining all the JSON files that
     * may be included in the root JSON file into one single object containing all the settings, but before default
     * values have been loaded or unit conversions have been applied.
     * \return The original JSON object.
     */
    nlohmann::json getOriginalJsonObject( ) const
    {
        return originalJsonObject_;
    }

    //! Get the JSON object.
    /*!
     * @copybrief getJsonObject
     * \remark The returned JSON object is the result of combining all the JSON files that
     * may be included in the root JSON file into one single object containing all the settings.
     * \return The JSON object.
     */
    nlohmann::json getJsonObject( ) const
    {
        return jsonObject_;
    }

    nlohmann::json at( const std::string& key ) const
    {
        return valueAt( jsonObject_, KeyPath( key ) );
    }

    nlohmann::json& operator[] ( const std::string& key )
    {
        return valueAt( jsonObject_, key, true );
    }

    //! Get simulation start epoch.
    TimeType getStartEpoch( ) const
    {
        if ( integratorSettings_ )
        {
            return integratorSettings_->initialTime_;
        }
        else
        {
            return getValue< TimeType >( jsonObject_, { Keys::initialEpoch,
                                                        Keys::integrator / Keys::Integrator::initialTime } );
        }
    }

    //! Get maximum simulation end epoch. Returns `TUDAT_NAN` if there is no time termination condition.
    TimeType getEndEpoch( ) const
    {
        TimeType endEpoch = getTerminationEpoch< TimeType >( propagatorSettings_->getTerminationSettings( ) );
        if ( ! isNaN( endEpoch ) )
        {
            return endEpoch;
        }
        else
        {
            return getValue< TimeType >( jsonObject_, Keys::finalEpoch, TUDAT_NAN );
        }
    }

    //! Get integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getIntegratorSettings( ) const
    {
        return integratorSettings_;
    }

    //! Get Spice settings (nullptr if Spice is not used).
    std::shared_ptr< SpiceSettings > getSpiceSettings( ) const
    {
        return spiceSettings_;
    }

    //! Get global frame origin.
    std::string getGlobalFrameOrigin( ) const
    {
        return globalFrameOrigin_;
    }

    //! Get global frame orientation.
    std::string getGlobalFrameOrientation( ) const
    {
        return globalFrameOrientation_;
    }

    //! Get map of body settings.
    simulation_setup::BodyListSettings getBodySettingsMap( ) const
    {
        return bodySettingsMap_;
    }

    //! Get system of bodies.
    simulation_setup::SystemOfBodies getBodyMap( ) const
    {
        return bodies_;
    }

    //! Add a body named \p bodyName.
    void addBody( const std::string& bodyName )
    {
        bodies_.createEmptyBody( bodyName );
    }

    //! Get body named \p bodyName.
    std::shared_ptr< simulation_setup::Body > getBody( const std::string& bodyName ) const
    {
        return bodies_.at( bodyName );
    }

    //! Get propagator settings.
    std::shared_ptr< propagators::MultiTypePropagatorSettings< StateScalarType > > getPropagatorSettings( ) const
    {
        return propagatorSettings_;
    }

    //! Get vector of export settings (each element corresponds to an output file).
    std::vector< std::shared_ptr< ExportSettings > > getExportSettingsVector( ) const
    {
        return exportSettingsVector_;
    }

    //! Get export settings at \p index.
    std::shared_ptr< ExportSettings > getExportSettings( const unsigned int index ) const
    {
        return exportSettingsVector_.at( index );
    }

    //! Get application options.
    std::shared_ptr< ApplicationOptions > getApplicationOptions( ) const
    {
        return applicationOptions_;
    }

    //! Get dynamics simulator.
    std::shared_ptr< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulator( ) const
    {
        return dynamicsSimulator_;
    }


protected:

    bool profiling = false;

    //! Reset integratorSettings_ from the current jsonObject_.
    /*!
     * @copybrief resetIntegratorSettings
     */
    virtual void resetIntegratorSettings( )
    {
        updateFromJSON( this->integratorSettings_, jsonObject_, Keys::integrator );

        if ( profiling )
        {
            std::cout << "resetIntegratorSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset spiceSettings_ from the current jsonObject_.
    /*!
     * @copybrief resetSpice
     * Loads the requested kernels in Tudat (if any).
     */
    virtual void resetSpice( )
    {
        spiceSettings_ = nullptr;
        updateFromJSONIfDefined( spiceSettings_, jsonObject_, Keys::spice );
        loadSpiceKernels( spiceSettings_ );

        if ( profiling )
        {
            std::cout << "resetSpice: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset bodySettingsMap_ and bodies_ from the current jsonObject_.
    /*!
     * @copybrief resetBodies
     */
    virtual void resetBodies( )
    {
        globalFrameOrigin_ = getValue< std::string >( jsonObject_, Keys::globalFrameOrigin, "SSB" );
        globalFrameOrientation_ = getValue< std::string >( jsonObject_, Keys::globalFrameOrientation, "ECLIPJ2000" );
        updateBodiesFromJSON( jsonObject_, bodies_, bodySettingsMap_, globalFrameOrigin_, globalFrameOrientation_,
                              spiceSettings_, integratorSettings_ );

        if ( profiling )
        {
            std::cout << "resetBodies: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset exportSettingsVector_ from the current jsonObject_.
    /*!
     * @copybrief resetExportSettings
     */
    virtual void resetExportSettings( )
    {
        exportSettingsVector_.clear( );
        updateFromJSONIfDefined( exportSettingsVector_, jsonObject_, Keys::xport );

        if ( profiling )
        {
            std::cout << "resetExportSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset propagatorSettings_ from the current jsonObject_
    /*!
     * @copybrief resetPropagatorSettings
     * Tries to infer the initial states from the body ephemeris if not provided.
     * Creates the integrated state models using bodies_
     */
    virtual void
    resetPropagatorSettings( )
    {
        // Update jsonObject_ by determining initial states if not provided directly to the propagator settings:
        // * By obtaining the initial states from body properties (and transforming to Cartesian if necessary)
        // * By infering initial states from body ephemeris
        determineInitialStates< TimeType, StateScalarType >( jsonObject_, bodies_, integratorSettings_ );

        if ( profiling )
        {
            std::cout << "resetPropagatorSettings@determineInitialStates: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }

        // Update propagatorSettings_ from jsonObject_
        updateFromJSON( propagatorSettings_, jsonObject_ );
        if ( profiling )
        {
            std::cout << "resetExportSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
        if ( profiling )
        {
            std::cout << "resetPropagatorSettings@updateFromJSON: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }

        // Create integrated state models (acceleration, mass-rate, rotational models)
        propagatorSettings_->resetIntegratedStateModels( bodies_ );

        if ( profiling )
        {
            std::cout << "resetPropagatorSettings@resetIntegratedStateModels: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }

        // Update dependent variables to save
        bool printOutputVariables = false;
        try
        {
            printOutputVariables = getValue< bool >( jsonObject_, "printVariableTypes", false );
        }
        catch( std::runtime_error const& ){ }
        resetDependentVariableSaveSettings< StateScalarType >( propagatorSettings_, exportSettingsVector_, printOutputVariables );

        if ( profiling )
        {
            std::cout << "resetPropagatorSettings@resetDependentVariableSaveSettings: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset applicationOptions_ from the current jsonObject_.
    /*!
     * @copybrief resetApplicationOptions
     */
    virtual void resetApplicationOptions( )
    {
        applicationOptions_ = std::make_shared< ApplicationOptions >( );
        updateFromJSONIfDefined( applicationOptions_, jsonObject_, Keys::options );

        if ( profiling )
        {
            std::cout << "resetApplicationOptions: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Reset dynamicsSimulator_ for the current bodies_, integratorSettings_ and propagatorSettings_.
    /*!
     * @copybrief resetDynamicsSimulator
     */
    virtual void resetDynamicsSimulator( )
    {
        dynamicsSimulator_ =
                std::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies_, integratorSettings_, propagatorSettings_, false, false, false, false, initialClockTime_ );

        if ( profiling )
        {
            std::cout << "resetDynamicsSimulator: " << std::chrono::duration_cast< std::chrono::milliseconds >(
                             std::chrono::steady_clock::now( ) - initialClockTime_ ).count( ) * 1.0e-3 << " s" << std::endl;
            initialClockTime_ = std::chrono::steady_clock::now( );
        }
    }

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Spice settings (nullptr if Spice is not used).
    std::shared_ptr< SpiceSettings > spiceSettings_;

    JsonSimulationTypes simulationType_;

    //! Global frame origin.
    std::string globalFrameOrigin_;

    //! Global frame orientation.
    std::string globalFrameOrientation_;

    //! Map of body settings.
    simulation_setup::BodyListSettings bodySettingsMap_;

    //! Body map.
    simulation_setup::SystemOfBodies bodies_;

    //! Propagation settings.
    std::shared_ptr< propagators::MultiTypePropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Vector of export settings (each element corresponds to an output file).
    std::vector< std::shared_ptr< ExportSettings > > exportSettingsVector_;

    //! Application options.
    std::shared_ptr< ApplicationOptions > applicationOptions_;

    //! Dynamics simulator.
    std::shared_ptr< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //! Update the JSON object with all the data from the current settings (objests).
    void updateJsonObjectFromSettings( )
    {
        jsonObject_ = std::make_shared< JsonSimulationManager< TimeType, StateScalarType > >( *this );
    }

    //! Absolute path to the input file.
    std::chrono::steady_clock::time_point initialClockTime_;

    //! Absolute path to the input file.
    boost::filesystem::path inputFilePath_;

    //! Original JSON object with all the settings read directly from the input file.
    nlohmann::json originalJsonObject_;

    //! JSON object with the current settings.
    nlohmann::json jsonObject_;

};

extern template class JsonSimulationManager< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class JsonSimulationManager< double, long double >;
//extern template class JsonSimulationManager< Time, double >;
//extern template class JsonSimulationManager< Time, long double >;
#endif

//! Function to create a `json` object from a Simulation object.
template< typename TimeType, typename StateScalarType >
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< JsonSimulationManager< TimeType, StateScalarType > >& jsonSimulationManager )
{
    if ( ! jsonSimulationManager )
    {
        return;
    }

    // integrator
    jsonObject[ Keys::integrator ] = jsonSimulationManager->getIntegratorSettings( );

    // spice
    assignIfNotnullptr( jsonObject, Keys::spice, jsonSimulationManager->getSpiceSettings( ) );

    // bodies
    jsonObject[ Keys::globalFrameOrigin ] = jsonSimulationManager->getGlobalFrameOrigin( );
    jsonObject[ Keys::globalFrameOrientation ] = jsonSimulationManager->getGlobalFrameOrientation( );
    jsonObject[ Keys::bodies ] = jsonSimulationManager->getBodySettingsMap( );

    // export
    assignIfNotEmpty( jsonObject, Keys::xport, jsonSimulationManager->getExportSettingsVector( ) );

    // options
    jsonObject[ Keys::options ] = jsonSimulationManager->getApplicationOptions( );

    // propagation + termination + options.printInterval
    propagators::to_json( jsonObject, std::dynamic_pointer_cast<
                          propagators::SingleArcPropagatorSettings< StateScalarType > >( jsonSimulationManager->getPropagatorSettings( ) ) );
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_H
