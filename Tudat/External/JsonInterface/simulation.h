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

//! Class for JSON-based simulations.
template< typename TimeType = double, typename StateScalarType = double >
class Simulation
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param inputFile Path to the root JSON input file. Can be absolute or relative (to the working directory).
     * \param initialClockTime Initial clock time from which the cummulative CPU time during the propagation will be
     * computed. Default is the moment at which the constructor was called.
     */
    Simulation( const std::string& inputFile,
                const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) )
        : initialClockTime_( initialClockTime )
    {
        setInputFile( inputFile );
    }

    //! Set the root JSON input file.
    /*!
     * @copybrief setInputFile
     * \param inputFile Path to the root JSON input file. Can be absolute or relative (to the working directory).
     */
    void setInputFile( const std::string& inputFile )
    {
        inputFilePath_ = getPathForJSONFile( inputFile );
        boost::filesystem::current_path( inputFilePath_.parent_path( ) );

        jsonObject_ = getParsedModularJSON( inputFilePath_ );
        originalJsonObject_ = jsonObject_;

        // std::cout << originalJsonObject.dump( 2 ) << std::endl;
        // throw;

        // Clear global variable keeping track of the keys that have been accessed
        clearAccessHistory( );

        updateSettingsFromJSONObject( );
    }

    //! Run the simulation.
    /*!
     * @copybrief run
     * <br/>
     * Before running the simulation, the JSON representation of `this` will be exported if requested in
     * applicationOptions_, and a message will be printed if requested in applicationOptions_.
     * <br/>
     * If some of the keys in jsonObject_ haven't been used, a message may be printed or an error may be thrown
     * depending on applicationOptions_.
     * <br/>
     * After running the simulation, a message will be printed if requested in applicationOptions_.
     */
    virtual void run( )
    {
        // Check if any keys in jsonObject_ haven't been used
        checkUnusedKeys( jsonObject_, applicationOptions_->unusedKey_ );

        // Export populated JSON file if requested
        if ( ! applicationOptions_->populatedFile_.empty( ) )
        {
            exportAsJSON( applicationOptions_->populatedFile_ );
        }

        // Print message on propagation start if requested
        if ( applicationOptions_->notifyOnPropagationStart_ )
        {
            std::cout << "Propagation of file " << inputFilePath_ << " started." << std::endl;
        }

        // Run simulation
        dynamicsSimulator_->integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );

        // Print message on propagation termination if requested
        if ( applicationOptions_->notifyOnPropagationTermination_ )
        {
            std::cout << "Propagation of file " << inputFilePath_ << " terminated with no errors.\n" << std::endl;
        }
    }

    //! Export the results of the dynamics simulation according to the export settings.
    /*!
     * @copybrief exportResults
     */
    virtual void exportResults( )
    {
        exportResultsOfDynamicsSimulator( dynamicsSimulator_, exportSettingsVector_ );
    }

    //! Export `this` as a `json` object.
    /*!
     * @copybrief getAsJSON
     * \return JSON representation of `this`.
     */
    json getAsJSON( )
    {
        updateJSONObjectFromSettings( );
        return jsonObject_;
    }

    //! Export `this` as a `json` object to the file \p exportPath.
    /*!
     * @copybrief exportAsJSON
     * \param exportPath Path to which the contents of the JSON object are to be exported.
     * \param tabSize Size of tabulations in the exported file (default = 2, i.e. 2 spaces). If set to 0, no
     * indentantion or line breaks will be used, so the exported file will contain just one line.
     */
    void exportAsJSON( const path& exportPath, const unsigned int tabSize = 2 )
    {
        std::ofstream outputFile( exportPath.string( ) );
        outputFile << getAsJSON( ).dump( tabSize );
        outputFile.close( );
    }

    //! Get original JSON object (defined at construction or last time setInputFile was called).
    /*!
     * @copybrief getOriginalJSONObject
     * \remark The returned JSON object is the result of combining all the JSON files that
     * may be included in the root JSON file into one single object containing all the settings, but before default
     * values have been loaded or unit conversions have been applied.
     * \return The original JSON object.
     */
    json getOriginalJSONObject( )
    {
        return originalJsonObject_;
    }

    //! Synchronize JSON object and class members.
    /*!
     * @copybrief sync
     * <br/>
     * Call this method after modifying any of the public members manually and before:
     * <ul>
     *   <li>Accessing any specific setting that may depend on the modified settings.</li>
     *   <li>Calling the method run.</li>
     *   <li>Calling the method exportResults.</li>
     *   <li>Converting `this` to `json` by calling getAsJSON, exportAsJSON or using the `json()` constructor.</li>
     * </ul>
     */
    void sync( )
    {
        updateJSONObjectFromSettings( );
        updateSettingsFromJSONObject( );
    }


    // Settings read from the JSON file:

    //! Integrator settings.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Spice settings (NULL if Spice is not used).
    boost::shared_ptr< simulation_setup::SpiceSettings > spiceSettings_;

    //! Map of body settings.
    std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap_;

    //! Propagation settings.
    boost::shared_ptr< propagators::MultiTypePropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Vector of export settings (each element corresponds to an output file).
    std::vector< boost::shared_ptr< simulation_setup::ExportSettings > > exportSettingsVector_;

    //! Application options.
    boost::shared_ptr< ApplicationOptions > applicationOptions_;


protected:

    //! Reset integratorSettings_ from the current jsonObject_.
    /*!
     * @copybrief resetIntegratorSettings
     */
    virtual void resetIntegratorSettings( )
    {
        updateFromJSON( integratorSettings_, jsonObject_, Keys::integrator );
    }

    //! Reset spiceSettings_ from the current jsonObject_.
    /*!
     * @copybrief resetSpice
     * Loads the requested kernels in Tudat (if any).
     */
    virtual void resetSpice( )
    {
        spiceSettings_ = NULL;
        updateFromJSONIfDefined( spiceSettings_, jsonObject_, Keys::spice );
        loadSpiceKernels( spiceSettings_ );
    }

    //! Reset bodySettingsMap_ and bodyMap_ from the current jsonObject_.
    /*!
     * @copybrief resetBodies
     */
    virtual void resetBodies( )
    {
        updateBodiesFromJSON( jsonObject_, bodyMap_, bodySettingsMap_, spiceSettings_, integratorSettings_ );
    }

    //! Reset propagatorSettings_ from the current jsonObject_
    /*!
     * @copybrief resetPropagatorSettings
     * Tries to infer the initial states from the body ephemeris if not provided.
     * Creates the integrated state models using bodyMap_
     */
    virtual void resetPropagatorSettings( )
    {
        updateFromJSON( propagatorSettings_, jsonObject_ );

        // Infer initial states from body ephemeris if not provided
        inferInitialStatesIfNecessary( propagatorSettings_, bodyMap_, integratorSettings_ );

        // Create integrated state models (acceleration, mass-rate, rotational models)
        propagatorSettings_->createIntegratedStateModels( bodyMap_ );
    }

    //! Reset exportSettingsVector_ from the current jsonObject_.
    /*!
     * @copybrief resetExportSettings
     */
    virtual void resetExportSettings( )
    {
        exportSettingsVector_.clear( );
        updateFromJSONIfDefined( exportSettingsVector_, jsonObject_, Keys::xport );
    }

    //! Reset applicationOptions_ from the current jsonObject_.
    /*!
     * @copybrief resetApplicationOptions
     */
    virtual void resetApplicationOptions( )
    {
        applicationOptions_ = boost::make_shared< ApplicationOptions >( );
        updateFromJSONIfDefined( applicationOptions_, jsonObject_, Keys::options );
    }

    //! Reset dynamicsSimulator_ for the current bodyMap_, integratorSettings_ and propagatorSettings_.
    /*!
     * @copybrief resetDynamicsSimulator
     */
    virtual void resetDynamicsSimulator( )
    {
        // FIXME: MultiArc

        dynamicsSimulator_ =
                boost::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodyMap_, integratorSettings_, propagatorSettings_, false, false, false, initialClockTime_ );
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
        // Integrator settings first, because then integratorSettings_->initialTime is used by other methods
        resetIntegratorSettings( );
        resetSpice( );
        resetBodies( );
        resetPropagatorSettings( );
        resetExportSettings( );
        resetApplicationOptions( );
        resetDynamicsSimulator( );
    }

    //! Absolute path to the input file.
    std::chrono::steady_clock::time_point initialClockTime_;

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
template< typename TimeType, typename StateScalarType >
void to_json( json& jsonObject, const Simulation< TimeType, StateScalarType >& simulation )
{
    jsonObject.clear( );

    // assignIfNot( jsonObject, Keys::simulationType, simulation.type, customSimulation );

    // integrator
    jsonObject[ Keys::integrator ] = simulation.integratorSettings_;

    // spice
    assignIfNotNull( jsonObject, Keys::spice, simulation.spiceSettings_ );

    // bodies
    jsonObject[ Keys::bodies ] = simulation.bodySettingsMap_;

    // export
    assignIfNotEmpty( jsonObject, Keys::xport, simulation.exportSettingsVector_ );

    // options
    jsonObject[ Keys::options ] = simulation.applicationOptions_;

    // propagation + termination + options.printInterval
    propagators::to_json( jsonObject, simulation.propagatorSettings_ );
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SIMULATION_H
