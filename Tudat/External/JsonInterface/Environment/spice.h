/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SPICE_H
#define TUDAT_JSONINTERFACE_SPICE_H

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

/// SimulationType

//! Frequently-used simulations.
enum SimulationType
{
    customSimulation,
    singlePerturbedBody
};

//! Map of `SimulationType` string representations.
static std::map< SimulationType, std::string > simulationTypes =
{
    { customSimulation, "custom" },
    { singlePerturbedBody, "singlePerturbedBody" }
};

//! `SimulationType` not supported by `json_interface`.
static std::vector< SimulationType > unsupportedSimulationTypes = {  };

//! Convert `SimulationType` to `json`.
inline void to_json( json& jsonObject, const SimulationType& simulationType )
{
    jsonObject = json_interface::stringFromEnum( simulationType, simulationTypes );
}

//! Convert `json` to `SimulationType`.
inline void from_json( const json& jsonObject, SimulationType& simulationType )
{
    simulationType = json_interface::enumFromString( jsonObject.get< std::string >( ), simulationTypes );
}


/// SpiceSettings

//! Get the set of spice kernels to be used for a SimulationType.
std::vector< boost::filesystem::path > getSpiceKernels( const SimulationType simulationType );

class SpiceSettings
{
public:
    //! Constructor.
    SpiceSettings( const SimulationType simulationType ) :
        kernels_( getSpiceKernels( simulationType ) ) { }

    //! Constructor.
    SpiceSettings( const std::vector< boost::filesystem::path >& kernels, const double preloadOffset = 300.0 ) :
        kernels_( kernels ), preloadOffsets_( { -preloadOffset, preloadOffset } ) { }

    //! Constructor.
    SpiceSettings( const std::vector< boost::filesystem::path >& kernels, const bool preloadKernels ) :
        kernels_( kernels ), preloadKernels_( preloadKernels ) { }

    //! Vector containing the paths to the spice kernel files to be used.
    std::vector< boost::filesystem::path > kernels_;

    //! Whether all the data from the Spice kernels should be preloaded before the simulation for the interval start
    //! epoch to end epoch (true), or whether the data from Spice should be accessed on request at every step (false).
    //! Preloading Spice data generally results in faster propagations, unless:
    //! * The simulation ends much earlier than the specified maximum simulation end epoch.
    //! * The integrator step-size is very large (in the order of several hours or days).
    bool preloadKernels_ = true;

    //! Offsets for the interval for which the spice kernels are to be preloaded (only used if preloadKernels_ = true).
    //! The kernels will be preloaded for the interval:
    //! [ startEpoch + preloadOffsets_.first, endEpoch + preloadOffsets_.second ]
    std::pair< double, double > preloadOffsets_ = { -300.0, 300.0 };
};

//! Create a `json` object from a shared pointer to a `SpiceSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< SpiceSettings >& spiceSettings );

//! Create a shared pointer to a `SpiceSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< SpiceSettings >& spiceSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SPICE_H
