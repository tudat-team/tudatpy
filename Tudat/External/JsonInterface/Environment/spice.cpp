/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "spice.h"

namespace tudat
{

namespace simulation_setup
{

/*
//! Get the set of spice kernels to be used for a SimulationType.
std::vector< boost::filesystem::path > getSpiceKernels( const SimulationType simulationType )
{
    std::vector< std::string > filenames;
    switch ( simulationType ) {
    case singlePerturbedBody:
        filenames = { "pck00009.tpc", "de-403-masses.tpc", "de421.bsp" };
    default:
        throw std::runtime_error( "Could not determine Spice kernels for the specified simulation type." );
    }

    std::vector< boost::filesystem::path > kernels;
    for ( const std::string filename : filenames )
    {
        kernels.push_back( boost::filesystem::path( input_output::getSpiceKernelPath( ) ) / filename );
    }
    return kernels;
}
*/


//! Create a `json` object from a shared pointer to a `SpiceSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< SpiceSettings >& spiceSettings )
{
    if ( ! spiceSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Spice;

    jsonObject[ K::kernels ] = spiceSettings->kernels_;
    jsonObject[ K::preloadKernels ] = spiceSettings->preloadKernels_;
    if ( spiceSettings->preloadKernels_ )
    {
        jsonObject[ K::preloadOffsets ] = spiceSettings->preloadOffsets_;
    }
}

//! Create a shared pointer to a `SpiceSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< SpiceSettings >& spiceSettings )
{
    using namespace json_interface;
    using K = Keys::Spice;

    /*
    const boost::shared_ptr< SimulationType > simulationType
            = getOptional< SimulationType >( jsonObject, SpecialKeys::root / Keys::simulationType );
    if ( simulationType )
    {
        spiceSettings = boost::make_shared< SpiceSettings >( *simulationType );
    }
    else
    {
    */
    spiceSettings = boost::make_shared< SpiceSettings >(
                getValue< std::vector< path > >( jsonObject, K::kernels ) );
    updateFromJSONIfDefined( spiceSettings->preloadKernels_, jsonObject, K::preloadKernels );
    if ( spiceSettings->preloadKernels_ )
    {
        spiceSettings->preloadOffsets_ = getValue( jsonObject, K::preloadOffsets, spiceSettings->preloadOffsets_ );
    }
    // }
}

} // namespace simulation_setup


namespace json_interface
{

//! Load in Tudat the Spice kernels specified in \p spiceSettings.
void loadSpiceKernels( const boost::shared_ptr< simulation_setup::SpiceSettings >& spiceSettings )
{
    spice_interface::clearSpiceKernels( );

    if ( spiceSettings )
    {
        for ( const path kernel : spiceSettings->kernels_ )
        {
            spice_interface::loadSpiceKernelInTudat( kernel.string( ) );
        }
    }
}

} // namespace json_interface

} // namespace tudat
