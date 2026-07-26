/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/io/readSp3File.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::shared_ptr< input_output::Sp3FileContents >& fileContents,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin,
                                                           const std::string& frameOrientation )
{
    if( fileContents == nullptr )
    {
        throw std::invalid_argument( "Cannot create SP3 ephemeris settings from null file contents." );
    }

    const auto satelliteIterator = fileContents->satelliteStates.find( satelliteId );
    if( satelliteIterator == fileContents->satelliteStates.end( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' is not present in the SP3 file." );
    }
    if( satelliteIterator->second.empty( ) )
    {
        throw std::invalid_argument( "Satellite '" + satelliteId + "' has no states in the SP3 file." );
    }

    std::map< double, Eigen::Vector6d > stateHistory;
    for( const auto& state : satelliteIterator->second )
    {
        if( state.second.size( ) != 6 )
        {
            throw std::runtime_error( "SP3 state for satellite '" + satelliteId + "' does not have six components." );
        }

        const Eigen::Vector6d fixedSizeState = state.second;
        if( !fixedSizeState.allFinite( ) )
        {
            throw std::runtime_error( "SP3 ephemeris settings require finite position and velocity records for satellite '" + satelliteId +
                                      "'. The selected SP3 file contains a missing value or does not include velocity records." );
        }
        stateHistory[ state.first ] = fixedSizeState;
    }

    const std::string resolvedFrameOrientation = frameOrientation.empty( ) ? fileContents->frameName : frameOrientation;
    if( resolvedFrameOrientation.empty( ) )
    {
        throw std::invalid_argument( "The SP3 file does not declare a reference frame; provide frameOrientation explicitly." );
    }

    return tabulatedEphemerisSettings( stateHistory, frameOrigin, resolvedFrameOrientation );
}

std::shared_ptr< EphemerisSettings > sp3EphemerisSettings( const std::string& fileName,
                                                           const std::string& satelliteId,
                                                           const std::string& frameOrigin,
                                                           const std::string& frameOrientation,
                                                           const double referenceJulianDay )
{
    return sp3EphemerisSettings( input_output::readSp3File( fileName, referenceJulianDay ), satelliteId, frameOrigin, frameOrientation );
}

}  // namespace simulation_setup

}  // namespace tudat
