/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MRO_DSN_OBSERVATION_MODEL_TEST_HELPERS_H
#define TUDAT_MRO_DSN_OBSERVATION_MODEL_TEST_HELPERS_H

#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/basics/timeType.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/interface/spice/spiceInterface.h"

namespace tudat
{
namespace unit_tests
{
namespace mro_dsn_test
{

using namespace tudat::observation_models;
using namespace tudat::simulation_setup;

static const double mroTransponderDelay = 1.4149e-6;
static const std::string mroDsnDataDirectory = tudat::paths::getTudatTestDataPath( ) + "mro_dsn_observation_model/";
// These CSVs are generated from mromagr2012_076_0840xmmmv1.tnf at the Tudat trk234 converter boundary:
// each row is data passed to SingleObservationSet creation after selecting a one-hour MRO tracking interval.
static const std::string trk234InputsDirectory = mroDsnDataDirectory + "trk234_single_observation_set_inputs/";
// These CSVs replace large MRO SPICE spacecraft kernels in the unit tests. They contain Tudat-generated
// SPICE samples of the MRO state, RotationalEphemeris::getRotationStateVector, and HGA offset over the test interval.
static const std::string tabulatedEphemeridesDirectory = mroDsnDataDirectory + "tabulated_ephemerides/";

inline std::vector< std::string > splitCsvLine( const std::string& line )
{
    std::vector< std::string > tokens;
    std::stringstream lineStream( line );
    std::string token;
    while( std::getline( lineStream, token, ',' ) )
    {
        tokens.push_back( token );
    }
    return tokens;
}

inline std::vector< std::vector< std::string > > readCsvRows( const std::string& fileName )
{
    std::ifstream inputFile( fileName );
    if( !inputFile.good( ) )
    {
        throw std::runtime_error( "Could not open MRO DSN CSV fixture: " + fileName );
    }

    std::vector< std::vector< std::string > > rows;
    std::string line;
    std::getline( inputFile, line );
    while( std::getline( inputFile, line ) )
    {
        if( !line.empty( ) )
        {
            rows.push_back( splitCsvLine( line ) );
        }
    }
    return rows;
}

inline FrequencyBands frequencyBandFromString( const std::string& frequencyBand )
{
    if( frequencyBand == "S" )
    {
        return s_band;
    }
    if( frequencyBand == "X" )
    {
        return x_band;
    }
    if( frequencyBand == "Ka" )
    {
        return ka_band;
    }
    throw std::runtime_error( "Unsupported DSN frequency band in MRO fixture: " + frequencyBand );
}

inline std::map< double, Eigen::Vector6d > readCartesianStateHistoryFromCsv( const std::string& fileName )
{
    std::map< double, Eigen::Vector6d > stateHistory;
    for( const std::vector< std::string >& row : readCsvRows( fileName ) )
    {
        Eigen::Vector6d state;
        for( int i = 0; i < 6; ++i )
        {
            state( i ) = std::stod( row.at( i + 1 ) );
        }
        stateHistory[ std::stod( row.at( 0 ) ) ] = state;
    }
    return stateHistory;
}

inline std::map< double, Eigen::Vector7d > readRotationalStateHistoryFromCsv( const std::string& fileName )
{
    std::map< double, Eigen::Vector7d > stateHistory;
    for( const std::vector< std::string >& row : readCsvRows( fileName ) )
    {
        Eigen::Vector7d rotationalState;
        for( int i = 0; i < 7; ++i )
        {
            rotationalState( i ) = std::stod( row.at( i + 1 ) );
        }
        // The stored seven-vector is Tudat's own rotational state format: target-to-base quaternion
        // followed by target-frame angular velocity, as consumed by TabulatedRotationalEphemeris.
        stateHistory[ std::stod( row.at( 0 ) ) ] = rotationalState;
    }
    return stateHistory;
}

inline std::map< double, Eigen::Vector6d > readAntennaStateHistoryFromCsv( const std::string& fileName )
{
    std::map< double, Eigen::Vector6d > antennaStateHistory;
    for( const std::vector< std::string >& row : readCsvRows( fileName ) )
    {
        Eigen::Vector6d antennaState = Eigen::Vector6d::Zero( );
        for( int i = 0; i < 3; ++i )
        {
            antennaState( i ) = std::stod( row.at( i + 1 ) );
        }
        antennaStateHistory[ std::stod( row.at( 0 ) ) ] = antennaState;
    }
    return antennaStateHistory;
}

inline LinkEndId linkEndIdFromFixtureString( const std::string& linkEnd )
{
    if( linkEnd.rfind( "DSS-", 0 ) == 0 )
    {
        return LinkEndId( "Earth", linkEnd );
    }
    return LinkEndId( linkEnd );
}

inline LinkEnds linkEndsFromFixtureStrings( const std::vector< std::string >& row )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = linkEndIdFromFixtureString( row.at( 3 ) );
    if( !row.at( 4 ).empty( ) )
    {
        linkEnds[ reflector1 ] = linkEndIdFromFixtureString( row.at( 4 ) );
    }
    linkEnds[ receiver ] = linkEndIdFromFixtureString( row.at( 5 ) );
    return linkEnds;
}

inline std::vector< double > linkDelaysFromFixtureRow( const std::vector< std::string >& row )
{
    return { std::stod( row.at( 8 ) ), std::stod( row.at( 9 ) ), std::stod( row.at( 10 ) ) };
}

inline LinkEndType referenceLinkEndFromString( const std::string& referenceLinkEnd )
{
    if( referenceLinkEnd == "receiver" )
    {
        return receiver;
    }
    if( referenceLinkEnd == "transmitter" )
    {
        return transmitter;
    }
    throw std::runtime_error( "Unsupported reference link end in MRO fixture: " + referenceLinkEnd );
}

inline std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettingsFromFixtureRow( const ObservableType observableType,
                                                                                                  const std::vector< std::string >& row )
{
    std::vector< FrequencyBands > frequencyBands = { frequencyBandFromString( row.at( 6 ) ), frequencyBandFromString( row.at( 7 ) ) };
    const std::vector< double > linkDelays = linkDelaysFromFixtureRow( row );
    if( observableType == dsn_n_way_averaged_doppler )
    {
        return getDsnNWayAveragedDopplerAncillarySettings(
                frequencyBands, frequencyBandFromString( row.at( 7 ) ), std::stod( row.at( 13 ) ), std::stod( row.at( 12 ) ), linkDelays );
    }
    else if( observableType == dsn_n_way_range )
    {
        return getDsnNWayRangeAncillarySettings( frequencyBands, std::stod( row.at( 12 ) ), linkDelays );
    }
    throw std::runtime_error( "Unsupported observable type for MRO DSN fixture." );
}

inline std::shared_ptr< ObservationCollection< HighPrecisionStateScalar, Time > > createObservationCollectionFromTrk234Csv(
        const std::string& fileName,
        const ObservableType observableType )
{
    std::map< int, std::vector< std::vector< std::string > > > rowsPerSet;
    for( const std::vector< std::string >& row : readCsvRows( fileName ) )
    {
        rowsPerSet[ std::stoi( row.at( 0 ) ) ].push_back( row );
    }

    std::vector< std::shared_ptr< SingleObservationSet< HighPrecisionStateScalar, Time > > > observationSets;
    for( auto& setRows : rowsPerSet )
    {
        std::vector< Eigen::Matrix< HighPrecisionStateScalar, Eigen::Dynamic, 1 > > observations;
        std::vector< Time > observationTimes;
        for( const std::vector< std::string >& row : setRows.second )
        {
            Eigen::Matrix< HighPrecisionStateScalar, Eigen::Dynamic, 1 > observation( 1 );
            observation( 0 ) = std::stold( row.at( 2 ) );
            observations.push_back( observation );
            observationTimes.push_back( Time( std::stod( row.at( 1 ) ) ) );
        }

        observationSets.push_back( std::make_shared< SingleObservationSet< HighPrecisionStateScalar, Time > >(
                observableType,
                LinkDefinition( linkEndsFromFixtureStrings( setRows.second.front( ) ) ),
                observations,
                observationTimes,
                referenceLinkEndFromString( setRows.second.front( ).at( 11 ) ),
                std::vector< Eigen::VectorXd >( ),
                nullptr,
                ancillarySettingsFromFixtureRow( observableType, setRows.second.front( ) ) ) );
    }

    return std::make_shared< ObservationCollection< HighPrecisionStateScalar, Time > >( observationSets );
}

inline void loadMroSpiceKernels( )
{
    spice_interface::loadStandardSpiceKernels( );
}

inline SystemOfBodies createMroSystemOfBodies( const Time& initialTime, const Time& finalTime )
{
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    const std::string globalFrameOrigin = "SSB";
    const std::string globalFrameOrientation = "J2000";

    BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTime - 3600.0, finalTime + 3600.0, globalFrameOrigin, globalFrameOrientation, 120.0 );

    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings =
            gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, globalFrameOrientation );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    bodySettings.addSettings( "MRO" );
    bodySettings.at( "MRO" )->ephemerisSettings =
            tabulatedEphemerisSettings( readCartesianStateHistoryFromCsv( tabulatedEphemeridesDirectory + "mro_state_mars_j2000.csv" ),
                                        "Mars",
                                        globalFrameOrientation );
    bodySettings.at( "MRO" )->rotationModelSettings = tabulatedRotationSettings(
            readRotationalStateHistoryFromCsv( tabulatedEphemeridesDirectory + "mro_rotation_state_j2000_to_spacecraft.csv" ),
            globalFrameOrientation,
            "MRO_SPACECRAFT" );

    SystemOfBodies bodies = createSystemOfBodies< HighPrecisionStateScalar, Time >( bodySettings );
    bodies.at( "MRO" )->getVehicleSystems( )->setDefaultTransponderTurnaroundRatio( );
    return bodies;
}

inline void setRampFrequencyInterpolatorsInBodies( SystemOfBodies& bodies )
{
    std::map< std::string, std::vector< Time > > startTimesPerStation;
    std::map< std::string, std::vector< Time > > endTimesPerStation;
    std::map< std::string, std::vector< double > > rampRatesPerStation;
    std::map< std::string, std::vector< double > > startFrequenciesPerStation;

    for( const std::vector< std::string >& row : readCsvRows( trk234InputsDirectory + "ramp_frequency_interpolator_inputs.csv" ) )
    {
        const std::string station = row.at( 0 );
        startTimesPerStation[ station ].push_back( Time( std::stod( row.at( 2 ) ) ) );
        endTimesPerStation[ station ].push_back( Time( std::stod( row.at( 3 ) ) ) );
        startFrequenciesPerStation[ station ].push_back( std::stod( row.at( 4 ) ) );
        rampRatesPerStation[ station ].push_back( std::stod( row.at( 5 ) ) );
    }

    for( const auto& stationAndTimes : startTimesPerStation )
    {
        bodies.at( "Earth" )
                ->getGroundStation( stationAndTimes.first )
                ->setTransmittingFrequencyCalculator( std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                        stationAndTimes.second,
                        endTimesPerStation.at( stationAndTimes.first ),
                        rampRatesPerStation.at( stationAndTimes.first ),
                        startFrequenciesPerStation.at( stationAndTimes.first ),
                        ground_stations::extrapolate_at_gaps ) );
    }
}

inline std::shared_ptr< ephemerides::Ephemeris > createMroAntennaEphemeris( const Time& initialTime, const Time& finalTime )
{
    return std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, double > >(
            interpolators::createOneDimensionalInterpolator(
                    readAntennaStateHistoryFromCsv( tabulatedEphemeridesDirectory + "mro_antenna_position_spacecraft.csv" ),
                    interpolators::linearInterpolation( ) ),
            "-74000",
            "MRO_SPACECRAFT" );
}

inline void applyMroNotebookObservationCollectionPostProcessing(
        const std::shared_ptr< ObservationCollection< HighPrecisionStateScalar, Time > >& observationCollection,
        SystemOfBodies& bodies )
{
    observationCollection->setTransponderDelay( "MRO", mroTransponderDelay );
    std::pair< Time, Time > timeBounds = observationCollection->getTimeBounds( );
    observationCollection->setReferencePoint(
            bodies, createMroAntennaEphemeris( timeBounds.first, timeBounds.second ), "Antenna", "MRO", reflector1 );
}

inline std::vector< std::shared_ptr< LightTimeCorrectionSettings > > getMroDsnLightTimeCorrections(
        const bool includeSolarCoronaCorrection )
{
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections = {
        firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ),
        tabulatedTroposphericCorrectionSettings( { mroDsnDataDirectory + "mromagr2012_061_2012_092.tro" } ),
        tabulatedIonosphericCorrectionSettings( { mroDsnDataDirectory + "mromagr2012_061_2012_092.ion" }, { { 74, "MRO" } } )
    };

    if( includeSolarCoronaCorrection )
    {
        lightTimeCorrections.push_back( inversePowerSeriesSolarCoronaCorrectionSettings( { 1.70e12 }, { 2.44 } ) );
    }
    return lightTimeCorrections;
}

inline std::vector< std::shared_ptr< ObservationModelSettings > > createMroObservationModelSettings(
        const std::shared_ptr< ObservationCollection< HighPrecisionStateScalar, Time > >& observationCollection,
        const bool includeSolarCoronaCorrection )
{
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings;
    const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections =
            getMroDsnLightTimeCorrections( includeSolarCoronaCorrection );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable = observationCollection->getLinkEndsPerObservableType( );
    for( const auto& observableAndLinkEnds : linkEndsPerObservable )
    {
        for( const LinkEnds& linkEnds : observableAndLinkEnds.second )
        {
            if( observableAndLinkEnds.first == dsn_n_way_averaged_doppler )
            {
                observationModelSettings.push_back( dsnNWayAveragedDopplerObservationSettings(
                        LinkDefinition( linkEnds ), lightTimeCorrections, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ) ) );
            }
            else if( observableAndLinkEnds.first == dsn_n_way_range )
            {
                observationModelSettings.push_back( dsnNWayRangeObservationSettings( LinkDefinition( linkEnds ), lightTimeCorrections ) );
            }
        }
    }
    return observationModelSettings;
}

inline Eigen::VectorXd simulateAndGetResiduals(
        const std::shared_ptr< ObservationCollection< HighPrecisionStateScalar, Time > >& observationCollection,
        SystemOfBodies& bodies,
        const bool includeSolarCoronaCorrection )
{
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings =
            createMroObservationModelSettings( observationCollection, includeSolarCoronaCorrection );
    std::vector< std::shared_ptr< ObservationSimulatorBase< HighPrecisionStateScalar, Time > > > observationSimulators =
            createObservationSimulators< HighPrecisionStateScalar, Time >( observationModelSettings, bodies );
    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< HighPrecisionStateScalar, Time >( observationCollection, bodies );
    std::shared_ptr< ObservationCollection< HighPrecisionStateScalar, Time > > simulatedObservationCollection =
            simulateObservations< HighPrecisionStateScalar, Time >( observationSimulationSettings, observationSimulators, bodies );

    return ( simulatedObservationCollection->getConcatenatedObservations( ) - observationCollection->getConcatenatedObservations( ) )
            .template cast< double >( );
}

}  // namespace mro_dsn_test
}  // namespace unit_tests
}  // namespace tudat

#endif  // TUDAT_MRO_DSN_OBSERVATION_MODEL_TEST_HELPERS_H
