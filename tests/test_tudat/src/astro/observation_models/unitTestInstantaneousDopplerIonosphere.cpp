/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved. See LICENSE for the terms of redistribution.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createDopplerPartials.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;
using namespace simulation_setup;

struct IonosphericDopplerFixture {
    const double transmittedFrequency = 7.153e9;
    double initialTime;
    SystemOfBodies bodies;
    LinkEnds linkEnds;
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > ionosphere;

    IonosphericDopplerFixture( const int linkCase = 1 )
    {
        spice_interface::loadStandardSpiceKernels( );
        initialTime = spice_interface::convertDateStringToEphemerisTime( "2016-08-27 06:00:00 UTC" );
        const std::string transmittingStation = ( linkCase == 0 ) ? "DSS-55" : "DSS-14";
        const std::string receivingStation = ( linkCase == 0 ) ? "DSS-55" : "DSS-25";
        const std::string complex = ( linkCase == 0 ) ? "C60" : "C10";
        BodyListSettings bodySettings = getDefaultBodySettings( { "Sun", "Earth", "Moon" }, "SSB", "J2000" );
        bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
        bodySettings.addSettings( "Probe" );
        bodySettings.at( "Probe" )->ephemerisSettings =
                keplerEphemerisSettings( ( Eigen::Vector6d( ) << 5.2 * 1.495978707e11, 0.05, 0.02, 0.3, 0.4, 1.0 ).finished( ),
                                         initialTime,
                                         1.32712440018e20,
                                         "Sun",
                                         "J2000" );
        bodies = createSystemOfBodies( bodySettings );
        bodies.at( "Earth" )
                ->getGroundStation( transmittingStation )
                ->setTransmittingFrequencyCalculator(
                        std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmittedFrequency ) );
        bodies.at( "Probe" )->getVehicleSystems( )->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );

        // Intentionally do not install a transmitted-frequency calculator on Probe:
        // the coherent downlink must inherit its frequency from the complete link.
        linkEnds = { { transmitter, LinkEndId( "Earth", transmittingStation ) },
                     { reflector1, LinkEndId( "Probe" ) },
                     { receiver, LinkEndId( "Earth", receivingStation ) } };

        // Parse the same synthetic TRK-2-23 CHPART card as the original Python test,
        // using the in-memory CSP constructor to avoid temporary files.
        const auto card = std::make_shared< input_output::CspRawFile >(
                std::vector< std::string >{ "ADJUST(DOPRNG)BY NRMPOW(100.0, 1000.0) "
                                            "MODEL(CHPART)FROM(16/08/27,00:00)TO(16/08/28,00:00)DSN(" +
                                            complex + ")SCID(61)." } );
        ionosphere = { std::make_shared< TabulatedIonosphericCorrectionSettings >(
                extractIonosphericCorrection( { card }, { { 61, "Probe" } } ) ) };
    }
};

BOOST_AUTO_TEST_SUITE( test_instantaneous_doppler_ionosphere )

BOOST_AUTO_TEST_CASE( testTabulatedIonosphericDoppler )
{
    const double referenceFrequency = 2295.0e6;
    const double speedOfLight = physical_constants::SPEED_OF_LIGHT;

    // Both links pass through the same DSN complex. X/Ka also checks that the
    // downlink correction uses the transponder frequency, rather than the uplink frequency.
    for( int linkCase = 0; linkCase < 2; ++linkCase )
    {
        const FrequencyBands downlinkBand = ( linkCase == 0 ) ? x_band : ka_band;
        const double turnaroundRatio = ( linkCase == 0 ) ? 880.0 / 749.0 : 3344.0 / 749.0;

        IonosphericDopplerFixture fixture( linkCase );
        const auto& bodies = fixture.bodies;
        const auto& linkEnds = fixture.linkEnds;
        const auto& ionosphere = fixture.ionosphere;
        const double initialTime = fixture.initialTime;
        const double transmittedFrequency = fixture.transmittedFrequency;

        // rho = 100 + 1000*x metres at the reference frequency, where
        // x = 2*(t-t_start)/86400 - 1. Its constant slope yields an independent
        // first-order Doppler prediction. A growing group delay advances the phase.
        const double rangeDelayRate = 2.0 * 1000.0 / 86400.0;
        const double downlinkFrequency = transmittedFrequency * turnaroundRatio;
        const double expectedFrequencyShift = downlinkFrequency / speedOfLight * rangeDelayRate *
                ( std::pow( referenceFrequency / transmittedFrequency, 2 ) + std::pow( referenceFrequency / downlinkFrequency, 2 ) );

        for( const ObservableType observableType : { doppler_measured_frequency, two_way_doppler } )
        {
            BOOST_TEST_CONTEXT( "link case " << linkCase << ", observable " << observableType )
            {
                const auto correctedModel = ObservationModelCreator< 1, double, Time >::createObservationModel(
                        std::make_shared< ObservationModelSettings >( observableType, linkEnds, ionosphere ), bodies );
                const auto vacuumModel = ObservationModelCreator< 1, double, Time >::createObservationModel(
                        std::make_shared< ObservationModelSettings >( observableType, linkEnds ), bodies );
                const auto ancillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
                ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, { x_band, static_cast< double >( downlinkBand ) } );

                // The plain two-way observable is c times the fractional Doppler.
                const double convertToHertz = ( observableType == two_way_doppler ) ? downlinkFrequency / speedOfLight : 1.0;
                double maximumError = 0.0;
                for( int hour = 0; hour < 12; ++hour )
                {
                    const Time time = initialTime + hour * 3600.0;
                    const double frequencyShift = convertToHertz *
                            ( correctedModel->computeObservations( time, receiver, ancillarySettings )( 0 ) -
                              vacuumModel->computeObservations( time, receiver, ancillarySettings )( 0 ) );
                    maximumError = std::max( maximumError, std::fabs( frequencyShift - expectedFrequencyShift ) );
                    BOOST_TEST_CONTEXT( "hour " << hour )
                    {
                        // 115.2 mHz for X/X and 266.6 mHz for X/Ka. The 0.5 mHz
                        // tolerance includes the moving-link terms omitted above.
                        BOOST_CHECK_SMALL( frequencyShift - expectedFrequencyShift, 0.5e-3 );
                    }
                }
                BOOST_TEST_MESSAGE( "Link case " << linkCase << ", observable " << observableType << ": expected ionospheric shift "
                                                 << expectedFrequencyShift << " Hz, maximum error " << maximumError << " Hz" );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testStandaloneDopplerFrequencySettings )
{
    IonosphericDopplerFixture fixture;
    const auto uplink = getUplinkFromTwoWayLinkEnds( fixture.linkEnds );
    for( const auto observable : { one_way_doppler, one_way_doppler_measured_frequency } )
    {
        const auto settings = [ & ]( const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& corrections )
                -> std::shared_ptr< ObservationModelSettings > {
            return ( observable == one_way_doppler ) ? std::make_shared< ObservationModelSettings >( observable, uplink, corrections )
                                                     : oneWayDopplerMeasuredFrequencySettings( uplink, corrections );
        };
        const auto corrected =
                ObservationModelCreator< 1, double, Time >::createObservationModel( settings( fixture.ionosphere ), fixture.bodies );
        const auto vacuum = ObservationModelCreator< 1, double, Time >::createObservationModel( settings( {} ), fixture.bodies );
        const double toHertz = ( observable == one_way_doppler ) ? fixture.transmittedFrequency / physical_constants::SPEED_OF_LIGHT : 1.0;
        const double expectedShift = fixture.transmittedFrequency / physical_constants::SPEED_OF_LIGHT * ( 2000.0 / 86400.0 ) *
                std::pow( 2295.0e6 / fixture.transmittedFrequency, 2 );
        for( const auto reference : { receiver, transmitter } )
        {
            if( observable == one_way_doppler_measured_frequency && reference == transmitter )
            {
                continue;  // Measured frequency is defined at reception only.
            }
            BOOST_TEST_CONTEXT( "observable " << observable << ", reference " << reference )
            {
                const double implicit = corrected->computeObservations( Time( fixture.initialTime ), reference )( 0 );
                const auto settings = std::make_shared< ObservationAncillarySimulationSettings >( );
                BOOST_CHECK_EQUAL( corrected->computeObservations( Time( fixture.initialTime ), reference, settings )( 0 ), implicit );
                BOOST_CHECK_EQUAL( settings->getIntermediateDoubleData( transmitter_frequency_intermediate ),
                                   fixture.transmittedFrequency );
                BOOST_CHECK_SMALL(
                        toHertz * ( implicit - vacuum->computeObservations( Time( fixture.initialTime ), reference )( 0 ) ) - expectedShift,
                        0.5e-3 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testIonosphericDopplerEstimationScaling )
{
    using namespace observation_partials;
    using namespace estimatable_parameters;
    IonosphericDopplerFixture fixture;
    const auto parameters = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >{},
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{},
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{
                    std::make_shared< InitialTranslationalStateParameter< double > >(
                            "Probe", Eigen::Vector6d::Zero( ), "SSB", "J2000" ) } );
    const auto ancillary = std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillary->setAncillaryDoubleVectorData( frequency_bands, { x_band, ka_band } );
    for( const auto observable : { two_way_doppler, doppler_measured_frequency } )
    {
        const auto model = ObservationModelCreator< 1, double, Time >::createObservationModel(
                std::make_shared< ObservationModelSettings >( observable, fixture.linkEnds, fixture.ionosphere ), fixture.bodies );
        std::vector< double > times;
        std::vector< Eigen::Vector6d > states;
        const auto observation =
                model->computeObservationsWithLinkEndData( Time( fixture.initialTime ), receiver, times, states, ancillary );
        const auto partials = createTwoWayDopplerPartials< double, Time >( model, fixture.bodies, parameters );
        partials.second->update( states, times, receiver, observation );
        BOOST_REQUIRE( !partials.first.empty( ) );
        for( const auto& partial : partials.first )
        {
            const auto values = partial.second->calculatePartial( states, times, receiver, ancillary, observation );
            BOOST_REQUIRE( !values.empty( ) );
            for( const auto& value : values )
            {
                BOOST_CHECK( value.first.allFinite( ) );
            }
        }
        // Building/updating approximate partials must not remove corrections from the observation.
        BOOST_CHECK_EQUAL( model->computeObservations( Time( fixture.initialTime ), receiver, ancillary )( 0 ), observation( 0 ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
