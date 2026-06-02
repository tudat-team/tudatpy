/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/observation_models/azimuthElevationObservationModel.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::unit_conversions;

BOOST_AUTO_TEST_SUITE( test_azimuth_elevation_model )

BOOST_AUTO_TEST_CASE( testAzimuthElevationModel )
{
    spice_interface::loadStandardSpiceKernels( );

    const double initialEphemerisTime = 0.0;
    const double finalEphemerisTime = 7.0 * 86400.0;
    const double buffer = 10.0 * 3600.0;

    for( bool useOblateEarthShape : { false, true } )
    {
        BOOST_TEST_CONTEXT( "Oblate Earth shape: " << useOblateEarthShape )
        {
            std::vector< std::string > bodiesToCreate = { "Earth", "Mars" };
            BodyListSettings defaultBodySettings =
                    getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
            if( useOblateEarthShape )
            {
                defaultBodySettings.at( "Earth" )->shapeModelSettings =
                        std::make_shared< OblateSphericalBodyShapeSettings >( 6378137.0, 1.0 / 298.257223563 );
            }
            else
            {
                defaultBodySettings.at( "Earth" )->shapeModelSettings =
                        std::make_shared< SphericalBodyShapeSettings >( spice_interface::getAverageRadius( "Earth" ) );
            }
            SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

            const std::string stationName = "AzElStation";
            createGroundStation( bodies.at( "Earth" ),
                                 stationName,
                                 ( Eigen::Vector3d( ) << 0.0, convertDegreesToRadians( 52.0 ), convertDegreesToRadians( 4.0 ) ).finished( ),
                                 coordinate_conversions::geodetic_position );

            LinkEnds linkEnds;
            linkEnds[ transmitter ] = LinkEndId( "Mars", "" );
            linkEnds[ receiver ] = LinkEndId( "Earth", stationName );

            Eigen::Vector2d absoluteBias = ( Eigen::Vector2d( ) << 3.2E-7, -1.5E-7 ).finished( );
            std::shared_ptr< ObservationModelSettings > idealObservableSettings =
                    azimuthElevationSettings( linkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), nullptr );
            std::shared_ptr< ObservationModelSettings > biasedObservableSettings =
                    azimuthElevationSettings( linkEnds,
                                              std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                              std::make_shared< ConstantObservationBiasSettings >( absoluteBias, true ) );

            std::shared_ptr< ObservationModel< 2, double, double > > idealObservationModel =
                    ObservationModelCreator< 2, double, double >::createObservationModel( idealObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 2, double, double > > biasedObservationModel =
                    ObservationModelCreator< 2, double, double >::createObservationModel( biasedObservableSettings, bodies );

            const double receiverObservationTime = 0.5 * ( finalEphemerisTime + initialEphemerisTime );
            std::vector< double > linkEndTimes;
            std::vector< Eigen::Vector6d > linkEndStates;
            Eigen::Vector2d idealObservation = idealObservationModel->computeIdealObservationsWithLinkEndData(
                    receiverObservationTime, receiver, linkEndTimes, linkEndStates );

            BOOST_CHECK_EQUAL( linkEndTimes.size( ), 2 );
            BOOST_CHECK_EQUAL( linkEndStates.size( ), 2 );

            std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
                    bodies.at( "Earth" )->getGroundStation( stationName )->getPointingAnglesCalculator( );
            Eigen::Vector3d vectorFromStationToTarget = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 );
            std::pair< double, double > referenceElevationAzimuth =
                    pointingAnglesCalculator->calculatePointingAngles( vectorFromStationToTarget, linkEndTimes.at( 1 ) );
            Eigen::Vector2d referenceObservation =
                    ( Eigen::Vector2d( ) << referenceElevationAzimuth.second, referenceElevationAzimuth.first ).finished( );

            const double comparisonTolerance = 10.0 * std::numeric_limits< double >::epsilon( );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( idealObservation, referenceObservation, comparisonTolerance );

            Eigen::Vector2d biasedObservation = biasedObservationModel->computeObservations( receiverObservationTime, receiver );
            Eigen::Vector2d expectedBiasedObservation = idealObservation + absoluteBias;
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( biasedObservation, expectedBiasedObservation, comparisonTolerance );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
