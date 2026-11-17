/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/support/celestialFrameReferenceValues.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;

namespace
{

// Observed minus computed right ascension (times cos(declination)) and declination [rad]
Eigen::Vector2d getAngularPositionResidual( const Eigen::Vector2d& observed, const Eigen::Vector2d& computed )
{
    const double rightAscensionDifference = observed( 0 ) - computed( 0 );
    return Eigen::Vector2d(
            std::cos( observed( 1 ) ) * std::atan2( std::sin( rightAscensionDifference ), std::cos( rightAscensionDifference ) ),
            observed( 1 ) - computed( 1 ) );
}

// Right ascension [rad] from hours, minutes and seconds
double convertHoursMinutesSecondsToRadians( const int hours, const int minutes, const double seconds )
{
    return unit_conversions::convertDegreesToRadians( 15.0 * ( hours + minutes / 60.0 + seconds / 3600.0 ) );
}

// Declination [rad] from sign, degrees, arcminutes and arcseconds
double convertDegreesMinutesSecondsToRadians( const bool isNegative, const int degrees, const int arcMinutes, const double arcSeconds )
{
    const double declination = unit_conversions::convertDegreesToRadians( degrees + arcMinutes / 60.0 + arcSeconds / 3600.0 );
    return isNegative ? -declination : declination;
}

// Split a block of records (one per line) into its non-empty lines
std::vector< std::string > splitRecords( const std::string& records )
{
    std::vector< std::string > lines;
    std::istringstream recordStream( records );
    std::string line;
    while( std::getline( recordStream, line ) )
    {
        if( !line.empty( ) )
        {
            lines.push_back( line );
        }
    }
    return lines;
}

// Decode a record of an IMCCE planetary observation database URSS file, which contains (separated by blanks): planet, date
// (year-day-month) and time (h:m:s), time scale, right ascension (h m s) and O-C [s], declination (sign, deg ' ") and O-C ["],
// observatory, observation type, position type (Ast: astrometric), equinox, centre (Geo: geocentric), catalogue, and references.
// The record must be an astrometric, geocentric FK4 B1950 position, in ET. ET is taken as TT, and converted to TDB with Tudat.
// Returns the epoch [TDB seconds since J2000], and the B1950 right ascension and declination [rad].
Eigen::Vector3d decodeImcceAstrometricB1950Record(
        const std::string& record,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& timeScaleConverter )
{
    char planet[ 16 ], month[ 4 ], timeScale[ 4 ], declinationSign, observationType[ 8 ], positionType[ 8 ], equinox[ 8 ], centre[ 8 ],
            catalogue[ 8 ];
    int year, day, hour, minute, rightAscensionHours, rightAscensionMinutes, declinationDegrees, declinationMinutes, observatory;
    double second, rightAscensionSeconds, rightAscensionResidual, declinationSeconds, declinationResidual;
    const int numberOfDecodedFields = std::sscanf( record.c_str( ),
                                                   "%15s %d-%d-%3s %d:%d:%lf %3s %d %d %lf %lf %c%d %d %lf %lf %d %7s %7s %7s %7s %7s",
                                                   planet,
                                                   &year,
                                                   &day,
                                                   month,
                                                   &hour,
                                                   &minute,
                                                   &second,
                                                   timeScale,
                                                   &rightAscensionHours,
                                                   &rightAscensionMinutes,
                                                   &rightAscensionSeconds,
                                                   &rightAscensionResidual,
                                                   &declinationSign,
                                                   &declinationDegrees,
                                                   &declinationMinutes,
                                                   &declinationSeconds,
                                                   &declinationResidual,
                                                   &observatory,
                                                   observationType,
                                                   positionType,
                                                   equinox,
                                                   centre,
                                                   catalogue );
    BOOST_REQUIRE_EQUAL( numberOfDecodedFields, 23 );
    BOOST_REQUIRE_EQUAL( std::string( timeScale ), "ET" );
    BOOST_REQUIRE_EQUAL( std::string( positionType ), "Ast" );
    BOOST_REQUIRE_EQUAL( std::string( equinox ), "B1950" );
    BOOST_REQUIRE_EQUAL( std::string( centre ), "Geo" );
    BOOST_REQUIRE_EQUAL( std::string( catalogue ), "FK4" );

    const std::string monthNames = "JanFebMarAprMayJunJulAugSepOctNovDec";
    const std::size_t monthIndex = monthNames.find( month );
    BOOST_REQUIRE( monthIndex != std::string::npos && monthIndex % 3 == 0 );
    const double ephemerisTime =
            basic_astrodynamics::DateTime(
                    year, static_cast< int >( monthIndex / 3 ) + 1, day, hour, minute, static_cast< long double >( second ) )
                    .epoch< double >( );

    return Eigen::Vector3d(
            timeScaleConverter->getCurrentTime( basic_astrodynamics::tt_scale, basic_astrodynamics::tdb_scale, ephemerisTime ),
            convertHoursMinutesSecondsToRadians( rightAscensionHours, rightAscensionMinutes, rightAscensionSeconds ),
            convertDegreesMinutesSecondsToRadians( declinationSign == '-', declinationDegrees, declinationMinutes, declinationSeconds ) );
}

// Decode a record of the JPL planetary optical observation files, in the format of https://ssd.jpl.nasa.gov/dat/planets/optical_format.txt:
// UT of observation as Julian date (columns 5-16, D12.5), catalogue code (26-28), type (29), right ascension (34-42, hhmmss^sss)
// and its equinox code (43), and declination (51-59, +/-ddmmss^ss) and its equinox code (60). The record must be a geocentric
// transit observation (type 1) w.r.t. the true equator and equinox of date (equinox codes 0), reduced with the given catalogue.
// Blanks in numeric fields are read as zeros. The UT is taken as UTC, and converted to TDB with Tudat. Returns the epoch [TDB
// seconds since J2000], and the right ascension and declination [rad].
Eigen::Vector3d decodeJplTransitRecord( const std::string& record,
                                        const std::string& catalogueCode,
                                        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& timeScaleConverter )
{
    BOOST_REQUIRE( record.size( ) >= 60 );
    BOOST_REQUIRE_EQUAL( record.substr( 25, 3 ), catalogueCode );
    BOOST_REQUIRE_EQUAL( record.at( 28 ), '1' );
    BOOST_REQUIRE_EQUAL( record.at( 42 ), '0' );
    BOOST_REQUIRE_EQUAL( record.at( 59 ), '0' );

    std::string numericRecord = record;
    std::replace( numericRecord.begin( ), numericRecord.end( ), ' ', '0' );
    const double universalTime = ( std::stod( numericRecord.substr( 4, 12 ) ) * 1.0E-5 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;
    const std::string rightAscension = numericRecord.substr( 33, 9 );
    const std::string declination = numericRecord.substr( 50, 9 );
    return Eigen::Vector3d(
            timeScaleConverter->getCurrentTime( basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, universalTime ),
            convertHoursMinutesSecondsToRadians( std::stoi( rightAscension.substr( 0, 2 ) ),
                                                 std::stoi( rightAscension.substr( 2, 2 ) ),
                                                 std::stod( rightAscension.substr( 4, 5 ) ) * 1.0E-3 ),
            convertDegreesMinutesSecondsToRadians( declination.at( 0 ) == '-',
                                                   std::stoi( declination.substr( 1, 2 ) ),
                                                   std::stoi( declination.substr( 3, 2 ) ),
                                                   std::stod( declination.substr( 5, 4 ) ) * 1.0E-2 ) );
}

// Check pre-fit residuals w.r.t. INPOP19a of geocentric right ascension/declination observations of a planet, each given as epoch
// [TDB seconds since J2000], right ascension and declination [rad]. The primary and comparison observations (at the same epochs)
// are compared with the predictions using the primary and comparison ancillary settings, respectively. For both primary residual
// components (right ascension times cos(declination), and declination), the absolute mean must be below residualBounds( 0 ), and
// the RMS between residualBounds( 1 ) and residualBounds( 2 ) [arcsec]. For the norm of the difference between the primary and
// comparison residuals, the RMS must be between correctionBounds( 0 ) and correctionBounds( 1 ), and the maximum between
// correctionBounds( 2 ) and correctionBounds( 3 ) [arcsec]. If requirePrimaryBelowComparison is true, the primary residual RMS
// must be the smaller one, for both components.
void checkAngularPositionPrefitResiduals( const std::string& observedBody,
                                          const std::vector< Eigen::Vector3d >& primaryObservations,
                                          const std::shared_ptr< ObservationAncillarySimulationSettings >& primarySettings,
                                          const std::vector< Eigen::Vector3d >& comparisonObservations,
                                          const std::shared_ptr< ObservationAncillarySimulationSettings >& comparisonSettings,
                                          const Eigen::Vector3d& residualBounds,
                                          const Eigen::Vector4d& correctionBounds,
                                          const bool requirePrimaryBelowComparison )
{
    spice_interface::clearSpiceKernels( );
    spice_interface::loadStandardSpiceKernels( );
    BodyListSettings bodySettings( "SSB", "J2000" );
    for( const std::string& bodyName : { std::string( "Earth" ), observedBody } )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), std::string( "" ) );
    linkEnds[ transmitter ] = std::make_pair( observedBody, std::string( "" ) );
    const std::shared_ptr< ObservationModel< 2 > > model = ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( angular_position, linkEnds ), bodies );

    // Residuals and corrections, in arcsec
    BOOST_REQUIRE_EQUAL( primaryObservations.size( ), comparisonObservations.size( ) );
    const int numberOfObservations = primaryObservations.size( );
    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    Eigen::Matrix2Xd primaryResiduals( 2, numberOfObservations ), comparisonResiduals( 2, numberOfObservations );
    Eigen::VectorXd corrections( numberOfObservations );
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    for( int i = 0; i < numberOfObservations; i++ )
    {
        BOOST_REQUIRE_EQUAL( primaryObservations.at( i )( 0 ), comparisonObservations.at( i )( 0 ) );
        const Eigen::Vector2d primaryPrediction = model->computeObservationsWithLinkEndData(
                primaryObservations.at( i )( 0 ), receiver, linkEndTimes, linkEndStates, primarySettings );
        const Eigen::Vector2d comparisonPrediction = model->computeObservationsWithLinkEndData(
                comparisonObservations.at( i )( 0 ), receiver, linkEndTimes, linkEndStates, comparisonSettings );
        primaryResiduals.col( i ) =
                radiansToArcseconds * getAngularPositionResidual( primaryObservations.at( i ).segment< 2 >( 1 ), primaryPrediction );
        comparisonResiduals.col( i ) =
                radiansToArcseconds * getAngularPositionResidual( comparisonObservations.at( i ).segment< 2 >( 1 ), comparisonPrediction );
        corrections( i ) = ( primaryResiduals.col( i ) - comparisonResiduals.col( i ) ).norm( );
    }
    const double correctionRms = corrections.norm( ) / std::sqrt( numberOfObservations );
    BOOST_TEST_MESSAGE( observedBody << " RA cos(Dec), Dec mean [arcsec]: " << primaryResiduals.row( 0 ).mean( ) << ", "
                                     << primaryResiduals.row( 1 ).mean( )
                                     << "; RMS [arcsec]: " << primaryResiduals.row( 0 ).norm( ) / std::sqrt( numberOfObservations ) << ", "
                                     << primaryResiduals.row( 1 ).norm( ) / std::sqrt( numberOfObservations )
                                     << "; correction RMS/max [arcsec]: " << correctionRms << " / " << corrections.maxCoeff( ) );

    BOOST_CHECK_GT( correctionRms, correctionBounds( 0 ) );
    BOOST_CHECK_LT( correctionRms, correctionBounds( 1 ) );
    BOOST_CHECK_GT( corrections.maxCoeff( ), correctionBounds( 2 ) );
    BOOST_CHECK_LT( corrections.maxCoeff( ), correctionBounds( 3 ) );
    for( int i = 0; i < 2; i++ )
    {
        const double primaryRms = primaryResiduals.row( i ).norm( ) / std::sqrt( numberOfObservations );
        BOOST_CHECK_LT( std::abs( primaryResiduals.row( i ).mean( ) ), residualBounds( 0 ) );
        BOOST_CHECK_GT( primaryRms, residualBounds( 1 ) );
        BOOST_CHECK_LT( primaryRms, residualBounds( 2 ) );
        if( requirePrimaryBelowComparison )
        {
            BOOST_CHECK_GT( comparisonResiduals.row( i ).norm( ) / std::sqrt( numberOfObservations ), primaryRms );
        }
    }
}

// Convert a SOFA rotation matrix to an Eigen matrix
Eigen::Matrix3d convertSofaMatrix( double matrix[ 3 ][ 3 ] )
{
    return Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( &matrix[ 0 ][ 0 ] );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_angular_position_model )

BOOST_AUTO_TEST_CASE( testAngularPositionModel )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    for( int test = 0; test < 2; test++ )
    {
        // Create observation settings
        std::shared_ptr< ObservationModelSettings > observableSettings = std::make_shared< ObservationModelSettings >(
                angular_position,
                linkEnds,
                lightTimeCorrectionSettings,
                std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 3.2E-9, -1.5E-8 ).finished( ), true ),
                std::make_shared< LightTimeConvergenceCriteria >( ),
                basic_astrodynamics::tdb_scale,
                static_cast< bool >( test ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 2, double, double > > observationModel =
                ObservationModelCreator< 2, double, double >::createObservationModel( observableSettings, bodies );
        std::shared_ptr< ObservationBias< 2 > > observationBias = observationModel->getObservationBiasCalculator( );

        // Compute observation separately with two functions.
        double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        Eigen::Vector2d observationFromReceptionTime = observationModel->computeObservations( receiverObservationTime, receiver );
        Eigen::Vector2d observationFromReceptionTime2 =
                observationModel->computeObservationsWithLinkEndData( receiverObservationTime, receiver, linkEndTimes, linkEndStates );
        BOOST_CHECK_EQUAL( linkEndTimes.size( ), 2 );
        BOOST_CHECK_EQUAL( linkEndStates.size( ), 2 );

        // Manually create and compute light time corrections
        std::shared_ptr< LightTimeCorrection > lightTimeCorrectionCalculator =
                createLightTimeCorrections( lightTimeCorrectionSettings.at( 0 ), bodies, linkEnds, transmitter, receiver );
        double lightTimeCorrection = lightTimeCorrectionCalculator->calculateLightTimeCorrection(
                linkEndStates.at( 0 ), linkEndStates.at( 1 ), linkEndTimes.at( 0 ), linkEndTimes.at( 1 ) );

        // Check equality of computed observations.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

        // Check consistency of link end states and time.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStates.at( 0 ),
                                           bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) ),
                                           std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStates.at( 1 ),
                                           bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 1 ) ),
                                           std::numeric_limits< double >::epsilon( ) );

        // Check that reception time is kept fixed.
        BOOST_CHECK_CLOSE_FRACTION(
                static_cast< double >( receiverObservationTime ), linkEndTimes[ 1 ], std::numeric_limits< double >::epsilon( ) );

        // Manually compute light time
        Eigen::Vector3d positionDifference = ( linkEndStates[ 0 ] - linkEndStates[ 1 ] ).segment( 0, 3 );
        BOOST_CHECK_CLOSE_FRACTION( positionDifference.norm( ) / physical_constants::SPEED_OF_LIGHT + lightTimeCorrection,
                                    linkEndTimes[ 1 ] - linkEndTimes[ 0 ],
                                    std::numeric_limits< double >::epsilon( ) * 1000.0 );
        // Poor tolerance due to rounding errors when subtracting times

        // Check computed right ascension/declination from link end states
        //    Eigen::Vector3d sphericalRelativeCoordinates = coordinate_conversions::convertCartesianToSpherical(
        //                positionDifference );
        double rightAscension = 2.0 *
                std::atan( positionDifference[ 1 ] /
                           ( std::sqrt( positionDifference[ 0 ] * positionDifference[ 0 ] +
                                        positionDifference[ 1 ] * positionDifference[ 1 ] ) +
                             positionDifference[ 0 ] ) );
        double declination = mathematical_constants::PI / 2.0 - std::acos( positionDifference[ 2 ] / positionDifference.norm( ) );
        if( test == 0 )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    rightAscension +
                            observationBias->getObservationBias( std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                    observationFromReceptionTime( 0 ),
                    std::numeric_limits< double >::epsilon( ) );
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    rightAscension * std::cos( declination ) +
                            observationBias->getObservationBias( std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).x( ),
                    observationFromReceptionTime( 0 ),
                    std::numeric_limits< double >::epsilon( ) );
        }
        BOOST_CHECK_CLOSE_FRACTION(
                declination + observationBias->getObservationBias( std::vector< double >( ), std::vector< Eigen::Vector6d >( ) ).y( ),
                observationFromReceptionTime( 1 ),
                std::numeric_limits< double >::epsilon( ) );

        // Compute transmission time from light time.
        double transmitterObservationTime = receiverObservationTime - ( linkEndTimes[ 1 ] - linkEndTimes[ 0 ] );

        // Compare computed against returned transmission time.
        BOOST_CHECK_CLOSE_FRACTION(
                static_cast< double >( transmitterObservationTime ), linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );

        // Recompute observable while fixing transmission time.
        std::vector< double > linkEndTimes2;
        std::vector< Eigen::Vector6d > linkEndStates2;
        Eigen::Vector2d observationFromTransmissionTime = observationModel->computeObservations( transmitterObservationTime, transmitter );
        Eigen::Vector2d observationFromTransmissionTime2 = observationModel->computeObservationsWithLinkEndData(
                transmitterObservationTime, transmitter, linkEndTimes2, linkEndStates2 );

        // Compare results against those obtained when keeping reception fixed.
        for( unsigned int i = 0; i < 2; i++ )
        {
            BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromTransmissionTime2( i ), 1.0E-15 );
            BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromReceptionTime( i ), 1.0E-15 );
            BOOST_CHECK_CLOSE_FRACTION( linkEndTimes2.at( i ), linkEndTimes2.at( i ), 1.0E-15 );
        }

        // Check consistency of link end states and time.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStates2.at( 0 ),
                                           bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes2.at( 0 ) ),
                                           std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStates2.at( 1 ),
                                           bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( linkEndTimes2.at( 1 ) ),
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

// The tests below use 16 observations, evenly spread in time, of series from the IMCCE planetary observation database
// (https://ftp.imcce.fr/pub/databases/PODB/, URSS files) and the JPL planetary optical observation files
// (https://ssd.jpl.nasa.gov/dat/planets/), retrieved on 2026-10-31. The records are copied verbatim from these files (trailing blanks
// removed), at the line numbers given per test, and are decoded and converted here with Tudat only:
//  - IMCCE records: Ephemeris Time (ET) is taken as TT;
//  - JPL records: the UT (as Julian date in columns 5-16) is taken as UTC (|UT1 - UTC| < 0.9 s). Column 17, which contains a further
//    digit not defined in the format description, is not used (< 0.9 s).
// All records are flagged as geocentric by their source (IMCCE: 'Geo'; JPL: type code 1), so that the receiver is the centre of the
// Earth. The planet system barycentre of INPOP19a is used (offset from the planet centre < 0.1 arcsec). Light deflection and the
// planetary phase are not modelled. The residual bounds allow for the scatter of these historical observations (0.2 to 0.8 arcsec).

//! Validate the FK4 B1950 to ICRF conversion, and the constant B1950 frame rotation, on Nikolaev photographic observations of Saturn
BOOST_AUTO_TEST_CASE( testSaturnB1950AngularPositionPrefitResiduals )
{
    // Lines 1, 21, 42, 62, 83, 103, 124, 144, 164, 184, 205, 225, 245, 265, 286 and 306 of
    // https://ftp.imcce.fr/pub/databases/PODB/URSS/SATU/PHOTO/nikol.pht.ftp (Nikolaev, MPC 089)
    const std::string records = R"(
 Saturn   1961- 6-Aug 20:42:25.6012 ET   19 48  8.25900 -0.76500  -21 18 15.0000 -0.7700   89 Photo Ast B1950 Geo FK4 GC GC XX  68- 68
 Saturn   1963- 3-Oct 18:24:34.6890 ET   21 17 29.82000 -0.82500  -17  4 38.5100  0.0900   89 Photo Ast B1950 Geo FK4 GC GC XX  68- 68
 Saturn   1964-25-Oct 17:43: 0.4800 ET   22  3 48.42400 -0.27000  -13 41 38.3400 -0.1600   89 Photo Ast B1950 Geo FK4 NS EW XX  77- 77
 Saturn   1966-27-Aug 23:32:13.3440 ET   23 55 28.40800 -0.01500  - 3  8 51.3400 -1.1700   89 Photo Ast B1950 Geo FK4 NS EW XX  77- 77
 Saturn   1966-19-Sep 22:33:57.6000 ET   23 49 16.91500 -0.60000  - 3 51 30.9700 -1.2700   89 Photo Ast B1950 Geo FK4 NS EW XX  77- 77
 Saturn   1966-27-Dec 15:52:22.3680 ET   23 39 47.89900 -1.15500  - 4 37 38.2100 -0.6900   89 Photo Ast B1950 Geo FK4 NS EW XX  77- 77
 Saturn   1968-26-Oct 21:23:31.9212 ET    1 21 28.14000 -0.63000  + 5 35 17.9900  0.4000   89 Photo Ast B1950 Geo FK4 GC GC XX  76- 76
 Saturn   1970-23-Sep  1:26:31.1842 ET    3 21  3.63400 -0.39000  +15 54 37.9500 -0.0300   89 Photo Ast B1950 Geo FK4 GC GC XX  76- 76
 Saturn   1972- 6-Oct  2:23:31.3783 ET    5 18  8.27600 -0.21000  +21 27 15.4400  0.1700   89 Photo Ast B1950 Geo FK4 GC GC XX  75- 75
 Saturn   1973-28-Oct  1: 2:10.9732 ET    6 18 36.33600 -0.37500  +22 17  1.4400  0.2500   89 Photo Ast B1950 Geo FK4 GC GC XX  75- 75
 Saturn   1975-17-Feb 19:21: 3.7955 ET    6 52 38.41200 -0.54000  +22 35  7.2400 -0.2800   89 Photo Ast B1950 Geo FK4 GC GC XX  75- 75
 Saturn   1976-14-Mar 17:58:41.7363 ET    7 51 28.74700 -0.04500  +21 21  2.8400  0.1500   89 Photo Ast B1950 Geo FK4 GC GC XX  75- 75
 Saturn   1978-15-Jan  0:16:24.0321 ET   10  6 59.15400 -0.31500  +13  9 58.2200 -0.1800   89 Photo Ast B1950 Geo FK4 NS EW XX  79- 79
 Saturn   1979-10-May 18:31:43.6019 ET   10 36 35.80400 -0.22500  +10 53 28.6400  0.3000   89 Photo Ast B1950 Geo FK4 NS EW XX  79- 79
 Saturn   1981-24-Mar 21:44: 3.0633 ET   12 26 52.81000 -0.79500  - 0  1  1.7100 -0.0300   89 Photo Ast B1950 Geo FK4 GC GC XX  80- 80
 Saturn   1982-22-May 19:13:33.0165 ET   13  1 42.48200 -1.35000  - 3 40 10.2200 -0.0400   89 Photo Ast B1950 Geo FK4 GC GC XX  80- 80
)";
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );
    std::vector< Eigen::Vector3d > b1950Observations, icrfObservations;
    for( const std::string& record : splitRecords( records ) )
    {
        const Eigen::Vector3d b1950Observation = decodeImcceAstrometricB1950Record( record, timeScaleConverter );
        b1950Observations.push_back( b1950Observation );
        icrfObservations.push_back(
                ( Eigen::Vector3d( ) << b1950Observation( 0 ),
                  convertFk4B1950ToIcrfRightAscensionDeclination( b1950Observation( 1 ), b1950Observation( 2 ), b1950Observation( 0 ) ) )
                        .finished( ) );
    }
    BOOST_REQUIRE_EQUAL( b1950Observations.size( ), 16 );

    // Primary: observations converted to the ICRF, compared with J2000 (ICRF) astrometric directions. Comparison: B1950 observations,
    // compared with directions rotated with the constant B1950 frame rotation. The difference consists of the E-terms of aberration
    // (up to 0.34 arcsec) and the epoch-dependent FK4 equinox and proper-motion errors (0.085 s of right ascension per century since
    // 1950, so up to about 0.4 arcsec for these observations), and is required to be between 0.05 and 1.0 arcsec RMS (maximum
    // between 0.1 and 1.5 arcsec). As it is of the order of the scatter, the primary residuals are not required to be the smaller.
    checkAngularPositionPrefitResiduals( "Saturn_barycenter",
                                         icrfObservations,
                                         nullptr,
                                         b1950Observations,
                                         getPositionAngleAncillarySettings( b1950_position_angle_reference_frame ),
                                         { 0.5, 0.1, 1.0 },
                                         { 0.05, 1.0, 0.1, 1.5 },
                                         false );
}

//! Validate the FK4 B1950 to ICRF conversion on Pulkovo observations of Uranus over 50 years, against aberrated directions
BOOST_AUTO_TEST_CASE( testUranusB1950AngularPositionPrefitResiduals )
{
    // Lines 1, 9, 15, 23, 31, 38, 46, 54, 61, 69, 77, 84, 93, 101, 108 and 116 of
    // https://ftp.imcce.fr/pub/databases/PODB/URSS/URAN/PHOTO/pulk.pht.ftp (Pulkovo, MPC 084)
    const std::string records = R"(
 Uranus   1919-21-Sep 20:11:25.1621 ET   22  5 50.28000  1.03500  -12 32 48.2000  1.5900   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1933-14-Nov 20: 2:17.0578 ET    1 32  2.21000  0.61500  + 8 58 35.7000  0.7200   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1934-21-Sep 21:54:46.0227 ET    1 55 44.71000  0.16500  +11 14 49.1000  0.4600   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1936-20-Feb 16:40:17.0472 ET    2  1 43.02000  0.01500  +11 53  2.3000  0.2300   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1939-12-Oct 22:54:13.7068 ET    3 15 50.83000  0.22500  +17 46 19.5000  0.8200   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1941- 1-Feb 17:58:26.5964 ET    3 19 43.45000  0.63000  +18  6  5.0000  0.8200   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1949-16-Feb 17:41:40.8132 ET    5 45 19.39000  0.04500  +23 36 52.2000  0.3000   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1951- 4-Feb 20:17:59.5133 ET    6 26  8.61000  0.33000  +23 38 19.4000  0.7200   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1953-18-Feb 19:41:17.0961 ET    7  4 35.97000  0.21000  +23  3 26.6000  0.5800   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1955- 4-Apr 18:34:28.8471 ET    7 41 54.42000  0.39000  +21 55  8.8000  0.4000   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1958-10-Feb 22:37:39.2419 ET    8 46 24.63000  0.39000  +18 40  2.9000  0.7100   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1960-25-Mar 19:33: 7.9304 ET    9 19 22.60000  0.13500  +16 23 32.4000  0.7300   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1963-18-Feb 21:53:59.0989 ET   10 21 50.59000  0.36000  +11  4 14.4000  0.0800   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1965-28-Jan  1: 4:52.6881 ET   11  2 11.20000  0.19500  + 7  3 10.0000  0.0000   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1967-10-Feb  1:35:45.9715 ET   11 37  6.42000  0.67500  + 3 20 56.1000  0.5800   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
 Uranus   1969-12-Mar 22:45:38.8215 ET   12  9  1.66000  1.05000  - 0  7 55.4000  0.8600   84 Photo Ast B1950 Geo FK4 GC GC XX 118-118
)";
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );
    std::vector< Eigen::Vector3d > icrfObservations;
    for( const std::string& record : splitRecords( records ) )
    {
        const Eigen::Vector3d b1950Observation = decodeImcceAstrometricB1950Record( record, timeScaleConverter );
        icrfObservations.push_back(
                ( Eigen::Vector3d( ) << b1950Observation( 0 ),
                  convertFk4B1950ToIcrfRightAscensionDeclination( b1950Observation( 1 ), b1950Observation( 2 ), b1950Observation( 0 ) ) )
                        .finished( ) );
    }
    BOOST_REQUIRE_EQUAL( icrfObservations.size( ), 16 );

    // Primary: astrometric J2000 (ICRF) directions. Comparison: directions including stellar aberration, which differ by the annual
    // aberration (at most 20.8 arcsec): required to be between 10 and 21 arcsec RMS (maximum between 15 and 21 arcsec).
    checkAngularPositionPrefitResiduals( "Uranus_barycenter",
                                         icrfObservations,
                                         getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, false ),
                                         icrfObservations,
                                         getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, true ),
                                         { 1.0, 0.1, 1.5 },
                                         { 10.0, 21.0, 15.0, 21.0 },
                                         true );
}

//! Validate apparent true-of-date (IAU 1976/1980) angular positions on FK5-based USNO transit observations of Neptune
BOOST_AUTO_TEST_CASE( testNeptuneTrueOfDateIau1976AngularPositionPrefitResiduals )
{
    // Lines 1762, 1798, 1834, 1870, 1906, 1942, 1978, 2014, 2051, 2087, 2123, 2159, 2195, 2231, 2267 and 2303 of
    // https://ssd.jpl.nasa.gov/dat/planets/trnstswash6c.txt
    // (USNO 6-inch transit circle, 1985-1993): apparent places w.r.t. the true equator and equinox of date, reduced with FK5 (catalogue
    // code 995), so that the IAU 1976/1980 precession-nutation applies
    const std::string records = R"(
P0082446159939509    78669951    1815377770       -221503260            1985USNO
P0082446280602070    78669951    1805263990       -221719570            1985USNO
P0082446596746910    78669951    1820262870       -221459300            1986USNO
P0082446673531789    78669951    1813237600       -221946010            1986USNO
P0082446968734670    78669951    1829237980       -221300530            1987USNO
P0082447037541830    78669951    1823001740       -221845540            1987USNO
P0082447334739269    78669951    1839017260       -220832120            1988USNO
P0082447411524329    78669951    1832147840       -221616460            1988USNO
P0082448076720339    78669951    1857056970       -215452020            1990USNO
P0082448378903380    78669951    1912042660       -213539420            1991USNO
P0082448469649100    78669951    1903400950       -214924490            1991USNO
P0082448531477580    78669951    1900271970       -215608110            1991USNO
P0082448824684520    78669951    1914255200       -213507180            1992USNO
P0082448901471150    78669951    1909548780       -214504740            1992USNO
P0082449173736820    78669951    1925549360       -211639120            1993USNO
P0082449273459280    78669951    1919258450       -213155470            1993USNO
)";
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );
    std::vector< Eigen::Vector3d > observations;
    for( const std::string& record : splitRecords( records ) )
    {
        observations.push_back( decodeJplTransitRecord( record, "995", timeScaleConverter ) );
    }
    BOOST_REQUIRE_EQUAL( observations.size( ), 16 );

    // Primary: apparent directions (with stellar aberration). Comparison: astrometric directions, which differ by the annual
    // aberration (at most 20.8 arcsec): required to be between 10 and 21 arcsec RMS (maximum between 15 and 21 arcsec).
    checkAngularPositionPrefitResiduals(
            "Neptune_barycenter",
            observations,
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, true ),
            observations,
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, false ),
            { 0.5, 0.05, 1.0 },
            { 10.0, 21.0, 15.0, 21.0 },
            true );
}

//! Validate apparent true-of-date (IAU 2006/2000A) angular positions on Carlsberg Meridian Circle observations of Uranus
BOOST_AUTO_TEST_CASE( testUranusTrueOfDateIau2006AngularPositionPrefitResiduals )
{
    // Lines 1392, 1451, 1510, 1568, 1627, 1686, 1745, 1804, 1862, 1921, 1980, 2039, 2098, 2156, 2215 and 2274 of
    // https://ssd.jpl.nasa.gov/dat/planets/trnstslplma.txt
    // (Carlsberg Automatic Meridian Circle, La Palma, 1984-1998): apparent places w.r.t. the true equator and equinox of date, referred
    // to the ICRF (catalogue code 888)
    const std::string records = R"(
P7002445833605450    95018881    1643338600       -221430570            1984LPA
P7002446161723630    95018881    1707224060       -225433900            1985LPA
P7002446282381899    95018881    1650592450       -223159670            1985LPA
P7002446618477532    95018881    1713469420       -230646540            1986LPA
P7002446946597299    95018881    1739530660       -233043140            1987LPA
P7002447318591637    95018881    1758210290       -233907370            1988LPA
P7002447697565925    95018881    1815277620       -233952460            1989LPA
P7002448054602802    95018881    1836125750       -233017030            1990LPA
P7002448382720708    95018881    1859367220       -230555400            1991LPA
P7002448458505641    95018881    1848422160       -232108050            1991LPA
P7002448772664355    95018881    1915502950       -224317690            1992LPA
P7002449192521749    95018881    1925478190       -222745000            1993LPA
P7002449585455823    95018881    1940014210       -215902000            1994LPA
P7002449975398930    95018881    1955277100       -212152410            1995LPA
P7002450351381261    95018881    2012206520       -203438540            1996LPA
P7002450747308952    95018881    2029107930       -194015860            1997LPA
)";
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );
    std::vector< Eigen::Vector3d > observations;
    for( const std::string& record : splitRecords( records ) )
    {
        observations.push_back( decodeJplTransitRecord( record, "888", timeScaleConverter ) );
    }
    BOOST_REQUIRE_EQUAL( observations.size( ), 16 );

    // Primary: apparent true-of-date directions. Comparison: apparent J2000 directions, which differ by the precession since J2000
    // (50.3 arcsec per year, 2.2 to 15.6 years) and nutation (at most 19 arcsec): required to be between 300 and 800 arcsec RMS
    // (maximum between 600 and 820 arcsec).
    checkAngularPositionPrefitResiduals(
            "Uranus_barycenter",
            observations,
            getPositionAngleAncillarySettings( true_of_date_iau_2006_2000a_position_angle_reference_frame, TUDAT_NAN, true ),
            observations,
            getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, true ),
            { 0.2, 0.05, 0.5 },
            { 300.0, 800.0, 600.0, 820.0 },
            true );
}

//! Compare the rotation matrices of the celestial frames with SOFA (frames of date) and SPICE (B1950), at several epochs
BOOST_AUTO_TEST_CASE( testAngularPositionReferenceFrameRotations )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );

    // Frame bias from ICRF to FK5 J2000, which precedes the IAU 1976 precession(-nutation)
    double frameBiasMatrix[ 3 ][ 3 ], precessionMatrix[ 3 ][ 3 ], biasPrecessionMatrix[ 3 ][ 3 ];
    iauBp00( basic_astrodynamics::JULIAN_DAY_ON_J2000, 0.0, frameBiasMatrix, precessionMatrix, biasPrecessionMatrix );
    const Eigen::Matrix3d frameBias = convertSofaMatrix( frameBiasMatrix );

    for( const double epoch : { -2.5E9, -9.5E8, 0.0, 7.5E8 } )
    {
        // SOFA rotation matrices from the ICRF at the epoch [TT days since J2000]
        const double daysSinceJ2000 =
                timeScaleConverter->getCurrentTime( basic_astrodynamics::tdb_scale, basic_astrodynamics::tt_scale, epoch ) /
                physical_constants::JULIAN_DAY;
        double sofaMatrix[ 3 ][ 3 ];
        std::map< PositionAngleReferenceFrame, Eigen::Matrix3d > sofaRotations;
        iauPmat76( basic_astrodynamics::JULIAN_DAY_ON_J2000, daysSinceJ2000, sofaMatrix );
        sofaRotations[ mean_of_date_iau_1976_position_angle_reference_frame ] = convertSofaMatrix( sofaMatrix ) * frameBias;
        iauPnm80( basic_astrodynamics::JULIAN_DAY_ON_J2000, daysSinceJ2000, sofaMatrix );
        sofaRotations[ true_of_date_iau_1976_1980_position_angle_reference_frame ] = convertSofaMatrix( sofaMatrix ) * frameBias;
        iauPmat06( basic_astrodynamics::JULIAN_DAY_ON_J2000, daysSinceJ2000, sofaMatrix );
        sofaRotations[ mean_of_date_iau_2006_position_angle_reference_frame ] = convertSofaMatrix( sofaMatrix );
        iauPnm06a( basic_astrodynamics::JULIAN_DAY_ON_J2000, daysSinceJ2000, sofaMatrix );
        sofaRotations[ true_of_date_iau_2006_2000a_position_angle_reference_frame ] = convertSofaMatrix( sofaMatrix );

        for( const auto& sofaRotation : sofaRotations )
        {
            BOOST_CHECK_SMALL(
                    ( getPositionAngleReferenceFrameRotationFromJ2000( epoch, getPositionAngleAncillarySettings( sofaRotation.first ) ) -
                      sofaRotation.second )
                            .cwiseAbs( )
                            .maxCoeff( ),
                    1.0E-15 );

            // A reference epoch in the ancillary settings replaces the observation epoch
            BOOST_CHECK_SMALL( ( getPositionAngleReferenceFrameRotationFromJ2000(
                                         epoch + 1.0E8, getPositionAngleAncillarySettings( sofaRotation.first, epoch ) ) -
                                 sofaRotation.second )
                                       .cwiseAbs( )
                                       .maxCoeff( ),
                               1.0E-15 );
        }

        BOOST_CHECK_SMALL( ( getPositionAngleReferenceFrameRotationFromJ2000(
                                     epoch, getPositionAngleAncillarySettings( b1950_position_angle_reference_frame ) ) -
                             spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "FK4", epoch ).toRotationMatrix( ) )
                                   .cwiseAbs( )
                                   .maxCoeff( ),
                           1.0E-14 );
        BOOST_CHECK_SMALL( ( getPositionAngleReferenceFrameRotationFromJ2000(
                                     epoch, getPositionAngleAncillarySettings( j2000_position_angle_reference_frame ) ) -
                             Eigen::Matrix3d::Identity( ) )
                                   .cwiseAbs( )
                                   .maxCoeff( ),
                           std::numeric_limits< double >::epsilon( ) );
    }
}

//! Compare angular positions with celestial frame settings with an independent computation (SOFA stellar aberration, rotation to
//! the celestial frame), for J2000 and ECLIPJ2000 global frame orientations
BOOST_AUTO_TEST_CASE( testAngularPositionCelestialFrames )
{
    spice_interface::loadStandardSpiceKernels( );
    BodyListSettings j2000BodySettings( "SSB", "J2000" );
    BodyListSettings eclipticBodySettings( "SSB", "ECLIPJ2000" );
    for( const std::string bodyName : { "Earth", "Mars" } )
    {
        j2000BodySettings.addSettings( bodyName );
        j2000BodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
        eclipticBodySettings.addSettings( bodyName );
        eclipticBodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "ECLIPJ2000", bodyName );
    }
    SystemOfBodies j2000Bodies = createSystemOfBodies( j2000BodySettings );
    SystemOfBodies eclipticBodies = createSystemOfBodies( eclipticBodySettings );

    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), std::string( "" ) );
    linkEnds[ transmitter ] = std::make_pair( std::string( "Mars" ), std::string( "" ) );
    const std::shared_ptr< ObservationModel< 2 > > j2000Model = ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( angular_position, linkEnds ), j2000Bodies );
    const std::shared_ptr< ObservationModel< 2 > > eclipticModel = ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( angular_position, linkEnds ), eclipticBodies );

    const double observationTime = 1.0E8;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    const std::vector< PositionAngleReferenceFrame > referenceFrames = { j2000_position_angle_reference_frame,
                                                                         b1950_position_angle_reference_frame,
                                                                         mean_of_date_iau_1976_position_angle_reference_frame,
                                                                         true_of_date_iau_1976_1980_position_angle_reference_frame,
                                                                         mean_of_date_iau_2006_position_angle_reference_frame,
                                                                         true_of_date_iau_2006_2000a_position_angle_reference_frame };
    for( const PositionAngleReferenceFrame referenceFrame : referenceFrames )
    {
        for( int useStellarAberration = 0; useStellarAberration < 2; useStellarAberration++ )
        {
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
                    getPositionAngleAncillarySettings( referenceFrame, TUDAT_NAN, static_cast< bool >( useStellarAberration ) );
            const Eigen::Vector2d j2000Observation = j2000Model->computeObservationsWithLinkEndData(
                    observationTime, receiver, linkEndTimes, linkEndStates, ancillarySettings );
            const Eigen::Vector2d eclipticObservation = eclipticModel->computeObservations( observationTime, receiver, ancillarySettings );

            // Independent computation of the direction from receiver to transmitter, with stellar aberration from SOFA (relativistic,
            // which differs from the classical model of Tudat by less than 2.5E-9 rad)
            Eigen::Vector3d direction = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment< 3 >( 0 ).normalized( );
            if( useStellarAberration )
            {
                const Eigen::Vector3d receiverVelocity = linkEndStates.at( 1 ).segment< 3 >( 3 ) / physical_constants::SPEED_OF_LIGHT;
                const Eigen::Vector3d sunPosition =
                        spice_interface::getBodyCartesianPositionAtEpoch( "Sun", "SSB", "J2000", "NONE", linkEndTimes.at( 1 ) );
                const double sunDistance =
                        ( linkEndStates.at( 1 ).segment< 3 >( 0 ) - sunPosition ).norm( ) / physical_constants::ASTRONOMICAL_UNIT;
                double naturalDirection[ 3 ] = { direction( 0 ), direction( 1 ), direction( 2 ) };
                double velocity[ 3 ] = { receiverVelocity( 0 ), receiverVelocity( 1 ), receiverVelocity( 2 ) };
                double properDirection[ 3 ];
                iauAb( naturalDirection, velocity, sunDistance, std::sqrt( 1.0 - receiverVelocity.squaredNorm( ) ), properDirection );
                direction = Eigen::Vector3d( properDirection[ 0 ], properDirection[ 1 ], properDirection[ 2 ] );
            }
            direction = getPositionAngleReferenceFrameRotationFromJ2000( linkEndTimes.at( 1 ), ancillarySettings ) * direction;
            const Eigen::Vector2d expectedObservation( std::atan2( direction( 1 ), direction( 0 ) ), std::asin( direction( 2 ) ) );

            const double tolerance = useStellarAberration ? 3.0E-9 : 1.0E-14;
            BOOST_CHECK_SMALL( getAngularPositionResidual( expectedObservation, j2000Observation ).norm( ), tolerance );
            BOOST_CHECK_SMALL( getAngularPositionResidual( j2000Observation, eclipticObservation ).norm( ), 1.0E-13 );
        }
    }

    // Without ancillary settings, the angles are w.r.t. the axes of the global frame orientation (ecliptic longitude and latitude
    // for ECLIPJ2000)
    const Eigen::Vector2d eclipticAngles = eclipticModel->computeObservations( observationTime, receiver );
    j2000Model->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimes, linkEndStates );
    const Eigen::Vector3d eclipticDirection = reference_frames::getJ2000toECLIPJ2000TransformationMatrix( ) *
            ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment< 3 >( 0 ).normalized( );
    const Eigen::Vector2d expectedEclipticAngles( std::atan2( eclipticDirection( 1 ), eclipticDirection( 0 ) ),
                                                  std::asin( eclipticDirection( 2 ) ) );
    BOOST_CHECK_SMALL( getAngularPositionResidual( expectedEclipticAngles, eclipticAngles ).norm( ), 1.0E-13 );

    // A custom pole defines no right ascension and declination
    BOOST_CHECK_THROW( j2000Model->computeObservations(
                               observationTime, receiver, getCustomPositionAngleAncillarySettings( Eigen::Vector3d::UnitZ( ) ) ),
                       std::runtime_error );
}

//! Compare angular positions of the Mars system barycentre (INPOP19a, Newtonian light time), for each celestial frame, with and
//! without stellar aberration, with independent values computed with skyfield, astropy and pyerfa (celestialFrameReferenceValues.h)
BOOST_AUTO_TEST_CASE( testAngularPositionIndependentReferenceValues )
{
    spice_interface::clearSpiceKernels( );
    spice_interface::loadStandardSpiceKernels( );
    BodyListSettings bodySettings( "SSB", "J2000" );
    for( const std::string bodyName : { "Earth", "Mars_barycenter" } )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), std::string( "" ) );
    linkEnds[ transmitter ] = std::make_pair( std::string( "Mars_barycenter" ), std::string( "" ) );
    const std::shared_ptr< ObservationModel< 2 > > model = ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( angular_position, linkEnds ), bodies );

    for( const CelestialFrameReferenceValue& referenceValue : getCelestialFrameReferenceValues( ) )
    {
        const double observationTime =
                ( referenceValue.julianDayTdb - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;
        const Eigen::Vector2d computedAngles = model->computeObservations(
                observationTime,
                receiver,
                getPositionAngleAncillarySettings( referenceValue.referenceFrame, TUDAT_NAN, referenceValue.useStellarAberration ) );

        // The IAU 2006/2000A true-of-date frame of skyfield differs from that of SOFA by up to 6E-11 rad
        BOOST_CHECK_SMALL(
                getAngularPositionResidual( Eigen::Vector2d( referenceValue.rightAscension, referenceValue.declination ), computedAngles )
                        .norm( ),
                1.0E-10 );
    }
}

//! Compare the FK4 B1950 to ICRF conversion with pyerfa (fk45z + fk5hz, the same algorithm), and with astropy, which converts FK5
//! to ICRS with the frame bias instead of the FK5 to Hipparcos rotation and spin, and differs by up to 0.06 arcsec for these epochs
BOOST_AUTO_TEST_CASE( testFk4B1950ToIcrfConversion )
{
    for( const Fk4ConversionReferenceValue& referenceValue : getFk4ConversionReferenceValues( ) )
    {
        const Eigen::Vector2d icrfAngles = convertFk4B1950ToIcrfRightAscensionDeclination(
                referenceValue.b1950RightAscension, referenceValue.b1950Declination, referenceValue.observationTime );
        BOOST_CHECK_SMALL( getAngularPositionResidual( Eigen::Vector2d( referenceValue.erfaRightAscension, referenceValue.erfaDeclination ),
                                                       icrfAngles )
                                   .norm( ),
                           1.0E-12 );
        BOOST_CHECK_SMALL( getAngularPositionResidual(
                                   Eigen::Vector2d( referenceValue.astropyRightAscension, referenceValue.astropyDeclination ), icrfAngles )
                                   .norm( ),
                           unit_conversions::convertArcSecondsToRadians( 0.06 ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
