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

#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

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
Eigen::Vector2d getAngularPositionResidual( const double rightAscension, const double declination, const Eigen::Vector2d& computed )
{
    const double rightAscensionDifference = rightAscension - computed( 0 );
    return Eigen::Vector2d(
            std::cos( declination ) * std::atan2( std::sin( rightAscensionDifference ), std::cos( rightAscensionDifference ) ),
            declination - computed( 1 ) );
}

// Check pre-fit residuals w.r.t. INPOP19a of geocentric right ascension/declination observations of a planet (TDB [s since J2000],
// right ascension and declination [deg] in the frame of the observations). For both residual components, the absolute mean must be
// below residualBounds( 0 ), and the RMS between residualBounds( 1 ) and residualBounds( 2 ) [arcsec]. For the angular distance between
// the predictions with the primary and comparison settings, the RMS must be between correctionBounds( 0 ) and correctionBounds( 1 ),
// and the maximum between correctionBounds( 2 ) and correctionBounds( 3 ) [arcsec]. The primary settings must give the smaller
// residual RMS. Right ascensions of date on the FK4 equinox are corrected to FK5.
void checkAngularPositionPrefitResiduals( const std::string& observedBody,
                                          const std::vector< Eigen::Vector3d >& observations,
                                          const std::shared_ptr< ObservationAncillarySimulationSettings >& primarySettings,
                                          const std::shared_ptr< ObservationAncillarySimulationSettings >& comparisonSettings,
                                          const Eigen::Vector3d& residualBounds,
                                          const Eigen::Vector4d& correctionBounds,
                                          const bool correctFk4Equinox = false )
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
    const int numberOfObservations = observations.size( );
    const double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;
    Eigen::Matrix2Xd primaryResiduals( 2, numberOfObservations ), comparisonResiduals( 2, numberOfObservations );
    Eigen::VectorXd corrections( numberOfObservations );
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    for( int i = 0; i < numberOfObservations; i++ )
    {
        const Eigen::Vector3d& observation = observations.at( i );

        // FK4 to FK5 equinox correction of 0.035 s + 0.085 s per Julian century since B1950.0 (Fricke, 1982)
        const double fk4EquinoxCorrection = ( 0.035 + 0.085 * ( observation( 0 ) / physical_constants::JULIAN_DAY + 18262.5 ) / 36525.0 ) *
                mathematical_constants::PI / 43200.0;
        const double rightAscension =
                unit_conversions::convertDegreesToRadians( observation( 1 ) ) + ( correctFk4Equinox ? fk4EquinoxCorrection : 0.0 );
        const double declination = unit_conversions::convertDegreesToRadians( observation( 2 ) );
        const Eigen::Vector2d primaryPrediction =
                model->computeObservationsWithLinkEndData( observation( 0 ), receiver, linkEndTimes, linkEndStates, primarySettings );
        const Eigen::Vector2d comparisonPrediction =
                model->computeObservationsWithLinkEndData( observation( 0 ), receiver, linkEndTimes, linkEndStates, comparisonSettings );
        primaryResiduals.col( i ) = radiansToArcseconds * getAngularPositionResidual( rightAscension, declination, primaryPrediction );
        comparisonResiduals.col( i ) =
                radiansToArcseconds * getAngularPositionResidual( rightAscension, declination, comparisonPrediction );
        corrections( i ) = ( primaryResiduals.col( i ) - comparisonResiduals.col( i ) ).norm( );
    }
    const double correctionRms = corrections.norm( ) / std::sqrt( numberOfObservations );
    BOOST_TEST_MESSAGE( observedBody << " RA cos(Dec), Dec RMS [arcsec]: "
                                     << primaryResiduals.row( 0 ).norm( ) / std::sqrt( numberOfObservations ) << ", "
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
        BOOST_CHECK_GT( comparisonResiduals.row( i ).norm( ) / std::sqrt( numberOfObservations ), primaryRms );
    }
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

// The angular position tests below use 16 observations, evenly spread in time, of series in the IMCCE planetary observation
// database (https://ftp.imcce.fr/pub/databases/PODB/, URSS files) and the JPL planetary optical observation files
// (https://ssd.jpl.nasa.gov/dat/planets/, transit files). All observations are geocentric, of the planet centre; the planet
// system barycentre of INPOP19a is used (offset < 0.1 arcsec). Epochs are converted from the published time scale (ET, UTC, or
// UT1 with the Espenak and Meeus Delta T model) to TDB. The stated scatter is the RMS of the residuals of these observations
// w.r.t. INPOP19a; the residual bounds allow a margin of a factor of about 2, and a mean offset up to about the scatter.

//! Validate astrometric B1950 angular positions on Nikolaev photographic observations of Saturn, against J2000.
BOOST_AUTO_TEST_CASE( testSaturnB1950AngularPositionPrefitResiduals )
{
    // IMCCE file URSS/SATU/PHOTO/nikol.pht.ftp (Nikolaev, MPC 089): astrometric, B1950 (FK4) photographic positions, in ET
    // Scatter of these observations about 0.34 arcsec (RA) and 0.53 arcsec (Dec). The B1950 and J2000 frames differ by a rotation of
    // 2515 arcsec (50 years of precession) about the ecliptic pole, which bounds the correction.
    const std::vector< Eigen::Vector3d > observations = {
        { -1211901454.399, 297.034412500, -21.304166667 }, { -1143826525.311, 319.374250000, -17.077363889 },
        { -1110305819.520, 330.951766667, -13.693983333 }, { -1052310466.656, 358.868366667, -3.147594444 },
        { -1050326762.400, 357.320479167, -3.858602778 },  { -1041797257.632, 354.949579167, -4.627280556 },
        { -983975788.079, 20.367250000, 5.588330556 },     { -923826808.816, 50.265141667, 15.910541667 },
        { -859541788.622, 79.534483333, 21.454288889 },    { -826109869.027, 94.651400000, 22.283733333 },
        { -784831136.204, 103.160050000, 22.585344444 },   { -751053678.264, 117.869779167, 21.350788889 },
        { -693056615.968, 151.746475000, 13.166172222 },   { -651518896.398, 159.149183333, 10.891288889 },
        { -592409756.937, 186.720041667, -0.017141667 },   { -555785186.984, 195.427008333, -3.669505556 }
    };
    checkAngularPositionPrefitResiduals( "Saturn_barycenter",
                                         observations,
                                         getPositionAngleAncillarySettings( b1950_position_angle_reference_frame ),
                                         nullptr,
                                         { 0.5, 0.1, 1.0 },
                                         { 1000.0, 2520.0, 2000.0, 2520.0 } );
}

//! Validate the astrometric B1950 frame on Pulkovo observations of Uranus over 50 years, against aberrated B1950 directions.
BOOST_AUTO_TEST_CASE( testUranusB1950AngularPositionPrefitResiduals )
{
    // IMCCE file URSS/URAN/PHOTO/pulk.pht.ftp (Pulkovo, MPC 084): astrometric, B1950 (FK4) photographic positions, in ET
    // Scatter of these observations about 0.39 arcsec (RA) and 0.77 arcsec (Dec), with a declination offset of about -0.7 arcsec. The
    // correction is the annual stellar aberration, of at most 20.8 arcsec.
    const std::vector< Eigen::Vector3d > observations = {
        { -2533391314.838, 331.459500000, -12.546722222 }, { -2086876662.942, 23.009208333, 8.976583333 },
        { -2059999513.977, 28.936291667, 11.246972222 },   { -2015349582.953, 30.429250000, 11.883972222 },
        { -1900415146.293, 48.961791667, 17.772083333 },   { -1859133693.404, 49.931041667, 18.101388889 },
        { -1605377899.187, 86.330791667, 23.614500000 },   { -1543333320.487, 96.535875000, 23.638722222 },
        { -1478967522.904, 106.149875000, 23.057388889 },  { -1412011531.153, 115.476750000, 21.919111111 },
        { -1321881740.758, 131.602625000, 18.667472222 },  { -1255019212.070, 139.844166667, 16.392333333 },
        { -1163426760.901, 155.460791667, 11.070666667 },  { -1102157707.312, 165.546666667, 7.052777778 },
        { -1037960654.028, 174.276750000, 3.348916667 },   { -972134061.179, 182.256916667, -0.132055556 }
    };
    checkAngularPositionPrefitResiduals( "Uranus_barycenter",
                                         observations,
                                         getPositionAngleAncillarySettings( b1950_position_angle_reference_frame, TUDAT_NAN, false ),
                                         getPositionAngleAncillarySettings( b1950_position_angle_reference_frame, TUDAT_NAN, true ),
                                         { 1.0, 0.1, 1.5 },
                                         { 10.0, 21.0, 15.0, 21.0 } );
}

//! Validate apparent true-of-date (IAU 1976/1980) angular positions on USNO transit observations of Neptune, against astrometric.
BOOST_AUTO_TEST_CASE( testNeptuneTrueOfDateIau1976AngularPositionPrefitResiduals )
{
    // JPL file trnstswash6a.txt (USNO 6-inch transit circle, MPC 786): apparent places w.r.t. the true equator and equinox of date,
    // with right ascensions on the FK4 equinox, in UT1 (UTC from 1972)
    // Scatter of these observations about 0.38 arcsec (RA) and 0.36 arcsec (Dec). The correction is the annual stellar aberration.
    const std::vector< Eigen::Vector3d > observations = {
        { -2333940291.851, 146.589616667, 13.844527778 },  { -2235884218.011, 152.180000000, 12.084552778 },
        { -2045372404.121, 164.765833333, 7.560694444 },   { -1947144055.575, 170.134375000, 5.493802778 },
        { -1820394720.984, 178.302204167, 2.203830556 },   { -1695454723.757, 186.923000000, -1.334416667 },
        { -1571893243.139, 196.003454167, -5.031447222 },  { -1470477042.826, 200.456958333, -6.761227778 },
        { -1343813840.306, 208.760154167, -9.871511111 },  { -1251788877.660, 215.911562500, -12.335652778 },
        { -1188026453.488, 220.071579167, -13.671427778 }, { -1092555115.130, 226.440833333, -15.554447222 },
        { -996222168.497, 232.671579167, -17.196297222 },  { -908591214.573, 241.328504167, -19.109608333 },
        { -776844388.407, 249.380379167, -20.466941667 },  { -738156584.209, 249.907091667, -20.608386111 }
    };
    checkAngularPositionPrefitResiduals(
            "Neptune_barycenter",
            observations,
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, true ),
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, false ),
            { 0.5, 0.1, 1.0 },
            { 10.0, 21.0, 15.0, 21.0 },
            true );
}

//! Validate apparent true-of-date (IAU 2006/2000A) angular positions on Carlsberg Meridian Circle observations of Uranus.
BOOST_AUTO_TEST_CASE( testUranusTrueOfDateIau2006AngularPositionPrefitResiduals )
{
    // JPL file trnstslplma.txt (Carlsberg Automatic Meridian Circle, La Palma, MPC 950): apparent places w.r.t. the true equator
    // and equinox of date, not on the FK4 equinox, in UTC
    // Scatter of these observations about 0.20 arcsec (RA) and 0.18 arcsec (Dec). The correction is the precession since J2000 (50.3 arcsec
    // per year, 2.2 to 15.6 years) and nutation (at most 19 arcsec).
    const std::vector< Eigen::Vector3d > observations = {
        { -493464434.935, 250.891083333, -22.241825000 }, { -465115024.182, 256.843358333, -22.909416667 },
        { -454690148.743, 252.746854167, -22.533241667 }, { -425651486.051, 258.445591667, -23.112927778 },
        { -397301938.181, 264.971108333, -23.511983333 }, { -365161626.378, 269.587620833, -23.652047222 },
        { -332418247.896, 273.865675000, -23.664572222 }, { -301570260.723, 279.052395833, -23.504730556 },
        { -273220872.643, 284.903008333, -23.098722222 }, { -266673054.434, 282.175900000, -23.352236111 },
        { -239529741.543, 288.959562500, -22.721580556 }, { -203254060.703, 291.449245833, -22.462500000 },
        { -169304555.710, 295.005920833, -21.983888889 }, { -135613471.266, 298.865458333, -21.364558333 },
        { -103128596.867, 303.086050000, -20.577372222 }, { -68920443.365, 307.294970833, -19.671072222 }
    };
    checkAngularPositionPrefitResiduals(
            "Uranus_barycenter",
            observations,
            getPositionAngleAncillarySettings( true_of_date_iau_2006_2000a_position_angle_reference_frame, TUDAT_NAN, true ),
            getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, true ),
            { 0.2, 0.05, 0.5 },
            { 300.0, 800.0, 600.0, 820.0 } );
}

//! Compare the celestial frames of date with ERFA (pmat76, pnm80, pmat06, pnm06a), and the B1950 frame with SPICE FK4, at one epoch.
BOOST_AUTO_TEST_CASE( testAngularPositionReferenceFrameRotations )
{
    // Rows of the rotation matrices from ICRF/J2000 at TDB -9.5E8 s since J2000 (TT from ERFA dtdb at the geocentre)
    const double epoch = -9.5E8;
    const std::map< PositionAngleReferenceFrame, std::vector< Eigen::Vector3d > > erfaRotations = {
        { mean_of_date_iau_1976_position_angle_reference_frame,
          { { 9.9997306732068258e-01, 6.7310278920048282e-03, 2.9253883131556209e-03 },
            { -6.7310278934951725e-03, 9.9997734632669155e-01, -9.8450582992428236e-06 },
            { -2.9253883097264846e-03, -9.8460771897104868e-06, 9.9999572099399114e-01 } } },
        { true_of_date_iau_1976_1980_position_angle_reference_frame,
          { { 9.9997315032192557e-01, 6.7206444601792103e-03, 2.9208857019559688e-03 },
            { -6.7205238961366752e-03, 9.9997741571917986e-01, -5.1089648512088931e-05 },
            { -2.9211630912162772e-03, 3.1458394613597844e-05, 9.9999573289917798e-01 } } },
        { mean_of_date_iau_2006_position_angle_reference_frame,
          { { 9.9997307076024111e-01, 6.7305567691246342e-03, 2.9252965510215131e-03 },
            { -6.7305567563067503e-03, 9.9997734949785166e-01, -9.8489515243382897e-06 },
            { -2.9252965805130508e-03, -9.8401881660215018e-06, 9.9999572126238956e-01 } } },
        { true_of_date_iau_2006_2000a_position_angle_reference_frame,
          { { 9.9997315385089758e-01, 6.7201614772488287e-03, 2.9207887990840100e-03 },
            { -6.7200408261038269e-03, 9.9997741896407877e-01, -5.1119805487442882e-05 },
            { -2.9210663779947808e-03, 3.1490613143170254e-05, 9.9999573318067536e-01 } } }
    };
    // The IAU 1976/1980 and IAU 2006/2000A matrices differ by up to 4.8E-7, so a tolerance of 1E-12 identifies each model
    for( const auto& erfaRotation : erfaRotations )
    {
        const Eigen::Matrix3d rotation =
                getPositionAngleReferenceFrameRotationFromJ2000( epoch, getPositionAngleAncillarySettings( erfaRotation.first ) );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( ( rotation.row( i ).transpose( ) - erfaRotation.second.at( i ) ).cwiseAbs( ).maxCoeff( ), 1.0E-12 );
        }
    }
    const Eigen::Matrix3d b1950Rotation = getPositionAngleReferenceFrameRotationFromJ2000(
            epoch, getPositionAngleAncillarySettings( b1950_position_angle_reference_frame ) );
    BOOST_CHECK_SMALL(
            ( b1950Rotation - spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "FK4", epoch ).toRotationMatrix( ) )
                    .cwiseAbs( )
                    .maxCoeff( ),
            1.0E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
