/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
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
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::mathematical_constants;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;

namespace
{

using PositionSeparationModel = std::shared_ptr< ObservationModel< 2 > >;
using PositionSeparationModels = std::map< int, PositionSeparationModel >;

constexpr double radiansToArcseconds = 180.0 / mathematical_constants::PI * 3600.0;

struct ReferenceObservation {
    double receiverTimeTdb;
    int modelCode;
    double positionAngleRadians;
    double separationRadians;
};

struct AngularStatistics {
    void addRadians( const double valueRadians )
    {
        count++;
        sumRadians += valueRadians;
        squaredSumRadians += valueRadians * valueRadians;
        maximumAbsoluteRadians = std::max( maximumAbsoluteRadians, std::abs( valueRadians ) );
    }

    double meanArcseconds( ) const
    {
        return sumRadians / static_cast< double >( count ) * radiansToArcseconds;
    }

    double rmsArcseconds( ) const
    {
        return std::sqrt( squaredSumRadians / static_cast< double >( count ) ) * radiansToArcseconds;
    }

    double maximumArcseconds( ) const
    {
        return maximumAbsoluteRadians * radiansToArcseconds;
    }

    double meanMilliArcseconds( ) const
    {
        return meanArcseconds( ) * 1000.0;
    }

    double rmsMilliArcseconds( ) const
    {
        return rmsArcseconds( ) * 1000.0;
    }

    int count = 0;
    double sumRadians = 0.0;
    double squaredSumRadians = 0.0;
    double maximumAbsoluteRadians = 0.0;
};

struct ResidualBounds {
    // All bounds are in arcseconds, including the milliarcsecond-scale HST case.
    double maximumAbsoluteMeanArcseconds;
    double minimumRmsArcseconds;
    double maximumRmsArcseconds;
};

struct CorrectionBounds {
    // Bounds on the RMS and largest absolute change between two model conventions.
    double minimumRmsArcseconds;
    double maximumRmsArcseconds;
    double minimumMaximumArcseconds;
    double maximumMaximumArcseconds;
};

struct ReferenceDataset {
    // Dataset-specific source, observation convention, and physically expected scales.
    std::string name;
    std::string dataFile;
    Eigen::Index expectedRows = 0;
    Eigen::Index expectedColumns = 0;
    bool allowAdditionalColumns = false;
    bool reportMilliArcseconds = false;
    PositionSeparationModels models;
    std::string primaryConvention;
    std::shared_ptr< ObservationAncillarySimulationSettings > primarySettings;
    std::string comparisonConvention;
    std::shared_ptr< ObservationAncillarySimulationSettings > comparisonSettings;
    ResidualBounds separationBounds;
    ResidualBounds transversePositionAngleBounds;
    std::optional< CorrectionBounds > separationCorrectionBounds;
    std::optional< CorrectionBounds > transversePositionAngleCorrectionBounds;
    bool requireSeparationInvariance = false;
    bool requirePrimaryTransverseRmsBelowComparison = false;
    std::function< ReferenceObservation( const Eigen::MatrixXd&, Eigen::Index ) > decodeRecord;
};

SystemOfBodies createJ2000SpiceBodies( const std::vector< std::string >& bodyNames, const bool useItrfEarthRotation )
{
    BodyListSettings bodySettings( "SSB", "J2000" );
    for( const std::string& bodyName : bodyNames )
    {
        bodySettings.addSettings( bodyName );
        bodySettings.at( bodyName )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "J2000", bodyName );
    }
    if( useItrfEarthRotation )
    {
        bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257223563 );
        bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    }
    return createSystemOfBodies( bodySettings );
}

PositionSeparationModel createPositionSeparationModel( SystemOfBodies& bodies,
                                                       const std::string& receiverStation,
                                                       const std::string& referenceBody,
                                                       const std::string& observedBody )
{
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), receiverStation );
    linkEnds[ transmitter ] = std::make_pair( referenceBody, std::string( "" ) );
    linkEnds[ transmitter2 ] = std::make_pair( observedBody, std::string( "" ) );
    return ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds ), bodies );
}

double convertUtcJulianDateToTdb( const double julianDate )
{
    std::ostringstream utcJulianDate;
    utcJulianDate << "JD " << std::setprecision( 16 ) << julianDate << " UTC";
    return spice_interface::convertDateStringToEphemerisTime( utcJulianDate.str( ) );
}

double convertFractionalUtcDateToTdb( const int year,
                                      const int month,
                                      const double decimalDay,
                                      const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& timeScaleConverter )
{
    const int day = static_cast< int >( std::floor( decimalDay ) );
    const double secondsOfDay = ( decimalDay - static_cast< double >( day ) ) * physical_constants::JULIAN_DAY;
    const int hour = static_cast< int >( secondsOfDay / 3600.0 );
    const int minute = static_cast< int >( ( secondsOfDay - 3600.0 * hour ) / 60.0 );
    const long double second = secondsOfDay - 3600.0 * hour - 60.0 * minute;
    const basic_astrodynamics::DateTime utcDateTime( year, month, day, hour, minute, second );
    return timeScaleConverter->getCurrentTime(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcDateTime.epoch< double >( ) );
}

double wrappedAngleDifferenceRadians( const double firstAngleRadians, const double secondAngleRadians )
{
    const double differenceRadians = firstAngleRadians - secondAngleRadians;
    return std::atan2( std::sin( differenceRadians ), std::cos( differenceRadians ) );
}

void checkResidualBounds( const AngularStatistics& statistics, const ResidualBounds& bounds )
{
    BOOST_CHECK_LT( std::abs( statistics.meanArcseconds( ) ), bounds.maximumAbsoluteMeanArcseconds );
    BOOST_CHECK_GT( statistics.rmsArcseconds( ), bounds.minimumRmsArcseconds );
    BOOST_CHECK_LT( statistics.rmsArcseconds( ), bounds.maximumRmsArcseconds );
}

void checkCorrectionBounds( const AngularStatistics& statistics, const CorrectionBounds& bounds )
{
    BOOST_CHECK_GT( statistics.rmsArcseconds( ), bounds.minimumRmsArcseconds );
    BOOST_CHECK_LT( statistics.rmsArcseconds( ), bounds.maximumRmsArcseconds );
    BOOST_CHECK_GT( statistics.maximumArcseconds( ), bounds.minimumMaximumArcseconds );
    BOOST_CHECK_LT( statistics.maximumArcseconds( ), bounds.maximumMaximumArcseconds );
}

void validateReferenceDataset( const ReferenceDataset& dataset )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    const Eigen::MatrixXd observations = input_output::readMatrixFromFile( validationDataDirectory + dataset.dataFile, " \t" );
    BOOST_REQUIRE_EQUAL( observations.rows( ), dataset.expectedRows );
    if( dataset.allowAdditionalColumns )
    {
        BOOST_REQUIRE( observations.cols( ) >= dataset.expectedColumns );
    }
    else
    {
        BOOST_REQUIRE_EQUAL( observations.cols( ), dataset.expectedColumns );
    }
    BOOST_REQUIRE( observations.rows( ) > 0 );

    AngularStatistics primarySeparationResiduals;
    AngularStatistics primaryTransversePositionAngleResiduals;
    AngularStatistics comparisonSeparationResiduals;
    AngularStatistics comparisonTransversePositionAngleResiduals;
    AngularStatistics separationCorrections;
    AngularStatistics transversePositionAngleCorrections;
    for( Eigen::Index observationIndex = 0; observationIndex < observations.rows( ); observationIndex++ )
    {
        const ReferenceObservation record = dataset.decodeRecord( observations, observationIndex );
        BOOST_REQUIRE( dataset.models.count( record.modelCode ) == 1 );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d primaryPrediction =
                dataset.models.at( record.modelCode )
                        ->computeObservationsWithLinkEndData(
                                record.receiverTimeTdb, receiver, linkEndTimes, linkEndStates, dataset.primarySettings );
        // Reduce every dataset in the same radial/transverse basis: observed
        // minus computed separation, and observed separation times wrapped PA error.
        const double primarySeparationResidualRadians = record.separationRadians - primaryPrediction( 1 );
        const double primaryTransversePositionAngleResidualRadians =
                record.separationRadians * wrappedAngleDifferenceRadians( record.positionAngleRadians, primaryPrediction( 0 ) );
        primarySeparationResiduals.addRadians( primarySeparationResidualRadians );
        primaryTransversePositionAngleResiduals.addRadians( primaryTransversePositionAngleResidualRadians );

        if( dataset.comparisonSettings != nullptr )
        {
            const Eigen::Vector2d comparisonPrediction =
                    dataset.models.at( record.modelCode )
                            ->computeObservationsWithLinkEndData(
                                    record.receiverTimeTdb, receiver, linkEndTimes, linkEndStates, dataset.comparisonSettings );
            // The comparison convention measures the frame or aberration signal
            // without changing the observations or the receiver/transmitter setup.
            comparisonSeparationResiduals.addRadians( record.separationRadians - comparisonPrediction( 1 ) );
            comparisonTransversePositionAngleResiduals.addRadians(
                    record.separationRadians * wrappedAngleDifferenceRadians( record.positionAngleRadians, comparisonPrediction( 0 ) ) );
            separationCorrections.addRadians( primaryPrediction( 1 ) - comparisonPrediction( 1 ) );
            transversePositionAngleCorrections.addRadians(
                    record.separationRadians * wrappedAngleDifferenceRadians( primaryPrediction( 0 ), comparisonPrediction( 0 ) ) );
        }
    }

    // Report both residual components for each convention; retain mas output for the precise HST data.
    if( dataset.reportMilliArcseconds )
    {
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.primaryConvention << " separation mean/RMS [mas]: " << std::setprecision( 15 )
                                         << primarySeparationResiduals.meanMilliArcseconds( ) << " / "
                                         << primarySeparationResiduals.rmsMilliArcseconds( ) );
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.primaryConvention << " transverse PA mean/RMS [mas]: "
                                         << primaryTransversePositionAngleResiduals.meanMilliArcseconds( ) << " / "
                                         << primaryTransversePositionAngleResiduals.rmsMilliArcseconds( ) );
    }
    else
    {
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.primaryConvention << " separation mean/RMS [arcsec]: " << std::setprecision( 15 )
                                         << primarySeparationResiduals.meanArcseconds( ) << " / "
                                         << primarySeparationResiduals.rmsArcseconds( ) );
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.primaryConvention << " transverse PA mean/RMS [arcsec]: "
                                         << primaryTransversePositionAngleResiduals.meanArcseconds( ) << " / "
                                         << primaryTransversePositionAngleResiduals.rmsArcseconds( ) );
    }
    if( dataset.comparisonSettings != nullptr )
    {
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.comparisonConvention
                                         << " separation mean/RMS [arcsec]: " << comparisonSeparationResiduals.meanArcseconds( ) << " / "
                                         << comparisonSeparationResiduals.rmsArcseconds( ) );
        BOOST_TEST_MESSAGE( dataset.name << " " << dataset.comparisonConvention << " transverse PA mean/RMS [arcsec]: "
                                         << comparisonTransversePositionAngleResiduals.meanArcseconds( ) << " / "
                                         << comparisonTransversePositionAngleResiduals.rmsArcseconds( ) );
        BOOST_TEST_MESSAGE( dataset.name << " separation correction RMS/max [arcsec]: " << separationCorrections.rmsArcseconds( ) << " / "
                                         << separationCorrections.maximumArcseconds( ) );
        BOOST_TEST_MESSAGE( dataset.name << " transverse PA correction RMS/max [arcsec]: "
                                         << transversePositionAngleCorrections.rmsArcseconds( ) << " / "
                                         << transversePositionAngleCorrections.maximumArcseconds( ) );
    }

    // Bounds are dataset-level expectations from the documented source precision, not exact outputs of the current model.
    checkResidualBounds( primarySeparationResiduals, dataset.separationBounds );
    checkResidualBounds( primaryTransversePositionAngleResiduals, dataset.transversePositionAngleBounds );
    if( dataset.comparisonSettings != nullptr )
    {
        // A changed celestial reference pole must not alter separation.  B1950
        // observations also need to improve transverse PA RMS over J2000.
        if( dataset.requireSeparationInvariance )
        {
            BOOST_CHECK_SMALL( separationCorrections.maximumAbsoluteRadians, 5.0 * std::numeric_limits< double >::epsilon( ) );
        }
        if( dataset.requirePrimaryTransverseRmsBelowComparison )
        {
            BOOST_CHECK_GT( comparisonTransversePositionAngleResiduals.rmsArcseconds( ),
                            primaryTransversePositionAngleResiduals.rmsArcseconds( ) );
        }
        if( dataset.separationCorrectionBounds )
        {
            checkCorrectionBounds( separationCorrections, *dataset.separationCorrectionBounds );
        }
        if( dataset.transversePositionAngleCorrectionBounds )
        {
            checkCorrectionBounds( transversePositionAngleCorrections, *dataset.transversePositionAngleCorrectionBounds );
        }
    }
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_position_angle_and_separation_observation_model )

//! Check north-through-east angle signs, right-angle separation, and geometries where position angle is undefined.
BOOST_AUTO_TEST_CASE( testPositionAngleConventionAndSingularities )
{
    const Eigen::Vector3d northPole = Eigen::Vector3d::UnitZ( );
    const Eigen::Vector3d xDirection = Eigen::Vector3d::UnitX( );
    const Eigen::Vector3d yDirection = Eigen::Vector3d::UnitY( );
    const Eigen::Vector3d negativeXDirection = -xDirection;
    const Eigen::Vector3d negativeYDirection = -yDirection;

    Eigen::Vector2d observation = calculatePositionAngleAndSeparation( xDirection, yDirection, northPole );
    // A source due east of the reference direction has PA +90 degrees and separation 90 degrees.
    BOOST_CHECK_SMALL( observation( 0 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    observation = calculatePositionAngleAndSeparation( xDirection, negativeYDirection, northPole );
    // Reversing east to west reverses the sign of PA.
    BOOST_CHECK_SMALL( observation( 0 ) + mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );

    // PA is undefined at the pole of its tangent plane and for coincident or antipodal sight lines.
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( northPole, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, xDirection, northPole ), std::runtime_error );
    BOOST_CHECK_THROW( calculatePositionAngleAndSeparation( xDirection, negativeXDirection, northPole ), std::runtime_error );

    // Separation remains well-defined when the first line of sight is at the reference pole.
    observation = calculatePositionAngleAndSeparation( northPole, xDirection, northPole, false );
    BOOST_CHECK_SMALL( observation( 1 ) - mathematical_constants::PI / 2.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Compare P-only, S-only, and combined models against each other and an independent RA/Dec construction.
BOOST_AUTO_TEST_CASE( testPositionAngleAndSeparationObservationModel )
{
    spice_interface::loadStandardSpiceKernels( { paths::getSpiceKernelPath( ) + "/de430_mar097_small.bsp" } );

    // Define bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Phobos" );

    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.addSettings( "Phobos" );
    bodySettings.at( "Phobos" )->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define link ends
    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );
    linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Phobos", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create observation settings for all three models
    std::shared_ptr< ObservationModelSettings > positionAngleSettings =
            std::make_shared< PositionAngleObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > separationSettings =
            std::make_shared< SeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > positionAngleAndSeparationSettings =
            std::make_shared< PositionAngleAndSeparationObservationModelSettings >( linkEnds, lightTimeCorrectionSettings );

    // Create observation models
    std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( positionAngleSettings, bodies );
    std::shared_ptr< ObservationModel< 1 > > separationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel( separationSettings, bodies );
    std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( positionAngleAndSeparationSettings, bodies );

    // Simulate all three observables at one common receiver epoch with identical link ends and light-time corrections.
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    std::vector< double > positionAngleLinkEndTimes;
    std::vector< double > separationLinkEndTimes;
    std::vector< double > positionAngleAndSeparationLinkEndTimes;
    std::vector< Eigen::Vector6d > positionAngleLinkEndStates;
    std::vector< Eigen::Vector6d > separationLinkEndStates;
    std::vector< Eigen::Vector6d > positionAngleAndSeparationLinkEndStates;

    Eigen::VectorXd positionAngleObservation = positionAngleModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, positionAngleLinkEndTimes, positionAngleLinkEndStates );
    Eigen::VectorXd separationObservation = separationModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, separationLinkEndTimes, separationLinkEndStates );
    Eigen::VectorXd positionAngleAndSeparationObservation = positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, positionAngleAndSeparationLinkEndTimes, positionAngleAndSeparationLinkEndStates );

    // The two components of the combined observable must equal the corresponding standalone observables.
    BOOST_CHECK_CLOSE_FRACTION(
            positionAngleObservation( 0 ), positionAngleAndSeparationObservation( 0 ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION(
            separationObservation( 0 ), positionAngleAndSeparationObservation( 1 ), std::numeric_limits< double >::epsilon( ) * 10.0 );

    // Each wrapper must return the same retarded transmitter and receiver states as the combined model.
    for( int i = 0; i < 3; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                positionAngleLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                separationLinkEndStates[ i ], positionAngleAndSeparationLinkEndStates[ i ], std::numeric_limits< double >::epsilon( ) );
    }

    // Independent geometric check: obtain retarded states from the relative-RA/Dec model, then use spherical trigonometry.
    std::shared_ptr< ObservationModelSettings > relativeAngularPositionSettings =
            std::make_shared< ObservationModelSettings >( relative_angular_position, linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModel< 2 > > relativeAngularPositionModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( relativeAngularPositionSettings, bodies );

    std::vector< double > relativeAngularPositionLinkEndTimes;
    std::vector< Eigen::Vector6d > relativeAngularPositionLinkEndStates;
    relativeAngularPositionModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, relativeAngularPositionLinkEndTimes, relativeAngularPositionLinkEndStates );

    // relativeAngularPositionLinkEndStates: [0]=Mars (transmitter), [1]=Phobos (transmitter2), [2]=Earth (receiver)
    const Eigen::Matrix3d globalToJ2000 = reference_frames::getECLIPJ2000toJ2000TransformationMatrix( );
    auto computeRightAscensionAndDeclination = [ globalToJ2000 ]( const Eigen::Vector6d& transmitterState,
                                                                  const Eigen::Vector6d& receiverState ) {
        const Eigen::Vector3d relativePosition = globalToJ2000 * ( transmitterState - receiverState ).segment( 0, 3 );
        double rightAscension = 2.0 *
                std::atan( relativePosition( 1 ) /
                           ( std::sqrt( relativePosition( 0 ) * relativePosition( 0 ) + relativePosition( 1 ) * relativePosition( 1 ) ) +
                             relativePosition( 0 ) ) );
        double declination = mathematical_constants::PI / 2.0 - std::acos( relativePosition( 2 ) / relativePosition.norm( ) );
        return std::make_pair( rightAscension, declination );
    };

    auto [ marsRightAscension, marsDeclination ] =
            computeRightAscensionAndDeclination( relativeAngularPositionLinkEndStates[ 0 ], relativeAngularPositionLinkEndStates[ 2 ] );
    auto [ phobosRightAscension, phobosDeclination ] =
            computeRightAscensionAndDeclination( relativeAngularPositionLinkEndStates[ 1 ], relativeAngularPositionLinkEndStates[ 2 ] );

    // Convert RA/DEC to position angle and separation
    double rightAscensionDifference = phobosRightAscension - marsRightAscension;
    double computedSeparationDistance =
            std::acos( std::sin( marsDeclination ) * std::sin( phobosDeclination ) +
                       std::cos( marsDeclination ) * std::cos( phobosDeclination ) * std::cos( rightAscensionDifference ) );
    double computedPositionAngle =
            std::atan2( std::cos( phobosDeclination ) * std::sin( rightAscensionDifference ),
                        std::sin( phobosDeclination ) * std::cos( marsDeclination ) -
                                std::cos( phobosDeclination ) * std::sin( marsDeclination ) * std::cos( rightAscensionDifference ) );

    // The independently reconstructed P/S values must agree with both dedicated and combined model outputs.
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - separationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedPositionAngle - positionAngleAndSeparationObservation( 0 ), 1.0e-11 );
    BOOST_CHECK_SMALL( computedSeparationDistance - positionAngleAndSeparationObservation( 1 ), 1.0e-11 );

    // P/S time tags are defined at reception; a transmitter-referenced epoch must be rejected.
    BOOST_CHECK_THROW( positionAngleModel->computeObservationsWithLinkEndData(
                               receiverObservationTime, transmitter, positionAngleLinkEndTimes, positionAngleLinkEndStates ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            positionAngleAndSeparationModel->computeObservationsWithLinkEndData(
                    receiverObservationTime, transmitter, positionAngleAndSeparationLinkEndTimes, positionAngleAndSeparationLinkEndStates ),
            std::runtime_error );
}

//! Validate default astrometric J2000 P/S against HST Pluto--Charon measurements without orbit fitting.
BOOST_AUTO_TEST_CASE( testPlutoCharonHstPrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "plu060_pm0001_subset.bsp" } );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", "Pluto", "Charon" }, false );

    ReferenceDataset dataset;
    dataset.name = "pm0001";
    dataset.dataFile = "pm0001.txt";
    dataset.expectedRows = 60;
    dataset.expectedColumns = 3;
    dataset.allowAdditionalColumns = true;
    dataset.reportMilliArcseconds = true;
    dataset.primaryConvention = "default astrometric J2000";
    dataset.primarySettings = nullptr;
    dataset.models[ 0 ] = createPositionSeparationModel( bodies, "", "Pluto", "Charon" );

    // Tholen and Buie (1997), NSDB pm0001: UTC Julian date at HST exposure
    // midpoint, separation [arcsec], and J2000 position angle [deg].  Exposure
    // midpoints agree with MAST to 6 ms; geocentric versus HST receiver changes
    // either observable by less than 0.002 mas for these records.
    dataset.decodeRecord = []( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        return ReferenceObservation{ convertUtcJulianDateToTdb( observations( row, 0 ) ),
                                     0,
                                     unit_conversions::convertDegreesToRadians( observations( row, 2 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 1 ) ) };
    };

    // Source-file fitted-orbit RMS is 2.93 mas radial and 2.70 mas transverse;
    // independent PLU060 pre-fit residuals should have the same order of magnitude.
    dataset.separationBounds = { 0.005, 0.001, 0.010 };
    dataset.transversePositionAngleBounds = { 0.005, 0.001, 0.010 };
    validateReferenceDataset( dataset );
}

//! Validate true-of-date PA on real Mars-satellite P/S data and isolate its effect from the invariant separation.
BOOST_AUTO_TEST_CASE( testMarsSatelliteTrueOfDatePrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "mar099_de441_mm0012_subset.bsp" } );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", "Mars", "Phobos", "Deimos" }, true );

    // Isaac Newton Group coordinates for the Jacobus Kapteyn Telescope.
    const Eigen::Vector3d jacobusKapteynTelescopeGeodeticPosition(
            2364.0,
            unit_conversions::convertDegreesToRadians( 28.0 + 45.0 / 60.0 + 40.1 / 3600.0 ),
            unit_conversions::convertDegreesToRadians( -( 17.0 + 52.0 / 60.0 + 41.2 / 3600.0 ) ) );
    createGroundStation( bodies.at( "Earth" ), "JKT", jacobusKapteynTelescopeGeodeticPosition, coordinate_conversions::geodetic_position );

    ReferenceDataset dataset;
    dataset.name = "mm0012";
    dataset.dataFile = "mm0012.txt";
    dataset.expectedRows = 166;
    dataset.expectedColumns = 18;
    dataset.primaryConvention = "true-of-date IAU 1976/1980";
    dataset.primarySettings = getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame );
    dataset.comparisonConvention = "J2000";
    dataset.comparisonSettings = getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    const std::map< int, std::pair< std::string, std::string > > observedBodyPairs = { { 1, { "Mars", "Phobos" } },
                                                                                       { 2, { "Mars", "Deimos" } },
                                                                                       { 3, { "Phobos", "Deimos" } } };
    for( const auto& observedBodyPair : observedBodyPairs )
    {
        dataset.models[ observedBodyPair.first ] =
                createPositionSeparationModel( bodies, "JKT", observedBodyPair.second.first, observedBodyPair.second.second );
    }

    // Jones, Sinclair, and Williams (1989), NSDB mm0012: UTC calendar date;
    // columns 14/15 give astrometric separation [arcsec] and PA [deg] referred
    // to the true equator and equinox of date.  Reject unexpected NSDB format flags.
    const auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );
    dataset.decodeRecord = [ timeScaleConverter ]( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( row, 0 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( row, 9 ) ), 0 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( row, 10 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( row, 11 ) ), 1 );
        BOOST_REQUIRE_EQUAL( static_cast< int >( observations( row, 12 ) ), 1 );
        return ReferenceObservation{ convertFractionalUtcDateToTdb( static_cast< int >( observations( row, 2 ) ),
                                                                    static_cast< int >( observations( row, 3 ) ),
                                                                    observations( row, 4 ),
                                                                    timeScaleConverter ),
                                     static_cast< int >( observations( row, 8 ) ),
                                     unit_conversions::convertDegreesToRadians( observations( row, 14 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 13 ) ) };
    };

    // Published photographic scatter is about 0.2 arcsec.  The frame signal is
    // independently reproduced with ERFA pnm80 and SPICE at about 6 mas RMS;
    // changing the pole must leave the geometric separation unchanged.
    dataset.separationBounds = { 0.05, 0.1, 0.5 };
    dataset.transversePositionAngleBounds = { 0.1, 0.1, 0.5 };
    dataset.transversePositionAngleCorrectionBounds = CorrectionBounds{ 0.001, 0.01, 0.005, 0.05 };
    dataset.requireSeparationInvariance = true;
    validateReferenceDataset( dataset );
}

//! Validate B1950 north on historical Nereid astrometry and resolve its effect on PA separately from residual noise.
BOOST_AUTO_TEST_CASE( testNereidB1950ReferencePolePrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "nep105_nm1007_subset.bsp" } );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", "Neptune", "Nereid" }, true );

    // Historical MPC parallax constants specify these sites in the Earth-fixed
    // frame without requiring a geodetic-datum interpretation.
    const double terrestrialEquatorialRadius = 6378137.0;
    const std::map< int, Eigen::Vector3d > stationParallaxConstants = {
        { 568, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 204.5278 ), 0.94171, 0.33725 ) },
        { 711, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 255.9785 ), 0.86114, 0.50731 ) },
        { 809, Eigen::Vector3d( unit_conversions::convertDegreesToRadians( 289.26626 ), 0.873440, -0.486052 ) }
    };

    ReferenceDataset dataset;
    dataset.name = "nm1007";
    dataset.dataFile = "nm1007.txt";
    dataset.expectedRows = 17;
    dataset.expectedColumns = 10;
    dataset.primaryConvention = "B1950";
    dataset.primarySettings = getPositionAngleAncillarySettings( b1950_position_angle_reference_frame );
    dataset.comparisonConvention = "J2000";
    dataset.comparisonSettings = getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    for( const auto& stationParallaxConstant : stationParallaxConstants )
    {
        const double longitude = stationParallaxConstant.second( 0 );
        const double distanceFromSpinAxis = terrestrialEquatorialRadius * stationParallaxConstant.second( 1 );
        const Eigen::Vector3d bodyFixedStationPosition( distanceFromSpinAxis * std::cos( longitude ),
                                                        distanceFromSpinAxis * std::sin( longitude ),
                                                        terrestrialEquatorialRadius * stationParallaxConstant.second( 2 ) );
        const std::string stationName = std::to_string( stationParallaxConstant.first );
        createGroundStation( bodies.at( "Earth" ), stationName, bodyFixedStationPosition, coordinate_conversions::cartesian_position );
        dataset.models[ stationParallaxConstant.first ] = createPositionSeparationModel( bodies, stationName, "Neptune", "Nereid" );
    }

    // Veillet (1982) and Veillet & Bois (1988), NSDB nm1007: UTC topocentric
    // astrometric B1950 offsets x/y [arcsec].  Convert those offsets to P/S;
    // undocumented atmospheric reduction is deliberately not modeled here.
    const auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );
    dataset.decodeRecord = [ timeScaleConverter ]( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        const double observedX = unit_conversions::convertArcSecondsToRadians( observations( row, 4 ) );
        const double observedY = unit_conversions::convertArcSecondsToRadians( observations( row, 5 ) );
        return ReferenceObservation{ convertFractionalUtcDateToTdb( static_cast< int >( observations( row, 1 ) ),
                                                                    static_cast< int >( observations( row, 2 ) ),
                                                                    observations( row, 3 ),
                                                                    timeScaleConverter ),
                                     static_cast< int >( observations( row, 6 ) ),
                                     std::atan2( observedX, observedY ),
                                     std::hypot( observedX, observedY ) };
    };

    // Yuan et al. (2021) report 0.17--0.22 arcsec coordinate RMS for these
    // historical series; the B1950 pole correction is larger than that noise.
    dataset.separationBounds = { 0.1, 0.05, 0.5 };
    dataset.transversePositionAngleBounds = { 0.2, 0.1, 0.5 };
    dataset.transversePositionAngleCorrectionBounds = CorrectionBounds{ 0.2, 0.5, 0.5, 1.0 };
    dataset.requireSeparationInvariance = true;
    dataset.requirePrimaryTransverseRmsBelowComparison = true;
    validateReferenceDataset( dataset );
}

//! Validate native B1950 P/S measurements of Uranian satellites and pole-invariance of separation.
BOOST_AUTO_TEST_CASE( testUranianSatelliteNativeB1950PrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "ura184_de441_um0027_subset.bsp" } );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", "Uranus", "Ariel", "Umbriel", "Titania", "Oberon" }, true );

    // MPC observatory code 371 (Tokyo-Okayama); Earth-fixed parallax constants
    // avoid imposing a modern geodetic datum on this historical observing site.
    const double terrestrialEquatorialRadius = 6378137.0;
    const double observatoryLongitude = unit_conversions::convertDegreesToRadians( 133.5965 );
    const double observatoryDistanceFromSpinAxis = terrestrialEquatorialRadius * 0.82433;
    const Eigen::Vector3d bodyFixedObservatoryPosition( observatoryDistanceFromSpinAxis * std::cos( observatoryLongitude ),
                                                        observatoryDistanceFromSpinAxis * std::sin( observatoryLongitude ),
                                                        terrestrialEquatorialRadius * 0.56431 );
    createGroundStation( bodies.at( "Earth" ), "Tokyo-Okayama", bodyFixedObservatoryPosition, coordinate_conversions::cartesian_position );

    ReferenceDataset dataset;
    dataset.name = "um0027";
    dataset.dataFile = "um0027.txt";
    dataset.expectedRows = 156;
    dataset.expectedColumns = 10;
    dataset.primaryConvention = "B1950";
    dataset.primarySettings = getPositionAngleAncillarySettings( b1950_position_angle_reference_frame );
    dataset.comparisonConvention = "J2000";
    dataset.comparisonSettings = getPositionAngleAncillarySettings( j2000_position_angle_reference_frame );
    const std::map< int, std::string > observedSatellites = { { 1, "Ariel" }, { 2, "Umbriel" }, { 3, "Titania" }, { 4, "Oberon" } };
    for( const auto& observedSatellite : observedSatellites )
    {
        dataset.models[ observedSatellite.first ] =
                createPositionSeparationModel( bodies, "Tokyo-Okayama", "Uranus", observedSatellite.second );
    }

    // Tomita and Soma (1979), NSDB um0027: native topocentric astrometric
    // B1950 P/S.  The file specifies ET calendar time, not UTC.
    dataset.decodeRecord = []( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        const double receiverTimeTdb = basic_astrodynamics::DateTime( static_cast< int >( observations( row, 0 ) ),
                                                                      static_cast< int >( observations( row, 1 ) ),
                                                                      static_cast< int >( observations( row, 2 ) ),
                                                                      static_cast< int >( observations( row, 3 ) ),
                                                                      static_cast< int >( observations( row, 4 ) ),
                                                                      static_cast< long double >( observations( row, 5 ) ) )
                                               .epoch< double >( );
        // Historical ET--TDB differences are negligible against the photographic precision.
        return ReferenceObservation{ receiverTimeTdb,
                                     static_cast< int >( observations( row, 7 ) ),
                                     unit_conversions::convertDegreesToRadians( observations( row, 9 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 8 ) ) };
    };

    // Jacobson (2014), Table 3 gives post-fit RMS 0.443 arcsec radial and
    // 0.363 arcsec transverse for retained records; this unfiltered set has
    // outliers.  B1950 still improves transverse PA over J2000.
    dataset.separationBounds = { 0.2, 0.2, 2.0 };
    dataset.transversePositionAngleBounds = { 0.2, 0.2, 2.0 };
    dataset.transversePositionAngleCorrectionBounds = CorrectionBounds{ 0.02, 0.1, 0.05, 0.2 };
    dataset.requireSeparationInvariance = true;
    dataset.requirePrimaryTransverseRmsBelowComparison = true;
    validateReferenceDataset( dataset );
}

//! Validate the apparent-direction option on Saturnian P/S data and resolve differential stellar-aberration effects.
BOOST_AUTO_TEST_CASE( testSaturnSatelliteApparentDirectionPrefitResiduals )
{
    const std::string validationDataDirectory = paths::getTudatTestDataPath( ) + "position_angle_and_separation/";
    spice_interface::loadStandardSpiceKernels( { validationDataDirectory + "sat441_de441_qiao1999_subset.bsp" } );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", "Enceladus", "Tethys", "Dione", "Rhea", "Titan" }, true );

    // Qiao et al. (1999) give these Sheshan coordinates, with longitude positive east.
    const Eigen::Vector3d sheshanGeodeticPosition( 97.0,
                                                   unit_conversions::convertDegreesToRadians( 31.0 + 5.0 / 60.0 + 46.1 / 3600.0 ),
                                                   unit_conversions::convertDegreesToRadians( 121.0 + 11.0 / 60.0 + 3.3 / 3600.0 ) );
    createGroundStation( bodies.at( "Earth" ), "Sheshan", sheshanGeodeticPosition, coordinate_conversions::geodetic_position );

    ReferenceDataset dataset;
    dataset.name = "Qiao 1999 Table 4";
    dataset.dataFile = "qiao1999_table4.txt";
    dataset.expectedRows = 41;
    dataset.expectedColumns = 4;
    dataset.primaryConvention = "aberrated true-of-date";
    dataset.primarySettings = getPositionAngleAncillarySettings(
            true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, aberrated_position_angle_direction );
    dataset.comparisonConvention = "astrometric true-of-date";
    dataset.comparisonSettings = getPositionAngleAncillarySettings(
            true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, astrometric_position_angle_direction );
    const std::map< int, std::string > observedSatellites = { { 602, "Enceladus" }, { 603, "Tethys" }, { 604, "Dione" }, { 605, "Rhea" } };
    for( const auto& observedSatellite : observedSatellites )
    {
        dataset.models[ observedSatellite.first ] = createPositionSeparationModel( bodies, "Sheshan", "Titan", observedSatellite.second );
    }

    // Qiao et al. (1999), Table 4: apparent topocentric P/S with differential
    // refraction removed, stellar aberration and parallax retained.  Columns:
    // UTC Julian date, measured-satellite NAIF ID, PA [deg], separation [arcsec].
    dataset.decodeRecord = []( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        return ReferenceObservation{ convertUtcJulianDateToTdb( observations( row, 0 ) ),
                                     static_cast< int >( observations( row, 1 ) ),
                                     unit_conversions::convertDegreesToRadians( observations( row, 2 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 3 ) ) };
    };

    // Reported coordinate RMS is roughly 0.08--0.20 arcsec, while independent
    // SpiceyPy calculations show a differential-aberration signal of a few mas.
    dataset.separationBounds = { 0.2, 0.1, 0.5 };
    dataset.transversePositionAngleBounds = { 0.1, 0.05, 0.2 };
    dataset.separationCorrectionBounds = CorrectionBounds{ 0.002, 0.01, 0.005, 0.02 };
    dataset.transversePositionAngleCorrectionBounds = CorrectionBounds{ 0.0001, 0.001, 0.0001, 0.001 };
    validateReferenceDataset( dataset );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
