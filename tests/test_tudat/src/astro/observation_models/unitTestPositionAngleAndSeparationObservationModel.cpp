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

// The tests below compare P/S observations from the literature with Tudat predictions computed from published ephemerides,
// without fitting any orbit (pre-fit residuals). The ResidualBounds of each data set are derived from the precision reported by the
// publication of the observations (or of a later orbit fit to the same data) as follows:
//  - The published values are post-fit RMS residuals of an orbit fitted to these (and other) data. The pre-fit residuals w.r.t. an
//    independent ephemeris are expected to be of the same order, so the RMS bounds bracket the published RMS by a factor of a few
//    (lower bound: residuals much smaller than the measurement precision would indicate a defective test setup; upper bound: margin
//    for ephemeris differences and outliers that were removed in the published fit).
//  - The mean residual is bounded by about the published RMS, as a mean offset larger than the scatter indicates a convention error
//    (e.g. in reference pole, time scale, or reference body).
// The CorrectionBounds are bounds on the difference between two conventions (reference pole, or aberrated vs. astrometric
// directions), and are derived from an independent computation of that difference (stated per data set), with a margin of a
// factor of a few. They verify that the convention options have the expected effect, independently of the observation noise.
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
            std::make_shared< ObservationModelSettings >( position_angle_and_separation, linkEnds ), bodies );
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

// Historical geocentric right ascension/declination observations of the giant planets, in celestial frames other than J2000,
// compared with angular position observables computed from INPOP19a (pre-fit residuals), using the same ResidualBounds
// (for both components) and CorrectionBounds (on the angular distance between the predictions of two conventions) as above.
struct AngularPositionDataset {
    std::string name;
    std::string observedBody;
    // TDB [s since J2000], right ascension [deg], declination [deg], as published, in the frame of the observations
    std::vector< Eigen::Vector3d > observations;
    // Right ascensions of date referred to the FK4 catalogue, corrected to the FK5 equinox
    bool correctFk4Equinox = false;
    std::string primaryConvention;
    std::shared_ptr< ObservationAncillarySimulationSettings > primarySettings;
    std::string comparisonConvention;
    std::shared_ptr< ObservationAncillarySimulationSettings > comparisonSettings;
    ResidualBounds residualBounds;
    CorrectionBounds correctionBounds;
    bool requirePrimaryRmsBelowComparison = true;
};

void validateAngularPositionDataset( const AngularPositionDataset& dataset )
{
    spice_interface::clearSpiceKernels( );
    spice_interface::loadStandardSpiceKernels( );
    SystemOfBodies bodies = createJ2000SpiceBodies( { "Earth", dataset.observedBody }, false );

    LinkDefinition linkEnds;
    linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), std::string( "" ) );
    linkEnds[ transmitter ] = std::make_pair( dataset.observedBody, std::string( "" ) );
    const std::shared_ptr< ObservationModel< 2 > > model = ObservationModelCreator< 2, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( angular_position, linkEnds ), bodies );

    AngularStatistics primaryRightAscensionResiduals, primaryDeclinationResiduals;
    AngularStatistics comparisonRightAscensionResiduals, comparisonDeclinationResiduals;
    AngularStatistics corrections;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    for( const Eigen::Vector3d& observation : dataset.observations )
    {
        // FK4 to FK5 equinox correction of 0.035 s + 0.085 s per Julian century since B1950.0 (Fricke, 1982)
        const double centuriesSinceB1950 = ( observation( 0 ) / physical_constants::JULIAN_DAY + 18262.5 ) / 36525.0;
        const double rightAscension = unit_conversions::convertDegreesToRadians( observation( 1 ) ) +
                ( dataset.correctFk4Equinox ? ( 0.035 + 0.085 * centuriesSinceB1950 ) * PI / 43200.0 : 0.0 );
        const double declination = unit_conversions::convertDegreesToRadians( observation( 2 ) );

        const Eigen::Vector2d primaryPrediction = model->computeObservationsWithLinkEndData(
                observation( 0 ), receiver, linkEndTimes, linkEndStates, dataset.primarySettings );
        const Eigen::Vector2d comparisonPrediction = model->computeObservationsWithLinkEndData(
                observation( 0 ), receiver, linkEndTimes, linkEndStates, dataset.comparisonSettings );
        primaryRightAscensionResiduals.addRadians( std::cos( declination ) *
                                                   wrappedAngleDifferenceRadians( rightAscension, primaryPrediction( 0 ) ) );
        primaryDeclinationResiduals.addRadians( declination - primaryPrediction( 1 ) );
        comparisonRightAscensionResiduals.addRadians( std::cos( declination ) *
                                                      wrappedAngleDifferenceRadians( rightAscension, comparisonPrediction( 0 ) ) );
        comparisonDeclinationResiduals.addRadians( declination - comparisonPrediction( 1 ) );
        corrections.addRadians(
                std::hypot( std::cos( declination ) * wrappedAngleDifferenceRadians( primaryPrediction( 0 ), comparisonPrediction( 0 ) ),
                            primaryPrediction( 1 ) - comparisonPrediction( 1 ) ) );
    }

    BOOST_TEST_MESSAGE( dataset.name << " " << dataset.primaryConvention
                                     << " RA cos(Dec) / Dec mean/RMS [arcsec]: " << primaryRightAscensionResiduals.meanArcseconds( )
                                     << " / " << primaryRightAscensionResiduals.rmsArcseconds( ) << ", "
                                     << primaryDeclinationResiduals.meanArcseconds( ) << " / "
                                     << primaryDeclinationResiduals.rmsArcseconds( ) );
    BOOST_TEST_MESSAGE( dataset.name << " " << dataset.comparisonConvention
                                     << " RA cos(Dec) / Dec mean/RMS [arcsec]: " << comparisonRightAscensionResiduals.meanArcseconds( )
                                     << " / " << comparisonRightAscensionResiduals.rmsArcseconds( ) << ", "
                                     << comparisonDeclinationResiduals.meanArcseconds( ) << " / "
                                     << comparisonDeclinationResiduals.rmsArcseconds( ) );
    BOOST_TEST_MESSAGE( dataset.name << " correction RMS/max [arcsec]: " << corrections.rmsArcseconds( ) << " / "
                                     << corrections.maximumArcseconds( ) );

    checkResidualBounds( primaryRightAscensionResiduals, dataset.residualBounds );
    checkResidualBounds( primaryDeclinationResiduals, dataset.residualBounds );
    checkCorrectionBounds( corrections, dataset.correctionBounds );
    if( dataset.requirePrimaryRmsBelowComparison )
    {
        BOOST_CHECK_GT( comparisonRightAscensionResiduals.rmsArcseconds( ), primaryRightAscensionResiduals.rmsArcseconds( ) );
        BOOST_CHECK_GT( comparisonDeclinationResiduals.rmsArcseconds( ), primaryDeclinationResiduals.rmsArcseconds( ) );
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
            std::make_shared< ObservationModelSettings >( position_angle, linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > separationSettings =
            std::make_shared< ObservationModelSettings >( separation_distance, linkEnds, lightTimeCorrectionSettings );
    std::shared_ptr< ObservationModelSettings > positionAngleAndSeparationSettings =
            std::make_shared< ObservationModelSettings >( position_angle_and_separation, linkEnds, lightTimeCorrectionSettings );

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

    // Source: Tholen, D.J. and Buie, M.W. (1997), The orbit of Charon, Icarus 125, 245-260, https://doi.org/10.1006/icar.1996.5494
    // (https://www2.boulder.swri.edu/~buie/biblio/pub024.pdf); data from NSDB record pm0001
    // (https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.html). Columns: Julian date at HST exposure midpoint, separation [arcsec],
    // and position angle [deg] w.r.t. J2000 north, astrometric (no stellar aberration). Although NSDB lists the time scale as TDT,
    // the dates agree with the UTC exposure midpoints in the MAST archive to 6 ms, so they are interpreted as UTC. The observations
    // are of the centres of light; the geocentric receiver used here (instead of HST) changes either observable by less than
    // 0.002 mas for these records.
    dataset.decodeRecord = []( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        return ReferenceObservation{ convertUtcJulianDateToTdb( observations( row, 0 ) ),
                                     0,
                                     unit_conversions::convertDegreesToRadians( observations( row, 2 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 1 ) ) };
    };

    // Tholen and Buie (1997) report a post-fit RMS of 2.93 mas (separation) and 2.70 mas (transverse PA). Bounds on pre-fit residuals
    // w.r.t. PLU060: |mean| < 5 mas and 1 mas < RMS < 10 mas, i.e. RMS within a factor of about 3 of the published values.
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

    // Source: Jones, D.H.P., Sinclair, A.T., and Williams, I.P. (1989), Secular acceleration of Phobos confirmed from positions
    // obtained on La Palma, MNRAS 237, 15P-19P, https://doi.org/10.1093/mnras/237.1.15P; data from NSDB record mm0012
    // (https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.html). UTC calendar date; columns 14/15 give astrometric separation [arcsec]
    // and PA [deg] w.r.t. the true equator and equinox of date (modelled here with IAU 1976/1980 precession-nutation).
    // Whether differential refraction was removed is not documented, and it is not modelled. Reject unexpected NSDB format flags.
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

    // Jones et al. (1989) report a scatter of about 0.2 arcsec. Residual bounds: |mean| < 0.05 arcsec (separation) and < 0.1 arcsec
    // (transverse PA), and 0.1 arcsec < RMS < 0.5 arcsec.
    // The difference between true-of-date and J2000 PA (times separation) was computed independently with ERFA pnm80 and SPICE as
    // about 6 mas RMS; it is required to be between 1 and 10 mas RMS, with a maximum between 5 and 50 mas. Changing the pole must
    // leave the separation unchanged (to within rounding error).
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

    // Source: Veillet, C. (1982), Orbital elements of Nereid from new observations, A&A 112, 277-280, and Veillet, C. and Bois, E.
    // (1988), 1982-1984 observations of Nereid - New orbital determination, A&A 189, 303-305; data from NSDB record nm1007
    // (https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm1007.html). UTC; topocentric astrometric offsets x/y [arcsec] w.r.t. the B1950
    // equator. These offsets are converted to P/S here, so that this data set provides a secondary check of the B1950 pole.
    // The atmospheric reduction is not documented, and is not modelled.
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

    // Yuan et al. (2021), A&A 645, A48, https://doi.org/10.1051/0004-6361/202038776, report 0.17--0.22 arcsec coordinate RMS for
    // these historical series. Residual bounds: |mean| < 0.1 arcsec (separation) and < 0.2 arcsec (transverse PA), and RMS between
    // 0.05 (separation) or 0.1 (PA) and 0.5 arcsec. The B1950 vs. J2000 PA difference (times separation), caused by the 0.28 deg
    // rotation of the pole, is required to be between 0.2 and 0.5 arcsec RMS, and between 0.5 and 1.0 arcsec maximum. As it exceeds
    // the noise, the B1950 pole must also give a smaller transverse PA RMS than J2000.
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

    // Source: Tomita, K. and Soma, M. (1979), Astrometric observations of satellites of Uranus, Tokyo Astronomical Bulletin, Second
    // Series 261, 2977-2981; data from NSDB record um0027 (https://nsdb.imcce.fr/obspos/OBS_COLL/U/um0027.html). Topocentric
    // astrometric P/S w.r.t. the B1950 equator, with ET calendar time (not UTC).
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

    // Jacobson, R.A. (2014), AJ 148, 76, https://doi.org/10.1088/0004-6256/148/5/76, Table 3, gives a post-fit RMS of 0.443 arcsec
    // (separation) and 0.363 arcsec (transverse) for the records retained in his fit; this unfiltered set contains outliers, so the
    // RMS upper bound is relaxed to 2 arcsec (lower bound 0.2 arcsec, |mean| < 0.2 arcsec). The B1950 vs. J2000 PA difference (times
    // separation) is required to be between 0.02 and 0.1 arcsec RMS (maximum between 0.05 and 0.2 arcsec); B1950 must still give a
    // smaller transverse PA RMS than J2000.
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
    dataset.primarySettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, true );
    dataset.comparisonConvention = "astrometric true-of-date";
    dataset.comparisonSettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, false );
    const std::map< int, std::string > observedSatellites = { { 602, "Enceladus" }, { 603, "Tethys" }, { 604, "Dione" }, { 605, "Rhea" } };
    for( const auto& observedSatellite : observedSatellites )
    {
        dataset.models[ observedSatellite.first ] = createPositionSeparationModel( bodies, "Sheshan", "Titan", observedSatellite.second );
    }

    // Source: Qiao, R.C., Shen, K.X., Liu, J.R., and Harper, D. (1999), 1994-1996 CCD astrometric observations of Saturn's
    // satellites and comparison with theories, A&AS 137, 1-5, https://doi.org/10.1051/aas:1999238, Table 4: apparent topocentric P/S
    // of satellites w.r.t. Titan, w.r.t. the true equator of date, with differential refraction removed, and stellar aberration and
    // parallax retained. Columns: UTC Julian date, measured-satellite NAIF ID, PA [deg], separation [arcsec].
    dataset.decodeRecord = []( const Eigen::MatrixXd& observations, const Eigen::Index row ) {
        return ReferenceObservation{ convertUtcJulianDateToTdb( observations( row, 0 ) ),
                                     static_cast< int >( observations( row, 1 ) ),
                                     unit_conversions::convertDegreesToRadians( observations( row, 2 ) ),
                                     unit_conversions::convertArcSecondsToRadians( observations( row, 3 ) ) };
    };

    // Qiao et al. (1999) report a coordinate RMS of roughly 0.08--0.20 arcsec. Residual bounds: |mean| < 0.2 arcsec and RMS between
    // 0.1 and 0.5 arcsec (separation), and |mean| < 0.1 arcsec and RMS between 0.05 and 0.2 arcsec (transverse PA).
    // The difference between aberrated and astrometric directions is only differential aberration (of order |v/c| ~ 1.0E-4 times
    // the separation), of a few mas as computed independently with SpiceyPy. It is required to be
    // between 2 and 10 mas RMS (maximum between 5 and 20 mas) in separation, and between 0.1 and 1 mas in transverse PA.
    dataset.separationBounds = { 0.2, 0.1, 0.5 };
    dataset.transversePositionAngleBounds = { 0.1, 0.05, 0.2 };
    dataset.separationCorrectionBounds = CorrectionBounds{ 0.002, 0.01, 0.005, 0.02 };
    dataset.transversePositionAngleCorrectionBounds = CorrectionBounds{ 0.0001, 0.001, 0.0001, 0.001 };
    validateReferenceDataset( dataset );
}

// The angular position tests below use 16 observations, evenly spread in time, of series in the IMCCE planetary observation
// database (https://ftp.imcce.fr/pub/databases/PODB/, URSS files) and the JPL planetary optical observation files
// (https://ssd.jpl.nasa.gov/dat/planets/, transit files). All observations are geocentric, of the planet centre; the planet
// system barycentre of INPOP19a is used (offset < 0.1 arcsec). Epochs are converted from the published time scale (ET, UTC, or
// UT1 with the Espenak and Meeus Delta T model) to TDB. The stated scatter is the RMS of the residuals of these observations
// w.r.t. INPOP19a; the residual bounds allow a margin of a factor of about 2, and a mean offset up to about the scatter.

//! Validate B1950 angular positions on Nikolaev photographic observations of Saturn, against J2000.
BOOST_AUTO_TEST_CASE( testSaturnB1950AngularPositionPrefitResiduals )
{
    AngularPositionDataset dataset;
    dataset.name = "Saturn Nikolaev 1961-1982";
    dataset.observedBody = "Saturn_barycenter";
    // IMCCE file URSS/SATU/PHOTO/nikol.pht.ftp (Nikolaev, MPC 089): astrometric, B1950 (FK4) photographic positions, in ET
    dataset.observations = { { -1211901454.399, 297.034412500, -21.304166667 }, { -1143826525.311, 319.374250000, -17.077363889 },
                             { -1110305819.520, 330.951766667, -13.693983333 }, { -1052310466.656, 358.868366667, -3.147594444 },
                             { -1050326762.400, 357.320479167, -3.858602778 },  { -1041797257.632, 354.949579167, -4.627280556 },
                             { -983975788.079, 20.367250000, 5.588330556 },     { -923826808.816, 50.265141667, 15.910541667 },
                             { -859541788.622, 79.534483333, 21.454288889 },    { -826109869.027, 94.651400000, 22.283733333 },
                             { -784831136.204, 103.160050000, 22.585344444 },   { -751053678.264, 117.869779167, 21.350788889 },
                             { -693056615.968, 151.746475000, 13.166172222 },   { -651518896.398, 159.149183333, 10.891288889 },
                             { -592409756.937, 186.720041667, -0.017141667 },   { -555785186.984, 195.427008333, -3.669505556 } };
    dataset.primaryConvention = "astrometric B1950";
    dataset.primarySettings = getPositionAngleAncillarySettings( b1950_position_angle_reference_frame );
    dataset.comparisonConvention = "astrometric J2000 (default)";
    dataset.comparisonSettings = nullptr;
    // Scatter of these observations about 0.34 arcsec (RA) and 0.53 arcsec (Dec). The B1950 and J2000 frames differ by a rotation of
    // 2515 arcsec (50 years of precession) about the ecliptic pole, which bounds the correction.
    dataset.residualBounds = { 0.5, 0.1, 1.0 };
    dataset.correctionBounds = CorrectionBounds{ 1000.0, 2520.0, 2000.0, 2520.0 };
    validateAngularPositionDataset( dataset );
}

//! Validate the astrometric B1950 frame on Pulkovo observations of Uranus over 50 years, against aberrated B1950 directions.
BOOST_AUTO_TEST_CASE( testUranusB1950AngularPositionPrefitResiduals )
{
    AngularPositionDataset dataset;
    dataset.name = "Uranus Pulkovo 1919-1969";
    dataset.observedBody = "Uranus_barycenter";
    // IMCCE file URSS/URAN/PHOTO/pulk.pht.ftp (Pulkovo, MPC 084): astrometric, B1950 (FK4) photographic positions, in ET
    dataset.observations = { { -2533391314.838, 331.459500000, -12.546722222 }, { -2086876662.942, 23.009208333, 8.976583333 },
                             { -2059999513.977, 28.936291667, 11.246972222 },   { -2015349582.953, 30.429250000, 11.883972222 },
                             { -1900415146.293, 48.961791667, 17.772083333 },   { -1859133693.404, 49.931041667, 18.101388889 },
                             { -1605377899.187, 86.330791667, 23.614500000 },   { -1543333320.487, 96.535875000, 23.638722222 },
                             { -1478967522.904, 106.149875000, 23.057388889 },  { -1412011531.153, 115.476750000, 21.919111111 },
                             { -1321881740.758, 131.602625000, 18.667472222 },  { -1255019212.070, 139.844166667, 16.392333333 },
                             { -1163426760.901, 155.460791667, 11.070666667 },  { -1102157707.312, 165.546666667, 7.052777778 },
                             { -1037960654.028, 174.276750000, 3.348916667 },   { -972134061.179, 182.256916667, -0.132055556 } };
    dataset.primaryConvention = "astrometric B1950";
    dataset.primarySettings = getPositionAngleAncillarySettings( b1950_position_angle_reference_frame, TUDAT_NAN, false );
    dataset.comparisonConvention = "aberrated B1950";
    dataset.comparisonSettings = getPositionAngleAncillarySettings( b1950_position_angle_reference_frame, TUDAT_NAN, true );
    // Scatter of these observations about 0.39 arcsec (RA) and 0.77 arcsec (Dec), with a declination offset of about -0.7 arcsec. The
    // correction is the annual stellar aberration, of at most 20.8 arcsec.
    dataset.residualBounds = { 1.0, 0.1, 1.5 };
    dataset.correctionBounds = CorrectionBounds{ 10.0, 21.0, 15.0, 21.0 };
    validateAngularPositionDataset( dataset );
}

//! Validate apparent true-of-date (IAU 1976/1980) angular positions on USNO transit observations of Neptune.
BOOST_AUTO_TEST_CASE( testNeptuneTrueOfDateIau1976AngularPositionPrefitResiduals )
{
    AngularPositionDataset dataset;
    dataset.name = "Neptune USNO 6-inch transit 1926-1976";
    dataset.observedBody = "Neptune_barycenter";
    // JPL file trnstswash6a.txt (USNO 6-inch transit circle, MPC 786): apparent places w.r.t. the true equator and equinox of date,
    // with right ascensions on the FK4 equinox, in UT1 (UTC from 1972)
    dataset.observations = { { -2333940291.851, 146.589616667, 13.844527778 },  { -2235884218.011, 152.180000000, 12.084552778 },
                             { -2045372404.121, 164.765833333, 7.560694444 },   { -1947144055.575, 170.134375000, 5.493802778 },
                             { -1820394720.984, 178.302204167, 2.203830556 },   { -1695454723.757, 186.923000000, -1.334416667 },
                             { -1571893243.139, 196.003454167, -5.031447222 },  { -1470477042.826, 200.456958333, -6.761227778 },
                             { -1343813840.306, 208.760154167, -9.871511111 },  { -1251788877.660, 215.911562500, -12.335652778 },
                             { -1188026453.488, 220.071579167, -13.671427778 }, { -1092555115.130, 226.440833333, -15.554447222 },
                             { -996222168.497, 232.671579167, -17.196297222 },  { -908591214.573, 241.328504167, -19.109608333 },
                             { -776844388.407, 249.380379167, -20.466941667 },  { -738156584.209, 249.907091667, -20.608386111 } };
    dataset.correctFk4Equinox = true;
    dataset.primaryConvention = "aberrated true-of-date IAU 1976/1980";
    dataset.primarySettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, true );
    dataset.comparisonConvention = "astrometric true-of-date IAU 1976/1980";
    dataset.comparisonSettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, false );
    // Scatter of these observations about 0.38 arcsec (RA) and 0.36 arcsec (Dec). The correction is the annual stellar aberration.
    dataset.residualBounds = { 0.5, 0.1, 1.0 };
    dataset.correctionBounds = CorrectionBounds{ 10.0, 21.0, 15.0, 21.0 };
    validateAngularPositionDataset( dataset );
}

//! Validate apparent true-of-date (IAU 2006/2000A) angular positions on Carlsberg Meridian Circle observations of Uranus.
BOOST_AUTO_TEST_CASE( testUranusTrueOfDateIau2006AngularPositionPrefitResiduals )
{
    AngularPositionDataset dataset;
    dataset.name = "Uranus Carlsberg Meridian Circle 1984-1997";
    dataset.observedBody = "Uranus_barycenter";
    // JPL file trnstslplma.txt (Carlsberg Automatic Meridian Circle, La Palma, MPC 950): apparent places w.r.t. the true equator
    // and equinox of date, not on the FK4 equinox, in UTC
    dataset.observations = { { -493464434.935, 250.891083333, -22.241825000 }, { -465115024.182, 256.843358333, -22.909416667 },
                             { -454690148.743, 252.746854167, -22.533241667 }, { -425651486.051, 258.445591667, -23.112927778 },
                             { -397301938.181, 264.971108333, -23.511983333 }, { -365161626.378, 269.587620833, -23.652047222 },
                             { -332418247.896, 273.865675000, -23.664572222 }, { -301570260.723, 279.052395833, -23.504730556 },
                             { -273220872.643, 284.903008333, -23.098722222 }, { -266673054.434, 282.175900000, -23.352236111 },
                             { -239529741.543, 288.959562500, -22.721580556 }, { -203254060.703, 291.449245833, -22.462500000 },
                             { -169304555.710, 295.005920833, -21.983888889 }, { -135613471.266, 298.865458333, -21.364558333 },
                             { -103128596.867, 303.086050000, -20.577372222 }, { -68920443.365, 307.294970833, -19.671072222 } };
    dataset.primaryConvention = "aberrated true-of-date IAU 2006/2000A";
    dataset.primarySettings =
            getPositionAngleAncillarySettings( true_of_date_iau_2006_2000a_position_angle_reference_frame, TUDAT_NAN, true );
    dataset.comparisonConvention = "aberrated J2000";
    dataset.comparisonSettings = getPositionAngleAncillarySettings( j2000_position_angle_reference_frame, TUDAT_NAN, true );
    // Scatter of these observations about 0.20 arcsec (RA) and 0.18 arcsec (Dec). The correction is the precession since J2000 (50.3 arcsec
    // per year, 2.2 to 15.6 years) and nutation (at most 19 arcsec).
    dataset.residualBounds = { 0.2, 0.05, 0.5 };
    dataset.correctionBounds = CorrectionBounds{ 300.0, 800.0, 600.0, 820.0 };
    validateAngularPositionDataset( dataset );
}

//! Compare the IAU 2006/2000A and IAU 1976/1980 true-of-date frames on Tokyo transit observations of Jupiter.
BOOST_AUTO_TEST_CASE( testJupiterTrueOfDateAngularPositionPrefitResiduals )
{
    AngularPositionDataset dataset;
    dataset.name = "Jupiter Tokyo transit 1963-1974";
    dataset.observedBody = "Jupiter_barycenter";
    // IMCCE file URSS/JUPI/TRANSIT/tok.trn.ftp (Tokyo, MPC 388): apparent places w.r.t. the true equator and equinox of date,
    // with right ascensions on the FK4 equinox, in UT1
    dataset.observations = { { -1147629415.887, 18.367933333, 6.144411111 },   { -1139101305.781, 9.446308333, 2.461691667 },
                             { -1108159143.168, 48.104729167, 16.627222222 },  { -1103679607.840, 43.937125000, 15.694888889 },
                             { -1074890446.138, 87.197300000, 23.012586111 },  { -1041622161.822, 124.565637500, 20.142552778 },
                             { -1033438447.544, 116.732433333, 21.772444444 }, { -1003617530.095, 150.745700000, 13.297125000 },
                             { -1000429993.860, 148.508475000, 14.028550000 }, { -975174999.223, 185.756525000, -0.896938889 },
                             { -939496662.623, 212.500791667, -11.549630556 }, { -905454487.707, 243.240225000, -20.129088889 },
                             { -864434277.905, 268.727487500, -23.361833333 }, { -839694187.769, 314.663150000, -17.678372222 },
                             { -830132229.041, 305.237745833, -20.308055556 }, { -824617006.884, 308.245458333, -19.576202778 } };
    dataset.correctFk4Equinox = true;
    dataset.primaryConvention = "aberrated true-of-date IAU 2006/2000A";
    dataset.primarySettings =
            getPositionAngleAncillarySettings( true_of_date_iau_2006_2000a_position_angle_reference_frame, TUDAT_NAN, true );
    dataset.comparisonConvention = "aberrated true-of-date IAU 1976/1980";
    dataset.comparisonSettings =
            getPositionAngleAncillarySettings( true_of_date_iau_1976_1980_position_angle_reference_frame, TUDAT_NAN, true );
    // Scatter of these observations about 0.64 arcsec (RA) and 0.72 arcsec (Dec). The two frames differ by the IAU 2000 precession-rate
    // correction (-0.30 arcsec per century in longitude, about 0.09 arcsec for these epochs) and the frame bias (about 0.02
    // arcsec), much less than the scatter, so no RMS improvement is required.
    dataset.residualBounds = { 0.5, 0.2, 1.5 };
    dataset.correctionBounds = CorrectionBounds{ 0.05, 0.2, 0.08, 0.2 };
    dataset.requirePrimaryRmsBelowComparison = false;
    validateAngularPositionDataset( dataset );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
