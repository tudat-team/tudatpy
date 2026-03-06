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

#include <array>
#include <filesystem>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/math/basic/coordinateConversions.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"

#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/io/readInpopEphemerisFile.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/interface/spice/spiceEphemeris.h"

#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"
#include "tudat/simulation/environment_setup/createMetric.h"

#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"
#include "tudat/math/basic/leastSquaresEstimation.h"



namespace tudat
{

namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::unit_conversions;
using namespace tudat::input_output;
using namespace tudat::basic_astrodynamics;

namespace
{

bool fileExists( const std::string& filePath )
{
    return std::filesystem::exists( std::filesystem::path( filePath ) );
}

bool areInpop19aResourcesAvailable( const std::string& spiceKernelsPath, const std::string& textKernelsPath )
{
    const std::array< std::string, 7 > requiredSpiceKernels = {
        "pck00010.tpc",
        "naif0012.tls",
        "inpop19a_TDB_m100_p100_spice.tpc",
        "inpop19a_TDB_m100_p100_spice.bsp",
        "inpop19a_TDB_m100_p100_spice.bpc",
        "codes_300ast_20100725.bsp",
        "codes_300ast_20100725.tf"
    };

    for( const std::string& kernelFile : requiredSpiceKernels )
    {
        if( !fileExists( spiceKernelsPath + "/" + kernelFile ) )
        {
            BOOST_TEST_MESSAGE( "Skipping INPOP test: missing SPICE kernel " + spiceKernelsPath + "/" + kernelFile );
            return false;
        }
    }

    if( !fileExists( textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_TCG.asc" ) )
    {
        BOOST_TEST_MESSAGE( "Skipping INPOP test: missing ASCII ephemeris directory " + textKernelsPath );
        return false;
    }

    return true;
}

void loadInpop19aSpiceKernels( const std::string& spiceKernelsPath )
{
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/pck00010.tpc" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/inpop19a_TDB_m100_p100_spice.tpc" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/inpop19a_TDB_m100_p100_spice.bsp" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/inpop19a_TDB_m100_p100_spice.bpc" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/codes_300ast_20100725.bsp" );
    spice_interface::loadSpiceKernelInTudat( spiceKernelsPath + "/codes_300ast_20100725.tf" );
}

}  // namespace


BOOST_AUTO_TEST_SUITE( test_RelativisticConversions )


BOOST_AUTO_TEST_CASE( test_tcb_to_tcg_conversion )
{

    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string textKernelsPath = spiceKernelsPath + "/inpop19a_TCB_m100_p100_asc";

    if( !areInpop19aResourcesAvailable( spiceKernelsPath, textKernelsPath ) )
    {
        return;
    }

    // Load SPICE kernels (INPOP19a + core kernels).
    loadInpop19aSpiceKernels( spiceKernelsPath );

    // Map from SPICE ID string to body name
    std::map< std::string, std::string > bodyIdToName = {
        { "1", "Mercury" },
        { "2", "Venus" },
        { "3", "Earth" },
        { "4", "Mars" },
        { "5", "Jupiter" },
        { "6", "Saturn" },
        { "7", "Uranus" },
        { "8", "Neptune" },
        { "9", "Pluto" },
        { "10", "Sun" },
        { "301", "Moon" },
        { "2000001", "Ceres" },
        { "2000002", "Pallas" },
        { "2000004", "Vesta" },
    };

    SystemOfBodies bodies;
    for ( const auto& idToNamePair : bodyIdToName )
    {
         const std::string& id = idToNamePair.first;
         const std::string& name = idToNamePair.second;
         std::shared_ptr< Body > body = std::make_shared< Body >();
         double gm = spice_interface::getBodyGravitationalParameter( id ) / ( 1.0 - physical_constants::LB_TIME_RATE_TERM );
         body->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( gm ) );
         bodies.addBody( body, name );
     }

    // Specify initial time
    double initialEphemerisTime = -365.25 * 86400.0 * 2.0;
    double finalEphemerisTime = 365.25 * 86400.0 * 2.0; 
    double maximumTimeStep = 3600.0;
    double numberOfTimeStepBuffer = 6.0;
    double buffer = numberOfTimeStepBuffer * maximumTimeStep;
    std::string centralBody = "Earth"; 

    bodies.at( "Sun" )->setEphemeris( createInpopEphemerisFromFiles(
                                        textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Sun.asc",
                                        textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Sun.asc" ) );    
    bodies.at( "Mercury" )->setEphemeris( createInpopEphemerisFromFiles(
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Mer.asc",
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Mer.asc" ) );
    bodies.at( "Venus" )->setEphemeris( createInpopEphemerisFromFiles(
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Ven.asc",
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Ven.asc" ) );
    bodies.at( "Earth" )->setEphemeris( createInpopEphemerisFromFiles(
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Ear.asc",
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Ear.asc" ) );
    bodies.at( "Moon" )->setEphemeris( createInpopEphemerisFromFiles(
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Moo.asc",
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Moo.asc",
                                         basic_astrodynamics::JULIAN_DAY_ON_J2000, 1 ) );
    bodies.at( "Mars" )->setEphemeris( createInpopEphemerisFromFiles(
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Mar.asc",
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Mar.asc" ) );
    bodies.at( "Jupiter" )->setEphemeris( createInpopEphemerisFromFiles(
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Jup.asc",
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Jup.asc" ) );
    bodies.at( "Saturn" )->setEphemeris( createInpopEphemerisFromFiles(
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Sat.asc",
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Sat.asc" ) );
    bodies.at( "Uranus" )->setEphemeris( createInpopEphemerisFromFiles(
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Ura.asc",
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Ura.asc" ) );
    bodies.at( "Neptune" )->setEphemeris( createInpopEphemerisFromFiles(
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Nep.asc",
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Nep.asc" ) );
    bodies.at( "Pluto" )->setEphemeris( createInpopEphemerisFromFiles(
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Plu.asc",
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Plu.asc" ) );

    bodies.at( "Ceres" )->setEphemeris( createTabulatedEphemerisFromSpice(
                                          "Ceres", initialEphemerisTime - buffer, finalEphemerisTime + buffer, 7200.0, "SSB", "ECLIPJ2000" ) );
    bodies.at( "Pallas" )->setEphemeris( createTabulatedEphemerisFromSpice(
                                          "Pallas", initialEphemerisTime - buffer, finalEphemerisTime + buffer, 7200.0, "SSB", "ECLIPJ2000" ) );
    bodies.at( "Vesta" )->setEphemeris( createTabulatedEphemerisFromSpice(
                                            "Vesta", initialEphemerisTime - buffer, finalEphemerisTime + buffer, 7200.0, "SSB", "ECLIPJ2000" ) );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    std::vector< std::string > externalBodies;
    for ( const auto& idToNamePair : bodyIdToName )
    {
        const std::string& name = idToNamePair.second;
        if ( name != centralBody )
        {
            externalBodies.push_back( name );
        }
    }

    double startTime = initialEphemerisTime;
    double endTime = finalEphemerisTime;
    double timeStep = 500; //6000.0;

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            numerical_integrators::rungeKutta4Settings( timeStep );
    integratorSettings->initialTimeDeprecated_ = startTime;
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( endTime );

    auto outputProcessingSettings = std::make_shared< SingleArcPropagatorProcessingSettings >(
                false,
                false,
                1,
                TUDAT_NAN,
                std::make_shared< PropagationPrintSettings >( ) );
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList{};

    std::shared_ptr< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings<double, double > > properTimeEquationSettings =
            std::make_shared< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings<double, double > >(  
                centralBody, externalBodies, startTime, integratorSettings, terminationSettings,
                ( std::map< std::string, std::pair< int, int > >( ) ),
                std::vector< std::string  >( ),
                []( const double inputTime ){ return inputTime; },
                1.0,
                dependentVariablesList
                //outputProcessingSettings 
            );
    properTimeEquationSettings->getOutputSettings( )->setIntegratedResult( true );

    SingleArcDynamicsSimulator< > timeEquationPropagator = SingleArcDynamicsSimulator< >( bodies, properTimeEquationSettings );

    std::string timeDifferenceFileName = textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_TCG.asc";

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, long double > > timeEphemerisInterpolator =
            input_output::createLongInpopTimeEphemerisInterpolator( timeDifferenceFileName );
    std::map< double, double > timeDifferences2;

    int counter = 0;
    long double initialDifference = 0.0L;
    long double rawDifference;

    std::shared_ptr< Body > earth = bodies.getBody( "Earth" );
    std::shared_ptr< TimeEphemeris > earthTimeEphemeris = bodies.getBody( "Earth" )->getTimeScaleConverter( );

    std::function< double( const double ) > timeDifferenceFunction =
            earthTimeEphemeris->getTimeDifferenceFunction( basic_astrodynamics::barycentric_coordinate_time_scale, basic_astrodynamics::body_centered_coordinate_time_scale, "" );
    
    double testTimeStep = 7100.0; //To prevent excessive resonance with integration step.
    double currentTime = initialEphemerisTime + 5.0 * timeStep;
    while( currentTime < finalEphemerisTime - 5.0 * timeStep )
    {
        rawDifference = timeEphemerisInterpolator->interpolate( currentTime )-
                static_cast< long double >( timeDifferenceFunction( currentTime ) );

        if( counter == 0 )
        {
            initialDifference = rawDifference;
            counter += 1;
        }

        timeDifferences2[ currentTime ] = static_cast< double >( rawDifference - initialDifference );
        currentTime += testTimeStep;
    }

    const std::string diagnosticsOutputDirectory = "/Users/michael.plumaris/Downloads/";
    if( std::filesystem::exists( diagnosticsOutputDirectory ) )
    {
        input_output::writeDataMapToTextFile(
                    timeDifferences2,
                    "tcgMinusTcbInpop2_tdb.dat",
                    diagnosticsOutputDirectory,
                    "",
                    16 );
    }

    Eigen::VectorXd timesVector = utilities::convertStlVectorToEigenVector(
                utilities::createVectorFromMapKeys( timeDifferences2 ) );
    Eigen::VectorXd resultDifferenceVector = utilities::convertStlVectorToEigenVector(
                utilities::createVectorFromMapValues( timeDifferences2 ) );
    
    std::vector<double> dummy = { 0.0, 1.0 };
    Eigen::VectorXd trendFit = linear_algebra::getLeastSquaresPolynomialFit( timesVector, resultDifferenceVector, dummy );

    BOOST_CHECK_SMALL( std::fabs( trendFit[ 1 ] ), 5.0E-18 );

    Eigen::VectorXd resultDifferenceWithoutTrend = resultDifferenceVector - (
        Eigen::VectorXd::Constant( resultDifferenceVector.rows( ), trendFit[ 0 ] ) + trendFit[ 1 ] * timesVector );

    double maximumDifference = resultDifferenceWithoutTrend.maxCoeff( );
    double minimumDifference = resultDifferenceWithoutTrend.minCoeff( );
    const double maxAbsDifference = std::max( std::fabs( maximumDifference ), std::fabs( minimumDifference ) );
    std::cout << "[test_tcb_to_tcg_conversion] max_abs_diff=" << maxAbsDifference << std::endl;

    BOOST_CHECK_SMALL( maximumDifference, 5.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( minimumDifference ), 5.0E-12 );

    std::cout<< "maximumDifference" << maximumDifference << std::endl;
    
}

BOOST_AUTO_TEST_CASE( test_concatenated_conversions )
{
    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string textKernelsPath = spiceKernelsPath + "/inpop19a_TCB_m100_p100_asc";
    if( !areInpop19aResourcesAvailable( spiceKernelsPath, textKernelsPath ) )
    {
        return;
    }
    loadInpop19aSpiceKernels( spiceKernelsPath );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );

    const std::map< std::string, std::string > gmBodyNameOverrides = {
        { "Jupiter", "Jupiter Barycenter" },
        { "Saturn", "Saturn Barycenter" }
    };

    // Specify initial time
    double initialEphemerisTime = -365.25 * 86400.0 * 1.0;
    double finalEphemerisTime = 365.25 * 86400.0 * 1.0;
    double maximumTimeStep = 3600.0;
    double numberOfTimeStepBuffer = 6.0;
    double buffer = numberOfTimeStepBuffer * maximumTimeStep;
    std::string centralBody = "Earth";


    SystemOfBodies bodies;
    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        if( bodyNames[ i ] != "Earth" )
        {
            std::shared_ptr< Body > body = std::make_shared< Body >( );
            const std::string gmBodyName = ( gmBodyNameOverrides.count( bodyNames[ i ] ) == 0 ) ?
                        bodyNames[ i ] : gmBodyNameOverrides.at( bodyNames[ i ] );
            body->setGravityFieldModel(
                        std::make_shared< gravitation::GravityFieldModel >(
                            spice_interface::getBodyGravitationalParameter( gmBodyName ) / ( 1.0 - physical_constants::LB_TIME_RATE_TERM ) ) );
            bodies.addBody( body, bodyNames[ i ] );
        }
    }

    
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    bodies.addBody( earth, "Earth" );
    earth->setShapeModel( createBodyShapeModel( getDefaultBodyShapeSettings( "Earth", initialEphemerisTime, finalEphemerisTime ), "Earth" ) );
    earth->setRotationalEphemeris( createRotationModel( getDefaultRotationModelSettings( "Earth", initialEphemerisTime, finalEphemerisTime ), "Earth" ) );

    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > earthGravityFieldSettings =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                getDefaultGravityFieldSettings( "Earth", initialEphemerisTime, finalEphemerisTime ) );
    earthGravityFieldSettings->resetAssociatedReferenceFrame( "IAU_Earth" );

    earth->setGravityFieldModel( createGravityFieldModel( earthGravityFieldSettings, "Earth", bodies ) );

    std::shared_ptr< Body > lro = std::make_shared< Body >( );

    Eigen::Vector6d lroInitialKeplerianElements;
    lroInitialKeplerianElements[ semiMajorAxisIndex ] = 2500.0E3;
    lroInitialKeplerianElements[ eccentricityIndex ] = 0.1;
    lroInitialKeplerianElements[ inclinationIndex ] = 0.375 * mathematical_constants::PI;
    lroInitialKeplerianElements[ argumentOfPeriapsisIndex ] = 0.0;
    lroInitialKeplerianElements[ longitudeOfAscendingNodeIndex ] = 0.0;
    lroInitialKeplerianElements[ trueAnomalyIndex ] = 0.0;

    lro->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                           lroInitialKeplerianElements, initialEphemerisTime, spice_interface::getBodyGravitationalParameter( "Moon" ), "Moon"  ) );

    bodies.addBody( lro, "LRO" );

    bodies.at( "Sun" )->setEphemeris( createInpopEphemerisFromFiles(
                                        textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Sun.asc",
                                        textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Sun.asc" ) );
    bodies.at( "Earth" )->setEphemeris( createInpopEphemerisFromFiles(
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Ear.asc",
                                          textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Ear.asc" ) );
    bodies.at( "Moon" )->setEphemeris( createInpopEphemerisFromFiles(
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Moo.asc",
                                         textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Moo.asc",
                                         basic_astrodynamics::JULIAN_DAY_ON_J2000, 1 ) );
    bodies.at( "Jupiter" )->setEphemeris( createInpopEphemerisFromFiles(
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Jup.asc",
                                            textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Jup.asc" ) );
    bodies.at( "Saturn" )->setEphemeris( createInpopEphemerisFromFiles(
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_pos_Sat.asc",
                                           textKernelsPath + "/inpop19a_TCB_m100_p100_asc_vel_Sat.asc" ) );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] = ( Eigen::Vector3d( ) << 4194511.7, 1162789.7, 4647362.5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Earth", "Yarragadee" ) ] = ( Eigen::Vector3d( ) << -2389008, 5043332, -3078526 ).finished( );

    createGroundStations( bodies, groundStationsToCreate );
    
    std::vector< std::string > externalBodies;
    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        if( bodyNames[ i ] != centralBody )
        {
            externalBodies.push_back( bodyNames[ i ] );
        }
    }

    double startTime = initialEphemerisTime;
    double endTime = finalEphemerisTime;
    double timeStep = 3000.0;

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            numerical_integrators::rungeKutta4Settings( timeStep );
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( endTime );

    std::vector< std::string > listOfPerturbingBodies{ "Earth",  "Moon",  "Sun", "Jupiter", "Saturn" };
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRelativisticTimeState = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 1 );

    auto outputProcessingSettings = std::make_shared< SingleArcPropagatorProcessingSettings >(
                false,
                false,
                1,
                TUDAT_NAN,
                std::make_shared< PropagationPrintSettings >( ) );
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList{};


    const std::vector< std::string > topocentricPerturbingBodies{ "Moon", "Sun", "Jupiter", "Saturn" };

    std::vector< std::shared_ptr< RelativisticTimeStatePropagatorSettings< double, double > > > bodyCentricToTopocentricConversionSettings;
    bodyCentricToTopocentricConversionSettings.push_back(
                std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                    std::make_pair( "Earth", "Graz" ), 0, 4, 0, topocentricPerturbingBodies,
                    initialRelativisticTimeState, initialEphemerisTime, integratorSettings, terminationSettings,
                    dependentVariablesList,
                    outputProcessingSettings
                ) );
    bodyCentricToTopocentricConversionSettings.push_back(
                std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                    std::make_pair( "Earth", "Yarragadee" ), 0, 4, 0, topocentricPerturbingBodies,
                    initialRelativisticTimeState, initialEphemerisTime, integratorSettings, terminationSettings,
                    dependentVariablesList,
                    outputProcessingSettings
                ) );

    std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings<> > > relativisticConverterSettings;
    relativisticConverterSettings[ "LRO" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
                std::make_shared< propagators::FirstOrderBodycentricRelativisticTimePropagatorSettings< double, double > >(
                    "LRO", listOfPerturbingBodies, initialEphemerisTime, integratorSettings, terminationSettings ),
                integratorSettings ); 

    std::vector< std::string > listOfPerturbingBodies2{ "Moon",  "Sun", "Jupiter", "Saturn" };
    relativisticConverterSettings[ "Earth" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
                std::make_shared< propagators::SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "Earth", listOfPerturbingBodies2, initialEphemerisTime, integratorSettings, terminationSettings ),
                integratorSettings,
                bodyCentricToTopocentricConversionSettings );

    setRelativisticTimeConverters( bodies, relativisticConverterSettings );

    std::shared_ptr< TimeEphemeris > earthTimeScaleConverter = earth->getTimeScaleConverter( );
    std::shared_ptr< TimeEphemeris > lroTimeScaleConverter = lro->getTimeScaleConverter( );

    BOOST_CHECK_EQUAL( ( lroTimeScaleConverter == nullptr ), 0 );
    BOOST_CHECK_EQUAL( ( earthTimeScaleConverter == nullptr ), 0 );

    BOOST_CHECK_SMALL( lroTimeScaleConverter->getTimeDifference(
                           body_centered_coordinate_time_scale, barycentric_coordinate_time_scale, initialEphemerisTime ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( earthTimeScaleConverter->getTimeDifference(
                           body_centered_coordinate_time_scale, barycentric_coordinate_time_scale, initialEphemerisTime ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( lroTimeScaleConverter->getTimeDifference(
                           barycentric_coordinate_time_scale, body_centered_coordinate_time_scale, initialEphemerisTime ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( earthTimeScaleConverter->getTimeDifference(
                           barycentric_coordinate_time_scale, body_centered_coordinate_time_scale, initialEphemerisTime ),
                       std::numeric_limits< double >::epsilon( ) );


    BOOST_CHECK_SMALL( earthTimeScaleConverter->getTimeDifference(
                           body_centered_coordinate_time_scale, local_proper_time_scale, initialEphemerisTime, "Graz" ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( earthTimeScaleConverter->getTimeDifference(
                           local_proper_time_scale, body_centered_coordinate_time_scale, initialEphemerisTime, "Graz" ),
                       std::numeric_limits< double >::epsilon( ) );

    std::shared_ptr< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > > directEarthTimeScaleConverter =
            std::make_shared< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                "Earth", listOfPerturbingBodies2, initialEphemerisTime, integratorSettings, terminationSettings );

    directEarthTimeScaleConverter->getOutputSettings( )->setIntegratedResult( true );

    // Get directly calculated map of tcg-tcb from tcb input (key)
    SingleArcDynamicsSimulator< > timeEquationPropagator = SingleArcDynamicsSimulator< >(
                bodies, directEarthTimeScaleConverter, true );

    std::map< double, Eigen::VectorXd > directTimeDifferencesVectors = timeEquationPropagator.getEquationsOfMotionNumericalSolution( );
    std::map< double, double > directTimeDifferences;
    for( std::map< double, Eigen::VectorXd >::iterator resultIterator = directTimeDifferencesVectors.begin( ); resultIterator !=
         directTimeDifferencesVectors.end( ); resultIterator++ )
    {
        directTimeDifferences[ resultIterator->first ] = resultIterator->second.x( );
    }

    // Create map of tcb-tcg from tcg input (key)
    std::map< double, double > directInverseTimeDifferences;
    for( std::map< double, double >::iterator differenceIterator = directTimeDifferences.begin( ); differenceIterator !=
         directTimeDifferences.end( ); differenceIterator++ )
    {
        directInverseTimeDifferences[ differenceIterator->first + differenceIterator->second ] = -differenceIterator->second;
    }

    // Get time difference functions from indirect calculator.
    std::function< double( const double ) > indirectDifferenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                barycentric_coordinate_time_scale, body_centered_coordinate_time_scale );
    std::function< double( const double ) > indirectInverseDifferenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                body_centered_coordinate_time_scale, barycentric_coordinate_time_scale );

    // Iterate over all directly calculated function values, and use indirect inverse function to check whether a zero difference results.
    Eigen::VectorXd forwardBackardTransformationResults =  Eigen::VectorXd( directTimeDifferences.size( ) );
    Eigen::VectorXd inverseForwardBackardTransformationResults =  Eigen::VectorXd( directTimeDifferences.size( ) );

    int counter = 0;
    double convertedValue = 0.0;
    for( std::map< double, double >::iterator differenceIterator = directTimeDifferences.begin( ); differenceIterator !=
         directTimeDifferences.end( ); differenceIterator++ )
    {
        convertedValue = indirectDifferenceFunction( differenceIterator->first );
        forwardBackardTransformationResults( counter ) = convertedValue - differenceIterator->second;
        counter++;
    }

    double maximumDifference = forwardBackardTransformationResults.maxCoeff( );
    double minimumDifference = forwardBackardTransformationResults.minCoeff( );

    BOOST_CHECK_SMALL( maximumDifference, 1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( minimumDifference ), 1.0E-9 );

    counter = 0;
    convertedValue = 0.0;
    forwardBackardTransformationResults.setZero( );
    for( std::map< double, double >::iterator differenceIterator = directInverseTimeDifferences.begin( ); differenceIterator !=
         directInverseTimeDifferences.end( ); differenceIterator++ )
    {
        convertedValue = indirectInverseDifferenceFunction( differenceIterator->first );
        forwardBackardTransformationResults( counter ) = convertedValue-differenceIterator->second;
        counter++;
    }

    maximumDifference = forwardBackardTransformationResults.maxCoeff( );
    minimumDifference = forwardBackardTransformationResults.minCoeff( );

    BOOST_CHECK_SMALL( maximumDifference, 1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( minimumDifference ), 1.0E-9 );

    std::vector< double > evaluationTimes = utilities::createVectorFromMapKeys( directTimeDifferences );

    std::function< double( const double ) > differenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                body_centered_coordinate_time_scale, local_proper_time_scale, "Graz" );
    std::function< double( const double ) > inverseDifferenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                local_proper_time_scale, body_centered_coordinate_time_scale, "Graz" );

    double convertedTime;
    forwardBackardTransformationResults.setZero( );
    for( unsigned int i = 0; i < evaluationTimes.size( ); i++ )
    {
        convertedTime = evaluationTimes[ i ] + differenceFunction( evaluationTimes[ i ] );
        forwardBackardTransformationResults( i ) = evaluationTimes[ i ] - ( convertedTime + inverseDifferenceFunction( convertedTime ) );
    }

    maximumDifference = forwardBackardTransformationResults.maxCoeff( );
    minimumDifference = forwardBackardTransformationResults.minCoeff( );

    BOOST_CHECK_SMALL( maximumDifference, 1.0E-9 );
    BOOST_CHECK_SMALL( std::fabs( minimumDifference ), 1.0E-9 );

    differenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                barycentric_coordinate_time_scale, local_proper_time_scale, "Graz" );
    inverseDifferenceFunction = earthTimeScaleConverter->getTimeDifferenceFunction(
                local_proper_time_scale, barycentric_coordinate_time_scale, "Graz" );


    forwardBackardTransformationResults.setZero( );
    for( unsigned int i = 0; i < evaluationTimes.size( ); i++ )
    {
        convertedTime = evaluationTimes[ i ] + differenceFunction( evaluationTimes[ i ] );
        forwardBackardTransformationResults( i ) = evaluationTimes[ i ] - ( convertedTime + inverseDifferenceFunction( convertedTime ) );
    }

    maximumDifference = forwardBackardTransformationResults.maxCoeff( );
    minimumDifference = forwardBackardTransformationResults.minCoeff( );
    const double maxAbsDifference = std::max( std::fabs( maximumDifference ), std::fabs( minimumDifference ) );
    std::cout << "[test_concatenated_conversions] max_abs_diff=" << maxAbsDifference << std::endl;

    BOOST_CHECK_SMALL( maximumDifference, 1.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( minimumDifference ), 1.0E-12 );
}

BOOST_AUTO_TEST_CASE( test_geoid_tt_tcg_sh_rotation_rate )
{
    spice_interface::loadStandardSpiceKernels( );

    const double initialEpoch = 0.0;
    const double finalEpoch = initialEpoch + 7.0 * physical_constants::JULIAN_DAY;
    const double integrationStep = 100.0;
    const double ephemerisBuffer = physical_constants::JULIAN_DAY;

    const std::string globalFrameOrigin = "Earth";
    const std::string globalFrameOrientation = "J2000";
    const std::vector< std::string > bodyNames{ "Sun", "Earth", "Moon" };
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodyNames, initialEpoch - ephemerisBuffer, finalEpoch + ephemerisBuffer,
                globalFrameOrigin, globalFrameOrientation );

    // WGS84 mean geoid geometry.
    const double wgs84Flattening = 1.0 / 298.257223563;
    const double wgs84EquatorialRadius = 6378137.0;
    bodySettings.at( "Earth" )->shapeModelSettings =
            std::make_shared< simulation_setup::OblateSphericalBodyShapeSettings >(
                wgs84EquatorialRadius, wgs84Flattening );

    // High-accuracy Earth rotation (GCRS<->ITRS).
    bodySettings.at( "Earth" )->rotationModelSettings =
            std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
                basic_astrodynamics::iau_2006, globalFrameOrientation );

    // Truncate Earth gravity field to degree/order 2 (SH-enabled PN integrand test).
    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > earthGravityFieldSettings =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                bodySettings.at( "Earth" )->gravityFieldSettings );
    BOOST_REQUIRE( earthGravityFieldSettings != nullptr );

    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients.block( 0, 0, 3, 3 ) =
            earthGravityFieldSettings->getCosineCoefficients( ).block( 0, 0, 3, 3 );
    sineCoefficients.block( 0, 0, 3, 3 ) =
            earthGravityFieldSettings->getSineCoefficients( ).block( 0, 0, 3, 3 );

    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                earthGravityFieldSettings->getGravitationalParameter( ),
                earthGravityFieldSettings->getReferenceRadius( ),
                cosineCoefficients,
                sineCoefficients,
                "ITRS",
                earthGravityFieldSettings->getScaledMeanMomentOfInertia( ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodies.getMap( ), globalFrameOrigin, globalFrameOrientation );

    for( const std::string& bodyName : bodyNames )
    {
        bodies.getBody( bodyName )->setStateFromEphemeris( initialEpoch );
    }
    bodies.getBody( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialEpoch );

    const std::string stationName = "Equator45";
    // Geodetic coordinates in Tudat are ellipsoidal; use a local undulation correction so h_geoid = 0
    // is represented on the WGS84 ellipsoid at this site.
    const double geoidHeight = 0.0;
    const double geoidUndulationAtStation = 15.0;
    const double ellipsoidalHeightForMeanGeoid = geoidHeight - geoidUndulationAtStation;
    const Eigen::Vector3d stationGeodeticCoordinates =
            ( Eigen::Vector3d( ) << convertDegreesToRadians( 0.0 ),
              convertDegreesToRadians( 45.0 ), ellipsoidalHeightForMeanGeoid ).finished( );
    createGroundStation(
                bodies.getBody( "Earth" ),
                stationName,
                stationGeodeticCoordinates,
                coordinate_conversions::geodetic_position );

    auto integratorSettings = numerical_integrators::rungeKuttaFixedStepSettings(
                integrationStep, numerical_integrators::rungeKutta87DormandPrince );
    auto terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( finalEpoch );

    const Eigen::VectorXd initialRelativisticState = Eigen::VectorXd::Zero( 1 );
    std::vector< std::shared_ptr< RelativisticTimeStatePropagatorSettings< double, double > > >
            bodyCentricToTopocentricSettings;
    const std::vector< std::string > topocentricPerturbingBodies{ "Sun", "Moon" };
    bodyCentricToTopocentricSettings.push_back(
                std::make_shared< BodycenteredToTopocentricTimePropagatorSettings< double, double > >(
                    std::make_pair( "Earth", stationName ),
                    false,
                    2,
                    true,
                    topocentricPerturbingBodies,
                    initialRelativisticState,
                    initialEpoch,
                    integratorSettings,
                    terminationSettings ) );

    const std::vector< std::string > earthPerturbingBodies{ "Sun", "Moon" };
    std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings<> > > converterSettings;
    converterSettings[ "Earth" ] = std::make_shared< DirectRelativisticTimeConverterSettings<> >(
                std::make_shared< SecondOrderBodyCenteredRelativisticTimeConverterSettings< double, double > >(
                    "Earth", earthPerturbingBodies, initialEpoch, integratorSettings, terminationSettings ),
                integratorSettings,
                bodyCentricToTopocentricSettings );

    setRelativisticTimeConverters( bodies, converterSettings );

    std::shared_ptr< TimeEphemeris > earthTimeScaleConverter = bodies.getBody( "Earth" )->getTimeScaleConverter( );
    BOOST_REQUIRE( earthTimeScaleConverter != nullptr );

    // Sample TCG-local_proper over one week and fit mean slope.
    std::vector< double > times;
    std::vector< double > tcgMinusProper;
    for( double epoch = initialEpoch; epoch <= finalEpoch + std::numeric_limits< double >::epsilon( ); epoch += integrationStep )
    {
        times.push_back( epoch - initialEpoch );
        tcgMinusProper.push_back(
                    earthTimeScaleConverter->getTimeDifference(
                        body_centered_coordinate_time_scale, local_proper_time_scale, epoch, stationName ) );
    }

    Eigen::VectorXd timeVector = utilities::convertStlVectorToEigenVector( times );
    Eigen::VectorXd resultVector = utilities::convertStlVectorToEigenVector( tcgMinusProper );

    const double meanTime = timeVector.mean( );
    const double meanResult = resultVector.mean( );
    const Eigen::VectorXd centeredTimes =
            timeVector - Eigen::VectorXd::Constant( timeVector.rows( ), meanTime );
    const Eigen::VectorXd centeredResults =
            resultVector - Eigen::VectorXd::Constant( resultVector.rows( ), meanResult );

    const double measuredRate = centeredTimes.dot( centeredResults ) / centeredTimes.squaredNorm( );
    const double fittedOffset = meanResult - measuredRate * meanTime;
    BOOST_CHECK_CLOSE_FRACTION( measuredRate, -physical_constants::LG_TIME_RATE_TERM, 1.0E-6 );

    const Eigen::VectorXd detrendedResults = resultVector - (
                Eigen::VectorXd::Constant( resultVector.rows( ), fittedOffset ) +
                measuredRate * timeVector );
    const double maximumResidualAmplitude = detrendedResults.cwiseAbs( ).maxCoeff( );
    BOOST_CHECK_SMALL( maximumResidualAmplitude, 1.0E-12 );
}



BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
