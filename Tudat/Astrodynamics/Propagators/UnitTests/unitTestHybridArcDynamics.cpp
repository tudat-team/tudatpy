#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{

namespace unit_tests
{

//! Using declarations.
using namespace tudat;
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace unit_conversions;

BOOST_AUTO_TEST_SUITE( test_hybrid_arc_dynamics )

//! Reset Mats ephemeris to tabulated ephemeris from default settings
void resetMarsEphemeris(
        const NamedBodyMap& bodyMap, const double ephemerisStartTime, const double ephemerisEndTime )
{
    std::shared_ptr< Ephemeris > marsEphemeris =
            createBodyEphemeris( getDefaultEphemerisSettings( "Mars", ephemerisStartTime, ephemerisEndTime ), "Mars" );

    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > defaultMarsStateInterpolator =
            std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( marsEphemeris )->getInterpolator( );

    std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
                bodyMap.at( "Mars" )->getEphemeris( ) )->resetInterpolator( defaultMarsStateInterpolator );
}

//! Test if hybrid-arc orbit propagation is done correctly (single-arc Mars w.r.t. SSB and multi-arc orbiter w.r.t. Mars)
BOOST_AUTO_TEST_CASE( testHybridArcDynamics )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.5E7;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create orbiter
    bodyMap[ "Orbiter" ] = std::make_shared< Body >( );
    bodyMap[ "Orbiter" ]->setConstantBodyMass( 5.0E3 );
    bodyMap[ "Orbiter" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );

    // Create and set radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );
    bodyMap[ "Orbiter" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Set accelerations for Mars
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );

    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );

    // Create single-arc propagator settings for Mars
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodyMap, initialEphemerisTime );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );

    // Set accelerations for Orbiter
    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );

    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );

    // Creater arc times for orbiter
    std::vector< double > integrationArcStarts, integrationArcEnds;
    double integrationStartTime = initialEphemerisTime + 1.0E4;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double arcDuration = 1.0E6;
    double arcOverlap = 1.0E4;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcEnds.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    // Create list of multi-arc Orbiter initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define (quasi-arbitrary) arc Orbiter initial states
    double marsGravitationalParameter =  bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
    {
        Eigen::Vector6d orbiterInitialStateInKeplerianElements;
        orbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 4000.0E3;
        orbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
        orbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        orbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 - j );
        orbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 + j );
        orbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 + j * 10.0 );

        // Convert state from Keplerian elements to Cartesian elements.
        multiArcSystemInitialStates[ j ]  = convertKeplerianToCartesianElements(
                    orbiterInitialStateInKeplerianElements,
                    marsGravitationalParameter );;
    }

    // Create multi-arc propagator settings for Orbiter
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( multiArcCentralBodies, multiArcAccelerationModelMap, multiArcBodiesToIntegrate,
                      multiArcSystemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
    }
    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList );

    // Perform separate single-arc propagation
    std::map< double, Eigen::VectorXd > singleArcSolution;
    {
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, initialEphemerisTime, 120.0 );
        SingleArcDynamicsSimulator< > singleArcDynamicsSimulator(
                    bodyMap, integratorSettings, singleArcPropagatorSettings, true, false, true );
        singleArcSolution = singleArcDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    }

    // Perform separate multi-arc propagation
    std::vector< std::map< double, Eigen::VectorXd > > multiArcSolution;
    {
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, initialEphemerisTime, 120.0 );
        MultiArcDynamicsSimulator< > multiArcDynamicsSimulator(
                    bodyMap, integratorSettings, multiArcPropagatorSettings, integrationArcStarts, true, false, false );
        multiArcSolution = multiArcDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    }

    // Perform hybrid arc propagation
    resetMarsEphemeris( bodyMap, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    std::map< double, Eigen::VectorXd > singleArcSolutionFromHybrid;
    std::vector< std::map< double, Eigen::VectorXd > > multiArcSolutionFromHybrid;
    {
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, initialEphemerisTime, 120.0 );
        HybridArcDynamicsSimulator< > hybridArcDynamicsSimulator(
                    bodyMap, integratorSettings, std::make_shared< HybridArcPropagatorSettings< > >(
                        singleArcPropagatorSettings, multiArcPropagatorSettings ), integrationArcStarts );

        singleArcSolutionFromHybrid = hybridArcDynamicsSimulator.getSingleArcDynamicsSimulator( )->
                getEquationsOfMotionNumericalSolution( );
        multiArcSolutionFromHybrid = hybridArcDynamicsSimulator.getMultiArcDynamicsSimulator( )->
                getEquationsOfMotionNumericalSolution( );
    }

    // Compare separate propagation results with hybrid arc
    BOOST_CHECK_EQUAL( singleArcSolutionFromHybrid.size( ), singleArcSolution.size( ) );

    std::map< double, Eigen::VectorXd >::const_iterator singleArcSolutionFromHybridIterator =
            singleArcSolutionFromHybrid.begin( );
    std::map< double, Eigen::VectorXd >::const_iterator singleArcSolutionIterator =
            singleArcSolution.begin( );

    for( unsigned int i = 0; i < singleArcSolution.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION(
                    singleArcSolutionFromHybridIterator->first, singleArcSolutionIterator->first,
                    4.0 * std::numeric_limits< double >::epsilon( ) );
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs(
                            singleArcSolutionFromHybridIterator->second( j ) -
                            singleArcSolutionIterator->second( j ) ), 1.0E-5 );
            BOOST_CHECK_SMALL(
                        std::fabs(
                            singleArcSolutionFromHybridIterator->second( j + 3 ) -
                            singleArcSolutionIterator->second( j + 3  ) ), 1.0E-11 );
        }
        singleArcSolutionFromHybridIterator++;
        singleArcSolutionIterator++;
    }

    BOOST_CHECK_EQUAL( multiArcSolutionFromHybrid.size( ), multiArcSolution.size( ) );
    for( unsigned int arc = 0; arc < multiArcSolutionFromHybrid.size( ); arc++ )
    {
        BOOST_CHECK_EQUAL( multiArcSolutionFromHybrid.at( arc ).size( ), multiArcSolution.at( arc ).size( ) );


        std::map< double, Eigen::VectorXd >::const_iterator currentArcSolutionFromHybridIterator =
                multiArcSolutionFromHybrid.at( arc ).begin( );
        std::map< double, Eigen::VectorXd >::const_iterator currentArcSolutionIterator =
                multiArcSolution.at( arc ).begin( );

        for( unsigned int i = 0; i < multiArcSolution.at( arc ).size( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                        currentArcSolutionFromHybridIterator->first, currentArcSolutionIterator->first,
                        4.0 * std::numeric_limits< double >::epsilon( ) );
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs(
                                currentArcSolutionFromHybridIterator->second( j ) -
                                currentArcSolutionIterator->second( j ) ), 1.0E-5 );
                BOOST_CHECK_SMALL(
                            std::fabs(
                                currentArcSolutionFromHybridIterator->second( j + 3 ) -
                                currentArcSolutionIterator->second( j + 3  ) ), 1.0E-11 );
            }
            currentArcSolutionFromHybridIterator++;
            currentArcSolutionIterator++;
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
