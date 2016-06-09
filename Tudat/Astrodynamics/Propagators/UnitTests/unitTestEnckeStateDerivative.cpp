#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <thread>
#include <omp.h>

#include <boost/make_shared.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_encke_propagator )

BOOST_AUTO_TEST_CASE( testEnckePopagator )
{

    for( unsigned int simulationCase = 0; simulationCase < 2; simulationCase++ )
    {
        //Using declarations.
        using namespace tudat::interpolators;
        using namespace tudat::numerical_integrators;
        using namespace tudat::spice_interface;
        using namespace tudat::simulation_setup;
        using namespace tudat::basic_astrodynamics;
        using namespace tudat::orbital_element_conversions;
        using namespace tudat::propagators;


        //Load spice kernels.
        std::string kernelsPath = input_output::getSpiceKernelPath( );
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

        unsigned int totalNumberOfBodies = 7;
        std::vector< std::string > bodyNames;
        bodyNames.resize( totalNumberOfBodies );
        bodyNames[ 0 ] = "Earth";
        bodyNames[ 1 ] = "Mars";
        bodyNames[ 2 ] = "Sun";
        bodyNames[ 3 ] = "Venus";
        bodyNames[ 4 ] = "Moon";
        bodyNames[ 5 ] = "Mercury";
        bodyNames[ 6 ] = "Jupiter";

        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = 2.0E7;
        double maximumTimeStep = 3600.0;
        double buffer = 5.0 * maximumTimeStep;

        // Create bodies needed in simulation
        NamedBodyMap bodyMap = createBodies(
                    getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer ) );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Earth" ] = accelerationsOfEarth;

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
        accelerationsOfMars[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Mars" ] = accelerationsOfMars;

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
        accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Moon" ] = accelerationsOfMoon;

        std::vector< std::string > bodiesToIntegrate;
        bodiesToIntegrate.push_back( "Earth" );
        bodiesToIntegrate.push_back( "Mars" );
        bodiesToIntegrate.push_back( "Moon" );

        unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

        std::vector< std::string > centralBodies;
        std::map< std::string, std::string > centralBodyMap;
        centralBodies.resize( numberOfNumericalBodies );

        for( int i = 0; i < 3; i++ )
        {
            if( i == 2 && simulationCase == 1 )
            {
                centralBodies[ i ] = "Earth";
            }
            else
            {
                centralBodies[ i ] = "Sun";
            }
            centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
        }



        Eigen::VectorXd systemInitialState = Eigen::VectorXd( bodiesToIntegrate.size( ) * 6 );
        for( unsigned int i = 0; i < numberOfNumericalBodies ; i++ )
        {
            systemInitialState.segment( i * 6 , 6 ) =
                    bodyMap[ bodiesToIntegrate[ i ] ]->getStateInBaseFrameFromEphemeris( initialEphemerisTime ) -
                    bodyMap[ centralBodies[ i ] ]->getStateInBaseFrameFromEphemeris( initialEphemerisTime );
        }

        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, centralBodyMap );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4,
                  initialEphemerisTime, finalEphemerisTime, 250.0 );

        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState );

        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
                    bodyMap, integratorSettings, propagatorSettings, true );

        double initialTestTime = initialEphemerisTime + 10.0 * maximumTimeStep;
        double finalTestTime = finalEphemerisTime - 10.0 * maximumTimeStep;
        double testTimeStep = 1.0E4;

        double currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 18, 1 > > enckeIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            enckeIntegrationResults[ currentTestTime ].segment( 0, 6 ) = bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            enckeIntegrationResults[ currentTestTime ].segment( 6, 6 ) = bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            enckeIntegrationResults[ currentTestTime ].segment( 12, 6 ) = bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( currentTestTime );

            currentTestTime += testTimeStep;
        }


        propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, encke );

        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true );

        currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 18, 1 > > cowellIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            cowellIntegrationResults[ currentTestTime ].segment( 0, 6 ) = bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            cowellIntegrationResults[ currentTestTime ].segment( 6, 6 ) = bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            cowellIntegrationResults[ currentTestTime ].segment( 12, 6 ) = bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            currentTestTime += testTimeStep;
        }

        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator cowellIterator = cowellIntegrationResults.begin( );
        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator enckeIterator = enckeIntegrationResults.begin( );

        for( unsigned int i = 0; i < cowellIntegrationResults.size( ); i++ )
        {
            for( int j= 0; j< 3; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 0.01 );
            }

            for( int j = 6; j < 9; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 0.01 );
            }

            for( int j = 12; j < 15; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 0.075 );
            }

            for( int j = 3; j < 6; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 1.0E-8 );

            }

            for( int j = 9; j < 12; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 1.0E-8 );

            }

            for( int j = 15; j < 18; j++ )
            {
                BOOST_CHECK_SMALL( ( cowellIterator->second - enckeIterator->second ).segment( j, 1 )( 0 ), 1.0E-6 );

            }
            cowellIterator++;
            enckeIterator++;
        }
    }
}
BOOST_AUTO_TEST_SUITE_END( )


}

}


