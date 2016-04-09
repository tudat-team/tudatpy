#define BOOST_TEST_MAIN

#include <string>
#include <thread>
#include <omp.h>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Propagators/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace unit_tests
{

//Using declarations.
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( testVariationalEquationCalculation )

template< typename TimeType = double , typename StateScalarType  = double >
        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
execute( const std::vector< std::string > centralBodies,
         const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference = Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ) )
{

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = initialEphemerisTime + 0.5E7;
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );


    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSun;
    //accelerationsOfSun[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    //accelerationsOfSun[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    //accelerationMap[ "Sun" ] = accelerationsOfSun;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    bodiesToIntegrate.push_back( "Earth" );
    //bodiesToIntegrate.push_back( "Sun" );

    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::map< std::string, std::string > centralBodyMap;

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, centralBodyMap );


    boost::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            boost::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime ),
              TimeType( finalEphemerisTime ), 300.0 );


    // Set initial states of bodies to integrate.
    TimeType initialIntegrationTime = initialEphemerisTime;

    boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings;

    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState;
        systemInitialState = getInitialStatesOfBodies< TimeType, StateScalarType >(
                    bodiesToIntegrate, centralBodies, bodyMap, initialIntegrationTime );
        systemInitialState += initialStateDifference;

        propagatorSettings =  boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState );
    }

    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >
                              ( "Sun", gravitational_parameter ) );
    {
        parameterNames.push_back( boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                      "Moon", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                                          "Moon", centralBodies[ 0 ], bodyMap, TimeType( initialEphemerisTime ) ) + initialStateDifference.segment( 0, 6 ),
                                  centralBodies[ 0 ] ) );
        parameterNames.push_back( boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                      "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                                          "Earth", centralBodies[ 1 ], bodyMap, TimeType( initialEphemerisTime ) ) + initialStateDifference.segment( 6, 6 ),
                                      centralBodies[ 1 ] ) );
    }

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

    {
        SingleArcVariationalEquationsSolver< StateScalarType, TimeType, double > dynamicsSimulator =
                SingleArcVariationalEquationsSolver< StateScalarType, TimeType, double >(
                    bodyMap, integratorSettings, propagatorSettings,
                    parametersToEstimate );

        double testEpoch = 1.4E7;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
        testStates.block( 0, 0, 6, 1 ) = bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( testEpoch );
        if( centralBodyMap[ "Moon" ] == "Earth" )
        {
            testStates.block( 0, 0, 6, 1 ) -= bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( testEpoch );
        }

        testStates.block( 6, 0, 6, 1 ) = bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( testEpoch );
        //testStates.block( 12, 0, 6, 1 ) = bodyMap[ "Sun" ]->getStateInBaseFrameFromEphemeris( testEpoch );

        results.first.push_back( dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinesStateTransitionAndSensitivityMatrix( testEpoch ) );
        results.second.push_back( testStates );
    }
    return results;
}

BOOST_AUTO_TEST_CASE( test_variational_equation_calculation )
{
    std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

    std::vector< std::vector< std::string > > centralBodiesSet;
    std::vector< std::string > centralBodies;
    centralBodies.resize( 2 );

    centralBodies[ 0 ] = "SSB";
    centralBodies[ 1 ] = "SSB";
    //    //centralBodies[ 2 ] = "SSB";
    centralBodiesSet.push_back( centralBodies );

    centralBodies[ 0 ] = "Earth";
    centralBodies[ 1 ] = "Sun";
    //centralBodies[ 2 ] = "SSB";
    centralBodiesSet.push_back( centralBodies );


    Eigen::Matrix< double, 12, 1>  perturbedState;

    Eigen::Matrix< double, 12, 1> statePerturbation;


    std::vector< Eigen::MatrixXd > manualMoonInitialStatePartials;

    for( unsigned int i = 0; i < centralBodiesSet.size( ); i++ )
    {
        currentOutput = execute< double, double >( centralBodiesSet[ i ] );
        Eigen::MatrixXd stateTransitionMatrixAtEpoch = currentOutput.first.at( 0 );
        Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 12, 12 );
        if( i == 0 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( )<< 100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1, 100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1 ).finished( );

        }
        else if( i == 1 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( )<< 100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1, 100000.0, 100000.0, 10000000.0, 0.1, 0.1, 10.0 ).finished( );
        }

        for( unsigned int j = 0; j < 12; j++ )
        {
            Eigen::VectorXd upPerturbedState, downPerturbedState;
            perturbedState.setZero( );
            perturbedState( j ) += statePerturbation( j );
            upPerturbedState = execute< double, double >( centralBodiesSet[ i ], perturbedState ).second.at( 0 );

            perturbedState.setZero( );
            perturbedState( j ) -= statePerturbation( j );
            downPerturbedState = execute< double, double >( centralBodiesSet[ i ], perturbedState ).second.at( 0 );
            manualPartial.block( 0, j, 12, 1 ) = ( upPerturbedState - downPerturbedState ) / ( 2.0 * statePerturbation( j ) );
        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateTransitionMatrixAtEpoch.block( 0, 0, 12, 12 ), manualPartial, 1.0E-3 );

        std::cout<<(stateTransitionMatrixAtEpoch.block( 0, 0, 12, 12 ) - manualPartial ).cwiseQuotient( manualPartial )<<std::endl;

    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

