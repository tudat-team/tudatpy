#ifndef ORBITDETERMINATIONTESTCASES_H
#define ORBITDETERMINATIONTESTCASES_H

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

Eigen::VectorXd getDefaultInitialParameterPerturbation( )
{
    Eigen::VectorXd parameterPerturbations = Eigen::VectorXd( 7 );
    for( int i = 0; i < 3; i++ )
    {
        parameterPerturbations( i ) = 1.0E3;
        parameterPerturbations( i + 3 ) = 1.0E-2;
    }
    parameterPerturbations( 6 ) = 5.0E6;

    return parameterPerturbations;
}

template< typename TimeType = double, typename StateScalarType  = double >
std::pair< boost::shared_ptr< PodOutput< StateScalarType > >, Eigen::VectorXd > executeParameterEstimation(
        const int observableType = 1,
        Eigen::VectorXd parameterPerturbation = getDefaultInitialParameterPerturbation( ),
        Eigen::MatrixXd inverseAPrioriCovariance  = Eigen::MatrixXd::Zero( 7, 7 ),
        const double weight = 1.0 )
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
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = TimeType( 3.0E7 );
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrigin( "Sun" );

    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Saturn" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationMap[ "Earth" ] = accelerationsOfEarth;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToEstimate;
    bodiesToEstimate.push_back( "Earth" );
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;

    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, centralBodyMap );

    // Set parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                  "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                                      "Earth", centralBodyMap[ "Earth" ], bodyMap, initialEphemerisTime ),
                              centralBodyMap[ "Earth" ] ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define integrator settings.
    boost::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            boost::make_shared< IntegratorSettings< TimeType > >(
                rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 900.0 );


    boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate,
              getInitialStateVectorOfBodiesToEstimate( parametersToEstimate ),
              TimeType( finalEphemerisTime + 4.0 * maximumTimeStep ),
              cowell, boost::shared_ptr< DependentVariableSaveSettings >( ), 600.0 );

    // Define light-time corrections.
    std::map< ObservableType, std::map< LinkEnds, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > > >
            observableCorrections;

    // Define link ends
    LinkEnds linkEnds;
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsMap;

    observation_models::ObservationSettingsMap observationSettingsMap;

    if(observableType == 0 )
    {
        linkEnds[ observed_body ] = std::make_pair( "Earth", "" );
        observationSettingsMap[ position_observable ][ linkEnds ] = boost::make_shared< ObservationSettings >(
                    position_observable );
    }
    else
    {
        linkEnds[ transmitter ] = std::make_pair( "Earth", "" );
        linkEnds[ receiver ] = std::make_pair( "Mars", "" );

        if( observableType == 1 )
        {
            observationSettingsMap[ one_way_range ][ linkEnds ] = boost::make_shared< ObservationSettings >(
                        one_way_range );
        }
        else if( observableType == 2 )
        {
            observationSettingsMap[ angular_position ][ linkEnds ] = boost::make_shared< ObservationSettings >(
                        angular_position );
        }
        else if( observableType == 3 )
        {
            observationSettingsMap[ one_way_doppler ][ linkEnds ] = boost::make_shared< ObservationSettings >(
                        one_way_doppler );
        }
    }



    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );


    // Define observation times.
    double observationTimeStep = 1000.0;
    TimeType observationTime = Time( initialEphemerisTime + 10.0E4 );
    int numberOfObservations = 18000;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservations );

    for( int i = 0; i < numberOfObservations; i++ )
    {
        initialObservationTimes[ i ] = observationTime;
        observationTime += observationTimeStep;
    }

    // Define observation simulation settings
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;
    initialObservationTimes = utilities::addScalarToVector( initialObservationTimes, 30.0 );
    if( observableType == 0 )
    {
        singleObservableSimulationInput[ linkEnds ] = std::make_pair( initialObservationTimes, observed_body );
        measurementSimulationInput[ position_observable ] = singleObservableSimulationInput;
    }
    else
    {
        singleObservableSimulationInput[ linkEnds ] = std::make_pair( initialObservationTimes, transmitter );

        if( observableType == 1 )
        {
            measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;
        }
        else if( observableType == 2 )
        {
            measurementSimulationInput[ angular_position ] = singleObservableSimulationInput;
        }
        else if( observableType == 3 )
        {
            measurementSimulationInput[ one_way_doppler ] = singleObservableSimulationInput;
        }
    }

    singleObservableSimulationInput.clear( );

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationManagers( ) );

    // Perturb parameter estimate
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    if( parameterPerturbation.rows( ) == 0 )
    {
        parameterPerturbation = Eigen::VectorXd::Zero( 7 );
    }
    for( unsigned int i = 0; i < 7; i++ )
    {
        initialParameterEstimate( i ) += parameterPerturbation( i );
    }

    // Define estimation input
    boost::shared_ptr< PodInput< StateScalarType, TimeType > > podInput =
            boost::make_shared< PodInput< StateScalarType, TimeType > >(
                observationsAndTimes, initialParameterEstimate.rows( ), inverseAPrioriCovariance,
                initialParameterEstimate - truthParameters );
    podInput->setConstantWeightsMatrix( weight );

    // Perform estimation
    boost::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, boost::make_shared< EstimationConvergenceChecker >( ), true, true, false, false );

    return std::make_pair( podOutput,
                           ( podOutput->parameterEstimate_.template cast< double >( ) -
                             truthParameters .template cast< double >( ) ) );
}

}

}
#endif // ORBITDETERMINATIONTESTCASES_H
