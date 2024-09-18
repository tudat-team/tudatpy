/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/estimation.h>

#include <tudat/io/applicationOutput.h>

//! Execute propagation of orbits of Vehicle and Obelix around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
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
    using namespace tudat::coordinate_conversions;
    using namespace tudat::ground_stations;
    using namespace tudat::observation_models;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );
    double finalEphemerisTime = initialEphemerisTime + 3.0 * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0,finalEphemerisTime + 3600.0 );
    setSimpleRotationSettingsFromSpice( bodySettings, "Earth", initialEphemerisTime );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    // Create and add aerodynamic coefficient interface
    double referenceArea = 4.0;
    double aerodynamicDragCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            constantAerodynamicCoefficientSettings(
                referenceArea, ( Eigen::Vector3d( ) << aerodynamicDragCoefficient, -0.01, 0.1 ).finished( ), 1, 1 );
    addAerodynamicCoefficientInterface(
                bodies, "Vehicle", aerodynamicCoefficientSettings );

    // Create and add radiation pressure interace
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Earth" };
    std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
            cannonBallRadiationPressureSettings(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );
    addRadiationPressureInterface(
                bodies, "Vehicle", vehicleRadiationPressureSettings );
    addEmptyTabulateEphemeris( bodies, "Vehicle", "Earth" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS               //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -1.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.8, 4.0 ).finished( ), geodetic_position );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationSettingsList;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Earth" ] = {
            sphericalHarmonicAcceleration( 32, 32 ),
            aerodynamicAcceleration( ) };
    accelerationsOfVehicle[ "Sun" ] = {
            pointMassGravityAcceleration( ),
            cannonBallRadiationPressureAcceleration( ) };
    accelerationsOfVehicle[ "Mars" ] = {
            pointMassGravityAcceleration( ) };
    accelerationsOfVehicle[ "Moon" ] = {
            pointMassGravityAcceleration( ) };
    accelerationSettingsList[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Earth" };

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationSettingsList, bodiesToIntegrate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Set initial state.
    double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ) );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                double( initialEphemerisTime ), 40.0,
                rungeKuttaFehlberg78,
                40.0, 40.0, 1.0, 1.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE OBSERVATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define settings for observable, no light-time corrections, and biases for selected 1-way range links
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           std::make_shared< ObservationSettings >( currentObservable ) ) );
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + 1000.0;

    // Define time between two observations
    double  observationInterval = 20.0;

    // Simulate observations for 3 days
    std::vector< double > baseTimeList;
    for( unsigned int i = 0; i < 3; i++ )
    {
        // Simulate 500 observations per day (observationInterval apart)
        for( unsigned int j = 0; j < 500; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    // Create measureement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        // Define observable type and link ends
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define observation times and reference link ends
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_pair( baseTimeList, receiver );
        }
    }

    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 10.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.0E-2 );
    parameterPerturbation( 6 ) = 0.01;
    parameterPerturbation( 7 ) = 0.01;
    initialParameterEstimate += parameterPerturbation;

    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, true, false, true );

    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 4 ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "EarthOrbiterStateEstimationExample/";

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;

    std::map< double, Eigen::VectorXd > propagatedErrors;

    propagateFormalErrors(
                propagatedErrors, podOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                60.0, initialEphemerisTime + 3600.0, finalEphemerisTime - 3600.0 );
    input_output::writeDataMapToTextFile( propagatedErrors,
                                          "earthOrbitBasicEstimationPropagatedErrors.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ),
                                     "earthOrbitBasicEstimationCovariance.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "earthOrbitBasicEstimationInformationMatrix.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "earthOrbitBasicEstimationWeightsDiagonal.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "earthOrbitBasicEstimationResiduals.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "earthOrbitBasicEstimationCorrelations.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "earthOrbitBasicResidualHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "earthOrbitBasicParameterHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "earthOrbitBasicObservationMeasurements.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "earthOrbitBasicObservationTimes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "earthOrbitBasicObservationLinkEnds.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "earthOrbitBasicObservationObservableTypes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( estimationError,
                                     "earthOrbitBasicObservationTrueEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "earthOrbitBasicObservationFormalEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
