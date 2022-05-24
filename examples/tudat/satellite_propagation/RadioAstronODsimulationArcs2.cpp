/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"


//! Execute propagation of orbits of Asterix and Obelix around the Earth.
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
    using namespace tudat::input_output;
    using namespace tudat::sofa_interface;
    using namespace tudat::statistics;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND RA       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0012.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter");



    // Reading Mikhail's data related to March 2016
    Eigen::MatrixXd MikhailData = readMatrixFromFile("/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/MikhailData.txt", " ", "");

    int rows = MikhailData.rows();
    int columns = MikhailData.cols();
    //rows = 212;

    std::cout<<"rows: "<<rows<<std::endl;
    std::cout<<"columns: "<<columns<<std::endl;

    std::cout<<"Mikhail data before: "<< MikhailData.block(0,5,1,6)<<std::endl;

    //Converting the units to the standard ones
    MikhailData.block( 0, 5, rows, 3 ) = MikhailData.block( 0, 5, rows, 3 ) * 1.0E6;
    MikhailData.block( 0, 8, rows, 3 ) = MikhailData.block( 0, 8, rows, 3 ) * 1.0E3;

    for ( unsigned int i = 0; i < rows; i ++ )
    {

        MikhailData.block( i, 11, 1, 3 ) = MikhailData.block( i, 11, 1, 3 ) * 1.0E12;
        MikhailData.block( i, 14, 1, 3 ) = MikhailData.block( i, 14, 1, 3 ) * 1.0E9;
        MikhailData.block( i, 17, 1, 2 ) = MikhailData.block( i, 17, 1, 2 ) * 1.0E12;
        MikhailData.block( i, 19, 1, 3 ) = MikhailData.block( i, 19, 1, 3 ) * 1.0E9;
        MikhailData.block( i, 22, 1, 1 ) = MikhailData.block( i, 22, 1, 1 ) * 1.0E12;
        MikhailData.block( i, 23, 1, 3 ) = MikhailData.block( i, 23, 1, 3 ) * 1.0E9;
        MikhailData.block( i, 26, 1, 6 ) = MikhailData.block( i, 26, 1, 6 ) * 1.0E6;


    }

    std::cout<<"Mikhail data after: "<< MikhailData.block(0,5,1,6)<<std::endl;

    int row1 = 0;
    int row_final = 696;
    
    

    

    int year = MikhailData(row1,0);
    int months = MikhailData(row1,1);
    int days = MikhailData(row1,2);
    int hours = MikhailData(row1,3);
    int minutes = MikhailData(row1,4);
    double seconds = 0;

    double InitialutcSecondsSinceJ2000 =
            convertCalendarDateToJulianDaysSinceEpoch< double >(
                year, months, days, hours, minutes, seconds,
                getJulianDayOnJ2000< double >( ) ) *
            physical_constants::getJulianDay< double >( );

    year = MikhailData(row_final - 1, 0);
    months = MikhailData(row_final - 1, 1);
    days = MikhailData(row_final - 1, 2);
    hours = MikhailData(row_final - 1, 3);
    minutes = MikhailData(row_final - 1, 4);
    seconds = 0;

    double FinalutcSecondsSinceJ2000 =
            convertCalendarDateToJulianDaysSinceEpoch< double >(
                year, months, days, hours, minutes, seconds,
                getJulianDayOnJ2000< double >( ) ) *
            physical_constants::getJulianDay< double >( );



    // Specify initial time
    double initialEphemerisTime = convertUTCtoTT( InitialutcSecondsSinceJ2000 );
    double finalEphemerisTime = convertUTCtoTT( FinalutcSecondsSinceJ2000 );
//    double finalEphemerisTime = initialEphemerisTime + 72*3600;

    std::cout<< "InitialEphemerisTime: "<< initialEphemerisTime<<std::endl;

    //sleep(10);

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 30000.0, finalEphemerisTime + 30000.0 );


    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        bodySettings[ bodyNames.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodyNames.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }


    SystemOfBodies bodies = createBodies( bodySettings );
    bodies[ "RA" ] = std::make_shared< Body >( );
    bodies[ "RA" ]->setConstantBodyMass( 3600.0 );

    // Create radiation pressure settings
    double referenceAreaRadiation = 100.0;
    double radiationPressureCoefficient = 1.2 + 0.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > RARadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies[ "RA" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    RARadiationPressureSettings, "RA", bodies ) );

    bodies[ "RA" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                       std::map< double, std::shared_ptr< Ephemeris > >( ), "Earth", "J2000" ) );

//    bodies[ "RA" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
//                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                            < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );


    setGlobalFrameBodyEphemerides( bodies, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS               //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1M" );
//    groundStationNames.push_back( "Station2Y" );
//    groundStationNames.push_back( "Station3H" );

    createGroundStation( bodies.at( "Earth" ), "Station1M",
                         ( Eigen::Vector3d( ) << 2892607.149, 1311813.079,  5512598.659 ).finished( ), cartesian_position );
//    createGroundStation( bodies.at( "Earth" ), "Station2Y",
//                         ( Eigen::Vector3d( ) << -2389008, 5043332,  -3078526 ).finished( ), cartesian_position );
//    createGroundStation( bodies.at( "Earth" ), "Station3H",
//                         ( Eigen::Vector3d( ) << 5085401.135, 2668330.108, -2768688.865 ).finished( ), cartesian_position );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on RA that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfRA;
    accelerationsOfRA[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfRA[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfRA[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
    accelerationsOfRA[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
    accelerationsOfRA[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfRA[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfRA[ "Earth" ].push_back( std::make_shared< EmpiricalAccelerationSettings >( ) );

    accelerationMap[ "RA" ] = accelerationsOfRA;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "RA" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Set Keplerian elements for RadioAstron.
//    Eigen::Vector6d asterixInitialStateInKeplerianElements;
//    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
//    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
//    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
//    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
//            = unit_conversions::convertDegreesToRadians( 235.7 );
//    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
//            = unit_conversions::convertDegreesToRadians( 23.4 );
//    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState = MikhailData.block( row1, 5, 1, 6 ).transpose();

    std::cout<< "Initial state:"<<systemInitialState.transpose()<<std::endl;

    Eigen::Matrix< double, 6, 1 > RAInitialKeplerianState = convertCartesianToKeplerianElements(
                systemInitialState, earthGravitationalParameter );


    std::cout<< "Initial state Keplerian:"<<RAInitialKeplerianState.transpose()<<std::endl;

    double orbitalPeriod = 2*mathematical_constants::PI * sqrt(pow(RAInitialKeplerianState[0],3)/earthGravitationalParameter);

    std::cout<< "Orbital Period (days): "<<orbitalPeriod/86400<<std::endl;


    // Defining Arcs
    std::vector<double> ArcInitialTimes;
    int row2, row3;
    row2 = 174;
    row3 = 348;
    std::vector<int> Rows = {row1, row2, 200, 220, 280};//  310, row3 //400, 420, 459, 480, 510, 550, 580, 610, 615};

    double arcDuration = 3600;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings < double > > > propagatorSettingsList;
    double currentTime;

    for ( unsigned int i = 0; i < Rows.size(); i ++ )
    {
        ArcInitialTimes.push_back( initialEphemerisTime + Rows[i] * 3600 );
    }

    Eigen::VectorXd SystemInitialState = Eigen::VectorXd( 6 * ArcInitialTimes.size( ) );


    for ( unsigned int i = 0; i < ArcInitialTimes.size(); i ++ )
    {

        currentTime = ArcInitialTimes[i];

        Eigen::Vector6d currentArcInitialState = MikhailData.block( Rows[i], 5, 1, 6 ).transpose();

        std::cout<<"Distance arc "<< i + 1<<": "<< sqrt( pow(currentArcInitialState(0),2) + pow(currentArcInitialState(1),2) + pow(currentArcInitialState(2),2))<<std::endl;

        SystemInitialState.segment( i * 6, 6 ) = MikhailData.block( Rows[i], 5, 1, 6 ).transpose();

        propagatorSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate, currentArcInitialState,
                      currentTime + arcDuration + 3600, cowell ) );


    }


//    // Create propagator settings
//    std::shared_ptr< PropagatorSettings< double > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
//            ( centralBodies, accelerationModelMap, bodiesToIntegrate, MikhailData.block( row1, 5, 1, 6 ).transpose( ),
//              initialEphemerisTime + arcDuration + 3600, cowell );// propagatorSettingsList.at(0) ;

    // Create propagator settings
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );



    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >
            ( rungeKuttaVariableStepSize, double( initialEphemerisTime ), 40.0,
              rungeKuttaFehlberg78,
              0.00001, 1.0E2, 1.0E-13, 1.0E-13);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > OneWayLinkEnds;
    std::vector< LinkEnds > TwoWayLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ reflector1 ] = std::make_pair( "RA", "" );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        TwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ transmitter ] = std::make_pair( "RA", "" );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        OneWayLinkEnds.push_back( linkEnds );
    }


    std::cout<<"check"<<std::endl;

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_doppler ].push_back( OneWayLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_doppler ].push_back( OneWayLinkEnds[ 1 ] );
//    linkEndsPerObservable[ one_way_doppler ].push_back( OneWayLinkEnds[ 2 ] );

    linkEndsPerObservable[ two_way_doppler ].push_back( TwoWayLinkEnds[ 0 ] );
//    linkEndsPerObservable[ two_way_doppler ].push_back( TwoWayLinkEnds[ 1 ] );
//    linkEndsPerObservable[ two_way_doppler ].push_back( TwoWayLinkEnds[ 2 ] );


    std::cout<<"check2"<<std::endl;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

    // Create concatenated list of arc initial states
//    Eigen::VectorXd SystemInitialState = Eigen::VectorXd( 6 * ArcInitialTimes.size( ) );
//    for( unsigned int i = 0; i < ArcInitialTimes.size( ); i++ )
//    {
//        SystemInitialState.segment( i * 6, 6 ) = propagatorSettingsList.at( i )->getInitialStates( );

//    }



    std::cout<< "System Initial States: "<<SystemInitialState.transpose()<<std::endl;

  //  sleep(10);

    parameterNames.push_back(
                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    "RA", SystemInitialState, ArcInitialTimes, "Earth" ) );
    parameterNames.push_back( std::make_shared< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >( "RA", ArcInitialTimes  ) );
//    parameterNames.push_back(
//                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                    "RA", propagatorSettings->getInitialStates( ), "Earth" ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "RA", radiation_pressure_coefficient ) );

    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", equivalence_principle_lpi_violation_parameter , "" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station1M" ) );

    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                  linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, false ) );

    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                  linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, false ) );



    std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > > empiricalAccelerationComponents;
    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( cosine_empirical );
    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( constant_empirical );
    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( cosine_empirical );
    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( sine_empirical );
    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
    empiricalAccelerationComponents[ radial_empirical_acceleration_component ].push_back( cosine_empirical );
    empiricalAccelerationComponents[ radial_empirical_acceleration_component ].push_back( sine_empirical );
    empiricalAccelerationComponents[ radial_empirical_acceleration_component ].push_back( constant_empirical );


    std::vector< double > empiricalAccelerationArcTimes;
    empiricalAccelerationArcTimes = ArcInitialTimes;

//    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
//                                 "RA", "Earth", empiricalAccelerationComponents, empiricalAccelerationArcTimes ) );



    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies, accelerationModelMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );
    
//    std::cout<<"T1"<<std::endl;
//    // Define covariance matrix
    
//    Eigen::MatrixXd Covariance;

//    Covariance.setZero( 7*Rows.size() + 6 , 7*Rows.size() + 6  );
    
    
//    std::map< int, Eigen::Matrix< double, 6, 6 > > StateCovMap;


//    for ( int i = 0; i < Rows.size(); i ++ )
//    {

//        int a = 11;

//        Eigen::Matrix< double, 6, 6 > StateCovArc;

//        for ( unsigned int k = 0; k < 6; k ++ )
//        {
//            for ( unsigned int j = k; j < 6 ; j ++ )
//            {
//                StateCovArc( k, j ) = MikhailData( Rows[i], a );
//                StateCovArc( j, k ) = StateCovArc( k, j );
//                a = a + 1;

//            }


//        }

//        StateCovMap[ i ] = StateCovArc;

//    }
    
//    std::cout<<"T2"<<std::endl;



//    for ( unsigned int i = 0; i < Rows.size(); i ++ )
//    {
////        for( int j = 0; j < 3; j++ )
////        {
////            Covariance( 6 * i + j, 6 * i + j ) = 1.0E6;
////            Covariance( 6 * i + 3 + j, 6 * i + 3 + j ) = 1.0E-2;
////        }
////        Covariance( 6 * Rows.size( ) + i,  6 * Rows.size( ) + i ) = 100.0;

//        Covariance.block( 6*i , 6*i , 6, 6 ) = StateCovMap[ i ];

//    }

//    std::cout<<"T3"<<std::endl;

//    Eigen::MatrixXd ParamCov;

//    ParamCov.setZero( Rows.size(), Rows.size() );

//    for ( unsigned int i = 0; i < Rows.size(); i ++ )
//    {

//        ParamCov( i, i ) = 100;

//    }

//    Covariance( 6*(Rows.size() ) , 6*(Rows.size() )  ) = 10E-6;

//    Covariance.block( 6*(Rows.size() ) + 1, 6*(Rows.size() ) + 1, Rows.size(), Rows.size()  ) = ParamCov;
//    std::cout<<"T4"<<std::endl;


//    Eigen::MatrixXd GroundStCov;

////    GroundStCov.setZero( 3, 3 );

////    for ( unsigned int i = 0; i < 9; i ++ )
////    {

////        GroundStCov( i, i ) = 100;

////    }

//    std::cout<<"T4"<<std::endl;

////    Covariance.block( 7*(Rows.size() ) + 1 , 7*(Rows.size() ) + 1 , 3, 3  ) = GroundStCov;

////    std::cout<<"T4"<<std::endl;

////    Covariance( 39, 39 ) = 10^-12;
////    Covariance( 40, 40 ) = 10^-12;


//    Eigen::VectorXd diagonal = Covariance.diagonal();

//    //std::cout<<diagonal<<std::endl;

//    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Covariance1 = diagonal.asDiagonal();


//    //std::cout<< Covariance<<std::endl;





    Eigen::MatrixXd InverseAprioriCov;// = (Covariance).inverse();

    Eigen::MatrixXd AprioriCov;// = InverseAprioriCov.inverse();
    
    
    std::cout<<"T5"<<std::endl;

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
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            std::shared_ptr< DopplerProperTimeRateSettings > properTimeRate = std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings > (
                        "Earth" );

             //Define biases for 1-way range observables
            if(  currentObservable == one_way_doppler  )
            {
                // Absolute bias for 1-way range
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Zero( ), false );

                std::shared_ptr< OneWayDopplerObservationSettings > OneWayDopplerSettings = std::make_shared< OneWayDopplerObservationSettings > (
                            std::shared_ptr< LightTimeCorrectionSettings >( ), properTimeRate, properTimeRate, biasSettings );

                observationSettingsMap.insert(
                            std::make_pair( currentLinkEndsList.at( i ), OneWayDopplerSettings ));



            }
            else if( currentObservable == two_way_doppler  )
            {
                // Only relative bias for 1-way doppler
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Zero( ), false );

                std::shared_ptr< OneWayDopplerObservationSettings > uplinkDopplerSettings = std::make_shared< OneWayDopplerObservationSettings >(
                            std::shared_ptr< LightTimeCorrectionSettings >( ), properTimeRate, properTimeRate, biasSettings );

                std::shared_ptr< OneWayDopplerObservationSettings > downlinkDopplerSettings = std::make_shared< OneWayDopplerObservationSettings >(
                            std::shared_ptr< LightTimeCorrectionSettings >( ), properTimeRate, properTimeRate, biasSettings );

                observationSettingsMap.insert(
                            std::make_pair( currentLinkEndsList.at( i ), std::make_shared< TwoWayDopplerObservationSettings >(
                                                uplinkDopplerSettings, downlinkDopplerSettings, biasSettings ) ) );


            }

        }
    }

    //std::cout<<"Parameters B: "<<parametersToEstimate->template getFullParameterValues< double >( ).transpose( )<<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout<<"Initial state: "<<propagatorSettings->getInitialStates( ).transpose( )<<std::endl;

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );

    input_output::writeDataMapToTextFile(
                orbitDeterminationManager.getVariationalEquationsSolver()->getDynamicsSimulatorBase()->getEquationsOfMotionNumericalSolutionBase().at( 0 ),
                "testStateOutput.dat" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::vector< double > observationTimeArcStart;

    // Define time of first observation
    for ( unsigned int i = 0; i < ArcInitialTimes.size(); i ++ )
    {
         observationTimeArcStart.push_back( ArcInitialTimes[i] + 600.0 );

    }


//    // Define time between two observations
    double  observationInterval = 1;

     //Simulate observations for 30 days
    std::vector< double > baseTimeList;
    for( unsigned int i = 0; i < observationTimeArcStart.size(); i++ )
    {
        // Simulate 500 observations per day (observationInterval apart)
        for( unsigned int j = 0; j < ( arcDuration - 600 ); j++ )
        {
            baseTimeList.push_back( observationTimeArcStart[i] + ( double ) j * observationInterval );
        }
    }

    // Create measureement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
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
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( receiver, baseTimeList );
        }
    }

    std::cout<<"Parameters C: "<<parametersToEstimate->template getFullParameterValues< double >( ).transpose( )<<std::endl;



    // Create observation viability settings and calculators
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair( "Earth", "" ), "",
                                                5.0 * mathematical_constants::PI / 180.0 ) );
    PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
                bodies, linkEndsPerObservable, observationViabilitySettings );



    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with reference
    // link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;


    // Define noise levels
    double dopplerNoise = 1.0E-12;

    // Create noise functions per observable
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;
    noiseFunctions[ one_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                         createBoostContinuousRandomVariableGeneratorFunction(
                             normal_boost_distribution,
                             std::assign::list_of( 0.0 )( dopplerNoise ), 0.0 ), _1 );


    noiseFunctions[ two_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                         createBoostContinuousRandomVariableGeneratorFunction(
                             normal_boost_distribution,
                             std::assign::list_of( 0.0 )( dopplerNoise ), 0.0 ), _1 );


    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions, viabilityCalculators );


//    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
//                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );

    std::cout<<"Parameters: "<<initialParameterEstimate.transpose( )<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
//    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
//    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
//    parameterPerturbation( 6 ) = 0.05;
//    parameterPerturbation( 7 ) = 0.05;
    initialParameterEstimate += parameterPerturbation;

//    // Define estimation input
//    std::shared_ptr< PodInput< double, double > > podInput =
//            std::make_shared< PodInput< double, double > >(
//                observationsAndTimes, initialParameterEstimate.rows( ),
//                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
//                initialParameterEstimate - truthParameters );



        std::shared_ptr< PodInput< double, double > > podInput =
                std::make_shared< PodInput< double, double > >(
                    observationsAndTimes, initialParameterEstimate.rows( ),
                    InverseAprioriCov,
                    initialParameterEstimate - truthParameters );



    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_doppler ] = 1.0 / ( dopplerNoise * dopplerNoise );
    weightPerObservable[ two_way_doppler ] = 1.0 / ( dopplerNoise * dopplerNoise );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->defineEstimationSettings( true, false, true, true, true );


    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 0 )); //true, true, false, true );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "RAStateEstimationExampleEEPArcs1sSpacingRelBias/";

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;

    std::cout<<"True estimation error is:   "<<std::endl<<( estimationError ).transpose( )<<std::endl;
    std::cout<<"Formal estimation error is: "<<std::endl<<podOutput->getFormalErrorVector( ).transpose( )<<std::endl;
    std::cout<<"True to form estimation error ratio is: "<<std::endl<<
               ( podOutput->getFormalErrorVector( ).cwiseQuotient( estimationError ) ).transpose( )<<std::endl;

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "RAEstimationInformationMatrix.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "RAEstimationInformationMatrixNormalization.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "RAEstimationWeightsDiagonal.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "RAEstimationResiduals.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "RAEstimationCorrelations.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getUnnormalizedInverseCovarianceMatrix( ),
                                     "RAEstimationInverseCovarianceMatrix.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "RAResidualHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "RAParameterHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "earthOrbitObservationMeasurements.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "RAObservationTimes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "RAObservationLinkEnds.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "RAObservationObservableTypes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "RAObservationMeasurements.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( estimationError,
                                     "RAObservationTrueEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "RAObservationFormalEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( AprioriCov,
                                     "RAaprioriCov.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );


    return EXIT_SUCCESS;
}
