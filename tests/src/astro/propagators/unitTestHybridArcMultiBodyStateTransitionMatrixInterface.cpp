/* git    Copyright (c) 2010-2019, Delft University of Technology
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

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{

    namespace unit_tests
    {

//Using declarations.
        using namespace tudat;
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

        BOOST_AUTO_TEST_SUITE( test_hybrid_arc_multi_body_state_transition_matrix_interface )

            static const std::vector< std::string > galileanSatelliteNames = { "Io", "Europa", "Ganymede", "Callisto" };
            static const std::map< std::string, double > satelliteOrbitalPeriods =
                    { { "Io", 1.769 * 86400.0 }, { "Europa", 3.551 * 86400.0 }, { "Ganymede", 7.155 * 86400.0 }, { "Callisto", 16.689 * 86400.0 } };
            static const double centralBodyRotationPeriod = 9.8 * 3600.0;
            static const double ganymedeEllipticalInsertionTime = 32.71 * physical_constants::JULIAN_YEAR;
            static const double ganymedeSphericalInsertionTime = 33.125 * physical_constants::JULIAN_YEAR;
            static const double ganymedeSphericalEndTime = 33.45 * physical_constants::JULIAN_YEAR;

            SystemOfBodies createBodies( const double initialEpoch,
                                         const double finalEpoch,
                                         const std::string globalFrameOrientation )
            {
                // Define bodies settings for simulation
                std::vector< std::string > bodiesToCreate;
                bodiesToCreate.push_back( "Io" );
                bodiesToCreate.push_back( "Europa" );
                bodiesToCreate.push_back( "Ganymede" );
                bodiesToCreate.push_back( "Callisto" );
                bodiesToCreate.push_back( "Jupiter" );
                bodiesToCreate.push_back( "Sun" );

                // Create body objects.
                BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEpoch, finalEpoch/*, "Jupiter"*/ );

                bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                        "Jupiter", globalFrameOrientation, "IAU_Io" );

                bodySettings.at( "Europa" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                        "Jupiter", globalFrameOrientation, "IAU_Europa" );

                bodySettings.at( "Ganymede" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                        "Jupiter", globalFrameOrientation, "IAU_Ganymede" );

                bodySettings.at( "Callisto" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                        "Jupiter", globalFrameOrientation, "IAU_Callisto" );

                // Define settings for JUICE spacecraft
                bodySettings.addSettings( "-28" );
                bodySettings.at( "-28" )->ephemerisSettings = std::make_shared< DirectSpiceEphemerisSettings >( "Jupiter", globalFrameOrientation );

                bodySettings.at( "Io" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
                bodySettings.at( "Europa" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
                bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
                bodySettings.at( "Callisto" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

                bodySettings.at( "Io" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
                bodySettings.at( "Europa" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
                bodySettings.at( "Ganymede" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
                bodySettings.at( "Callisto" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );

                // Create body map
                SystemOfBodies bodies = createSystemOfBodies( bodySettings );

                // Add JUICE to list of bodies
                bodies.createEmptyBody( "JUICE" );
                std::shared_ptr< ephemerides::Ephemeris > juiceEphemeris = std::make_shared< ephemerides::MultiArcEphemeris >(
                        std::map< double, std::shared_ptr< ephemerides::Ephemeris > >( ), "Jupiter", globalFrameOrientation );
                bodies.at( "JUICE" )->setEphemeris( juiceEphemeris );

                return bodies;
            }

            static const std::vector< std::string > flybysBodies = { "Ganymede", "Ganymede", "Ganymede", "Ganymede", "Europa", "Europa", "Callisto", "Ganymede", "Callisto",
                                                                     "Callisto", "Callisto", "Callisto", "Callisto", "Callisto", "Callisto", "Callisto", "Callisto", "Ganymede",
                                                                     "Ganymede", "Ganymede", "Callisto", "Callisto", "Ganymede" };
            static const std::vector< double > flybysTimes = { 957842194.921875, 962785494.140625, 966094442.578125, 967948681.640625, 969142826.953125, 970373478.515625,
                                                               971388256.640625, 973245765.234375, 976662144.140625, 978104000.390625, 979545181.640625, 980986391.015625,
                                                               988195320.703125, 989636516.015625, 991077999.609375, 992519462.109375, 993959258.203125, 995547638.671875,
                                                               998637697.265625, 1000051238.671875, 1001522351.953125, 1006640363.671875, 1018340314.453125 };

// for the first 10 flybys only
            static const std::vector< Eigen::Vector6d > juiceArcWiseInitialStates = {
                    ( Eigen::Vector6d( ) << 269010886.6440233, 109168378.2819441, 44289110.33410433, -6779.580191369188, -2799.333519625582, -1039.517630079758 ).finished( ),
                    ( Eigen::Vector6d( ) << 259598377.0870934, 144339073.2350956, 26715483.5516398, -6609.120943712102, -3572.453422312712, -695.9893371167046 ).finished( ),
                    ( Eigen::Vector6d( ) << 293969429.8000965, 29394060.46455703, 3480581.478900564, -7405.511587264501, -39.53749302205063, -76.72352611412057 ).finished( ),
                    ( Eigen::Vector6d( ) << 299003743.3363333, -6312246.895748794, 24661.48360344982, -7664.812160858565, 779.9926998155253, -5.742191691451843 ).finished( ),
                    ( Eigen::Vector6d( ) << 157165588.7250499, 7777707.313425752, 1294202.497287864, -3988.249166765505, -1221.234587619946, -110.0910244056724 ).finished( ),
                    ( Eigen::Vector6d( ) << 154071000.5605309, 21684830.9608587, -19419226.27090531, -3859.090030510767, -1559.112118690796, 356.4029133965716 ).finished( ),
                    ( Eigen::Vector6d( ) << 209171502.4842756, -65154470.31047933, 1036446.041218763, -4889.40811845498, 1520.87827695955, -23.96290337449799 ).finished( ),
                    ( Eigen::Vector6d( ) << 224970514.1768301, 3797278.908960023, 977664.86709076, -5712.485378054969, -56.12413597720118, -162.3267157603887 ).finished( ),
                    ( Eigen::Vector6d( ) << 210605526.4677481, 10704436.71676197, 22436862.81605911, -4983.707602658618, -297.2719633275694, -449.235406706365 ).finished( ),
                    ( Eigen::Vector6d( ) << 201898490.2017055, -2446115.307207319, 63585109.59455179, -4788.549101122165, 13.37353581103535, -1383.55446795813 ).finished( ) };

            double findLocalMinimumOfTargetDistance(
                    const double lowerBound,
                    const double upperBound,
                    const double tolerance,
                    const std::function< Eigen::Vector3d( const double )> positionFunction )
            {
                double currentUpperBound = upperBound;
                double currentLowerBound = lowerBound;
                double newTestValue = ( currentUpperBound + currentLowerBound ) / 2.0;
                double currentTestValue;

                double currentUpperDistance, currentLowerDistance, currentTestDistance ;
                int currentUpperDerivative, currentLowerDerivative, currentTestDerivative;

                int counter = 0;
                do
                {
                    currentTestValue = newTestValue;

                    currentUpperDistance = positionFunction( currentUpperBound ).norm( );
                    currentLowerDistance  = positionFunction( currentLowerBound ).norm( );
                    currentTestDistance = positionFunction( currentTestValue ).norm( );

                    currentUpperDerivative = utilities::sgn< double >( ( currentUpperDistance - positionFunction( currentUpperBound - tolerance ).norm( ) ) / tolerance );
                    currentLowerDerivative  = utilities::sgn< double >( ( currentLowerDistance - positionFunction( currentLowerBound - tolerance ).norm( ) ) / tolerance );
                    currentTestDerivative = utilities::sgn< double >( ( currentTestDistance - positionFunction( currentTestValue - tolerance ).norm( ) ) / tolerance );

                    if(  currentUpperDerivative > 0 && currentTestDerivative < 0 )
                    {

                        newTestValue = ( currentUpperBound + currentTestValue ) / 2.0;
                        currentLowerBound = currentTestValue;
                    }
                    else if( currentLowerDerivative < 0 && currentTestDerivative > 0 )
                    {
                        newTestValue = ( currentTestValue + currentLowerBound ) / 2.0;
                        currentUpperBound = currentTestValue;

                    }
                    else
                    {

                    }

                    counter ++;
                    if( counter > 1E5 )
                    {
                        std::cerr<<"Error 2 when finding minium target distance"<<std::endl;
                    }

                } while( std::fabs( newTestValue - currentTestValue ) > tolerance );

                return newTestValue;
            }

            std::map< double, double > findLocalMinimaOfTargetDistance(
                    const double lowerBound,
                    const double upperBound,
                    const double threshold,
                    const double tolerance,
                    const double initialSearchTimeStep,
                    const std::function< Eigen::Vector3d( const double )> positionFunction )
            {
                std::map< double, double > minima;

                double upperTime = lowerBound + 2.0 * initialSearchTimeStep,
                        middleTime = lowerBound + initialSearchTimeStep,
                        lowerTime = lowerBound;

                double upperValue, middleValue, lowerValue;

                double candidateTime;

                while( upperTime <= upperBound )
                {
                    upperValue = ( positionFunction( upperTime ) ).norm( );
                    middleValue = ( positionFunction( middleTime ) ).norm( );
                    lowerValue = ( positionFunction( lowerTime ) ).norm( );

                    if( utilities::sgn( upperValue - middleValue ) > 0 &&
                        utilities::sgn( middleValue - lowerValue ) < 0 )
                    {
                        candidateTime = findLocalMinimumOfTargetDistance( lowerTime, upperTime, tolerance, positionFunction );

                        if( positionFunction( candidateTime ).norm( ) <= threshold )
                        {
                            minima[ candidateTime ] = positionFunction( candidateTime ).norm( );
                        }
                    }

                    upperTime += initialSearchTimeStep;
                    middleTime += initialSearchTimeStep;
                    lowerTime += initialSearchTimeStep;


                }
                return minima;
            }

            std::map< double, double > findLocalMinimaOfTargetDistance(
                    const double lowerBound,
                    const double upperBound,
                    const double threshold,
                    const double tolerance,
                    const double initialSearchTimeStep,
                    const std::shared_ptr< Ephemeris > ephemeris )
            {
                return findLocalMinimaOfTargetDistance(
                        lowerBound, upperBound, threshold, tolerance, initialSearchTimeStep,
                        std::bind( &Ephemeris::getCartesianPosition, ephemeris, std::placeholders::_1 ) );
            }

            void getCloseApproachTimes(
                    const double initialTime, const double finalTime, const double approachThreshold,
                    std::vector< double >& flyByTimeVector,
                    std::vector< std::string >& flyByBodyVector )
            {
                std::map< std::string, std::map< double, double > >  closeApproachTimes;

                for( unsigned int i = 0; i < 4; i++ )
                {
                    std::map< double, double > localDistanceMinima = findLocalMinimaOfTargetDistance(
                            initialTime, finalTime, approachThreshold, 5.0, 1800.0,
                            std::bind( spice_interface::getBodyCartesianPositionAtEpoch, "-28", galileanSatelliteNames.at( i ),
                                       "ECLIPJ2000", "None", std::placeholders::_1 ) );
                    if( localDistanceMinima.size( ) > 0 )
                    {
                        closeApproachTimes[ galileanSatelliteNames.at( i ) ] = localDistanceMinima;
                    }
                }

                std::map< double, std::string > timeOrderedFlybyTimes;
                for( auto bodyIterator : closeApproachTimes )
                {
                    for( auto timeIterator : bodyIterator.second )
                    {
                        timeOrderedFlybyTimes[ timeIterator.first ] = bodyIterator.first;

                        std::cout<<"Flyby: "<<bodyIterator.first<<", "<<timeIterator.first <<" "<<timeIterator.second<<std::endl;
                    }
                }

                flyByTimeVector = utilities::createVectorFromMapKeys( timeOrderedFlybyTimes );
                flyByBodyVector = utilities::createVectorFromMapValues( timeOrderedFlybyTimes );
            }

            basic_astrodynamics::AccelerationMap getMultiArcAccelerationModelMap(
                    const SystemOfBodies& bodies,
                    const std::string& multiArcBodyToPropagate,
                    const std::string& multiArcCentralBody )
            {
                using namespace tudat::basic_astrodynamics;
                using namespace tudat::ephemerides;

                std::pair< int, int > maximumEstimatedDegreeAndOrder = std::make_pair< int, int >( 2, 2 );
                int maximumEstimatedDegree = maximumEstimatedDegreeAndOrder.first;
                int maximumEstimatedOrder = maximumEstimatedDegreeAndOrder.second;

                SelectedAccelerationMap accelerationSettingsJuice;
                accelerationSettingsJuice[ "JUICE" ][ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                accelerationSettingsJuice[ "JUICE" ][ multiArcCentralBody ].push_back(
                        std::make_shared< SphericalHarmonicAccelerationSettings >( maximumEstimatedDegree, maximumEstimatedOrder ) );

    for( unsigned int j = 0; j < galileanSatelliteNames.size( ); j++ )
    {
        if( galileanSatelliteNames.at( j ) != multiArcCentralBody )
        {
            accelerationSettingsJuice[ "JUICE" ][ galileanSatelliteNames.at( j ) ].push_back(
                    std::make_shared< AccelerationSettings >( central_gravity ) );
        }
    }

                accelerationSettingsJuice[ "JUICE" ][ "Sun" ].push_back(
                        std::make_shared< AccelerationSettings >( central_gravity ) );

                bool estimateAccelerometerCalibrationsPerArc = true;
                if( estimateAccelerometerCalibrationsPerArc )
                {
                    accelerationSettingsJuice[ "JUICE" ][ multiArcCentralBody ].push_back(
                            std::make_shared< EmpiricalAccelerationSettings >( ) );
                }
                return createAccelerationModelsMap(
                        bodies, accelerationSettingsJuice, { multiArcBodyToPropagate }, { multiArcCentralBody } );

            }





            basic_astrodynamics::AccelerationMap getMoonsAccelerationMap(
                    const SystemOfBodies& bodies,
                    const std::vector< std::string >& bodiesToPropagate,
                    const std::vector< std::string >& centralBodies )
            {
                // Retrieve and process acceleration settings
                int maximumSatelliteDegree = 2;
                int maximumSatelliteOrder = 2;
                int maximumJupiterDegree = 8;
                int maximumJupiterOrder = 0;

                SelectedAccelerationMap accelerationSettingsMoons;

                for( unsigned int i = 0; i < bodiesToPropagate.size( ) ; i++ )
                {
                    accelerationSettingsMoons[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back(
                            std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                                    maximumJupiterDegree, maximumJupiterOrder,
                                    maximumSatelliteDegree, maximumSatelliteOrder ) );
                    for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
                    {
                        if ( i != j )
                        {
                            accelerationSettingsMoons[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back(
                                    std::make_shared< AccelerationSettings >( central_gravity ) );
                        }
                    }
                }

                basic_astrodynamics::AccelerationMap accelerationMap = createAccelerationModelsMap(
                        bodies, accelerationSettingsMoons, bodiesToPropagate, centralBodies );

                return accelerationMap;
            }

            void getMultiArcInitialAndFinalConditions(
                    const double initialTime,
                    const double finalTime,
                    const std::string globalFrameOrientation,
                    std::vector< std::string >& multiArcCentralBodies,
                    std::vector< double >& flybyTimes,
                    std::vector< double >& multiArcStartTimes,
                    std::vector< double >& multiArcEndTimes,
                    std::vector< Eigen::VectorXd >& multiArcJuiceInitialStates,
                    const double arcDuration,
                    const SystemOfBodies& bodies,
                    const std::shared_ptr< numerical_integrators::IntegratorSettings< > > multiArcIntegratorSettings,
                    std::map< int, basic_astrodynamics::AccelerationMap >& multiArcJuiceAccelerationMap )
            {
                multiArcCentralBodies.clear( );
                flybyTimes.clear( );
                multiArcJuiceAccelerationMap.clear( );

                std::vector< Eigen::VectorXd > multiArcSystemInitialStates;

                std::vector< double > juiceArcEndTimes;
                {
                    std::vector< double > allFlybyTimes;
                    std::vector< std::string > allMultiArcCentralBodies;
                    getCloseApproachTimes(
                            initialTime, ( ganymedeEllipticalInsertionTime < finalTime ) ?
                                         ganymedeEllipticalInsertionTime : finalTime, 2.0E7, allFlybyTimes, allMultiArcCentralBodies );

                    for( unsigned int i = 0; i < allFlybyTimes.size( ); i++ )
                    {
                        flybyTimes.push_back( allFlybyTimes.at( i ) );
                        multiArcCentralBodies.push_back( allMultiArcCentralBodies.at( i ) );
                    }
                    for( unsigned int i = 0; i < flybyTimes.size( ); i++ )
                    {
                        juiceArcEndTimes.push_back(flybyTimes.at(i) - arcDuration / 2.0);
                        multiArcSystemInitialStates.push_back( spice_interface::getBodyCartesianStateAtEpoch(
                                "-28", multiArcCentralBodies.at(i), globalFrameOrientation,
                                "None", flybyTimes.at(i) ) );
                    }
                }

                std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
                for( unsigned int i = 0; i < juiceArcEndTimes.size( ); i++ )
                {
                    multiArcJuiceAccelerationMap[ i ] = getMultiArcAccelerationModelMap(
                            bodies, "JUICE", multiArcCentralBodies.at( i ) );
                    arcPropagationSettingsList.push_back(
                            std::make_shared<TranslationalStatePropagatorSettings<double> >
                                    (std::vector<std::string>({multiArcCentralBodies.at(i)}), multiArcJuiceAccelerationMap.at( i ),
                                     std::vector<std::string>({"JUICE"}),
                                     multiArcSystemInitialStates.at(i), juiceArcEndTimes.at(i)));
                }


                std::shared_ptr< propagators::MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
                        std::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList );

                MultiArcDynamicsSimulator< > backwardsFlybyMultiArcDynamicsSimulator =
                        MultiArcDynamicsSimulator< >(
                                bodies, multiArcIntegratorSettings, multiArcPropagatorSettings, flybyTimes,
                                true, false, false );

                std::vector< std::map< double, Eigen::VectorXd > > backwardsFlybyMultiArcStates =
                        backwardsFlybyMultiArcDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                for( unsigned int i = 0; i < flybyTimes.size( ); i++ )
                {
                    multiArcStartTimes.push_back( backwardsFlybyMultiArcStates.at( i ).begin( )->first );
                    multiArcEndTimes.push_back( backwardsFlybyMultiArcStates.at( i ).begin( )->first + arcDuration );
                    multiArcJuiceInitialStates.push_back( backwardsFlybyMultiArcStates.at( i ).begin( )->second );
                }
            }



            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > getParametersToEstimate(
                    const std::shared_ptr< propagators::HybridArcPropagatorSettings< > > hybridArcPropagatorSettings,
                    const SystemOfBodies& bodies,
                    const std::vector< std::string >& multiArcBodiesToPropagate, const std::map< std::string, Eigen::VectorXd > multiArcInitialStates,
                    const std::vector< double >& arcStartTimes, const std::map< std::string, std::vector< double > >& arcStartTimesPerBody,
                    const std::map< std::string, std::vector< std::string > >& multiArcCentralBodies,
                    const std::vector< std::string >& singleArcBodiesToPropagate, const std::vector< std::string >& singleArcCentralBodies,
                    const Eigen::VectorXd& singleArcInitialStates )
            {
                using namespace tudat::propagators;
                using namespace tudat::estimatable_parameters;
                using namespace tudat::basic_astrodynamics;
                using namespace tudat::observation_models;

                std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

                for ( unsigned int i = 0 ; i < singleArcBodiesToPropagate.size( ) ; i++ )
                {
                    parameterNames.push_back(
                            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                    singleArcBodiesToPropagate.at( i ), singleArcInitialStates.segment( i * 6, 6 ), singleArcCentralBodies.at( i ) ) );
                }

                for ( unsigned int i = 0 ; i < multiArcBodiesToPropagate.size( ) ; i++ ) {
                    parameterNames.push_back(
                            std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                    multiArcBodiesToPropagate[ i ], multiArcInitialStates.at( multiArcBodiesToPropagate[ i ] ),
                                    arcStartTimesPerBody.at( multiArcBodiesToPropagate[ i ] ), multiArcCentralBodies.at( multiArcBodiesToPropagate[ i ] ) ) );
                }

                for( unsigned int i = 0; i < galileanSatelliteNames.size( ); i++ )
                {
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                            galileanSatelliteNames.at( i ), gravitational_parameter ) );

                    {
                        std::pair< int, int > maximumDegreeAndOrder = std::make_pair< int, int >( 2, 2 );

                        int maximumEstimatedDegree = maximumDegreeAndOrder.first;
                        int maximumEstimatedOrder = maximumDegreeAndOrder.second;

                        if( maximumEstimatedDegree >= 2 )
                        {
                            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                    2, 0,  2, 2, galileanSatelliteNames.at( i ),
                                    spherical_harmonics_cosine_coefficient_block ) );
                            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                    2, 1, 2, 2,
                                    galileanSatelliteNames.at( i ), spherical_harmonics_sine_coefficient_block ) );
                        }
                    }
                }

                std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
                        std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > empiricalComponentsToEstimate;
                empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( constant_empirical );
                empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
                empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );

                parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                        "JUICE", "Europa", empiricalComponentsToEstimate, arcStartTimes ) );
                parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                        "JUICE", "Ganymede", empiricalComponentsToEstimate, arcStartTimes ) );
//    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
//            "JUICE", "Callisto", empiricalComponentsToEstimate, arcStartTimes ) );

                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                        createParametersToEstimate< double >( parameterNames, bodies, hybridArcPropagatorSettings );
                printEstimatableParameterEntries( parametersToEstimate );

                return parametersToEstimate;

            }

            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > getSingleArcParametersToEstimate(
                    const std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings,
                    const SystemOfBodies& bodies,
                    const std::vector< std::string >& bodiesToPropagate, const Eigen::VectorXd initialStates, const std::vector< double >& arcStartTimes,
                    const std::vector< std::string > centralBodies,
                    const std::string centralBodyJuice )
            {
                using namespace tudat::propagators;
                using namespace tudat::estimatable_parameters;
                using namespace tudat::basic_astrodynamics;
                using namespace tudat::observation_models;

                std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

//    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ ) {
//        parameterNames.push_back(
//                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                        bodiesToPropagate.at( i ), initialStates.segment( i * 6, 6 ), centralBodies.at( i ) ) );
//    }

                for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ ) {
                    parameterNames.push_back(
                            std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                    bodiesToPropagate[ i ], initialStates.segment( 6*i, 6 ),
                                    arcStartTimes, centralBodies[ i ] ) );
                }

                for( unsigned int i = 0; i < galileanSatelliteNames.size( ); i++ )
                {
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                            galileanSatelliteNames.at( i ), gravitational_parameter ) );

                    {
                        std::pair< int, int > maximumDegreeAndOrder = std::make_pair< int, int >( 2, 2 );

                        int maximumEstimatedDegree = maximumDegreeAndOrder.first;
                        int maximumEstimatedOrder = maximumDegreeAndOrder.second;

                        if( maximumEstimatedDegree >= 2 )
                        {
                            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                    2, 0,  2, 2, galileanSatelliteNames.at( i ),
                                    spherical_harmonics_cosine_coefficient_block ) );
                            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                    2, 1, 2, 2,
                                    galileanSatelliteNames.at( i ), spherical_harmonics_sine_coefficient_block ) );
                        }
                    }
                }

                std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
                        std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > empiricalComponentsToEstimate;
                empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( constant_empirical );
                empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
                empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );

                parameterNames.push_back( std::make_shared< EmpiricalAccelerationEstimatableParameterSettings >(
                        "JUICE", centralBodyJuice, empiricalComponentsToEstimate ) );

                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                        createParametersToEstimate< double >( parameterNames, bodies, singleArcPropagatorSettings );
                printEstimatableParameterEntries( parametersToEstimate );

                return parametersToEstimate;

            }

            BOOST_AUTO_TEST_CASE( testHybridArcMultiBodyVariationalEquationCalculation1 )
            {
//    std::string outputFolder = "/home/mfayolle/Documents/PHD/MultiArcMultiBody/";

                //Load spice kernels
                spice_interface::loadStandardSpiceKernels( );

                // Simulation parameters
                double initialEpoch = 946728000.0;
                double finalEpoch = 1200000000.0;

                int numberArcs = 5;

                double propagationTimeStep = 1800.0;
                double flybyDuration = 8.0 * 3600.0;

                std::string globalFrameOrientation = "ECLIPJ2000";
                std::string globalFrameOrigin = "Jupiter";

                // Create body map
                SystemOfBodies bodies = createBodies( initialEpoch, finalEpoch, globalFrameOrientation );

                // Create integrator settings
                std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        std::make_shared< IntegratorSettings< double > >( rungeKutta4, initialEpoch, propagationTimeStep );


                // Compute flybys times and associated central bodies
                std::vector< std::string > multiArcCentralBodies;
                std::vector< double > flybyTimes;
                std::vector< double > multiArcStartTimes, multiArcEndTimes;
                std::vector< Eigen::VectorXd > multiArcJuiceInitialStates;
                std::map< int, basic_astrodynamics::AccelerationMap > multiArcJuiceAccelerationMap;
                getMultiArcInitialAndFinalConditions(
                        initialEpoch, finalEpoch, globalFrameOrientation, multiArcCentralBodies, flybyTimes, multiArcStartTimes, multiArcEndTimes, multiArcJuiceInitialStates,
                        flybyDuration, bodies, integratorSettings, multiArcJuiceAccelerationMap );

                // Different test cases
                // testCase = 0 : same bodies for each arc
                // testCase = 1 : different bodies / nb of bodies for each arc
                for ( unsigned int testCase = 1 ; testCase < 2 ; testCase++ ) {

                    // Define bodies to propagate & associated central bodies
                    std::map<int, std::vector<std::string> > bodiesToPropagatePerArc, centralBodiesPerArc;
                    if ( testCase == 0 )
                    {
                        for (int i = 0; i < numberArcs; i++) {
                            bodiesToPropagatePerArc[i] = {"Io", "Europa", "Ganymede", "Callisto", "JUICE"};
                            centralBodiesPerArc[i] = {"Jupiter", "Jupiter", "Jupiter", "Jupiter", multiArcCentralBodies.at(i)};
                        }
                    }
                    else if ( testCase == 1 )
                    {
                        bodiesToPropagatePerArc[ 0 ] = { "Io", "Europa", "Ganymede", "Callisto", "JUICE" };
                        centralBodiesPerArc[ 0 ] = { "Jupiter", "Jupiter", "Jupiter", "Jupiter", multiArcCentralBodies.at( 0 ) };

                        bodiesToPropagatePerArc[ 1 ] = { "Ganymede", "JUICE" };
                        centralBodiesPerArc[ 1 ] = { "Jupiter", multiArcCentralBodies.at( 1 ) };

                        bodiesToPropagatePerArc[ 2 ] = { "Europa", "Ganymede", "JUICE" };
                        centralBodiesPerArc[ 2 ] = { "Jupiter", "Jupiter", multiArcCentralBodies.at( 2 ) };

                        bodiesToPropagatePerArc[ 3 ] = { "Io", "Ganymede", "Callisto", "JUICE" };
                        centralBodiesPerArc[ 3 ] = { "Jupiter", "Jupiter", "Jupiter", multiArcCentralBodies.at( 3 ) };

                        bodiesToPropagatePerArc[ 4 ] = { "Europa", "Callisto", "JUICE" };
                        centralBodiesPerArc[ 4 ] = { "Jupiter", "Jupiter", multiArcCentralBodies.at( 4 ) };
                    }

                    // Set accelerations map for the moons.
                    std::vector<std::string> moonsToPropagate, centralBodiesForMoons;
                    std::map<int, AccelerationMap> multiArcMoonsAccelerationMap;
                    for (int i = 0; i < numberArcs; i++) {
                        std::vector<std::string> arcWiseMoonsToPropagate, arcWiseMoonsCentralBodies;
                        for (int j = 0; j < bodiesToPropagatePerArc.at( i ).size() - 1; j++) {
                            arcWiseMoonsToPropagate.push_back( bodiesToPropagatePerArc.at( i ).at(j ) );
                            arcWiseMoonsCentralBodies.push_back( centralBodiesPerArc.at( i ).at( j ) );
                        }
                        multiArcMoonsAccelerationMap[ i ] = getMoonsAccelerationMap(bodies, arcWiseMoonsToPropagate, arcWiseMoonsCentralBodies);
                    }

                    // Merge arc-wise acceleration maps
                    std::map<int, AccelerationMap> multiArcCompleteAccelerationMaps = multiArcMoonsAccelerationMap;
                    for (unsigned int k = 0; k < numberArcs; k++) {
                        AccelerationMap arcWiseAccelerationMap = multiArcCompleteAccelerationMaps.at( k );
                        arcWiseAccelerationMap[ "JUICE" ] = multiArcJuiceAccelerationMap.at( k ).at("JUICE");
                        multiArcCompleteAccelerationMaps.at( k ) = arcWiseAccelerationMap;
                    }



                    // Create multi-arc propagator settings
                    std::vector<std::shared_ptr<SingleArcPropagatorSettings<> > > propagatorSettingsList;
                    std::map<std::string, std::vector<std::string> > multiArcCentralBodiesPerBody;
                    std::map<std::string, Eigen::VectorXd> multiArcInitialStatesPerBody;
                    std::map<std::string, std::vector<std::pair<int, int> > > counterArcsPerBody;
                    std::map<int, Eigen::VectorXd> multiArcInitialStates;
                    std::vector<double> arcStartTimes;
                    std::map< std::string, std::vector< double > > arcStartTimesPerBody;
                    std::vector< std::string > listBodiesToPropagate;

                    for (int arc = 0; arc < numberArcs; arc++) {
                        arcStartTimes.push_back(multiArcStartTimes.at(arc));

                        Eigen::VectorXd arcWiseConcatenatedStates;
                        arcWiseConcatenatedStates.resize(6 * bodiesToPropagatePerArc.at(arc).size());
                        for (unsigned int i = 0; i < bodiesToPropagatePerArc.at(arc).size(); i++) {

                            counterArcsPerBody[bodiesToPropagatePerArc.at(arc).at(i)].push_back(std::make_pair(arc, i));

                            arcStartTimesPerBody[ bodiesToPropagatePerArc.at(arc).at(i) ].push_back( multiArcStartTimes.at( arc ) );

                            bool bodyAlreadyIncluded = false;
                            for ( unsigned int k = 0 ; k < listBodiesToPropagate.size( ) ; k++ )
                            {
                                if ( listBodiesToPropagate[ k ] == bodiesToPropagatePerArc.at(arc).at(i) )
                                {
                                    bodyAlreadyIncluded = true;
                                }
                            }
                            if ( !bodyAlreadyIncluded )
                            {
                                listBodiesToPropagate.push_back( bodiesToPropagatePerArc.at(arc).at(i) );
                            }



                            if (bodiesToPropagatePerArc.at(arc).at(i) == "JUICE") {
                                arcWiseConcatenatedStates.segment(i * 6, 6) = multiArcJuiceInitialStates.at(arc).segment(0, 6);
                            } else {
                                arcWiseConcatenatedStates.segment(i * 6, 6) = spice_interface::getBodyCartesianStateAtEpoch(
                                        bodiesToPropagatePerArc.at(arc).at(i), centralBodiesPerArc.at(arc).at(i), globalFrameOrientation,
                                        "None", multiArcStartTimes.at(arc));
                            }

                            multiArcCentralBodiesPerBody[bodiesToPropagatePerArc.at(arc).at(i)].push_back(centralBodiesPerArc.at(arc).at(i));
                        }

                        multiArcInitialStates[arc] = arcWiseConcatenatedStates;
                        propagatorSettingsList.push_back(std::make_shared<TranslationalStatePropagatorSettings<> >
                                                                 (centralBodiesPerArc.at(arc), multiArcCompleteAccelerationMaps.at(arc), bodiesToPropagatePerArc.at(arc),
                                                                  arcWiseConcatenatedStates, multiArcEndTimes.at(arc)));
                    }
                    std::cout << "create multi-arc propagator settings" << "\n\n";
                    std::shared_ptr<MultiArcPropagatorSettings<> > multiArcPropagatorSettings =
                            std::make_shared<MultiArcPropagatorSettings<> >(propagatorSettingsList);

                    // Define single-arc propagator settings for Jupiter
                    std::vector< std::string > singleArcCentralBody = { "Sun" };
                    std::vector< std::string > singleArcPropagatedBody = { "Jupiter" };
                    SelectedAccelerationMap accelerationSettingsJupiter;
                    accelerationSettingsJupiter[ "Jupiter" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                    basic_astrodynamics::AccelerationMap jupiterAccelerationModelMap = createAccelerationModelsMap(
                            bodies, accelerationSettingsJupiter, singleArcPropagatedBody, singleArcCentralBody );

                    Eigen::VectorXd singleArcInitialStates = propagators::getInitialStatesOfBodies(
                            singleArcPropagatedBody, singleArcCentralBody, bodies, initialEpoch );

                    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
                            std::make_shared< TranslationalStatePropagatorSettings< > >(
                                    singleArcCentralBody, jupiterAccelerationModelMap, singleArcPropagatedBody, singleArcInitialStates, finalEpoch );

                    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
                            std::make_shared< HybridArcPropagatorSettings< > >( singleArcPropagatorSettings, multiArcPropagatorSettings );


                    for (auto itr: multiArcCentralBodiesPerBody) {
                        Eigen::VectorXd arcWiseStatesCurrentBody;
                        arcWiseStatesCurrentBody.resize(6 * itr.second.size());
                        std::vector<std::pair<int, int> > currentBodyArcIndices = counterArcsPerBody.at(itr.first);
                        for (unsigned int j = 0; j < itr.second.size(); j++) {
                            arcWiseStatesCurrentBody.segment(6 * j, 6) =
                                    multiArcInitialStates.at(currentBodyArcIndices.at(j).first).segment(currentBodyArcIndices.at(j).second * 6, 6);
                        }
                        multiArcInitialStatesPerBody[itr.first] = arcWiseStatesCurrentBody;
                    }

                    // Create parameters to estimate
                    std::shared_ptr<estimatable_parameters::EstimatableParameterSet<double> > parametersToEstimate = getParametersToEstimate(
                            hybridArcPropagatorSettings, bodies, listBodiesToPropagate, multiArcInitialStatesPerBody,
                            arcStartTimes, arcStartTimesPerBody, multiArcCentralBodiesPerBody, singleArcPropagatedBody, singleArcCentralBody, singleArcInitialStates );



//    SingleArcVariationalEquationsSolver< double, double > singleArcVariationalEquations =
//            SingleArcVariationalEquationsSolver< double, double >(
//                    bodies, integratorSettings, std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< > >( propagatorSettingsList.at( 0 ) ),
//                    singleArcParametersToEstimate/*, true, nullptr, false, true, false*/ );

//        MultiArcVariationalEquationsSolver<double, double> multiArcVariationalEquations =
//                MultiArcVariationalEquationsSolver<double, double>(
//                        bodies, integratorSettings, multiArcPropagatorSettings, parametersToEstimate, arcStartTimes, true,
//                        std::shared_ptr<numerical_integrators::IntegratorSettings<double> >(), false, true, true);

//        SystemOfBodies bodies2 = createBodies( initialEpoch, finalEpoch, globalFrameOrientation );

                    std::vector< std::shared_ptr< EstimatableParameterSettings > > singleArcParameterNames;
                    singleArcParameterNames.push_back(
                            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                    "Jupiter", singleArcInitialStates.segment( 0, 6 ), "Sun" ) );
                    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > singleArcParametersToEstimate =
                            createParametersToEstimate< double >( singleArcParameterNames, bodies, singleArcPropagatorSettings );
                    printEstimatableParameterEntries( singleArcParametersToEstimate );


                    std::vector< std::shared_ptr< EstimatableParameterSettings > > multiArcParameterNames;
                    for ( unsigned int i = 0 ; i < listBodiesToPropagate.size( ) ; i++ ) {
                        multiArcParameterNames.push_back(
                                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                        listBodiesToPropagate[ i ], multiArcInitialStatesPerBody.at( listBodiesToPropagate[ i ] ),
                                        arcStartTimesPerBody.at( listBodiesToPropagate[ i ] ), multiArcCentralBodiesPerBody.at( listBodiesToPropagate[ i ] ) ) );
                    }
                    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > multiArcParametersToEstimate =
                            createParametersToEstimate< double >( multiArcParameterNames, bodies, multiArcPropagatorSettings );
                    printEstimatableParameterEntries( multiArcParametersToEstimate );


                    integratorSettings->initialTime_ = initialEpoch;
                    HybridArcVariationalEquationsSolver< double, double > hybridArcVariationalEquationsSolver =
                            HybridArcVariationalEquationsSolver< double, double >( bodies, integratorSettings, hybridArcPropagatorSettings, parametersToEstimate,
                                                                                   arcStartTimes, true, false, true );

                    std::cout << "full parameter vector size from hybrid arc state transition matrix interface: " << "\n\n";
                    std::cout << hybridArcVariationalEquationsSolver.getStateTransitionMatrixInterface( )->getFullParameterVectorSize( ) << "\n\n";
                    Eigen::MatrixXd combinedMatrix = hybridArcVariationalEquationsSolver.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix(
                            ( arcStartTimes[ 0 ] + multiArcEndTimes[ 0 ] ) / 2.0 );
                    std::cout << "combinedMatrix: " << "\n\n";
                    std::cout << combinedMatrix << "\n\n";
//                    std::cout << "test combined matrix arc 0: " << "\n\n";
//                    std::cout << hybridArcVariationalEquationsSolver.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix(
//                            ( arcStartTimes[ 0 ] + multiArcEndTimes[ 0 ] ) / 2.0 ) << "\n\n";
//                    std::cout << "single-arc solution" << "\n\n";
//                    std::cout << hybridArcVariationalEquationsSolver.getSingleArcSolver()->getStateTransitionMatrixInterface()->getCombinedStateTransitionAndSensitivityMatrix(
//                            ( arcStartTimes[ 0 ] + multiArcEndTimes[ 0 ] ) / 2.0 ) << "\n\n";
//                    std::cout << "multi-arc solution: " << "\n\n";
//                    std::cout << hybridArcVariationalEquationsSolver.getMultiArcSolver()->getStateTransitionMatrixInterface()->getCombinedStateTransitionAndSensitivityMatrix(
//                            ( arcStartTimes[ 0 ] + multiArcEndTimes[ 0 ] ) / 2.0 ) << "\n\n";

                    std::vector< std::map< double, Eigen::MatrixXd > > singleArcVariationalSolution =
                            hybridArcVariationalEquationsSolver.getSingleArcSolver( )->getNumericalVariationalEquationsSolution( );
                    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > fullMultiArcVariationalSolution =
                            hybridArcVariationalEquationsSolver.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );

                    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > > singleArcStateTransitionMatrixInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( singleArcVariationalSolution[ 0 ] ),
                                    utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( singleArcVariationalSolution[ 0 ] ), 4 );
                    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > > singleArcSensitivityMatrixInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( singleArcVariationalSolution[ 1 ] ),
                                    utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( singleArcVariationalSolution[ 1 ] ), 4 );

                    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > > multiArc0StateTransitionMatrixInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( fullMultiArcVariationalSolution[ 0 ][ 0 ] ),
                                    utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( fullMultiArcVariationalSolution[ 0 ][ 0 ] ), 4 );
                    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > > multiArc0SensitivityMatrixInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::MatrixXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::MatrixXd, double >( fullMultiArcVariationalSolution[ 0 ][ 1 ] ),
                                    utilities::createVectorFromMapValues< Eigen::MatrixXd, double >( fullMultiArcVariationalSolution[ 0 ][ 1 ] ), 4 );

//                    std::cout << "size STM arc 0: " << fullMultiArcVariationalSolution[ 0 ][ 0 ].begin( )->second.rows( ) << " & " <<
//                              fullMultiArcVariationalSolution[ 0 ][ 0 ].begin( )->second.cols( ) << "\n\n";
//                    std::cout << "size STM arc 1: " << fullMultiArcVariationalSolution[ 1 ][ 0 ].begin( )->second.rows( ) << " & " <<
//                              fullMultiArcVariationalSolution[ 1 ][ 0 ].begin( )->second.cols( ) << "\n\n";

                    double interpolationTime = ( arcStartTimes[ 1 ] + multiArcEndTimes[ 1 ] ) / 2.0;
                    Eigen::MatrixXd singleArcStateTransitionMatrix = singleArcStateTransitionMatrixInterpolator->interpolate( interpolationTime );
                    Eigen::MatrixXd singleArcSensitivityMatrix = singleArcSensitivityMatrixInterpolator->interpolate( interpolationTime );

                    Eigen::MatrixXd multiArcStateTransitionMatrix = multiArc0StateTransitionMatrixInterpolator->interpolate( interpolationTime );
                    Eigen::MatrixXd multiArcSensitivityMatrix = multiArc0SensitivityMatrixInterpolator->interpolate( interpolationTime );
                    std::cout << "multiArcStateTransitionMatrix: " << multiArcStateTransitionMatrix << "\n\n";

                    Eigen::MatrixXd reconstructedHybridSolution =
                            Eigen::MatrixXd::Zero( fullMultiArcVariationalSolution[ 0 ][ 0 ].begin( )->second.rows( ),
                                                   fullMultiArcVariationalSolution[ 0 ][ 0 ].begin( )->second.cols( ) + fullMultiArcVariationalSolution[ 0 ][ 1 ].begin( )->second.cols( ) );
                    std::cout << "size reconstructed hybrid solution: " << reconstructedHybridSolution.rows( ) << " & " << reconstructedHybridSolution.cols( ) << "\n\n";



                    Eigen::MatrixXd fullCombinedMatrix = hybridArcVariationalEquationsSolver.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix(
                            ( arcStartTimes[ 0 ] + multiArcEndTimes[ 0 ] ) / 2.0 );
                    std::cout << "full combined matrix: " << "\n\n";
                    std::cout << fullCombinedMatrix << "\n\n";






//        integratorSettings->initialTime_ = initialEpoch;
//        SingleArcVariationalEquationsSolver< double, double > singleArcVariationalEquationsSolverTest =
//                SingleArcVariationalEquationsSolver< double, double >( bodies, integratorSettings, singleArcPropagatorSettings, singleArcParametersToEstimate,
//                                                                       true, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, true );
//
//        MultiArcVariationalEquationsSolver< double, double > multiArcVariationalEquationsSolverTest =
//                MultiArcVariationalEquationsSolver< double, double >( bodies, integratorSettings, multiArcPropagatorSettings, multiArcParametersToEstimate, arcStartTimes,
//                                                                      true, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, false );
//
//        std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > singleArcVariationalEquationsSolver =
//                hybridArcVariationalEquationsSolver.getSingleArcSolver( );
//
//        std::map< double, Eigen::MatrixXd > singleArcSTM = singleArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 0 ];


//        std::shared_ptr< MultiArcVariationalEquationsSolver< double, double > > multiArcVariationalEquationsSolver =
//                hybridArcVariationalEquationsSolver.getMultiArcSolver( );
//
//
//
//        std::map< double, Eigen::MatrixXd > singleArcStmTest = singleArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ 0 ];
//        std::cout << "last element single-arc STM history: " << "\n\n";
//        std::cout << singleArcSTM.rbegin( )->second << "\n\n";
//        std::cout << "last element single-arc STM history TEST: " << "\n\n";
//        std::cout << singleArcStmTest.rbegin( )->second << "\n\n";
//
//        std::map< double, Eigen::MatrixXd > multiArcStmArc1 = multiArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 0 ][ 0 ];
//        std::map< double, Eigen::MatrixXd > multiArcStmArc3 = multiArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 2 ][ 0 ];
//        std::cout << "last element multi-arc (arc 1) STM history: " << "\n\n";
//        std::cout << multiArcStmArc1.rbegin( )->second << "\n\n";

//        std::map< double, Eigen::MatrixXd > multiArcStmArc1Test = multiArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ 0 ][ 0 ];
//        std::map< double, Eigen::MatrixXd > multiArcStmArc3Test = multiArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ 2 ][ 0 ];
//        std::cout << "last element multi-arc (arc 1) STM history: " << "\n\n";
//        std::cout << multiArcStmArc1Test.rbegin( )->second << "\n\n";
//
//        std::cout << "last element multi-arc (arc 3) STM history: " << "\n\n";
//        std::cout << multiArcStmArc3.rbegin( )->second << "\n\n";
//        std::cout << "last element multi-arc (arc 3) STM history: " << "\n\n";
//        std::cout << multiArcStmArc3Test.rbegin( )->second << "\n\n";
//
//        // Test consistency of single-arc variational equations solutions.
//        for ( auto itr : singleArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 0 ] )
//        {
//            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                    itr.second, singleArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ 0 ].at( itr.first ), 1.0E-12 );
//            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                    singleArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 1 ].at( itr.first ),
//                    singleArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ 1 ].at( itr.first ), 1.0E-12 );
//        }

//        // Test consistency of multi-arc variational equations solutions.
//        for ( unsigned int k = 0 ; k < numberArcs ; k++ )
//        {
//            for ( auto itr : multiArcVariationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ k ][ 0 ] )
//            {
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        itr.second.block( 6, 6, itr.second.rows( ) - 6, itr.second.cols( ) - 6 ),
//                        multiArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ k ][ 0 ].at( itr.first ), 1.0E-12 );
//
//                Eigen::MatrixXd sensitivityMatrixFromHybridSolution =
//                        multiArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ k ][ 1 ].at( itr.first );
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        sensitivityMatrixFromHybridSolution.block( 6, 6, sensitivityMatrixFromHybridSolution.rows( ) - 6, sensitivityMatrixFromHybridSolution.cols( ) - 6 ),
//                        multiArcVariationalEquationsSolverTest.getNumericalVariationalEquationsSolution( )[ k ][ 1 ].at( itr.first ), 1.0E-12 );
//
//
//            }
//        }

                }

            }


        BOOST_AUTO_TEST_SUITE_END( )

    }

}