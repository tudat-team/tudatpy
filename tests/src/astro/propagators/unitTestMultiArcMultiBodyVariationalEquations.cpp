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
#include "tudat/simulation/estimation.h"

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
using namespace tudat::observation_models;

BOOST_AUTO_TEST_SUITE( test_multi_arc_multi_body_variational_equations_calculation )

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
    bodiesToCreate.push_back( "Earth" );

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

//    for( unsigned int j = 0; j < galileanSatelliteNames.size( ); j++ )
//    {
//        if( galileanSatelliteNames.at( j ) != multiArcCentralBody )
//        {
//            accelerationSettingsJuice[ "JUICE" ][ galileanSatelliteNames.at( j ) ].push_back(
//                    std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        }
//    }

    accelerationSettingsJuice[ "JUICE" ][ "Sun" ].push_back(
            std::make_shared< AccelerationSettings >( point_mass_gravity ) );

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
        const std::shared_ptr< propagators::MultiArcPropagatorSettings< > > multiArcPropagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< std::string >& multiArcBodiesToPropagate, const std::map< std::string, Eigen::VectorXd > multiArcInitialStates,
        const std::vector< double >& arcStartTimes, const std::map< std::string, std::vector< double > >& arcStartTimesPerBody,
        const std::map< std::string, std::vector< std::string > >& multiArcCentralBodies,
        const LinkEnds linkEnds )
{
    using namespace tudat::propagators;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::observation_models;

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

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

//    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
//            linkEnds, position_observable, arcStartTimes, observed_body, true ) );

    std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > empiricalComponentsToEstimate;
    empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( constant_empirical );
    empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
    empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );

    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(//    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
//            linkEnds, position_observable, arcStartTimes, observed_body, true ) );
            "JUICE", "Europa", empiricalComponentsToEstimate, arcStartTimes ) );
    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
            "JUICE", "Ganymede", empiricalComponentsToEstimate, arcStartTimes ) );
//    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
//            "JUICE", "Callisto", empiricalComponentsToEstimate, arcStartTimes ) );

//    parameterNames.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
//            "global_metric", ppn_parameter_gamma ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, multiArcPropagatorSettings );
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

BOOST_AUTO_TEST_CASE( testMultiArcMultiBodyVariationalEquationCalculation1 )
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
            for (unsigned int j = 0; j < bodiesToPropagatePerArc.at( i ).size() - 1; j++) {
                arcWiseMoonsToPropagate.push_back( bodiesToPropagatePerArc.at( i ).at(j ) );
                arcWiseMoonsCentralBodies.push_back( centralBodiesPerArc.at( i ).at( j ) );
            }
            multiArcMoonsAccelerationMap[ i ] = getMoonsAccelerationMap(bodies, arcWiseMoonsToPropagate, arcWiseMoonsCentralBodies);
        }

        // Merge arc-wise acceleration maps
        std::map<int, AccelerationMap> multiArcCompleteAccelerationMaps = multiArcMoonsAccelerationMap;
        for (int k = 0; k < numberArcs; k++) {
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
        std::shared_ptr<MultiArcPropagatorSettings<> > multiArcPropagatorSettings =
                std::make_shared<MultiArcPropagatorSettings<> >(propagatorSettingsList);

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

//    SingleArcDynamicsSimulator::SingleArcDynamicsSimulator< double, double >(  )

        // Define links in simulation.
        LinkEnds linkEndsJuice;
        linkEndsJuice[ observed_body ] = std::make_pair< std::string, std::string >( "JUICE", "" );
        LinkEnds linkEndsGanymede;
        linkEndsGanymede[ transmitter ] = std::make_pair< std::string, std::string >( "Ganymede", "" );
        linkEndsGanymede[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "" );

        // Create parameters to estimate
        std::shared_ptr<estimatable_parameters::EstimatableParameterSet<double> > parametersToEstimate = getParametersToEstimate(
                multiArcPropagatorSettings, bodies, listBodiesToPropagate, multiArcInitialStatesPerBody,
                arcStartTimes, arcStartTimesPerBody, multiArcCentralBodiesPerBody, linkEndsJuice );

        std::cout << "parameters values: " << parametersToEstimate->getFullParameterValues< double >( ).transpose( ) << "\n\n";



//    SingleArcVariationalEquationsSolver< double, double > singleArcVariationalEquations =
//            SingleArcVariationalEquationsSolver< double, double >(
//                    bodies, integratorSettings, std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< > >( propagatorSettingsList.at( 0 ) ),
//                    singleArcParametersToEstimate/*, true, nullptr, false, true, false*/ );

        MultiArcVariationalEquationsSolver<double, double> multiArcVariationalEquations =
                MultiArcVariationalEquationsSolver<double, double>(
                        bodies, integratorSettings, multiArcPropagatorSettings, parametersToEstimate, arcStartTimes, true,
                        std::shared_ptr<numerical_integrators::IntegratorSettings<double> >(), false, true, true);

        std::vector<std::map<double, Eigen::VectorXd> > multiArcStateHistory =
                multiArcVariationalEquations.getDynamicsSimulator()->getEquationsOfMotionNumericalSolution();
        std::vector<std::vector<std::map<double, Eigen::MatrixXd> > > variationalEquationsSolution =
                multiArcVariationalEquations.getNumericalVariationalEquationsSolution();

        std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface =
                multiArcVariationalEquations.getStateTransitionMatrixInterface( );
        std::shared_ptr< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface > multiArcStateTransitionInterface =
                std::dynamic_pointer_cast< MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >( stateTransitionMatrixInterface );

        for ( int arc = 0 ; arc < numberArcs ; arc++ )
        {
            std::cout << "ARC " << arc << "\n\n";
            double finalArcEpoch = multiArcStateHistory.at( arc ).rbegin( )->first;
            Eigen::VectorXd finalStates = multiArcStateHistory.at( arc ).rbegin( )->second;
            for ( unsigned int i = 0 ; i < bodiesToPropagatePerArc.at( arc ).size( ) ; i++ )
            {
                std::cout << "BODY: " << bodiesToPropagatePerArc.at( arc ).at( i ) << " - CENTRAL BODY: " << centralBodiesPerArc.at( arc ).at( i ) << "\n\n";
                Eigen::Vector6d currentSolutionFromEphemeris = bodies.at( bodiesToPropagatePerArc.at( arc ).at( i ) )->getEphemeris( )->getCartesianState( finalArcEpoch );
                Eigen::Vector6d stateCentralBody = bodies.at( centralBodiesPerArc.at( arc ).at( i ) )->getEphemeris( )->getCartesianState( finalArcEpoch );
                Eigen::Vector6d differenceStateHistoryWrtEphemeris = currentSolutionFromEphemeris - finalStates.segment( i * 6, 6 );
//                if ( centralBodiesPerArc.at( arc ).at( i )  == bodies.at( bodiesToPropagatePerArc.at( arc ).at( i ) )->getEphemeris( )->getReferenceFrameOrigin( ) )
//                {
//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                            differenceStateHistoryWrtEphemeris, Eigen::Vector6d::Zero( ), 1.0E-16);
//                }
//                else
//                {
//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                            differenceStateHistoryWrtEphemeris, stateCentralBody, 1.0E-12);
//                }


            }

//            if ( arc == 0 ) {
                Eigen::MatrixXd combinedStateTransitionMatrix = multiArcStateTransitionInterface->getCombinedStateTransitionAndSensitivityMatrix(
                        (arcStartTimes[arc] + multiArcEndTimes[arc]) / 2.0, true, std::vector< std::string >( ) );
                Eigen::MatrixXd fullCombinedStateTransitionMatrix = multiArcStateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix(
                        (arcStartTimes[arc] + multiArcEndTimes[arc]) / 2.0, true, std::vector< std::string >( ) );
                std::cout << "arc: " << arc << " - size combinedStateTransitionMatrix: " << combinedStateTransitionMatrix.rows() <<
                          " & " << combinedStateTransitionMatrix.cols() << "\n\n";
                std::cout << "arc: " << arc << " - size fullCombinedStateTransitionMatrix: " << fullCombinedStateTransitionMatrix.rows() <<
                          " & " << fullCombinedStateTransitionMatrix.cols() << "\n\n";

//                std::cout << "combined state transition matrix: " << "\n\n";
//                std::cout << combinedStateTransitionMatrix << "\n\n";
//                std::cout << "full combined state transition matrix: " << "\n\n";
//                std::cout << fullCombinedStateTransitionMatrix << "\n\n";
//            }



        }

        std::vector< Eigen::VectorXd > arcBiases;
        for ( unsigned int k = 0 ; k < arcStartTimes.size( ) ; k++ )
        {
            arcBiases.push_back( ( Eigen::Vector3d( ) <<  0.0, 0.0, 0.0 ).finished( ) );
        }
        std::shared_ptr< ArcWiseConstantObservationBiasSettings > biasSettings =
                std::make_shared< ArcWiseConstantObservationBiasSettings >( arcStartTimes, arcBiases, observed_body, true );

        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections;
        std::vector< std::string > relativisticPerturbingBodies = { "Sun" };
        lightTimeCorrections.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                relativisticPerturbingBodies ) );

        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                position_observable, linkEndsJuice, std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                one_way_range, linkEndsGanymede, lightTimeCorrections /*std::shared_ptr< LightTimeCorrectionSettings >( )*/ ) );

        OrbitDeterminationManager< double, double > orbitDeterminationManager =
                OrbitDeterminationManager< double, double >(
                        bodies, parametersToEstimate,
                        observationSettingsList, integratorSettings, multiArcPropagatorSettings );

        // Compute observation times
        std::vector< double > observationTimes;
        double vlbiObservationCadence = 600.0;
        for ( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
        {
            double currentTime = arcStartTimes[ i ] + 1000.0;
            while ( currentTime < multiArcEndTimes[ i ] - 1000.0 ) {
                observationTimes.push_back( currentTime );
                currentTime += vlbiObservationCadence;
            }
        }

        // Simulate JUICE observations
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                        position_observable, linkEndsJuice, observationTimes, observed_body ) );
        measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                        one_way_range, linkEndsGanymede, observationTimes, receiver ) );

        std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        // Set observations weights.
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightPerObservationParser;
        weightPerObservationParser[ observationParser( position_observable ) ] = 1.0 / ( 1.0 * 1.0 );
        weightPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
        observationsAndTimes->setConstantWeightPerObservable( weightPerObservationParser );

        // Define POD input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >(
                        observationsAndTimes );

        estimationInput->defineEstimationSettings( false, false, true, true, true, false );
        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        std::cout << "from OD manager, state transition matrix size: " << orbitDeterminationManager.
        getStateTransitionAndSensitivityMatrixInterface()->getStateTransitionMatrixSize( ) << "\n\n";
        std::cout << "from OD manager, sensitivity matrix size: " << orbitDeterminationManager.
                getStateTransitionAndSensitivityMatrixInterface()->getSensitivityMatrixSize( ) << "\n\n";

//        ///////////////////////////////////////////
//        /// Test resetting parameters
//        ///////////////////////////////////////////
//        Eigen::VectorXd fullParametersValues = parametersToEstimate->getFullParameterValues< double >( );
//        std::cout << "parameters values: " << fullParametersValues.transpose( ) << "\n\n";
//
//        unsigned int numberStates = 0;
//        for ( unsigned int i = 0 ; i < numberArcs ; i++ )
//        {
//            numberStates += bodiesToPropagatePerArc.at( i ).size( );
//        }
//        Eigen::VectorXd concatenatedMultiArcStates( numberStates * 6 );
//        concatenatedMultiArcStates << multiArcInitialStatesPerBody.at( "Io" ), multiArcInitialStatesPerBody.at( "Europa" ),
//                multiArcInitialStatesPerBody.at( "Ganymede" ), multiArcInitialStatesPerBody.at( "Callisto" ), multiArcInitialStatesPerBody.at( "JUICE" );
//        std::cout << "multi-arc initial states: " << concatenatedMultiArcStates.transpose( ) << "\n\n";
//
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                fullParametersValues.segment( 0, 6 * numberStates ), concatenatedMultiArcStates, 1.0E-16);
//
//        Eigen::VectorXd dummyNewParametersValues( fullParametersValues.size( ) );
//        for ( unsigned int i = 0 ; i < numberStates ; i++ )
//        {
//            dummyNewParametersValues.segment( 6 * i, 6 ) = ( i + 1 ) * ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
//        }
//        for ( unsigned int i = 6 * numberStates ; i < fullParametersValues.size( ) ; i++ )
//        {
//            dummyNewParametersValues[ i ] = - 1.0;
//        }
////        std::cout << "dummyNewParametersValues: " << dummyNewParametersValues.transpose( ) << "\n\n";
//        parametersToEstimate->resetParameterValues( dummyNewParametersValues );
////        std::cout << "new parameters values after reset: " << parametersToEstimate->getFullParameterValues< double >( ).transpose( ) << "\n\n";
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                parametersToEstimate->getFullParameterValues< double >( ), dummyNewParametersValues, 1.0E-16);
//
//        std::cout << "orbit determination manager parameter estimate before reset: " << "\n\n";
//        std::cout << orbitDeterminationManager.getCurrentParameterEstimate( ).transpose( ) << "\n\n";
//        orbitDeterminationManager.resetParameterEstimate( fullParametersValues /*dummyNewParametersValues*/, false );
//        std::cout << "orbit determination manager parameter estimated after reset: " << "\n\n";
//        std::cout << orbitDeterminationManager.getCurrentParameterEstimate( ).transpose( ) << "\n\n";
//
//
//        parametersToEstimate->resetParameterValues( dummyNewParametersValues );






        // Test arc per arc
        for ( int arc = 0; arc < numberArcs; arc++) {

            std::cout << "TEST - ARC : " << arc << "\n\n";
            std::vector<std::shared_ptr<SingleArcPropagatorSettings<> > > perArcPropagatorSettingsList;
            perArcPropagatorSettingsList.push_back( propagatorSettingsList.at( arc ) );
            std::shared_ptr<MultiArcPropagatorSettings<> > perArcPropagatorSettings =
                    std::make_shared<MultiArcPropagatorSettings<> >( perArcPropagatorSettingsList );

            // Create single-arc parameters to estimate
            std::shared_ptr<estimatable_parameters::EstimatableParameterSet<double> > singleArcParametersToEstimate = getSingleArcParametersToEstimate(
                    std::dynamic_pointer_cast<TranslationalStatePropagatorSettings<> >( propagatorSettingsList.at( arc ) ), bodies, bodiesToPropagatePerArc.at( arc ),
                    multiArcInitialStates.at(arc), { arcStartTimes.at( arc ) }, centralBodiesPerArc.at( arc ), multiArcCentralBodies.at( arc ) );

            MultiArcVariationalEquationsSolver<double, double> perArcVariationalEquations =
                    MultiArcVariationalEquationsSolver<double, double>(
                            bodies, integratorSettings, perArcPropagatorSettings,
                            singleArcParametersToEstimate, {arcStartTimes.at(arc)}, true,
                            std::shared_ptr<numerical_integrators::IntegratorSettings<double> >(), false, true, true);

            // Comparison - state histories
            std::vector< std::map< double, Eigen::VectorXd > > perArcStateHistory =
                    perArcVariationalEquations.getDynamicsSimulator()->getEquationsOfMotionNumericalSolution();

            // Comparison - variational equations solutions
            std::vector<std::vector<std::map<double, Eigen::MatrixXd> > > perArcVariationalEquationsSolution =
                    perArcVariationalEquations.getNumericalVariationalEquationsSolution();

            std::map< double, std::vector< std::string > > bodiesToEstimatePerArc;
            std::vector< int > multiArcStateParametersSizePerArc;
            std::map< double, std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > >
            perArcMultiArcParametersToEstimate = getMultiArcDynamicalStateToEstimatePerArc(
                            parametersToEstimate->getEstimatedMultiArcInitialStateParameters( ),
                            bodiesToEstimatePerArc, multiArcStateParametersSizePerArc );

            std::vector< double > arcStartingTimes = parametersToEstimate->getArcStartingTimes( );
            for ( int i = 0 ; i < numberArcs ; i++ )
            {
                std::cout << "size parameters arc " << i << " = " << getSingleArcParameterSetSize( parametersToEstimate, i ) << "\n\n";
                std::cout << "size dynamical parameters arc " << i << " = " << getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate, i ) << "\n\n";
            }

            for ( unsigned int i = 0 ; i < multiArcStateParametersSizePerArc.size( ) ; i++ )
            {
                std::cout << "arc " << i+1 << " - size multi-arc states: " << multiArcStateParametersSizePerArc.at( i ) << "\n\n";
            }
            for ( auto itr : bodiesToEstimatePerArc )
            {
                std::cout << "arc:" << itr.first+1 <<  "\n\n";
                for ( unsigned int i = 0 ; i < itr.second.size( ); i++ )
                {
                    std::cout << itr.second.at( i ) << " & " ;
                }
                std::cout << "\n\n";
            }


            std::cout << "TEST CENTRAL BODY DEPENDENCIES IN INTEGRATED MULTI-ARC STATES PER ARC - arc " << arc << "\n\n";
            double finalArcEpoch = perArcStateHistory.at( 0 ).rbegin( )->first;
            Eigen::VectorXd finalStates = perArcStateHistory.at( 0 ).rbegin( )->second;
            for ( unsigned int i = 0 ; i < bodiesToPropagatePerArc.at( arc ).size( ) ; i++ )
            {
                std::cout << "BODY: " << bodiesToPropagatePerArc.at( arc ).at( i ) << " - CENTRAL BODY: " << centralBodiesPerArc.at( arc ).at( i ) << "\n\n";
                Eigen::Vector6d currentSolutionFromEphemeris = bodies.at( bodiesToPropagatePerArc.at( arc ).at( i ) )->getEphemeris( )->getCartesianState( finalArcEpoch );
                Eigen::Vector6d stateCentralBody = bodies.at( centralBodiesPerArc.at( arc ).at( i ) )->getEphemeris( )->getCartesianState( finalArcEpoch );
                Eigen::Vector6d differenceStateHistoryWrtEphemeris = currentSolutionFromEphemeris - finalStates.segment( i * 6, 6 );
//                if ( centralBodiesPerArc.at( arc ).at( i )  == bodies.at( bodiesToPropagatePerArc.at( arc ).at( i ) )->getEphemeris( )->getReferenceFrameOrigin( ) )
//                {
//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                            differenceStateHistoryWrtEphemeris, Eigen::Vector6d::Zero( ), 1.0E-16);
//                }
//                else
//                {
//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                            differenceStateHistoryWrtEphemeris, stateCentralBody, 1.0E-12);
//                }
            }





            std::shared_ptr< ArcWiseInitialTranslationalStateParameter< double > > arcWiseTranslationalStateParameter =
                    std::dynamic_pointer_cast< ArcWiseInitialTranslationalStateParameter< double > >(
                            parametersToEstimate->getEstimatedMultiArcInitialStateParameters( ).at( 0 ) );
            std::vector< double > currentArcTimes = arcWiseTranslationalStateParameter->getArcStartTimes( );

            for ( unsigned int i = 0 ; i < currentArcTimes.size( ) ; i++ ) {
                // Multi-arc state parameter limited to current arc.
                std::cout << "arc: " << i+1 << "\n\n";
                std::cout << "test: " << arcWiseTranslationalStateParameter->getParameterName().second.first << "\n\n";
                std::cout << "state: " << arcWiseTranslationalStateParameter->getParameterValue( ).segment( i * 6, 6).transpose( ) << "\n\n";
                std::vector< double > bla1 = { currentArcTimes.at(i) };
                std::vector< std::string > bla2 = {arcWiseTranslationalStateParameter->getCentralBodies().at(i)};
                 std::shared_ptr<ArcWiseInitialTranslationalStateParameter<> > currentArcTranslationalStateParameter =
                        std::make_shared<ArcWiseInitialTranslationalStateParameter< double > >(
                                arcWiseTranslationalStateParameter->getParameterName().second.first,
                                bla1, arcWiseTranslationalStateParameter->getParameterValue().segment( i * 6, 6),
                                bla2, arcWiseTranslationalStateParameter->getFrameOrientation());
            }


//            for (auto itr: multiArcStateHistory.at(arc)) {
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        itr.second, perArcStateHistory.at(0).at(itr.first), 1.0E-16);
//            }
//
//            for (auto itr: variationalEquationsSolution.at(arc).at(0)) {
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        itr.second, perArcVariationalEquationsSolution.at(0).at(0).at(itr.first), 1.0E-16);
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        variationalEquationsSolution.at(arc).at(1).at(itr.first).block(0, 0, 6 * bodiesToPropagatePerArc.at( arc ).size( ), 24),
//                        perArcVariationalEquationsSolution.at(0).at(1).at(itr.first).block(0, 0, 6 * bodiesToPropagatePerArc.at( arc ).size( ), 24), 1.0E-16);
//            }
        }


        ///////////////////////////////////////////////////////////////////////////////////////////
        //// Test retrieved state transition and sensitivity matrices outside arc bounds
        ///////////////////////////////////////////////////////////////////////////////////////////

        Eigen::MatrixXd undefinedCombinedStateTransitionMatrix = multiArcStateTransitionInterface->getCombinedStateTransitionAndSensitivityMatrix(
                ( arcStartTimes[ 1 ] + multiArcEndTimes[ 0 ] ) / 2.0, true, std::vector< std::string >( ) );
        Eigen::MatrixXd undefinedFullCombinedStateTransitionMatrix = multiArcStateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix(
                ( arcStartTimes[ 1 ] + multiArcEndTimes[ 0 ] ) / 2.0, true, std::vector< std::string >( ) );
        std::cout << "undefined combined STM & SEM: " << undefinedCombinedStateTransitionMatrix << "\n\n";
        std::cout << "undefined full combined STM & SEM: " << undefinedFullCombinedStateTransitionMatrix << "\n\n";


    }

}


BOOST_AUTO_TEST_SUITE_END( )

}

}

