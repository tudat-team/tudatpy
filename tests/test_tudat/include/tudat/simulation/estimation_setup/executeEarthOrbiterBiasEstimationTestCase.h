/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef EXECUTEEARTHORBITERBIASESTIMATIONTESTCASE_H
#define EXECUTEEARTHORBITERBIASESTIMATIONTESTCASE_H

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCaseUtilities.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{

using namespace basic_astrodynamics;
using namespace simulation_setup;
using namespace observation_models;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace estimatable_parameters;
using namespace propagators;
using namespace numerical_integrators;
using namespace ephemerides;
using namespace coordinate_conversions;

namespace unit_tests
{

using namespace ephemerides;
using namespace basic_astrodynamics;
using namespace simulation_setup;
using namespace orbital_element_conversions;
using namespace coordinate_conversions;
using namespace physical_constants;

template< typename TimeType = double, typename StateScalarType = double >
std::pair< Eigen::VectorXd, bool > executeEarthOrbiterBiasEstimation( const bool estimateRangeBiases = true,
                                                                      const bool estimateTwoWayBiases = false,
                                                                      const bool useSingleBiasModel = true,
                                                                      const bool estimateAbsoluteBiases = true,
                                                                      const bool omitRangeData = false,
                                                                      const bool useMultiArcBiases = false,
                                                                      const bool estimateTimeBiases = false )
{
    const int numberOfDaysOfData = 1;

    int numberOfIterations = 3;
    if( estimateRangeBiases && estimateTimeBiases )
    {
        numberOfIterations = 6;
    }

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    // Specify initial time
    TimeType initialEphemerisTime = 1.0E7;
    TimeType finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    std::vector< double > biasArcs;
    biasArcs.push_back( initialEphemerisTime );
    biasArcs.push_back( initialEphemerisTime + 4.0 * 3600.0 );
    biasArcs.push_back( initialEphemerisTime + 12.0 * 3600.0 );

    std::vector< Eigen::VectorXd > biasPerArc;
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Matrix< StateScalarType, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, TimeType( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >(
                    TimeType( initialEphemerisTime ), 120.0, CoefficientSets::rungeKuttaFehlberg78, 120.0, 120.0, 1.0, 1.0 );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::vector< LinkEnds > stationTwoWayLinkEnds;
    std::vector< LinkEnds > stationTwoWayInverseLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ reflector1 ] = LinkEndId( "Vehicle", "" );
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        stationTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        linkEnds[ reflector1 ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        stationTwoWayInverseLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    if( estimateRangeBiases )
    {
        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 1 ] );

        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayLinkEnds[ 0 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayInverseLinkEnds[ 0 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayLinkEnds[ 1 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayInverseLinkEnds[ 1 ] );
    }

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 1 ] );

    //    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayLinkEnds[ 0 ] );
    //    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayInverseLinkEnds[ 0 ] );
    //    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayLinkEnds[ 1 ] );
    //    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayInverseLinkEnds[ 1 ] );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
            "Vehicle", systemInitialState, "Earth" ) );

    if( useMultiArcBiases )
    {
        if( estimateRangeBiases )
        {
            if( !estimateTimeBiases )
            {
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            }
            else
            {
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, biasArcs, transmitter, biasArcs ) );
            }

            if( estimateTwoWayBiases )
            {
                if( !estimateTimeBiases )
                {
                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                }
                else
                {
                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, biasArcs, transmitter, biasArcs ) );
                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, biasArcs, transmitter, biasArcs ) );
                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, biasArcs, transmitter, biasArcs ) );
                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, biasArcs, transmitter, biasArcs ) );
                }
            }
        }
        else
        {
            if( !estimateTimeBiases )
            {
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 0 ),
                        one_way_doppler,
                        biasArcs,
                        transmitter,
                        estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 1 ),
                        one_way_doppler,
                        biasArcs,
                        transmitter,
                        estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 2 ),
                        one_way_doppler,
                        biasArcs,
                        transmitter,
                        estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 3 ),
                        one_way_doppler,
                        biasArcs,
                        transmitter,
                        estimateAbsoluteBiases ) );
            }
            else
            {
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 1 ), one_way_doppler, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 2 ), one_way_doppler, biasArcs, transmitter, biasArcs ) );
                parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 3 ), one_way_doppler, biasArcs, transmitter, biasArcs ) );
            }

            if( estimateTwoWayBiases )
            {
                //                if ( !estimateTimeBiases )
                //                {
                //                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings
                //                    >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, biasArcs, transmitter,
                //                            estimateAbsoluteBiases ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings
                //                    >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, biasArcs, transmitter,
                //                            estimateAbsoluteBiases ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings
                //                    >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, biasArcs, transmitter,
                //                            estimateAbsoluteBiases ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings
                //                    >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, biasArcs, transmitter,
                //                            estimateAbsoluteBiases ) );
                //                }
                //                else
                //                {
                //                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, biasArcs, transmitter,
                //                            biasArcs ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, biasArcs, transmitter,
                //                            biasArcs ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, biasArcs, transmitter,
                //                            biasArcs ) );
                //                    parameterNames.push_back( std::make_shared< ArcWiseTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, biasArcs, transmitter,
                //                            biasArcs ) );
                //                }
            }
        }
    }
    else
    {
        if( estimateRangeBiases )
        {
            if( !estimateTimeBiases )
            {
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, estimateAbsoluteBiases ) );
            }
            else
            {
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, transmitter, initialEphemerisTime ) );
            }

            if( estimateTwoWayBiases )
            {
                if( !estimateTimeBiases )
                {
                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, estimateAbsoluteBiases ) );
                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, estimateAbsoluteBiases ) );
                }
                else
                {
                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, transmitter, initialEphemerisTime ) );
                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, transmitter, initialEphemerisTime ) );
                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, transmitter, initialEphemerisTime ) );
                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                            linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, transmitter, initialEphemerisTime ) );
                }
            }
        }
        else
        {
            if( !estimateTimeBiases )
            {
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 1 ), one_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 2 ), one_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 3 ), one_way_doppler, estimateAbsoluteBiases ) );
            }
            else
            {
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 1 ), one_way_doppler, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 2 ), one_way_doppler, transmitter, initialEphemerisTime ) );
                parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                        linkEndsPerObservable.at( one_way_doppler ).at( 3 ), one_way_doppler, transmitter, initialEphemerisTime ) );
            }

            if( estimateTwoWayBiases )
            {
                //                if ( !estimateTimeBiases )
                //                {
                //                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, estimateAbsoluteBiases )
                //                            );
                //                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, estimateAbsoluteBiases )
                //                            );
                //                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, estimateAbsoluteBiases )
                //                            );
                //                    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, estimateAbsoluteBiases )
                //                            );
                //                }
                //                else
                //                {
                //                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, transmitter,
                //                            initialEphemerisTime ) );
                //                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, transmitter,
                //                            initialEphemerisTime ) );
                //                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, transmitter,
                //                            initialEphemerisTime ) );
                //                    parameterNames.push_back( std::make_shared< ConstantTimeDriftBiasEstimatableParameterSettings >(
                //                            linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, transmitter,
                //                            initialEphemerisTime ) );
                //                }
            }
        }
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( useMultiArcBiases )
            {
                if( useSingleBiasModel )
                {
                    if( estimateTimeBiases )
                    {
                        biasSettings = std::make_shared< ArcWiseTimeDriftBiasSettings >( biasArcs, biasPerArc, transmitter, biasArcs );
                    }
                    else
                    {
                        if( estimateAbsoluteBiases )
                        {
                            biasSettings =
                                    std::make_shared< ArcWiseConstantObservationBiasSettings >( biasArcs, biasPerArc, transmitter, true );
                        }
                        else
                        {
                            biasSettings =
                                    std::make_shared< ArcWiseConstantObservationBiasSettings >( biasArcs, biasPerArc, transmitter, false );
                        }
                    }
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                    biasSettingsList.push_back(
                            std::make_shared< ArcWiseConstantObservationBiasSettings >( biasArcs, biasPerArc, transmitter, true ) );
                    biasSettingsList.push_back(
                            std::make_shared< ArcWiseConstantObservationBiasSettings >( biasArcs, biasPerArc, transmitter, false ) );
                    biasSettingsList.push_back(
                            std::make_shared< ArcWiseTimeDriftBiasSettings >( biasArcs, biasPerArc, transmitter, biasArcs ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
            }
            else
            {
                if( useSingleBiasModel )
                {
                    if( estimateTimeBiases )
                    {
                        biasSettings = std::make_shared< ConstantTimeDriftBiasSettings >(
                                Eigen::Vector1d::Zero( ), transmitter, initialEphemerisTime );
                    }
                    else
                    {
                        if( estimateAbsoluteBiases )
                        {
                            biasSettings = std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), true );
                        }
                        else
                        {
                            biasSettings = std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), false );
                        }
                    }
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;

                    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), true ) );
                    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), false ) );
                    biasSettingsList.push_back( std::make_shared< ConstantTimeDriftBiasSettings >(
                            Eigen::Vector1d::Zero( ), transmitter, initialEphemerisTime ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
            }
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                    bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< TimeType > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 600.0;
    double observationInterval = 600.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 100; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::cout << "Final time: " << baseTimeList.at( baseTimeList.size( ) - 1 ) << std::endl;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        if( !( omitRangeData && ( ( currentObservable == one_way_range ) || ( currentObservable == n_way_range ) ) ) )
        {
            std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                        currentObservable, currentLinkEndsList.at( i ), baseTimeList, receiver ) );
            }
        }
    }

    // Simulate observations
    std::shared_ptr< ObservationCollection< StateScalarType, TimeType > > simulatedObservations =
            simulateObservations< StateScalarType, TimeType >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
    weightsPerObservationParser[ observationParser( n_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
    weightsPerObservationParser[ observationParser( one_way_doppler ) ] =
            1.0 / ( 1.0E-12 * 1.0E-12 * physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    // Perturb parameter estimate
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    if( numberOfIterations > 0 )
    {
        parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

        for( unsigned int i = 6; i < initialParameterEstimate.rows( ); i++ )
        {
            if( estimateTimeBiases )
            {
                parameterPerturbation( i ) = 1.0e-7;
            }
            else
            {
                if( estimateAbsoluteBiases )
                {
                    if( estimateRangeBiases )
                    {
                        parameterPerturbation( i ) = 1.0;
                    }
                    else
                    {
                        parameterPerturbation( i ) = 1.0E-8;
                    }
                }
                else
                {
                    parameterPerturbation( i ) = 1.0E-6;
                }
            }
        }
        initialParameterEstimate += parameterPerturbation;
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define estimation input
    std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType > >( simulatedObservations );

    estimationInput->defineEstimationSettings( true, false, false, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< StateScalarType > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );

    Eigen::VectorXd estimationError =
            estimationOutput->parameterHistory_.at( estimationOutput->parameterHistory_.size( ) - 1 ) - truthParameters;

    std::cout << "initial error: " << ( parameterPerturbation ).transpose( ) << std::endl << std::endl;
    std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl << std::endl;

    return std::make_pair(
            estimationError,
            ( estimationOutput->exceptionDuringInversion_ ||
              !( estimationOutput->getUnnormalizedCovarianceMatrix( ) == estimationOutput->getUnnormalizedCovarianceMatrix( ) ) ) );
}

}  // namespace unit_tests

}  // namespace tudat

#endif  // EXECUTEEARTHORBITERBIASESTIMATIONTESTCASE_H
