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

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{

namespace unit_tests
{

// Using declarations.
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

BOOST_AUTO_TEST_SUITE( test_hybrid_arc_multi_body_variational_equations_calculation )

static const std::vector< std::string > galileanSatelliteNames = { "Io", "Europa", "Ganymede", "Callisto" };
static const std::vector< std::string > galileanSatelliteCentralNames = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

static const std::map< std::string, double > satelliteOrbitalPeriods = { { "Io", 1.769 * 86400.0 },
                                                                        { "Europa", 3.551 * 86400.0 },
                                                                        { "Ganymede", 7.155 * 86400.0 },
                                                                        { "Callisto", 16.689 * 86400.0 } };
static const double centralBodyRotationPeriod = 9.8 * 3600.0;
static const double ganymedeEllipticalInsertionTime = 32.71 * physical_constants::JULIAN_YEAR;
static const double ganymedeSphericalInsertionTime = 33.125 * physical_constants::JULIAN_YEAR;
static const double ganymedeSphericalEndTime = 33.45 * physical_constants::JULIAN_YEAR;

SystemOfBodies createBodies(  )
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
   BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );

   bodySettings.at( "Io" )->rotationModelSettings =
           std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "IAU_Io" );

   bodySettings.at( "Europa" )->rotationModelSettings =
           std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "IAU_Europa" );

   bodySettings.at( "Ganymede" )->rotationModelSettings =
           std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "IAU_Ganymede" );

   bodySettings.at( "Callisto" )->rotationModelSettings =
           std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "IAU_Callisto" );

   // Define settings for JUICE spacecraft
   bodySettings.addSettings( "JUICE" );
   bodySettings.addSettings( "JUICE2" );

   bodySettings.at( "Io" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
   bodySettings.at( "Europa" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
   bodySettings.at( "Ganymede" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );
   bodySettings.at( "Callisto" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );

   // Create body map
   SystemOfBodies bodies = createSystemOfBodies( bodySettings );

   return bodies;
}

static const std::vector< std::string > flybysBodies = { "Ganymede", "Ganymede", "Ganymede", "Ganymede", "Europa",   "Europa",
                                                        "Callisto", "Ganymede", "Callisto", "Callisto", "Callisto", "Callisto",
                                                        "Callisto", "Callisto", "Callisto", "Callisto", "Callisto", "Ganymede",
                                                        "Ganymede", "Ganymede", "Callisto", "Callisto", "Ganymede" };
static const std::vector< double > flybysTimes = { 957842194.921875,  962785494.140625,  966094442.578125, 967948681.640625,
                                                  969142826.953125,  970373478.515625,  971388256.640625, 973245765.234375,
                                                  976662144.140625,  978104000.390625,  979545181.640625, 980986391.015625,
                                                  988195320.703125,  989636516.015625,  991077999.609375, 992519462.109375,
                                                  993959258.203125,  995547638.671875,  998637697.265625, 1000051238.671875,
                                                  1001522351.953125, 1006640363.671875, 1018340314.453125 };


basic_astrodynamics::AccelerationMap getMultiArcAccelerationModelMap( const SystemOfBodies& bodies,
                                                                      const std::vector< std::string >& bodyNames,
                                                                      const std::vector< std::string >& centralBodyNames )
{
   using namespace tudat::basic_astrodynamics;
   using namespace tudat::ephemerides;

   SelectedAccelerationMap accelerationSettingsJuice;
   for( unsigned int i = 0; i < bodyNames.size( ); i++ )
   {
       accelerationSettingsJuice[ bodyNames.at( i ) ][ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
       accelerationSettingsJuice[ bodyNames.at( i ) ][ centralBodyNames.at( i ) ].push_back(
               std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
       accelerationSettingsJuice[ bodyNames.at( i ) ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
   }
   
   return createAccelerationModelsMap( bodies, accelerationSettingsJuice, bodyNames, centralBodyNames);
}

basic_astrodynamics::AccelerationMap getMoonsAccelerationMap( const SystemOfBodies& bodies,
                                                             const std::vector< std::string >& bodiesToPropagate,
                                                             const std::vector< std::string >& centralBodies )
{
   SelectedAccelerationMap accelerationSettingsMoons;
   for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
   {
       accelerationSettingsMoons[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back(
               std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                       2, 2, 2, 2 ) );
   }

   basic_astrodynamics::AccelerationMap accelerationMap =
           createAccelerationModelsMap( bodies, accelerationSettingsMoons, bodiesToPropagate, centralBodies );

   return accelerationMap;
}

std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > getParametersToEstimate(
       const std::shared_ptr< propagators::HybridArcPropagatorSettings<> > hybridArcPropagatorSettings,
       const SystemOfBodies& bodies )
{
   using namespace tudat::propagators;
   using namespace tudat::estimatable_parameters;
   using namespace tudat::basic_astrodynamics;
   using namespace tudat::observation_models;

   std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = getInitialStateParameterSettings< double, double >(
           hybridArcPropagatorSettings, bodies );

   for( unsigned int i = 0; i < galileanSatelliteNames.size( ); i++ )
   {
       parameterNames.push_back(
               std::make_shared< EstimatableParameterSettings >( galileanSatelliteNames.at( i ), gravitational_parameter ) );
       parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
               2, 0, 2, 2, galileanSatelliteNames.at( i ), spherical_harmonics_cosine_coefficient_block ) );
       parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
               2, 1, 2, 2, galileanSatelliteNames.at( i ), spherical_harmonics_sine_coefficient_block ) );
   }

   std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
           createParametersToEstimate< double >( parameterNames, bodies, hybridArcPropagatorSettings );
   printEstimatableParameterEntries( parametersToEstimate );

   return parametersToEstimate;
}


BOOST_AUTO_TEST_CASE( testPropagatorParameterConsistency )
{
   // Load spice kernels
   spice_interface::loadStandardSpiceKernels( );
   spice_interface::loadSpiceKernelInTudat( tudat::paths::getTudatTestDataPath( ) + "juice_mat_crema_5_1_150lb_v01.bsp" );

   // Simulation parameters
   double initialEpoch = 946728000.0;
   double finalEpoch = 1200000000.0;

   double propagationTimeStep = 1800.0;


   // Create body map
   SystemOfBodies bodies = createBodies( );


   std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagationSettingsList;
   for( unsigned int i = 0; i < flybysBodies.size( ); i++ )
   {
       std::vector< std::string > propagatedBodies;
       std::vector< std::string > centralBodies;

       propagatedBodies.push_back( "JUICE" );
       centralBodies.push_back( flybysBodies.at( i ) );

       if( i % 2 == 0 )
       {
           propagatedBodies.push_back( "JUICE2" );
           centralBodies.push_back( flybysBodies.at( i % 2 ) );
       }

       Eigen::VectorXd initialStates = Eigen::VectorXd::Zero( 6 * propagatedBodies.size( ) );

       for( unsigned int j = 0; j < propagatedBodies.size( ); j++ )
       {
           initialStates.segment( 6 * j, 6 ) =
                   spice_interface::getBodyCartesianStateAtEpoch( "JUICE", centralBodies.at( j ), "ECLIPJ2000", "NONE", flybysTimes.at( i ) ) +
                   static_cast< double >( j ) * Eigen::Vector6d::Ones( ) * 12.5;
       }

       propagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
               centralBodies,
               getMultiArcAccelerationModelMap( bodies, propagatedBodies, centralBodies ),
               propagatedBodies,
               initialStates,
               flybysTimes.at( i ),
               rungeKuttaFixedStepSettings( 50.0, numerical_integrators::CoefficientSets::rungeKutta4Classic ),
               std::make_shared< PropagationTimeTerminationSettings >( flybysTimes.at( i ) + 3600.0 ) ) );
   }

   Eigen::VectorXd initialStates = Eigen::VectorXd::Zero( 24 );
   for( int i = 0; i < 4; i++ )
   {
       initialStates.segment( 6 * i, 6 ) =
               spice_interface::getBodyCartesianStateAtEpoch( galileanSatelliteNames.at( i ), "Jupiter", "ECLIPJ2000", "NONE", flybysTimes.at( 0 ) );

   }
   std::shared_ptr< SingleArcPropagatorSettings< double > > singleArcPropagatorSettings =
           std::make_shared< TranslationalStatePropagatorSettings< double > >(
                   galileanSatelliteCentralNames,
                   getMoonsAccelerationMap( bodies, galileanSatelliteNames, galileanSatelliteCentralNames ),
                   galileanSatelliteNames,
                   initialStates,
                   flybysTimes.at( 0 ),
                   rungeKuttaFixedStepSettings( 50.0, numerical_integrators::CoefficientSets::rungeKutta4Classic ),
                   std::make_shared< PropagationTimeTerminationSettings >( flybysTimes.at( 22 ) ) );

   std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
           std::make_shared< MultiArcPropagatorSettings< double > >( propagationSettingsList );
   std::shared_ptr< HybridArcPropagatorSettings< double > > hybridArcPropagatorSettings =
           std::make_shared< HybridArcPropagatorSettings< double > >(
                   singleArcPropagatorSettings, multiArcPropagatorSettings );


   std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > hybridArcParameters = getParametersToEstimate(
           hybridArcPropagatorSettings, bodies );
   std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > singleArcParameters =
           hybridArcParameters->getInitialSingleArcStateParameters( );
   std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > multiArcParameters =
           hybridArcParameters->getInitialMultiArcStateParameters( );
   printEstimatableParameterEntries( hybridArcParameters );

   Eigen::Vector6d statePerturbation = Eigen::Vector6d::Zero( );
   statePerturbation( 0 ) = 10.0E3;
   statePerturbation( 4 ) = -0.53E3;
   
   for( int singleArcTest = 0; singleArcTest < 3; singleArcTest++ )
   {
       Eigen::VectorXd originalInitialStates = singleArcPropagatorSettings->getInitialStates( );
       for( int singleArcBody = 0; singleArcBody < 4; singleArcBody++ )
       {
            int startIndex = 6 * singleArcBody;
            if( singleArcTest == 0 )
            {
                Eigen::VectorXd nominalValue = singleArcParameters[ startIndex ]->getParameterValue( );
                Eigen::VectorXd perturbedValue = nominalValue + statePerturbation;
                singleArcParameters[ startIndex ]->setParameterValue( perturbedValue );
            }
            else if( singleArcTest == 1 )
            {
                Eigen::VectorXd nominalValue = hybridArcPropagatorSettings->getInitialStates( );
                Eigen::VectorXd perturbedValue = nominalValue;
                perturbedValue.segment( startIndex, 6 ) = perturbedValue.segment( startIndex, 6 ) + statePerturbation;
                hybridArcPropagatorSettings->resetInitialStates( perturbedValue );
            }
            else if( singleArcTest == 2 )
            {
                Eigen::VectorXd nominalValue = singleArcPropagatorSettings->getInitialStates( );
                Eigen::VectorXd perturbedValue = nominalValue;
                perturbedValue.segment( startIndex, 6 ) = perturbedValue.segment( startIndex, 6 ) + statePerturbation;
                singleArcPropagatorSettings->resetInitialStates( perturbedValue );
                hybridArcPropagatorSettings->updateInitialState( );
            }

            Eigen::Vector6d stateDifference = singleArcPropagatorSettings->getInitialStates( ).segment( startIndex, 6 ) - originalInitialStates.segment( startIndex, 6 ) ;
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    stateDifference,
                    statePerturbation,
                    ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( singleArcPropagatorSettings->getInitialStates( ).segment( 0, 24 ) ),
                    ( hybridArcPropagatorSettings->getInitialStates( ).segment( 0, 24 ) ),
                    ( std::numeric_limits< double >::epsilon( ) ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( singleArcPropagatorSettings->getInitialStates( ).segment( 0, 24 ) ),
                    ( hybridArcParameters->getFullParameterValues< double >( ).segment( 0, 24 ) ),
                    ( std::numeric_limits< double >::epsilon( ) ) );
       }
   }

   std::vector< int > arcStartIndex;
   int currentIndex = hybridArcPropagatorSettings->getSingleArcPropagatorSettings( )->getConventionalStateSize( );
   for( unsigned int i = 0; i < hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).size( ); i++ )
   {
       arcStartIndex.push_back( currentIndex );
       currentIndex += hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->getConventionalStateSize( );
   }

   for( int multiArcTest = 0; multiArcTest < 4; multiArcTest++ )
   {
       Eigen::VectorXd originalInitialStates = multiArcPropagatorSettings->getInitialStates( );
       for( int bodyIndex = 0; bodyIndex < 2; bodyIndex++ )
       {
           auto it = multiArcParameters.begin( );
           if( bodyIndex == 1 )
           {
               it++;
           }
           std::shared_ptr< EstimatableParameter< Eigen::Matrix< double, Eigen::Dynamic, 1 > > > currentParameter = it->second;
           for( int arcIndex = 0; arcIndex < 4; arcIndex++ )
           {
               std::cout<<multiArcTest<<" "<<bodyIndex<<" "<<arcIndex<<std::endl;
               int parameterStartIndex = arcIndex * 6;
               int propagatorArcIndex = ( bodyIndex == 0 ? 1 : 2 ) * arcIndex;
               int singleArcPropagatorStartIndex = 6 * bodyIndex;
               int hybridArcPropagatorStartIndex = arcStartIndex.at( propagatorArcIndex ) + singleArcPropagatorStartIndex;

               Eigen::VectorXd originalValueFromParameter = currentParameter->getParameterValue( ).segment( parameterStartIndex, 6 );
               Eigen::VectorXd originalValueFromSingleArc = hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( propagatorArcIndex )->getInitialStates( ).segment( singleArcPropagatorStartIndex, 6 );
               Eigen::VectorXd originalValueFromHybridArc = hybridArcPropagatorSettings->getInitialStates( ).segment( hybridArcPropagatorStartIndex, 6 );
               Eigen::VectorXd originalValueFromMultiArc = multiArcPropagatorSettings->getInitialStates( ).segment( hybridArcPropagatorStartIndex - 24, 6 );

              TUDAT_CHECK_MATRIX_CLOSE_FRACTION( originalValueFromParameter, originalValueFromSingleArc,( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
              TUDAT_CHECK_MATRIX_CLOSE_FRACTION( originalValueFromParameter, originalValueFromMultiArc,( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
              TUDAT_CHECK_MATRIX_CLOSE_FRACTION( originalValueFromParameter, originalValueFromHybridArc,( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

                if( multiArcTest == 0 )
               {
                   Eigen::VectorXd perturbedParameterValue = currentParameter->getParameterValue( );
                   perturbedParameterValue.segment( parameterStartIndex, 6 ) += statePerturbation;
                   currentParameter->setParameterValue( perturbedParameterValue );
               }
               else if( multiArcTest == 1 )
               {
                   Eigen::VectorXd perturbedSingleArcValue = hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( propagatorArcIndex )->getInitialStates( );
                   perturbedSingleArcValue.segment( singleArcPropagatorStartIndex, 6 ) += statePerturbation;
                   hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( propagatorArcIndex )->resetInitialStates( perturbedSingleArcValue );
               }
               else if( multiArcTest == 2 )
               {
                   Eigen::VectorXd perturbedMultiArcValue = hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getInitialStates( );
                   perturbedMultiArcValue.segment( hybridArcPropagatorStartIndex - 24, 6 ) += statePerturbation;
                   hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->resetInitialStates( perturbedMultiArcValue );
               }
               else if( multiArcTest == 3 )
               {
                   Eigen::VectorXd perturbedHybridArcValue = hybridArcPropagatorSettings->getInitialStates( );
                   perturbedHybridArcValue.segment( hybridArcPropagatorStartIndex, 6 ) += statePerturbation;
                   hybridArcPropagatorSettings->resetInitialStates( perturbedHybridArcValue );
               }

               Eigen::VectorXd perturbedValueFromParameter = currentParameter->getParameterValue( ).segment( parameterStartIndex, 6 );
               Eigen::VectorXd perturbedValueFromSingleArc = hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( propagatorArcIndex )->getInitialStates( ).segment( singleArcPropagatorStartIndex, 6 );
               Eigen::VectorXd perturbedValueFromHybridArc = hybridArcPropagatorSettings->getInitialStates( ).segment( hybridArcPropagatorStartIndex, 6 );
               Eigen::VectorXd perturbedValueFromMultiArc = multiArcPropagatorSettings->getInitialStates( ).segment( hybridArcPropagatorStartIndex - 24, 6 );

               std::cout<<( perturbedValueFromParameter - originalValueFromParameter ).transpose( )<<std::endl;
               std::cout<<( perturbedValueFromSingleArc - originalValueFromSingleArc ).transpose( )<<std::endl;
               std::cout<<( perturbedValueFromHybridArc - originalValueFromHybridArc ).transpose( )<<std::endl;
               std::cout<<( perturbedValueFromMultiArc - originalValueFromMultiArc ).transpose( )<<std::endl<<std::endl<<std::endl;

               Eigen::Vector6d stateDifference = perturbedValueFromParameter - originalValueFromParameter;
               TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                       stateDifference,
                       statePerturbation,
                       ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
               TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                       perturbedValueFromParameter,
                       perturbedValueFromSingleArc,
                       ( std::numeric_limits< double >::epsilon( ) ) );
               TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                       perturbedValueFromParameter,
                       perturbedValueFromMultiArc,
                       ( std::numeric_limits< double >::epsilon( ) ) );
               TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                       perturbedValueFromParameter,
                       perturbedValueFromHybridArc,
                       ( std::numeric_limits< double >::epsilon( ) ) );
           }
       }
   }


}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
