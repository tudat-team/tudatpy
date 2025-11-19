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

//   bodySettings.at( "JUICE" )->ephemerisSettings = directSpiceEphemerisSettings( "Jupiter", "ECLIPJ2000", "-28" );
//   bodySettings.at( "JUICE" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

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
                                                                     const std::string& multiArcBodyToPropagate,
                                                                     const std::string& multiArcCentralBody )
{
   using namespace tudat::basic_astrodynamics;
   using namespace tudat::ephemerides;

   SelectedAccelerationMap accelerationSettingsJuice;
   accelerationSettingsJuice[ "JUICE" ][ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
   accelerationSettingsJuice[ "JUICE" ][ multiArcCentralBody ].push_back(
           std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
   accelerationSettingsJuice[ "JUICE" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
   
   return createAccelerationModelsMap( bodies, accelerationSettingsJuice, { multiArcBodyToPropagate }, { multiArcCentralBody } );
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
   spice_interface::loadSpiceKernelInTudat( "/home/dominic/Downloads/juice_mat_crema_5_1_150lb_v01.bsp" );

   // Simulation parameters
   double initialEpoch = 946728000.0;
   double finalEpoch = 1200000000.0;

   double propagationTimeStep = 1800.0;


   // Create body map
   SystemOfBodies bodies = createBodies( );


   std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagationSettingsList;
   for( unsigned int i = 0; i < flybysBodies.size( ); i++ )
   {
       propagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
               std::vector< std::string >( { flybysBodies.at( i ) } ),
               getMultiArcAccelerationModelMap( bodies, "JUICE", flybysBodies.at( i ) ),
               std::vector< std::string >( { "JUICE" } ),
               spice_interface::getBodyCartesianStateAtEpoch( "JUICE", flybysBodies.at( i ), "ECLIPJ2000", "NONE", flybysTimes.at( i ) ),
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
                   galileanSatelliteNames,
                   getMoonsAccelerationMap( bodies, galileanSatelliteNames, galileanSatelliteCentralNames ),
                   galileanSatelliteCentralNames,
                   initialStates,
                   flybysTimes.at( 0 ),
                   rungeKuttaFixedStepSettings( 50.0, numerical_integrators::CoefficientSets::rungeKutta4Classic ),
                   std::make_shared< PropagationTimeTerminationSettings >( flybysTimes.at( 22 ) ) );

   std::shared_ptr< HybridArcPropagatorSettings< double > > hybridArcPropagatorSettings =
           std::make_shared< HybridArcPropagatorSettings< double > >(
                   singleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< double > >(
                                                                           propagationSettingsList ) );


   std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > hybridArcParameters = getParametersToEstimate(
           hybridArcPropagatorSettings, bodies );

   std::cout<<"Propagator settings "<<hybridArcPropagatorSettings->getInitialStates( ).transpose( )<<std::endl<<std::endl;
   std::cout<<"Parameters "<<hybridArcParameters->getFullParameterValues< double >( ).transpose( )<<std::endl<<std::endl;


}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
