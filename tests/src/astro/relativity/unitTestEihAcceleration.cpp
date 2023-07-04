/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//
//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/simulation.h"

//namespace tudat
//{
//
//namespace unit_tests
//{

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::basic_astrodynamics;
using namespace tudat::unit_conversions;

//BOOST_AUTO_TEST_SUITE( test_eih_accelerations )

//BOOST_AUTO_TEST_CASE( testEihPropagation )
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0 * tudat::physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = 25.0 * tudat::physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector<std::string> bodiesToCreate =
        { "Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
    BodyListSettings bodySettings =
        getDefaultBodySettings( bodiesToCreate );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings.at( bodiesToCreate.at( i ) )->gravityFieldSettings = centralGravityFromSpiceSettings( );
    }

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    for( unsigned int test = 0; test < 2; test++ )
    {
        auto accelerationType = ( test == 0 ? einstein_infeld_hoffmann_acceleration : point_mass_gravity );
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector<std::string> bodiesToPropagate = { "Venus" };
        std::vector<std::string> centralBodies = { "SSB" };

        // Define propagation settings.
        for ( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfCurrentBody;
            for ( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
            {
                if ( bodiesToPropagate.at( i ) != bodiesToCreate.at( j ))
                {
                    accelerationsOfCurrentBody[ bodiesToCreate.at( j ) ].push_back(
                        std::make_shared<AccelerationSettings>(
                            accelerationType ));
                }

            }
            accelerationMap[ bodiesToPropagate.at( i ) ] = accelerationsOfCurrentBody;
        }

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToPropagate, centralBodies );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
            bodiesToPropagate, centralBodies, bodies, simulationStartEpoch );

        std::shared_ptr<IntegratorSettings<> >
            integratorSettings = std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 86400.0, rungeKuttaFehlberg78 );
//              std::make_shared<PerElementIntegratorStepSizeControlSettings<double> >( 1.0E-12, 1.0E-12 ),
//              std::make_shared<IntegratorStepSizeValidationSettings>( std::numeric_limits<double>::min( ),
//                                                                      std::numeric_limits<double>::max( ),
//                                                                      set_to_minimum_step_silently ));

        std::shared_ptr<TranslationalStatePropagatorSettings<double> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<double> >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationStartEpoch,
                  integratorSettings,
                  std::make_shared<PropagationTimeTerminationSettings>( simulationEndEpoch ),
                  cowell );
        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 20 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator<> dynamicsSimulator(
            bodies, propagatorSettings );
        std::map<double, Eigen::VectorXd>
            integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        for ( auto it: integrationResult )
        {
            std::cout << it.first << " " << ( it.second - getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodies, it.first ) )
            .transpose( ) << std::endl;
        }
    }


//BOOST_AUTO_TEST_SUITE_END( )

}


