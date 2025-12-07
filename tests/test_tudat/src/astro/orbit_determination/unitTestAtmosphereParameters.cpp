/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"

#include <tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h>

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_atmosphere_parameters )

// Using declarations.
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


void updateFlightConditionsWithPerturbedState( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions,
                                                const double timeToUpdate )
{
    flightConditions->resetCurrentTime( );
    flightConditions->updateConditions( timeToUpdate );
}

//! Unit test to check if analytical atmosphere parameters (global and arc-wise) are computed correctly
BOOST_AUTO_TEST_CASE( test_ExponentialAtmosphereParameters )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    using namespace tudat;
    // Create Titan object
    BodyListSettings defaultBodySettings = getDefaultBodySettings( { "Titan" } );
    defaultBodySettings.at( "Titan" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );


    bodies.at( "Titan" )->setAtmosphereModel(
        createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( 7.58e04, 175, 3.0553447e-04 ), "Titan" ));
    std::shared_ptr< aerodynamics::ExponentialAtmosphere > retrievedAtmosphere =
        std::dynamic_pointer_cast< aerodynamics::ExponentialAtmosphere >(bodies.at( "Titan" )->getAtmosphereModel( ) );


    // Create vehicle objects.
    double vehicleMass = 2.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );


    // Create aerodynamic coefficient interface settings.
    double referenceArea = 22.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    referenceArea,
                    aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.01, 0.1 ).finished( ),
                    negative_aerodynamic_frame_coefficients );

    bodies.at( "Vehicle" )
            ->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CREATE ACCELERATION, DEFINE TEST STATES            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel =
            simulation_setup::createAerodynamicAcceleratioModel( bodies.at( "Vehicle" ), bodies.at( "Titan" ), "Vehicle", "Titan" );


    std::vector< double > testTimes;
    std::vector< Eigen::Vector6d > testStates;

    // State 0: well outside the atmosphere (expecting zero partials)
    testTimes.push_back( 2.205694765899772644e+08 );
    testStates.push_back( Eigen::Vector6d{-6.566106067480740137e+06, -6.894041365884457715e+06, 1.398789723666163161e+07, 3.329983483136099039e+03, 2.298663138859752053e+03, -3.904458109255926956e+03} );

    // State 1: inside the atmosphere (equivalent to conditions of T022 Cassini during peak dynamic pressure --> expecting non-zero partials)
    testTimes.push_back( 2.205723865899772644e+08 );
    testStates.push_back( Eigen::Vector6d{3.142204805038403720e+06, -6.317878925963304937e+04, 2.261319995938756503e+06, 3.190984327192136789e+03, 2.437742656628870009e+03, -4.367135216475068773e+03} );

    // State 2: inside the atmosphere (equivalent to conditions of T068 Cassini during peak dynamic pressure --> expecting non-zero partials)
    testTimes.push_back( 3.275979261739959717e+08 );
    testStates.push_back( Eigen::Vector6d{-1.691630151542660315e+06, 1.972704502105712891e+06, -3.004776720845708158e+06, 4.710368198265113278e+03, 3.497139953413709463e+03, -3.547699344472313783e+02} );


    std::shared_ptr< acceleration_partials::AccelerationPartial > aerodynamicAccelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel,
                                                 std::make_pair( "Vehicle", bodies.at( "Vehicle" ) ),
                                                 std::make_pair( "Titan", bodies.at( "Titan" ) ),
                                                 bodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CREATE PARTIAL, PARAMETER OBJECTS (simple)        ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // atmosphere base density parameter
    std::shared_ptr< EstimatableParameter< double > > baseDensityAtmosphericParameter =
            std::make_shared< ExponentialAtmosphereParameter >( std::dynamic_pointer_cast< aerodynamics::ExponentialAtmosphere >(
                                                                        bodies.at( "Titan" )->getAtmosphereModel( ) ),
                                                                        estimatable_parameters::exponential_atmosphere_base_density,
                                                                        "Vehicle");


    // atmosphere scale height parameter
    std::shared_ptr< EstimatableParameter< double > > scaleHeightAtmosphericParameter =
            std::make_shared< ExponentialAtmosphereParameter >( std::dynamic_pointer_cast< aerodynamics::ExponentialAtmosphere >(
                                            bodies.at( "Titan" )->getAtmosphereModel( ) ),
                                                                 estimatable_parameters::exponential_atmosphere_scale_height,
                                                                 "Vehicle");


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       EVALUATE PARTIALS                                  //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< Eigen::Vector3d > partialsWrtBaseDensity;
    std::vector< Eigen::Vector3d > partialsWrtScaleHeight;

    std::vector< Eigen::Vector3d > testPartialsWrtBaseDensity;
    std::vector< Eigen::Vector3d > testPartialsWrtScaleHeight;


    std::function< void( ) > environmentUpdateFunction =
        std::bind( &updateFlightConditionsWithPerturbedState, bodies.at( "Vehicle" )->getFlightConditions( ), 0.0 );


    for (unsigned int i = 0; i < testTimes.size( ); i++)
    {

        bodies.at( "Titan" )->setStateFromEphemeris( testTimes.at(i) );
        bodies.at( "Titan" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes.at(i) );

        bodies.at( "Vehicle" )->setState( testStates.at(i) );
        bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( testTimes.at(i) );

        accelerationModel->updateMembers( testTimes.at(i) );

        // Analytically
        aerodynamicAccelerationPartial->update( testTimes.at(i) );
        partialsWrtBaseDensity.push_back( aerodynamicAccelerationPartial->wrtParameter( baseDensityAtmosphericParameter ) );
        partialsWrtScaleHeight.push_back( aerodynamicAccelerationPartial->wrtParameter( scaleHeightAtmosphericParameter ) );

        // Numerically
        testPartialsWrtBaseDensity.push_back(
                acceleration_partials::calculateAccelerationWrtParameterPartials( baseDensityAtmosphericParameter, accelerationModel, 1.0E-8, environmentUpdateFunction ) );
        testPartialsWrtScaleHeight.push_back(
                acceleration_partials::calculateAccelerationWrtParameterPartials( scaleHeightAtmosphericParameter, accelerationModel, 10, environmentUpdateFunction ) );


        std::cout << "Base Density:" << std::endl << "\t Analytical: " << std::endl << partialsWrtBaseDensity.at( i ).transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialsWrtBaseDensity.at( i ).transpose(  ) << std::endl;
        std::cout << "Scale Height:" << std::endl << "\t Analytical: " << std::endl << partialsWrtScaleHeight.at( i ).transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialsWrtScaleHeight.at( i ).transpose(  ) << std::endl;

    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       TEST RESULTS OF SIMPLE PARAMETERS                    ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (unsigned int i = 0; i < testTimes.size( ); i++)
    {
        if( i == 0)
        {
            BOOST_CHECK_SMALL((partialsWrtBaseDensity.at( i ) - testPartialsWrtBaseDensity.at( i )).norm(), 1e-12);
            BOOST_CHECK_SMALL((partialsWrtScaleHeight.at( i ) - testPartialsWrtScaleHeight.at( i )).norm(), 1e-12);
        }

        else
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialsWrtBaseDensity.at( i ), testPartialsWrtBaseDensity.at( i ), 1.0E-8 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialsWrtScaleHeight.at( i ), testPartialsWrtScaleHeight.at( i ), 1.0E-6 );
        }

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CREATE PARAMETER OBJECTS (arc-wise)               ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< double > arcStartTimes;
    for (unsigned int i = 1; i < testTimes.size( ); i++)
    {
        arcStartTimes.push_back( testTimes.at( i ) - 60*60);
    }

    // atmosphere base density parameter
    std::shared_ptr< ArcWiseExponentialAtmosphereParameter > arcwiseBaseDensityAtmosphericParameter =
            std::make_shared< ArcWiseExponentialAtmosphereParameter >( std::dynamic_pointer_cast< aerodynamics::ExponentialAtmosphere >(
                                                                        bodies.at( "Titan" )->getAtmosphereModel( ) ),
                                                                        estimatable_parameters::arc_wise_exponential_atmosphere_base_density,
                                                                        arcStartTimes,
                                                                        "Vehicle");

    // atmosphere scale height parameter
    std::shared_ptr< ArcWiseExponentialAtmosphereParameter > arcwiseScaleHeightAtmosphericParameter =
            std::make_shared< ArcWiseExponentialAtmosphereParameter >( std::dynamic_pointer_cast< aerodynamics::ExponentialAtmosphere >(
                                                                        bodies.at( "Titan" )->getAtmosphereModel( ) ),
                                                                        estimatable_parameters::arc_wise_exponential_atmosphere_scale_height,
                                                                        arcStartTimes,
                                                                        "Vehicle");



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       EVALUATE PARTIALS                                  //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< Eigen::MatrixXd > partialsWrtArcWiseBaseDensity;
    std::vector< Eigen::MatrixXd > partialsWrtArcWiseScaleHeight;

    for (unsigned int i = 1; i < testTimes.size( ); i++)
    {

        bodies.at( "Titan" )->setStateFromEphemeris( testTimes.at(i) );
        bodies.at( "Titan" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes.at(i) );

        bodies.at( "Vehicle" )->setState( testStates.at(i) );
        bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( testTimes.at(i) );

        accelerationModel->updateMembers( testTimes.at(i) );

        // Analytically
        aerodynamicAccelerationPartial->update( testTimes.at(i) );
        partialsWrtArcWiseBaseDensity.push_back( aerodynamicAccelerationPartial->wrtParameter( arcwiseBaseDensityAtmosphericParameter ) );
        partialsWrtArcWiseScaleHeight.push_back( aerodynamicAccelerationPartial->wrtParameter( arcwiseScaleHeightAtmosphericParameter ) );


        std::cout << "Base Density:" << std::endl << "\t Analytical: " << std::endl << partialsWrtArcWiseBaseDensity.at( i-1 ) << std::endl << "\t Numerically: " << std::endl << testPartialsWrtBaseDensity.at( i ) << std::endl;
        std::cout << "Scale Height:" << std::endl << "\t Analytical: " << std::endl << partialsWrtArcWiseScaleHeight.at( i-1 ) << std::endl << "\t Numerically: " << std::endl << testPartialsWrtScaleHeight.at( i ) << std::endl;

    }

    for (unsigned int i = 0; i < arcStartTimes.size( ); i++)
    {

        for (int k = 0; k < 3; k++) {
            BOOST_CHECK_CLOSE_FRACTION(partialsWrtArcWiseBaseDensity.at( i ).col( i )(k), testPartialsWrtBaseDensity.at( i+1 )(k), 1.0E-8);
            BOOST_CHECK_CLOSE_FRACTION(partialsWrtArcWiseScaleHeight.at( i ).col( i )(k), testPartialsWrtScaleHeight.at( i+1 )(k), 1.0E-6);

        }

    }

}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
