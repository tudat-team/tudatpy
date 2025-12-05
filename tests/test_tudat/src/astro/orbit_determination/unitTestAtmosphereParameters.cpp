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

//! Unit test to check if tidal time lag parameters are estimated correctly
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
    // TODO: Isse here with the aero frame - analytical partial does not respond
    bodies.at( "Vehicle" )
            ->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CREATE ACCELERATION, PARTIAL, PARAMETERS           //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel =
            simulation_setup::createAerodynamicAcceleratioModel( bodies.at( "Vehicle" ), bodies.at( "Titan" ), "Vehicle", "Titan" );


    double testTime0 = 2.205694765899772644e+08;
    Eigen::Vector6d vehicleInitialState0;
    vehicleInitialState0 << -6.566106067480740137e+06, -6.894041365884457715e+06, 1.398789723666163161e+07, 3.329983483136099039e+03, 2.298663138859752053e+03, -3.904458109255926956e+03;

    double testTime1 = 2.205723865899772644e+08;
    Eigen::Vector6d vehicleInitialState1;
    vehicleInitialState1 << 3.142204805038403720e+06, -6.317878925963304937e+04, 2.261319995938756503e+06, 3.190984327192136789e+03, 2.437742656628870009e+03, -4.367135216475068773e+03;


    std::shared_ptr< acceleration_partials::AccelerationPartial > aerodynamicAccelerationPartial =
            createAnalyticalAccelerationPartial( accelerationModel,
                                                 std::make_pair( "Vehicle", bodies.at( "Vehicle" ) ),
                                                 std::make_pair( "Titan", bodies.at( "Titan" ) ),
                                                 bodies );

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


    bodies.at( "Titan" )->setStateFromEphemeris( testTime0 );
    bodies.at( "Titan" )->setCurrentRotationToLocalFrameFromEphemeris( testTime0 );

    bodies.at( "Vehicle" )->setState( vehicleInitialState0 );
    bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( testTime0 );

    accelerationModel->updateMembers( testTime0 );
    std::cout << "acceleration: " << accelerationModel->getUnscaledAcceleration(  ) << std::endl << std::flush;

    std::shared_ptr< AtmosphericFlightConditions > currentAtmFlightConditions = std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Vehicle" )->getFlightConditions( ) );
    std::cout << "Current Altitude: " << currentAtmFlightConditions->getCurrentAltitude( ) << std::endl << std::flush;
    std::cout << "Current Density: " << currentAtmFlightConditions->getCurrentDensity( ) << std::endl << std::flush;

    // Analytically
    Eigen::Vector3d partialWrtBaseDensity = aerodynamicAccelerationPartial->wrtParameter( baseDensityAtmosphericParameter );
    Eigen::Vector3d partialWrtScaleHeight = aerodynamicAccelerationPartial->wrtParameter( scaleHeightAtmosphericParameter );

    // Numerically
    Eigen::Vector3d testPartialWrtBaseDensity =
            acceleration_partials::calculateAccelerationWrtParameterPartials( baseDensityAtmosphericParameter, accelerationModel, 1.0E-8 );
    Eigen::Vector3d testPartialWrtScaleHeight =
            acceleration_partials::calculateAccelerationWrtParameterPartials( scaleHeightAtmosphericParameter, accelerationModel, 10 );


    std::cout << "Base Density:" << std::endl << "\t Analytical: " << std::endl << partialWrtBaseDensity.transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialWrtBaseDensity.transpose(  ) << std::endl;
    std::cout << "Scale Height:" << std::endl << "\t Analytical: " << std::endl << partialWrtScaleHeight.transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialWrtScaleHeight.transpose(  ) << std::endl;


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtBaseDensity, testPartialWrtBaseDensity, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtScaleHeight, testPartialWrtScaleHeight, 1.0E-8 );



    bodies.at( "Titan" )->setStateFromEphemeris( testTime1 );
    bodies.at( "Titan" )->setCurrentRotationToLocalFrameFromEphemeris( testTime1 );

    bodies.at( "Vehicle" )->setState( vehicleInitialState1 );
    bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( testTime1 );

    currentAtmFlightConditions = std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Vehicle" )->getFlightConditions( ) );
    std::cout << "Current Altitude: " << currentAtmFlightConditions->getCurrentAltitude( ) << std::endl << std::flush;
    std::cout << "Current Density: " << currentAtmFlightConditions->getCurrentDensity( ) << std::endl << std::flush;

    accelerationModel->updateMembers( testTime1 );
    accelerationModel->getUnscaledAcceleration(  );
    std::cout << "acceleration: " << accelerationModel->getUnscaledAcceleration(  ) << std::endl << std::flush;


    // Analytically
    partialWrtBaseDensity = aerodynamicAccelerationPartial->wrtParameter( baseDensityAtmosphericParameter );
    partialWrtScaleHeight = aerodynamicAccelerationPartial->wrtParameter( scaleHeightAtmosphericParameter );

    std::function< void( ) > environmentUpdateFunction =
        std::bind( &updateFlightConditionsWithPerturbedState, bodies.at( "Vehicle" )->getFlightConditions( ), 0.0 );

    // Numerically
    testPartialWrtBaseDensity =
            acceleration_partials::calculateAccelerationWrtParameterPartials( baseDensityAtmosphericParameter, accelerationModel, 1.0E-8, environmentUpdateFunction);
    testPartialWrtScaleHeight =
            acceleration_partials::calculateAccelerationWrtParameterPartials( scaleHeightAtmosphericParameter, accelerationModel, 10, environmentUpdateFunction );

    std::cout << "Base Density:" << std::endl << "\t Analytical: " << std::endl << partialWrtBaseDensity.transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialWrtBaseDensity.transpose(  ) << std::endl;
    std::cout << "Scale Height:" << std::endl << "\t Analytical: " << std::endl << partialWrtScaleHeight.transpose(  ) << std::endl << "\t Numerically: " << std::endl << testPartialWrtScaleHeight.transpose(  ) << std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtBaseDensity, testPartialWrtBaseDensity, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtScaleHeight, testPartialWrtScaleHeight, 1.0E-6 );

}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
