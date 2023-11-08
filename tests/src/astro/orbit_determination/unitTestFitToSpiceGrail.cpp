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

#include <string>
#include <thread>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
//
//
//namespace tudat
//{
//namespace unit_tests
//{
//BOOST_AUTO_TEST_SUITE( test_fit_to_spice )

//Using declarations.
using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
//
//BOOST_AUTO_TEST_CASE( test_FitToSpice )
//{
int main( )
{

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_de440_200625.tf");
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_pa_de440_200625.bpc");

    double initialTimeEnvironment = DateTime( 2012, 03, 02, 0, 0, 0.0 ).epoch< double >();
    double finalTimeEnvironment = DateTime( 2012, 05, 29, 0, 0, 0.0 ).epoch< double >();

    DateTime initialDateTime = basic_astrodynamics::getCalendarDateFromTime< double >( initialTimeEnvironment );
    DateTime finalDateTime = basic_astrodynamics::getCalendarDateFromTime< double >( finalTimeEnvironment );

    double initialPropagationTime = DateTime( 2012, 04, 07, 0, 0, 0.0 ).epoch< double >();
    double finalPropagtionTime = initialPropagationTime + 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );

    BodyListSettings bodySettings =
        getDefaultBodySettings( bodyNames, "Mars" );
    bodySettings.at( "Moon" )->rotationModelSettings = spiceRotationModelSettings(
        bodySettings.getFrameOrientation( ), "MOON_PA_DE440", "MOON_PA_DE440" );
    bodySettings.at( "Moon" )->gravityFieldSettings =
        std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( gggrx1200, 500 );
    std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
        bodySettings.at( "Moon" )->gravityFieldSettings )->resetAssociatedReferenceFrame( "MOON_PA_DE440" );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back(
        fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Earth", 0.02405, 2 ) );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back(
        fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Sun", 0.02405, 2 ) );

    bodySettings.addSettings( "GRAIL-A" );
    bodySettings.at( "GRAIL-A" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
        initialTimeEnvironment + 3600.0, finalTimeEnvironment - 3600.0, 30.0, "Moon" );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create radiation pressure settings
    double referenceAreaRadiation = 16.0;
    double radiationPressureCoefficient = 1.5;
    std::vector<std::string> occultingBodies;
    occultingBodies.push_back( "Moon" );
    std::shared_ptr<RadiationPressureInterfaceSettings> grailRadiationPressureSettings =
        std::make_shared<CannonBallRadiationPressureInterfaceSettings>(
            "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "GRAIL-A" )->setRadiationPressureInterface(
        "Sun", createRadiationPressureInterface(
            grailRadiationPressureSettings, "GRAIL-A", bodies ));
    bodies.at( "GRAIL-A" )->setConstantBodyMass( 2000.0 );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfSpacecraft;
    accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( cannon_ball_radiation_pressure ));
    accelerationsOfSpacecraft[ "Earth" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<SphericalHarmonicAccelerationSettings>( 256, 256 ));
    accelerationsOfSpacecraft[ "Moon" ].push_back( empiricalAcceleration( ));
    accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Jupiter" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Saturn" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));

    accelerationMap[ "GRAIL-A" ] =accelerationsOfSpacecraft;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector<std::string> bodiesToEstimate = { "GRAIL-A" };
    std::vector<std::string> centralBodies = { "Moon" };

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToEstimate, centralBodies );

    std::vector<std::shared_ptr<estimatable_parameters::EstimatableParameterSettings> > additionalParameterNames;

    additionalParameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ));

    std::map<basic_astrodynamics::EmpiricalAccelerationComponents,
        std::vector<basic_astrodynamics::EmpiricalAccelerationFunctionalShapes> > empiricalComponentsToEstimate;
//    empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( constant_empirical );
//    empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( sine_empirical );
//    empiricalComponentsToEstimate[ radial_empirical_acceleration_component ].push_back( cosine_empirical );
//
//    empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
//    empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( sine_empirical );
//    empiricalComponentsToEstimate[ along_track_empirical_acceleration_component ].push_back( cosine_empirical );
//
//    empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );
//    empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
//    empiricalComponentsToEstimate[ across_track_empirical_acceleration_component ].push_back( cosine_empirical );

    additionalParameterNames.push_back( std::make_shared<EmpiricalAccelerationEstimatableParameterSettings>(
        "GRAIL-A", "Moon", empiricalComponentsToEstimate ));


    std::shared_ptr<EstimationOutput<> > estimationOutput =
        createBestFitToCurrentEphemeris( bodies, accelerationModelMap, bodiesToEstimate, centralBodies,
                                         numerical_integrators::rungeKuttaFixedStepSettings( 60.0,
                                                                                             numerical_integrators::rungeKuttaFehlberg78 ),
                                         initialPropagationTime, finalPropagtionTime, 60.0, additionalParameterNames );




}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//}
//
//}
