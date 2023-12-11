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


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_full_planetary_rotational_parameters_estimation )

//Using declarations.
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
using namespace tudat;


BOOST_AUTO_TEST_CASE( test_CannonballPartials )
{
    std::string initialTimeString = "2012-04-06 01:44:04.686";
    std::string finalTimeString = "2012-04-07 21:25:39.686";

    double initialTime = basic_astrodynamics::timeFromIsoString< double >( initialTimeString );
    double finalTime = basic_astrodynamics::timeFromIsoString< double >( finalTimeString );

    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_v07.tf" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_sclkscet_00013.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_sclkscet_00014.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/gra_rec_120402_120408.bc" );

    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    std::string globalFrameOrigin = "Earth";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodiesToCreate, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Add Moon radiation pressure models
    std::vector< std::shared_ptr< PanelRadiosityModelSettings > > panelRadiosityModels;
    panelRadiosityModels.push_back(angleBasedThermalPanelRadiosityModelSettings( 95.0, 385.0, 0.95, "Sun" ) );
    panelRadiosityModels.push_back(albedoPanelRadiosityModelSettings( SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun" ) );
    std::map< std::string, std::vector< std::string > > originalSourceToSourceOccultingBodies;
    originalSourceToSourceOccultingBodies[ "Sun" ].push_back( "Earth" );
    bodySettings.at( "Moon" )->radiationSourceModelSettings =
        extendedRadiationSourceModelSettingsWithOccultationMap(
            panelRadiosityModels, { 4, 8  }, originalSourceToSourceOccultingBodies );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
        std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialTime - 3600.0, finalTime + 3600.0, 10.0, spacecraftCentralBody, globalFrameOrientation );


    bodySettings.at( spacecraftName )->constantMass = 150.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 5.0;
    double radiationPressureCoefficient = 1.5;
    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies;
    sourceToTargetOccultingBodies[ "Sun" ].push_back( "Moon" );
    bodySettings.at( spacecraftName )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
        referenceAreaRadiation, radiationPressureCoefficient, sourceToTargetOccultingBodies );
    bodySettings.at( spacecraftName )->rotationModelSettings = spiceRotationModelSettings(
        globalFrameOrientation, spacecraftName + "_SPACECRAFT", "" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

     // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfSpacecraft;
    accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
    accelerationsOfSpacecraft[ "Earth" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
    accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<SphericalHarmonicAccelerationSettings>( 8, 8 ));
    accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
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

    // Get initial state from Spice kernel
    Eigen::VectorXd initialState = getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialTime );

    // Define propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< long double, Time > >
            ( centralBodies, accelerationModelMap, bodiesToEstimate, initialState.template cast< long double >( ), Time( initialTime ),
              numerical_integrators::rungeKuttaFixedStepSettings< Time >( 30.0, numerical_integrators::rungeKuttaFehlberg78 ),
              std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    // Create parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
        getInitialStateParameterSettings< long double, Time >( propagatorSettings, bodies );
    std::vector<std::shared_ptr<estimatable_parameters::EstimatableParameterSettings> > additionalParameterNames;
    parameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ));
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
        createParametersToEstimate< long double, Time >( parameterNames, bodies, propagatorSettings );

    SingleArcVariationalEquationsSolver< long double, Time > variationalEquationsSolver(
        bodies, propagatorSettings, parametersToEstimate );
    auto sensitivityResults = variationalEquationsSolver.getSensitivityMatrixSolution( );

    {
        double parameterPerturbation = 0.1;
        auto nominalParameters = parametersToEstimate->getFullParameterValues< long double >( );
        auto perturbedParameters = nominalParameters;
        perturbedParameters( 6 ) += parameterPerturbation;
        parametersToEstimate->resetParameterValues( perturbedParameters );
        SingleArcDynamicsSimulator< long double, Time > upperturbedDynamics(
            bodies, propagatorSettings );
        auto upperturbedResults = upperturbedDynamics.getEquationsOfMotionNumericalSolution( );

        perturbedParameters = nominalParameters;
        perturbedParameters( 6 ) -= parameterPerturbation;
        parametersToEstimate->resetParameterValues( perturbedParameters );
        SingleArcDynamicsSimulator< long double, Time > downperturbedDynamics(
            bodies, propagatorSettings );
        auto downperturbedResults = downperturbedDynamics.getEquationsOfMotionNumericalSolution( );

        parametersToEstimate->resetParameterValues( nominalParameters );

        for( auto it : sensitivityResults )
        {
            std::cout<<it.first - initialTime<<std::endl;
            std::cout<<it.second.transpose( )<<std::endl;
            std::cout<<( upperturbedResults.at( it.first ) - downperturbedResults.at( it.first ) ).transpose( ) / ( 2.0 * parameterPerturbation )<<std::endl<<std::endl;
        }


    }


}

BOOST_AUTO_TEST_SUITE_END( )

}

}
