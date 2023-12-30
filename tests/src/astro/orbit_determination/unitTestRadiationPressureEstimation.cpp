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

    double initialTime = basic_astrodynamics::timeFromIsoString<double>( initialTimeString );
    double finalTime = initialTime + 7200.0;

    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    // Create settings for default bodies
    std::vector<std::string>
        bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    std::string globalFrameOrigin = "Earth";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodiesToCreate, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Add Moon radiation pressure models
    std::vector<std::shared_ptr<PanelRadiosityModelSettings> > panelRadiosityModels;
    panelRadiosityModels.push_back( angleBasedThermalPanelRadiosityModelSettings( 95.0, 385.0, 0.95, "Sun" ));
    panelRadiosityModels.push_back(
        albedoPanelRadiosityModelSettings( SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun" ));
    std::map<std::string, std::vector<std::string> > originalSourceToSourceOccultingBodies;
    originalSourceToSourceOccultingBodies[ "Sun" ].push_back( "Earth" );
    bodySettings.at( "Moon" )->radiationSourceModelSettings =
        extendedRadiationSourceModelSettingsWithOccultationMap(
            panelRadiosityModels, { 4, 8 }, originalSourceToSourceOccultingBodies );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );

    bodySettings.at( spacecraftName )->constantMass = 150.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 5.0;
    double radiationPressureCoefficient = 1.5;
    std::map<std::string, std::vector<std::string> > sourceToTargetOccultingBodies;
    sourceToTargetOccultingBodies[ "Sun" ].push_back( "Moon" );
    bodySettings.at( spacecraftName )->radiationPressureTargetModelSettings =
        cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
            referenceAreaRadiation, radiationPressureCoefficient, sourceToTargetOccultingBodies );

    for ( unsigned int test = 0; test < 3; test++ )
    {
        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies<long double, Time>( bodySettings );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfSpacecraft;
        accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        if( test == 0 || test == 2 )
        {
            accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
        }
        accelerationsOfSpacecraft[ "Earth" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<SphericalHarmonicAccelerationSettings>( 8, 8 ));
        if( test == 1 || test == 2 )
        {
            accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared<AccelerationSettings>( radiation_pressure ));
        }
        accelerationsOfSpacecraft[ "Moon" ].push_back( empiricalAcceleration( ));
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Jupiter" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));
        accelerationsOfSpacecraft[ "Saturn" ].push_back( std::make_shared<AccelerationSettings>( point_mass_gravity ));

        accelerationMap[ "GRAIL-A" ] = accelerationsOfSpacecraft;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector<std::string> bodiesToEstimate = { "GRAIL-A" };
        std::vector<std::string> centralBodies = { "Moon" };
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToEstimate, centralBodies );

        // Get initial state from Spice kernel
        Eigen::VectorXd initialState =
            ( Eigen::VectorXd( 6 )
                << -652685.231403348, 648916.108348616, 1549361.00176057, 495.505044473055, -1364.12825227884, 774.676881961036 ).finished( );
        // Define propagator settings
        std::shared_ptr<TranslationalStatePropagatorSettings<long double, Time> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<long double, Time> >
                ( centralBodies, accelerationModelMap, bodiesToEstimate, initialState.template cast<long double>( ),
                  Time( initialTime ),
                  numerical_integrators::rungeKuttaFixedStepSettings<Time>( 120.0,
                                                                            numerical_integrators::rungeKuttaFehlberg78 ),
                  std::make_shared<PropagationTimeTerminationSettings>( finalTime ));

        // Create parameters
        std::vector<std::shared_ptr<EstimatableParameterSettings> > parameterNames =
            getInitialStateParameterSettings<long double, Time>( propagatorSettings, bodies );
        std::vector<std::shared_ptr<estimatable_parameters::EstimatableParameterSettings> > additionalParameterNames;
        parameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ));
        std::shared_ptr<estimatable_parameters::EstimatableParameterSet<long double> > parametersToEstimate =
            createParametersToEstimate<long double, Time>( parameterNames, bodies, propagatorSettings );

        SingleArcVariationalEquationsSolver<long double, Time> variationalEquationsSolver(
            bodies, propagatorSettings, parametersToEstimate );
        auto sensitivityResults = variationalEquationsSolver.getSensitivityMatrixSolution( );

        {
            double parameterPerturbation = 1.0;
            auto nominalParameters = parametersToEstimate->getFullParameterValues<long double>( );
            auto perturbedParameters = nominalParameters;
            perturbedParameters( 6 ) += parameterPerturbation;
            parametersToEstimate->resetParameterValues( perturbedParameters );
            SingleArcDynamicsSimulator<long double, Time> upperturbedDynamics(
                bodies, propagatorSettings );
            auto upperturbedResults = upperturbedDynamics.getEquationsOfMotionNumericalSolution( );

            perturbedParameters = nominalParameters;
            perturbedParameters( 6 ) -= parameterPerturbation;
            parametersToEstimate->resetParameterValues( perturbedParameters );
            SingleArcDynamicsSimulator<long double, Time> downperturbedDynamics(
                bodies, propagatorSettings );
            auto downperturbedResults = downperturbedDynamics.getEquationsOfMotionNumericalSolution( );

            parametersToEstimate->resetParameterValues( nominalParameters );

            for ( auto it: sensitivityResults )
            {
                Eigen::VectorXd analyticalValue = it.second;
                Eigen::VectorXd numericalValue =
                    (( upperturbedResults.at( it.first ) - downperturbedResults.at( it.first )).transpose( ) /
                     ( 2.0 * parameterPerturbation )).cast<double>( );
                for ( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( analyticalValue( i ) - numericalValue( i ))),
                                       1.0E-4 * analyticalValue.block( 0, 0, 3, 1 ).norm( ));
                    BOOST_CHECK_SMALL(
                        std::fabs( static_cast< double >( analyticalValue( i + 3 ) - numericalValue( i + 3 ))),
                        1.0E-4 * analyticalValue.block( 0, 0, 3, 1 ).norm( ));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
