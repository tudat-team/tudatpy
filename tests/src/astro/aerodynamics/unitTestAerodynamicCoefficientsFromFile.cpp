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

#include <limits>


#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <tudat/basics/testMacros.h>

#include <tudat/simulation/simulation.h>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_coefficients_from_file )

BOOST_AUTO_TEST_CASE( testAerodynamicCoefficientsFromFile )
{

    using namespace tudat;
    using namespace ephemerides;
    using namespace interpolators;
    using namespace numerical_integrators;
    using namespace spice_interface;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace propagators;
    using namespace aerodynamics;
    using namespace basic_mathematics;
    using namespace input_output;
    using namespace unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( std::vector<std::string>() );

    for( unsigned int i = 0; i < 3; i++ )
    {
        // Set simulation start epoch.
        const double simulationStartEpoch = 0.0;

        // Set simulation end epoch.
        const double simulationEndEpoch = 300.0;

        // Set numerical integration fixed step size.
        const double fixedStepSize = 1.0;


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define simulation body settings.
        BodyListSettings bodySettings =
                getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                        simulationEndEpoch + 10.0 * fixedStepSize, "SSB", "J2000" );
        bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ), "SSB", "J2000" );

        // Create Earth object
        simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Create vehicle objects.
        bodies.createEmptyBody( "SpacePlane" );
        bodies.at( "SpacePlane" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
        bodies.at( "SpacePlane" )->getVehicleSystems( )->setCurrentControlSurfaceDeflection( "TestSurface", 0.1 );

        // Create vehicle  coefficients
        std::map< int, std::string > forceCoefficientFiles;
        forceCoefficientFiles[ 0 ] = tudat::paths::getTudatTestDataPath( )
                + "/aurora_CD.txt";
        forceCoefficientFiles[ 2 ] = tudat::paths::getTudatTestDataPath( )
                + "/aurora_CL.txt";
        std::map< int, std::string > momentCoefficientFiles;
        momentCoefficientFiles[ 1 ] = tudat::paths::getTudatTestDataPath( )
                + "/aurora_Cm.txt";
        std::map< int, std::string > controlSurfaceForceCoefficientFiles;
        controlSurfaceForceCoefficientFiles[ 0 ] = tudat::paths::getTudatTestDataPath( )
                + "/dCDwTest.txt";
        std::map< int, std::string > controlSurfaceMomentCoefficientFiles;
        controlSurfaceMomentCoefficientFiles[ 1 ] = tudat::paths::getTudatTestDataPath( )
                + "/dCDwTest.txt";

        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles, momentCoefficientFiles, 60.734, 600.0, 60.734, Eigen::Vector3d::Zero( ),
        { aerodynamics::mach_number_dependent, aerodynamics::angle_of_attack_dependent },
                    true, true );
        if( i == 1 )
        {
            aerodynamicCoefficientSettings->setControlSurfaceSettings(
                        simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                            controlSurfaceForceCoefficientFiles, controlSurfaceMomentCoefficientFiles,
            { aerodynamics::mach_number_dependent, aerodynamics::angle_of_attack_dependent, aerodynamics::control_surface_deflection_dependent } ),
                        "TestSurface" );
        }
        else if( i == 2 )
        {
            aerodynamicCoefficientSettings->setControlSurfaceSettings(
                        simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                            controlSurfaceForceCoefficientFiles,
            { aerodynamics::mach_number_dependent, aerodynamics::angle_of_attack_dependent, aerodynamics::control_surface_deflection_dependent } ),
                        "TestSurface" );
        }

        bodies.at( "SpacePlane" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "SpacePlane" ) );

        bodies.at( "SpacePlane" )->setConstantBodyMass( 50.0E3 );

        std::shared_ptr< ephemerides::RotationalEphemeris > vehicleRotationModel =
                createRotationModel(
                    std::make_shared< PitchTrimRotationSettings >( "Earth", "J2000", "VehicleFixed" ),
                    "SpacePlane", bodies );

        bodies.at( "SpacePlane" )->setRotationalEphemeris( vehicleRotationModel  );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        ////////////////////////////////////// Variable thrust calculations /////////////////////////////////////////////
        // Define thrustdependent variables
        std::vector< propulsion::ThrustIndependentVariables > thrustDependencies;
        thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
        thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

        // Read tables with thrust and specific impulse as function of dependent variables
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/Tmax_test.txt" );
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > specificImpulseValues =
                MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/Isp_test.txt" );

        // Use multilinear interpolation to find Thrust and specific impulse from tables at current values
        // of thrustdpendent variables
        std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    thrustValues.second, thrustValues.first );
        std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    specificImpulseValues.second, specificImpulseValues.first );

        addEngineModel(
                    "SpacePlane", "MainEngine",
                    std::make_shared< ParameterizedThrustMagnitudeSettings >(
                        thrustMagnitudeInterpolator, thrustDependencies,
                        specificImpulseInterpolator, thrustDependencies ), bodies );

        //////////////////////////////////////// End Variable thrust calculations /////////////////////////////////////////////


        // Define acceleration model settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacePlane;
        accelerationsOfSpacePlane[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
        accelerationsOfSpacePlane[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

        accelerationsOfSpacePlane[ "SpacePlane" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                                  "MainEngine" ) );

        accelerationMap[ "SpacePlane" ] = accelerationsOfSpacePlane;

        bodiesToPropagate.push_back( "SpacePlane" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set spherical elements for SpacePlane.
        Eigen::Vector6d SpacePlaneInitialState;
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::radiusIndex ) = spice_interface::getAverageRadius( "Earth" ) + 2.0E3;
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::latitudeIndex ) = convertDegreesToRadians( 0.0 );
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::longitudeIndex ) = convertDegreesToRadians( 0.0 );
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::speedIndex ) = 170;
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::flightPathIndex ) = convertDegreesToRadians( 20.0 );
        SpacePlaneInitialState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = convertDegreesToRadians( 90.0 );

        // Convert SpacePlane state from spherical elements to Cartesian elements.
        const Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                    SpacePlaneInitialState );


        // Define list of dependent variables to save.
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "SpacePlane" ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "SpacePlane", reference_frames::angle_of_attack ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "SpacePlane", reference_frames::angle_of_sideslip ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "SpacePlane", reference_frames::bank_angle ) );
        if( i > 0 )
        {
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            control_surface_deflection_dependent_variable, "SpacePlane", "TestSurface" ) );
        }
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_moment_coefficients_dependent_variable, "SpacePlane" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_force_coefficients_dependent_variable, "SpacePlane" ) );

        // Create total propagatortermination settings.
        std::shared_ptr< PropagationTerminationSettings > propagationTerminationSettings;

        std::vector< std::shared_ptr< PropagationTerminationSettings > >  propagationTerminationSettingsList;
        propagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationDependentVariableTerminationSettings >(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            altitude_dependent_variable, "SpacePlane", "Earth" ), 120.0E3, false ) );
        propagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationDependentVariableTerminationSettings >(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            altitude_dependent_variable, "SpacePlane", "Earth" ), 1.0E3, true ) );
        propagationTerminationSettingsList.push_back(
                    std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

        \
        propagationTerminationSettings = std::make_shared< PropagationHybridTerminationSettings >(
                    propagationTerminationSettingsList, true );

        // Create propagation settings.
        std::shared_ptr< TranslationalStatePropagatorSettings < double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, propagationTerminationSettings,
                  cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableSolution =
                dynamicsSimulator.getDependentVariableHistory( );

        // Iterate over results for dependent variables, and check against manually retrieved values.
        std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > apolloCoefficientInterface =
                bodies.at( "SpacePlane" )->getAerodynamicCoefficientInterface( );


        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > cdMuliArray =
                tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/aurora_CD.txt" );
        std::shared_ptr< interpolators::Interpolator< double, double > > cdInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    cdMuliArray.second, cdMuliArray.first );

        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > clMuliArray =
                tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/aurora_CL.txt" );
        std::shared_ptr< interpolators::Interpolator< double, double > > clInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    clMuliArray.second, clMuliArray.first );

        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > >  cmMuliArray =
                tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/aurora_Cm.txt" );
        std::shared_ptr< interpolators::Interpolator< double, double > > cmInterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    cmMuliArray.second, cmMuliArray.first );

        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > >  controlSurfaceIncrementMuliArray =
                tudat::input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables(
                    tudat::paths::getTudatTestDataPath( ) + "/dCDwTest.txt" );
        std::shared_ptr< interpolators::Interpolator< double, double > > controlSurfaceIncrementinterpolator =
                std::make_shared< interpolators::MultiLinearInterpolator< double, double, 3 > >(
                    controlSurfaceIncrementMuliArray.second, controlSurfaceIncrementMuliArray.first );

        int parameterAddition = 0;
        if( i > 0 )
        {
            parameterAddition++;
        }
        for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
             variableIterator != dependentVariableSolution.end( ); variableIterator++ )
        {
            double machNumber = variableIterator->second( 0 );
            double angleOfAttack = variableIterator->second( 1 );
            double sideslipAngle = variableIterator->second( 2 );
            double bankAngle = variableIterator->second( 3 );
            std::cout<<"AoA "<<angleOfAttack<<std::endl;

            double controlSurfaceDeflection  = 0.0;
            if( i > 0 )
            {
                controlSurfaceDeflection = variableIterator->second( 4 );
                std::cout<<"CSD "<<controlSurfaceDeflection<<std::endl;

            }

            Eigen::Vector3d momentCoefficients = variableIterator->second.segment( 4 + parameterAddition, 3 );
            Eigen::Vector3d forceCoefficients = variableIterator->second.segment( 7 + parameterAddition, 3 );

            std::vector< double > aerodynamicCoefficientInput;
            aerodynamicCoefficientInput.push_back( machNumber );
            aerodynamicCoefficientInput.push_back( angleOfAttack );

            std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput;
            if( i > 0 )
            {
                controlSurfaceCoefficientInput[ "TestSurface" ] = aerodynamicCoefficientInput;
                controlSurfaceCoefficientInput[ "TestSurface" ].push_back( controlSurfaceDeflection );
            }

            apolloCoefficientInterface->updateFullCurrentCoefficients(
                        aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
            Eigen::Vector3d computedForceCoefficients = apolloCoefficientInterface->getCurrentForceCoefficients( );
            Eigen::Vector3d computedMomentCoefficients = apolloCoefficientInterface->getCurrentMomentCoefficients( );

            BOOST_CHECK_SMALL( std::fabs( computedForceCoefficients( 1 ) ), 1.0E-14 );

            // Check force and moment coefficients from output to direct computation from interface
            for( unsigned int j = 0; j < 3; j ++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( computedForceCoefficients( j ) - forceCoefficients( j ) ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL(
                            std::fabs( computedMomentCoefficients( j ) - momentCoefficients( j ) ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
            }

            // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
            BOOST_CHECK_EQUAL( momentCoefficients( 0 ), 0.0 );
            BOOST_CHECK_SMALL( std::fabs( momentCoefficients( 1 ) ), 1.0E-16 );
            BOOST_CHECK_EQUAL( momentCoefficients( 2 ), 0.0 );

            // Check if coefficients are correctly parsed from files
            double directCd = cdInterpolator->interpolate( aerodynamicCoefficientInput );
            if( i > 0 )
            {
                directCd += controlSurfaceIncrementinterpolator->interpolate( controlSurfaceCoefficientInput[ "TestSurface" ] );
            }
            double directCl = clInterpolator->interpolate( aerodynamicCoefficientInput );
            double directCm = cmInterpolator->interpolate( aerodynamicCoefficientInput );
            if( i == 1 )
            {
                directCm += controlSurfaceIncrementinterpolator->interpolate( controlSurfaceCoefficientInput[ "TestSurface" ] );
            }
            BOOST_CHECK_SMALL( std::fabs( directCd - forceCoefficients( 0 ) ) , 1.0E-16 );
            BOOST_CHECK_SMALL( std::fabs( directCm - momentCoefficients( 1 ) ) , 1.0E-16 );
            BOOST_CHECK_SMALL( std::fabs( directCl - forceCoefficients( 2 ) ) , 1.0E-16 );
            BOOST_CHECK_EQUAL( std::fabs( forceCoefficients( 1 ) ), 0.0 );

            if( i == 1 )
            {
                controlSurfaceCoefficientInput[ "TestSurface" ][ 2 ] = 0.0;
                apolloCoefficientInterface->updateFullCurrentCoefficients(
                            aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
                Eigen::Vector3d computedPerturbedMomentCoefficients =
                        apolloCoefficientInterface->getCurrentMomentCoefficients( );
                BOOST_CHECK_EQUAL( ( std::fabs( computedPerturbedMomentCoefficients( 1 ) ) > 1.0E-8 ), true );
                controlSurfaceCoefficientInput[ "TestSurface" ][ 2 ] = controlSurfaceDeflection;

            }

            // Check zero angles.
            BOOST_CHECK_EQUAL( sideslipAngle, 0.0 );
            BOOST_CHECK_EQUAL( bankAngle, 0.0 );


            // Check if trim is not trivial
            std::vector< double > perturbedAerodynamicCoefficientInput;
            perturbedAerodynamicCoefficientInput.push_back( machNumber + 1.0 );
            perturbedAerodynamicCoefficientInput.push_back( angleOfAttack - 0.1 );
            apolloCoefficientInterface->updateFullCurrentCoefficients(
                        perturbedAerodynamicCoefficientInput, controlSurfaceCoefficientInput );
            Eigen::Vector3d computedPerturbedMomentCoefficients = apolloCoefficientInterface->getCurrentMomentCoefficients( );
            BOOST_CHECK_EQUAL( ( std::fabs( computedPerturbedMomentCoefficients( 1 ) ) > 1.0E-8 ), true );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

