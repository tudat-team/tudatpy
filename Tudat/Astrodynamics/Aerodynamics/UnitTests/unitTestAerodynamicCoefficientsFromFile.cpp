

#define BOOST_TEST_MAIN

#include <limits>


#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

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
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

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
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create vehicle objects.
    bodyMap[ "SpacePlane" ] = boost::make_shared< simulation_setup::Body >( );

    // Create vehicle aerodynamic coefficients
    std::map< int, std::string > forceCoefficientFiles;
    forceCoefficientFiles[ 0 ] = tudat::input_output::getTudatRootPath( )
            + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt";
    forceCoefficientFiles[ 2 ] = tudat::input_output::getTudatRootPath( )
            + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt";
    std::map< int, std::string > momentCoefficientFiles;
    momentCoefficientFiles[ 1 ] = tudat::input_output::getTudatRootPath( )
            + "/Astrodynamics/Aerodynamics/UnitTests/aurora_Cm.txt";

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles, momentCoefficientFiles, 60.734, 600.0, 60.734, Eigen::Vector3d::Zero( ),
                boost::assign::list_of( aerodynamics::mach_number_dependent )( aerodynamics::angle_of_attack_dependent ),
                true, true );
    bodyMap[ "SpacePlane" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "SpacePlane" ) );


    bodyMap[ "SpacePlane" ]->setConstantBodyMass( 50.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    ////////////////////////////////////// Variable thrust calculations /////////////////////////////////////////////
    // Define thrustdependent variables
    std::vector< propulsion::ThrustDependentVariables > thrustDependencies;
    thrustDependencies.push_back( propulsion::mach_number_dependent_thrust );
    thrustDependencies.push_back( propulsion::dynamic_pressure_dependent_thrust );

    // Read tables with thrust and specific impulse as function of dependent variables
    std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustValues =
            MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Tmax_test.txt" );
    std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > specificImpulseValues =
            MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Propulsion/UnitTests/Isp_test.txt" );

    // Use multilinear interpolation to find Thrust and specific impulse from tables at current values
    // of thrustdpendent variables
    boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                thrustValues.second, thrustValues.first );
    boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                specificImpulseValues.second, specificImpulseValues.first );

    //////////////////////////////////////// End Variable thrust calculations /////////////////////////////////////////////


    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSpacePlane;
    accelerationsOfSpacePlane[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
    accelerationsOfSpacePlane[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    accelerationsOfSpacePlane[ "SpacePlane" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                             boost::make_shared< ThrustDirectionGuidanceSettings >(
                                                                 thrust_direction_from_existing_body_orientation, "Earth" ),
                                                             boost::make_shared< ParameterizedThrustMagnitudeSettings >(
                                                                 thrustMagnitudeInterpolator, thrustDependencies,
                                                                 specificImpulseInterpolator, thrustDependencies ) ) );





    accelerationMap[ "SpacePlane" ] = accelerationsOfSpacePlane;

    bodiesToPropagate.push_back( "SpacePlane" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    setTrimmedConditions( bodyMap.at( "SpacePlane" ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set spherical elements for SpacePlane.
    Vector6d SpacePlaneInitialState;
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::radiusIndex ) = spice_interface::getAverageRadius( "Earth" ) + 2.0E3;
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::latitudeIndex ) = convertDegreesToRadians( 0.0 );
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::longitudeIndex ) = convertDegreesToRadians( 0.0 );
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::speedIndex ) = 170;
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::flightPathIndex ) = convertDegreesToRadians( 20.0 );
    SpacePlaneInitialState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = convertDegreesToRadians( 90.0 );

    // Convert SpacePlane state from spherical elements to Cartesian elements.
    const Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                SpacePlaneInitialState );


    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "SpacePlane" ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "SpacePlane", reference_frames::angle_of_attack ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "SpacePlane", reference_frames::angle_of_sideslip ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "SpacePlane", reference_frames::bank_angle ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_moment_coefficients_dependent_variable, "SpacePlane" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "SpacePlane" ) );

    // Create total propagatortermination settings.
    boost::shared_ptr< PropagationTerminationSettings > propagationTerminationSettings;

    std::vector< boost::shared_ptr< PropagationTerminationSettings > >  propagationTerminationSettingsList;
    propagationTerminationSettingsList.push_back( boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                                      boost::make_shared< SingleDependentVariableSaveSettings >(
                                                          altitude_dependent_variable, "SpacePlane", "Earth" ), 120.0E3, false ) );
    propagationTerminationSettingsList.push_back( boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                                      boost::make_shared< SingleDependentVariableSaveSettings >(
                                                          altitude_dependent_variable, "SpacePlane", "Earth" ), 1.0E3, true ) );
    propagationTerminationSettingsList.push_back( boost::make_shared< PropagationTimeTerminationSettings >( 100.0 ) );

    \
    propagationTerminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                propagationTerminationSettingsList, true );

    // Create propagation settings.
    boost::shared_ptr< TranslationalStatePropagatorSettings < double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, propagationTerminationSettings,
              cowell, boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator.getDependentVariableHistory( );

    // Iterate over results for dependent variables, and check against manually retrieved values.
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > apolloCoefficientInterface =
            bodyMap.at( "SpacePlane" )->getAerodynamicCoefficientInterface( );


    std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > cdMuliArray =
            tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt" );
    boost::shared_ptr< interpolators::Interpolator< double, double > > cdInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                cdMuliArray.second, cdMuliArray.first );

    std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > clMuliArray =
            tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt" );
    boost::shared_ptr< interpolators::Interpolator< double, double > > clInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                clMuliArray.second, clMuliArray.first );

    std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > >  cmMuliArray =
            tudat::input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                tudat::input_output::getTudatRootPath( ) + "/Astrodynamics/Aerodynamics/UnitTests/aurora_Cm.txt" );
    boost::shared_ptr< interpolators::Interpolator< double, double > > cmInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                cmMuliArray.second, cmMuliArray.first );


    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        double machNumber = variableIterator->second( 0 );
        double angleOfAttack = variableIterator->second( 1 );
        double sideslipAngle = variableIterator->second( 2 );
        double bankAngle = variableIterator->second( 3 );
        Eigen::Vector3d momentCoefficients = variableIterator->second.segment( 4, 3 );
        Eigen::Vector3d forceCoefficients = variableIterator->second.segment( 7, 3 );

        std::vector< double > aerodynamicCoefficientInput;
        aerodynamicCoefficientInput.push_back( machNumber );
        aerodynamicCoefficientInput.push_back( angleOfAttack );
        apolloCoefficientInterface->updateCurrentCoefficients( aerodynamicCoefficientInput );
        Eigen::Vector3d computedForceCoefficients = apolloCoefficientInterface->getCurrentForceCoefficients( );
        Eigen::Vector3d computedMomentCoefficients = apolloCoefficientInterface->getCurrentMomentCoefficients( );

        BOOST_CHECK_SMALL( std::fabs( computedForceCoefficients( 1 ) ), 1.0E-14 );

        // Check force and moment coefficients from output to direct computation from interface
        for( unsigned int i = 0; i < 3; i ++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( computedForceCoefficients( i ) - forceCoefficients( i ) ),
                        6.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( computedMomentCoefficients( i ) - momentCoefficients( i ) ),
                        6.0 * std::numeric_limits< double >::epsilon( ) );
        }

        // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
        BOOST_CHECK_EQUAL( momentCoefficients( 0 ), 0.0 );
        BOOST_CHECK_SMALL( std::fabs( momentCoefficients( 1 ) ), 1.0E-16 );
        BOOST_CHECK_EQUAL( momentCoefficients( 2 ), 0.0 );

        // Checkc if coefficients are correctly parsed from files
        double directCd = cdInterpolator->interpolate( aerodynamicCoefficientInput );
        double directCl = clInterpolator->interpolate( aerodynamicCoefficientInput );
        double directCm = cmInterpolator->interpolate( aerodynamicCoefficientInput );
        BOOST_CHECK_SMALL( std::fabs( directCd - forceCoefficients( 0 ) ) , 1.0E-16 );
        BOOST_CHECK_SMALL( std::fabs( directCm - momentCoefficients( 1 ) ) , 1.0E-16 );
        BOOST_CHECK_SMALL( std::fabs( directCl - forceCoefficients( 2 ) ) , 1.0E-16 );
        BOOST_CHECK_EQUAL( std::fabs( forceCoefficients( 1 ) ), 0.0 );

        // Check zero angles.
        BOOST_CHECK_EQUAL( sideslipAngle, 0.0 );
        BOOST_CHECK_EQUAL( bankAngle, 0.0 );


        // Check if trim is not trivial
        std::vector< double > perturbedAerodynamicCoefficientInput;
        perturbedAerodynamicCoefficientInput.push_back( machNumber + 1.0 );
        perturbedAerodynamicCoefficientInput.push_back( angleOfAttack - 0.1 );
        apolloCoefficientInterface->updateCurrentCoefficients( perturbedAerodynamicCoefficientInput );
        Eigen::Vector3d computedPerturbedMomentCoefficients = apolloCoefficientInterface->getCurrentMomentCoefficients( );
        BOOST_CHECK_EQUAL( ( std::fabs( computedPerturbedMomentCoefficients( 1 ) ) > 1.0E-8 ), true );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

