/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 *    Notes:
 *      Test tolerance was set at 5.0e-15 (or 5.0e-7 for floats) instead of epsilon due to
 *      rounding errors in Eigen types with entries over a number of orders of magnitude,
 *      presumably causing the observed larger than epsilon relative differences between
 *      expected and computed values.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>

#include <memory>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

// #include "tudat/astro/basic_astro/gravityDeformationModel.h"
// #include "tudat/astro/basic_astro/testBody.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/propagation_setup/createGravityDeformationModels.h"

// #include "tudat/basics/testMacros.h"

// #include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
// #include "tudat/math/basic/sphericalHarmonics.h"

namespace tudat
{
namespace unit_tests
{

// using basic_astrodynamics::GravityDeformationModel;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;

BOOST_AUTO_TEST_CASE( test_gravityDeformationModel )
{
    //     std::cout.precision( 20 );

    //    // Load spice kernels.
    //     spice_interface::loadStandardSpiceKernels( );

    //     // Specify initial time
    //     double initialTime = 0.0;
    //     double finalTime = 600.0; //1.0 * physical_constants::JULIAN_DAY;

    //     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io", "Europa", "Ganymede", "Callisto" };

    //     // Get body settings.
    //     BodyListSettings bodySettings =
    //             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, "Jupiter", "J2000" );

    //     double effectiveMu = 126692494120023072.0;
    //     Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
    //     std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >(
    //     initialKeplerianState, 0.0, effectiveMu, "Jupiter", "J2000" ); bodySettings.at( "Io" )->ephemerisSettings =
    //     keplerEphemerisSettings;

    //     bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000",
    //     "Io_IAU" );

    //     // Create bodies needed in simulation
    //     SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
    //     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
    //     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );

    //     Eigen::Quaterniond testRotation = bodies.at( "Io" )->getRotationalEphemeris( )->getRotationToBaseFrame( 540.0 );
    //     Eigen::Vector3d testRelativeInertialPosition = - bodies.at( "Io" )->getPositionInBaseFrameFromEphemeris( 540.0 );
    //     Eigen::Vector3d testRelativePosition = testRotation.inverse( ) * ( testRelativeInertialPosition );
    //     Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical( testRelativePosition );
    //     double testLongitude = testSphericalPosition[ 2 ];
    //     std::cout << "testRotation " << testRotation.toRotationMatrix( ) << std::endl;
    //     std::cout << "test relativeInertialPosition " << testRelativeInertialPosition.transpose( ) << std::endl;
    //     std::cout << "test longitude " << testLongitude << std::endl;

    //     const double maxwellRelaxationTime = 1000.0;
    //     const double globalRelaxationTime = 2000.0;
    //     const double loveNumber = 0.4;
    //     const double rotationRate =( 2.0 * mathematical_constants::PI ) / ( 1.77 * 86400.0 );
    //     const int maximumDegree = 2;
    //     const int maximumOrder = 2;

    //     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >(
    //         maxwellRelaxationTime, globalRelaxationTime, loveNumber, rotationRate, maximumDegree, maximumOrder, "Jupiter" );

    //     std::shared_ptr< Body > deformingBody = bodies.at( "Io" );
    //     std::shared_ptr< Body > perturbingBody = bodies.at( "Jupiter" );
    //     std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel =
    //         createMaxwellGravityFieldDeformationModel( deformingBody, perturbingBody, "Io", "Jupiter", maxwellDeformationSettings );

    //     std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
    //     gravityDeformationModels[ "Io" ] = maxwellDeformationModel;

    //     double timeStep = 60.0;
    //     // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >
    //     //             ( initialTime, timeStep, rungeKuttaFehlberg78, timeStep, timeStep );
    //     std::shared_ptr< IntegratorSettings< > > integratorSettings =
    //                     std::make_shared< IntegratorSettings< > >
    //                     ( rungeKutta4, initialTime, timeStep );

    //     std::vector< std::string > bodiesToPropagate = { "Io" };
    //     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialBodyGravity = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 3, 1 );

    //     Eigen::Vector3d computedEquilibriumCoefficients = Eigen::Vector3d::Zero( );
    //     std::shared_ptr< SphericalHarmonicsGravityField > shModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
    //     bodies.at( "Io" )->getGravityFieldModel( ) ); double distance = 4.2e8; double radius = shModel->getReferenceRadius( ); double
    //     muIo = shModel->getGravitationalParameter( ); double muJup = bodies.at( "Jupiter" )->getGravitationalParameter( ); double
    //     ratioDistanceRadiusPowerThree = radius * radius * radius / ( distance * distance * distance ); double muRatio = muJup / muIo;
    //     computedEquilibriumCoefficients[ 0 ] = - loveNumber * ( rotationRate * rotationRate * radius * radius * radius / ( 3.0 * muIo ) +
    //     0.5 * muRatio * ratioDistanceRadiusPowerThree ); computedEquilibriumCoefficients[ 1 ] = loveNumber / 4 * muRatio *
    //     ratioDistanceRadiusPowerThree; computedEquilibriumCoefficients[ 2 ] = 0.0; std::cout << "computedEquilibriumCoefficients " <<
    //     computedEquilibriumCoefficients.transpose( ) << std::endl;

    //     Eigen::MatrixXd originalCosineMatrix = shModel->getCosineCoefficients();
    //     Eigen::MatrixXd originalSineMatrix = shModel->getSineCoefficients();
    //     Eigen::Vector3d computedInitialCoefficients = Eigen::Vector3d::Zero( );
    //     computedInitialCoefficients[ 0 ] = originalCosineMatrix( 2, 0 );
    //     computedInitialCoefficients[ 1 ] = originalCosineMatrix( 2, 2 );
    //     computedInitialCoefficients[ 2 ] = originalSineMatrix( 2, 2 );
    //     std::cout << "computedInitialCoefficients " << computedInitialCoefficients.transpose( ) << std::endl;

    //     Eigen::Vector3d computedInitialTransientCoefficients = Eigen::Vector3d::Zero( );
    //     computedInitialTransientCoefficients = ( 1.0 / ( globalRelaxationTime - maxwellRelaxationTime ) ) *
    //         ( globalRelaxationTime * computedInitialCoefficients - maxwellRelaxationTime * computedEquilibriumCoefficients );
    //     std::cout << "computedInitialTransientCoefficients: " << computedInitialTransientCoefficients.transpose( ) << std::endl;

    //     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings =
    //         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels,
    //         computedInitialCoefficients, integratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    //     SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, gravityPropagatorSettings );
    //     std::map< double, Eigen::VectorXd > results = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    //     std::cout << "size results: " << results.size( ) << std::endl;

    //     std::map< double, Eigen::Vector3d > computedSolution;
    //     for ( double time = initialTime ; time <= finalTime ; time += timeStep )
    //     {
    //         Eigen::Vector3d currentTransientCoefficients =
    //             ( computedInitialTransientCoefficients - computedEquilibriumCoefficients ) * std::exp( - time / globalRelaxationTime )
    //             + computedEquilibriumCoefficients;
    //         std::cout << "currentTransientCoefficients " << currentTransientCoefficients.transpose( ) << std::endl;

    //         Eigen::Vector3d currentNominalCoefficients =
    //             ( 1.0 - maxwellRelaxationTime / globalRelaxationTime ) * currentTransientCoefficients
    //             + maxwellRelaxationTime / globalRelaxationTime * computedEquilibriumCoefficients;
    //         std::cout << "currentNominalCoefficients " << currentNominalCoefficients.transpose( ) << std::endl;
    //     }

    // -----------------------------------------------------------

    std::cout.precision( 20 );

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 600.0 / 10.0;  // 1.0 * physical_constants::JULIAN_DAY;

    std::vector< std::string > bodiesToCreate = { "Jupiter", "Io", "Europa", "Ganymede", "Callisto" };

    // Get body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, "Jupiter", "J2000" );

    // double muJupiter = 1.26683750521837e17;
    double muIo = 5959924010272.5136719;
    double muJupiter = 126686534196012800.0;
    double muEffective = 126692494120023072.0;
    Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
    std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings =
            std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
    // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

    // orbital period
    double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
    double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );
    std::cout << "rotationRateIo " << rotationRateIo << std::endl;

    Eigen::Matrix3d initialOrientation = Eigen::Matrix3d::Identity( );
    initialOrientation( 0, 0 ) = -1.0;
    initialOrientation( 1, 1 ) = -1.0;
    // bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );
    bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< simulation_setup::SimpleRotationModelSettings >(
            "J2000", "IAU_Io", Eigen::Quaterniond( initialOrientation ), initialTime, rotationRateIo );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
    bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
    bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
    // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

    double scaledMeanMomentOfInertiaIo = 0.37685;
    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )
            ->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertiaIo );
    std::cout << "initial inertia tensor " << std::endl;
    std::cout << bodies.at( "Io" )->getGravityFieldModel( )->getInertiaTensor( ) << std::endl;

    const double maxwellRelaxationTime = 1000.0;
    const double globalRelaxationTime = 2000.0;
    const double loveNumber = 0.4;
    // const double rotationRate = ( 2.0 * mathematical_constants::PI ) / ( 1.77 * 86400.0 );
    const int maximumDegree = 2;
    const int maximumOrder = 2;

    // std::cout << "initial rotation rate " << rotationRate << std::endl;

    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings =
            std::make_shared< MaxwellDeformationSettings >( maxwellRelaxationTime,
                                                            globalRelaxationTime,
                                                            loveNumber,
                                                            /*rotationRateIo,*/ maximumDegree,
                                                            maximumOrder,
                                                            std::vector< std::string >( { "Jupiter" } ) );

    // std::shared_ptr< Body > deformingBody = bodies.at( "Io" );
    // std::shared_ptr< Body > perturbingBody = bodies.at( "Jupiter" );
    // std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel =
    //     createMaxwellGravityFieldDeformationModel( deformingBody, perturbingBody, "Io", "Jupiter", maxwellDeformationSettings );

    // std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
    // gravityDeformationModels[ "Io" ] = maxwellDeformationModel;

    std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;
    gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

    basic_astrodynamics::GravityDeformationModelMap deformationModels =
            createGravityDeformationModelsMap( bodies, gravityDeformationModelMap );

    std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
    gravityDeformationModels[ "Io" ] = deformationModels.at( "Io" )[ 0 ];

    double timeStep = 60.0;
    // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >
    // ( initialTime, timeStep, rungeKuttaFehlberg78, timeStep, timeStep );
    std::shared_ptr< IntegratorSettings<> > integratorSettings =
            std::make_shared< IntegratorSettings<> >( rungeKutta4, initialTime, timeStep );

    std::vector< std::string > bodiesToPropagate = { "Io" };
    // Eigen::Matrix< double, Eigen::Dynamic, 1 > initialBodyGravity = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 3, 1 );

    Eigen::Vector3d computedEquilibriumCoefficients = Eigen::Vector3d::Zero( );
    std::shared_ptr< SphericalHarmonicsGravityField > shModel =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
    double distance = 4.2e8;
    double radius = shModel->getReferenceRadius( );
    double muJup = bodies.at( "Jupiter" )->getGravitationalParameter( );
    double ratioDistanceRadiusPowerThree = radius * radius * radius / ( distance * distance * distance );
    double muRatio = muJup / muIo;
    computedEquilibriumCoefficients[ 0 ] = -loveNumber *
            ( rotationRateIo * rotationRateIo * radius * radius * radius / ( 3.0 * muIo ) + 0.5 * muRatio * ratioDistanceRadiusPowerThree );
    computedEquilibriumCoefficients[ 1 ] = loveNumber / 4 * muRatio * ratioDistanceRadiusPowerThree;
    computedEquilibriumCoefficients[ 2 ] = 0.0;
    std::cout << "computedEquilibriumCoefficients " << computedEquilibriumCoefficients.transpose( ) << std::endl;

    Eigen::MatrixXd originalCosineMatrix = shModel->getCosineCoefficients( );
    Eigen::MatrixXd originalSineMatrix = shModel->getSineCoefficients( );
    Eigen::Vector3d computedInitialCoefficients = Eigen::Vector3d::Zero( );
    computedInitialCoefficients[ 0 ] = originalCosineMatrix( 2, 0 );
    computedInitialCoefficients[ 1 ] = originalCosineMatrix( 2, 2 );
    computedInitialCoefficients[ 2 ] = originalSineMatrix( 2, 2 );
    std::cout << "computedInitialCoefficients " << computedInitialCoefficients.transpose( ) << std::endl;

    Eigen::Vector3d computedInitialTransientCoefficients = Eigen::Vector3d::Zero( );
    computedInitialTransientCoefficients = ( 1.0 / ( globalRelaxationTime - maxwellRelaxationTime ) ) *
            ( globalRelaxationTime * computedInitialCoefficients - maxwellRelaxationTime * computedEquilibriumCoefficients );
    // std::cout << "computedInitialTransientCoefficients: " << computedInitialTransientCoefficients.transpose( ) << std::endl;

    std::shared_ptr< GravityDeformationPropagatorSettings<> > gravityPropagatorSettings =
            std::make_shared< GravityDeformationPropagatorSettings<> >(
                    bodiesToPropagate,
                    gravityDeformationModels,
                    computedInitialCoefficients,
                    integratorSettings,
                    std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    // Translational dynamics propagator
    std::vector< std::string > centralBodies = { "Jupiter" };
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings = std::make_shared< SingleArcPropagatorProcessingSettings >( );
    outputSettings->setIntegratedResult( false );

    AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
    Eigen::Vector6d initialState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
    std::shared_ptr< TranslationalStatePropagatorSettings<> > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings<> >( centralBodies,
                                                                        accelerationsMap,
                                                                        bodiesToPropagate,
                                                                        initialState,
                                                                        initialTime,
                                                                        integratorSettings,
                                                                        std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
                                                                        cowell );

    // Create torque models
    SelectedTorqueMap torqueSettings;
    torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000", bodies, initialTime );
    // std::cout << "initialRotationState " << initialRotationState << std::endl;

    // Create propagator settings for rotational dynamics
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >(
                    torqueModelMap,
                    bodiesToPropagate,
                    initialRotationState,
                    initialTime,
                    integratorSettings,
                    std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    propagatorSettingsList.push_back( gravityPropagatorSettings );
    std::shared_ptr< MultiTypePropagatorSettings<> > fullPropagatorSettings =
            std::make_shared< MultiTypePropagatorSettings<> >( propagatorSettingsList,
                                                               integratorSettings,
                                                               initialTime,
                                                               std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
                                                               std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
                                                               outputSettings );

    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, fullPropagatorSettings );
    std::map< double, Eigen::VectorXd > results = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    std::cout << "initialState " << initialState.segment( 0, 3 ).norm( ) << std::endl;

    std::cout << "size results: " << results.size( ) << std::endl;
    int stateSize = results.begin( )->second.size( );
    for( auto it : results )
    {
        std::cout << it.first << " " << it.second.segment( 0, 3 ).norm( ) << " " << it.second.segment( stateSize - 3, 3 ).transpose( )
                  << std::endl;
    }

    std::cout << "computedInitialCoefficients " << computedInitialCoefficients.transpose( ) << std::endl;

    std::map< double, Eigen::Vector3d > computedSolution;
    for( double time = initialTime; time <= finalTime; time += timeStep )
    {
        Eigen::Vector3d currentTransientCoefficients =
                ( computedInitialTransientCoefficients - computedEquilibriumCoefficients ) * std::exp( -time / globalRelaxationTime ) +
                computedEquilibriumCoefficients;
        // std::cout << "currentTransientCoefficients " << currentTransientCoefficients.transpose( ) << std::endl;

        Eigen::Vector3d currentNominalCoefficients = ( 1.0 - maxwellRelaxationTime / globalRelaxationTime ) * currentTransientCoefficients +
                maxwellRelaxationTime / globalRelaxationTime * computedEquilibriumCoefficients;
        std::cout << "currentNominalCoefficients " << time << " "
                  << ( currentNominalCoefficients - results.at( time ).segment( stateSize - 3, 3 ) ).transpose( ) << std::endl;
    }
}

}  // namespace unit_tests
}  // namespace tudat
