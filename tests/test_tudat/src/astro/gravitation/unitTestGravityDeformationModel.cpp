/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::numerical_integrators;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;

namespace
{

const double initialTime = 0.0;
const double timeStep = 3600.0;
// Keep the test long enough to exercise the coupled deformation update over a meaningful interval.
// With twelve propagation cases, 1,080 one-hour RK4 steps gives a total test runtime of approximately 14 seconds.
const int numberOfPropagationSteps = 1080;
const double finalTime = initialTime + timeStep * numberOfPropagationSteps;
const double moonLoveNumber = 0.024059;
// The 45-day propagation spans three global-relaxation e-folding times. The elastic relaxation
// time is shorter by a factor of six; it does not affect the static-equilibrium solution because
// the equilibrium-coefficient derivative is zero there.
const double moonRelaxationTime = 15.0 * physical_constants::JULIAN_DAY;
const double moonElasticRelaxationTime = 2.5 * physical_constants::JULIAN_DAY;

enum MoonEnvironmentCase { static_moon_environment, circular_synchronous_moon_environment, spice_moon_environment };

struct GravityDeformationPropagationResults {
    std::map< double, Eigen::VectorXd > stateHistory_;
    std::map< double, Eigen::VectorXd > dependentVariableHistory_;
    Eigen::VectorXd nominalCoefficients_;
    double moonReferenceRadius_;
    double moonGravitationalParameter_;
    double earthGravitationalParameter_;
};

Eigen::Matrix3d getRotationToGlobalFrameWithAxisTowardsDirection( const Eigen::Vector3d& direction, const int bodyFixedAxis )
{
    const Eigen::Vector3d normalizedDirection = direction.normalized( );
    const Eigen::Vector3d perpendicularDirection = normalizedDirection.unitOrthogonal( );
    Eigen::Matrix3d rotationToGlobalFrame;

    switch( bodyFixedAxis )
    {
        case 0:
            rotationToGlobalFrame.col( 0 ) = normalizedDirection;
            rotationToGlobalFrame.col( 1 ) = perpendicularDirection;
            rotationToGlobalFrame.col( 2 ) = normalizedDirection.cross( perpendicularDirection );
            break;
        case 1:
            rotationToGlobalFrame.col( 1 ) = normalizedDirection;
            rotationToGlobalFrame.col( 2 ) = perpendicularDirection;
            rotationToGlobalFrame.col( 0 ) = normalizedDirection.cross( perpendicularDirection );
            break;
        case 2:
            rotationToGlobalFrame.col( 2 ) = normalizedDirection;
            rotationToGlobalFrame.col( 0 ) = perpendicularDirection;
            rotationToGlobalFrame.col( 1 ) = normalizedDirection.cross( perpendicularDirection );
            break;
        default:
            throw std::runtime_error( "Error when creating Moon orientation, requested body-fixed axis is invalid." );
    }
    return rotationToGlobalFrame;
}

Eigen::VectorXd getUnnormalizedMoonDegreeTwoCoefficients(
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
        const Eigen::VectorXd& staticCoefficients )
{
    Eigen::VectorXd nominalCoefficients = Eigen::VectorXd::Zero( 5 );
    const Eigen::MatrixXd cosineCoefficients = gravityField->getCosineCoefficients( );
    const Eigen::MatrixXd sineCoefficients = gravityField->getSineCoefficients( );

    nominalCoefficients << cosineCoefficients( 2, 0 ), cosineCoefficients( 2, 1 ), cosineCoefficients( 2, 2 ), sineCoefficients( 2, 1 ),
            sineCoefficients( 2, 2 );
    const Eigen::VectorXd normalizationFactors =
            ( Eigen::VectorXd( 5 ) << basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
              basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ),
              basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ),
              basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ),
              basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) )
                    .finished( );

    // The deformation settings accept normalized coefficients, whereas the Maxwell model works internally with
    // unnormalized coefficients. Apply the same conversion before removing the static contribution.
    return nominalCoefficients.cwiseProduct( normalizationFactors ) - staticCoefficients.cwiseProduct( normalizationFactors );
}

Eigen::VectorXd calculateEquilibriumCoefficients( const Eigen::Vector3d& relativePositionOfEarthWrtMoonInBodyFixedFrame,
                                                  const double moonReferenceRadius,
                                                  const double moonGravitationalParameter,
                                                  const double earthGravitationalParameter )
{
    const double distance = relativePositionOfEarthWrtMoonInBodyFixedFrame.norm( );
    const double sineLatitude = relativePositionOfEarthWrtMoonInBodyFixedFrame.z( ) / distance;
    const double cosineLatitude = std::sqrt( 1.0 - sineLatitude * sineLatitude );
    const double longitude =
            std::atan2( relativePositionOfEarthWrtMoonInBodyFixedFrame.y( ), relativePositionOfEarthWrtMoonInBodyFixedFrame.x( ) );
    const double tidalCoefficient =
            moonLoveNumber * earthGravitationalParameter / moonGravitationalParameter * std::pow( moonReferenceRadius / distance, 3 );

    Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero( 5 );
    equilibriumCoefficients( 0 ) = 0.5 * tidalCoefficient * ( 3.0 * sineLatitude * sineLatitude - 1.0 );
    equilibriumCoefficients( 1 ) = tidalCoefficient * cosineLatitude * sineLatitude * std::cos( longitude );
    equilibriumCoefficients( 2 ) = 0.25 * tidalCoefficient * cosineLatitude * cosineLatitude * std::cos( 2.0 * longitude );
    equilibriumCoefficients( 3 ) = tidalCoefficient * cosineLatitude * sineLatitude * std::sin( longitude );
    equilibriumCoefficients( 4 ) = 0.25 * tidalCoefficient * cosineLatitude * cosineLatitude * std::sin( 2.0 * longitude );
    return equilibriumCoefficients;
}

//! Compare vectors fractionally, except for coefficients whose expected scale is numerically zero.
void checkVectorCloseFractionOrAbsolute( const Eigen::VectorXd& actualValues,
                                         const Eigen::VectorXd& expectedValues,
                                         const double fractionalTolerance,
                                         const double absoluteTolerance )
{
    BOOST_REQUIRE_EQUAL( actualValues.size( ), expectedValues.size( ) );
    for( int i = 0; i < actualValues.size( ); ++i )
    {
        if( std::max( std::abs( actualValues( i ) ), std::abs( expectedValues( i ) ) ) < absoluteTolerance )
        {
            BOOST_CHECK_SMALL( actualValues( i ) - expectedValues( i ), absoluteTolerance );
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION( actualValues( i ), expectedValues( i ), fractionalTolerance );
        }
    }
}

//! Compare vectors using a combined absolute-plus-relative tolerance.
void checkVectorCloseCombined( const Eigen::VectorXd& actualValues,
                               const Eigen::VectorXd& expectedValues,
                               const double fractionalTolerance,
                               const double absoluteTolerance )
{
    BOOST_REQUIRE_EQUAL( actualValues.size( ), expectedValues.size( ) );
    for( int i = 0; i < actualValues.size( ); ++i )
    {
        const double scale = std::max( std::abs( actualValues( i ) ), std::abs( expectedValues( i ) ) );
        BOOST_CHECK_LE( std::abs( actualValues( i ) - expectedValues( i ) ), absoluteTolerance + fractionalTolerance * scale );
    }
}

GravityDeformationPropagationResults propagateMoonGravityDeformation(
        const MoonEnvironmentCase environmentCase,
        const double relaxationTime,
        const double elasticRelaxationTime,
        const Eigen::VectorXd& staticCoefficients,
        const Eigen::Matrix3d& staticRotationToGlobalFrame = Eigen::Matrix3d::Identity( ) )
{
    std::vector< std::string > bodyNames = { "Earth", "Moon" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialTime - timeStep, finalTime + timeStep, "SSB", "J2000" );

    if( environmentCase == static_moon_environment )
    {
        Eigen::Vector6d moonState = spice_interface::getBodyCartesianStateAtEpoch( "Moon", "Earth", "J2000", "NONE", initialTime );
        moonState.segment( 3, 3 ).setZero( );
        bodySettings.at( "Earth" )->ephemerisSettings =
                std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
        bodySettings.at( "Moon" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( moonState, "SSB", "J2000" );
        bodySettings.at( "Moon" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                // SimpleRotationModelSettings takes the rotation from the base frame to the target frame. The
                // supplied matrix has body-fixed axes as columns in the global frame, so its inverse is required.
                "J2000",
                "MoonFixed",
                Eigen::Quaterniond( staticRotationToGlobalFrame.transpose( ) ),
                initialTime,
                0.0 );
    }
    else if( environmentCase == circular_synchronous_moon_environment )
    {
        const double orbitalRadius = 384400000.0;
        const double meanMotion = std::sqrt( physical_constants::GRAVITATIONAL_CONSTANT * 5.97219e24 / std::pow( orbitalRadius, 3 ) );
        Eigen::Vector6d initialKeplerianState = Eigen::Vector6d::Zero( );
        initialKeplerianState( 0 ) = orbitalRadius;

        bodySettings.at( "Earth" )->ephemerisSettings =
                std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
        bodySettings.at( "Moon" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                initialKeplerianState, initialTime, physical_constants::GRAVITATIONAL_CONSTANT * 5.97219e24, "Earth", "J2000" );
        bodySettings.at( "Moon" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000",
                "MoonFixed",
                Eigen::Quaterniond( ( Eigen::Matrix3d( ) << -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0 ).finished( ) ),
                initialTime,
                meanMotion );
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > moonGravityField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodies.at( "Moon" )->getGravityFieldModel( ) );
    BOOST_REQUIRE( moonGravityField != nullptr );
    const Eigen::VectorXd initialNominalCoefficients = getUnnormalizedMoonDegreeTwoCoefficients( moonGravityField, staticCoefficients );

    SelectedGravityDeformationModelMap deformationSettings;
    deformationSettings[ "Moon" ] = { maxwellDeformationSettings(
            elasticRelaxationTime, relaxationTime, moonLoveNumber, 2, 2, "Earth", staticCoefficients ) };
    const basic_astrodynamics::GravityDeformationModelMap deformationModels =
            createGravityDeformationModelsMap( bodies, deformationSettings );

    const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables = {
        gravityDeformationStateDerivativeDependentVariable( basic_astrodynamics::maxwell_deformation, "Moon" ),
        maxwellGravityDeformationEquilibriumCoefficientsDependentVariable( "Moon" ),
        maxwellGravityDeformationEquilibriumCoefficientDerivativeDependentVariable( "Moon" ),
        std::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable, "Moon", "Earth" ),
        std::make_shared< SingleDependentVariableSaveSettings >( inertial_to_body_fixed_rotation_matrix_variable, "Moon" )
    };

    const std::shared_ptr< IntegratorSettings<> > integratorSettings =
            std::make_shared< IntegratorSettings<> >( rungeKutta4, initialTime, timeStep );
    const std::shared_ptr< GravityDeformationPropagatorSettings<> > propagatorSettings =
            std::make_shared< GravityDeformationPropagatorSettings<> >( std::vector< std::string >( { "Moon" } ),
                                                                        deformationModels,
                                                                        Eigen::VectorXd::Zero( 5 ),
                                                                        integratorSettings,
                                                                        std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
                                                                        dependentVariables );

    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );
    return { dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
             dynamicsSimulator.getDependentVariableHistory( ),
             initialNominalCoefficients,
             moonGravityField->getReferenceRadius( ),
             moonGravityField->getGravitationalParameter( ),
             bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) };
}

void checkSavedDependentVariables( const GravityDeformationPropagationResults& propagationResults )
{
    BOOST_REQUIRE_EQUAL( propagationResults.stateHistory_.size( ), numberOfPropagationSteps + 1 );
    BOOST_REQUIRE_EQUAL( propagationResults.dependentVariableHistory_.size( ), numberOfPropagationSteps + 1 );

    for( const auto& stateEntry : propagationResults.stateHistory_ )
    {
        const Eigen::VectorXd& dependentVariables = propagationResults.dependentVariableHistory_.at( stateEntry.first );
        const Eigen::VectorXd savedEquilibriumCoefficients = dependentVariables.segment( 5, 5 );
        const Eigen::Vector3d savedMoonPositionWrtEarth = dependentVariables.segment( 15, 3 );
        const Eigen::Matrix3d rotationToMoonFixedFrame = getMatrixFromVectorRotationRepresentation( dependentVariables.segment( 18, 9 ) );

        const Eigen::VectorXd directlyCalculatedEquilibriumCoefficients =
                calculateEquilibriumCoefficients( -rotationToMoonFixedFrame * savedMoonPositionWrtEarth,
                                                  propagationResults.moonReferenceRadius_,
                                                  propagationResults.moonGravitationalParameter_,
                                                  propagationResults.earthGravitationalParameter_ );
        // A pure relative comparison is ill-conditioned when a coefficient crosses zero. The reconstruction differs
        // only at round-off level, so use a combined tolerance for this independently calculated quantity.
        checkVectorCloseCombined( directlyCalculatedEquilibriumCoefficients, savedEquilibriumCoefficients, 5.0e-14, 1.0e-20 );
    }
}

//! Check the exact solution of the Maxwell deformation equations for constant equilibrium coefficients.
/*!
 * With \f$\dot{\Delta C}^{eq}=\dot{\Delta S}^{eq}=0\f$, Eqs. (74)-(75) of Fayolle et al. (2026) give
 * \f$\dot{x}=(x^{eq}-x)/\tau\f$. Hence the state derivative decays exponentially and the state approaches
 * equilibrium exponentially; neither is linear in time.
 */
void checkConstantEquilibriumRelaxationSolution( const GravityDeformationPropagationResults& propagationResults,
                                                 const double globalRelaxationTime )
{
    const Eigen::VectorXd initialState = propagationResults.stateHistory_.begin( )->second;
    const Eigen::VectorXd initialStateDerivative = propagationResults.dependentVariableHistory_.begin( )->second.segment( 0, 5 );

    for( const auto& stateEntry : propagationResults.stateHistory_ )
    {
        const double elapsedTime = stateEntry.first - propagationResults.stateHistory_.begin( )->first;
        const double decayFactor = std::exp( -elapsedTime / globalRelaxationTime );
        const Eigen::VectorXd expectedStateDerivative = decayFactor * initialStateDerivative;
        const Eigen::VectorXd expectedState = initialState + globalRelaxationTime * ( 1.0 - decayFactor ) * initialStateDerivative;
        const Eigen::VectorXd& currentDependentVariables = propagationResults.dependentVariableHistory_.at( stateEntry.first );

        // The derivative inherits the accumulated RK4 state error. For a one-hour step and a 15-day relaxation
        // time, its expected relative truncation error over three e-folds is of order 1e-12.
        checkVectorCloseCombined( currentDependentVariables.segment( 0, 5 ), expectedStateDerivative, 5.0e-12, 1.0e-28 );
        // The RK4 propagation of this coupled state agrees with the exact exponential solution to better than 1e-12.
        checkVectorCloseFractionOrAbsolute( stateEntry.second, expectedState, 1.0e-12, 1.0e-20 );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testStaticMoonGravityDeformationForThreeOrientations )
{
    spice_interface::loadStandardSpiceKernels( );
    const Eigen::Vector3d moonToJupiter =
            spice_interface::getBodyCartesianPositionAtEpoch( "Jupiter", "Moon", "J2000", "NONE", initialTime );

    for( int bodyFixedAxis = 0; bodyFixedAxis < 3; bodyFixedAxis++ )
    {
        const Eigen::Matrix3d rotationToGlobalFrame = getRotationToGlobalFrameWithAxisTowardsDirection( moonToJupiter, bodyFixedAxis );
        const GravityDeformationPropagationResults propagationResults = propagateMoonGravityDeformation(
                static_moon_environment, moonRelaxationTime, moonElasticRelaxationTime, Eigen::VectorXd::Zero( 5 ), rotationToGlobalFrame );
        checkSavedDependentVariables( propagationResults );

        const Eigen::VectorXd firstDependentVariables = propagationResults.dependentVariableHistory_.begin( )->second;
        const Eigen::Matrix3d savedRotationToMoonFixedFrame =
                getMatrixFromVectorRotationRepresentation( firstDependentVariables.segment( 18, 9 ) );
        BOOST_CHECK_CLOSE_FRACTION(
                savedRotationToMoonFixedFrame.transpose( ).col( bodyFixedAxis ).dot( moonToJupiter.normalized( ) ), 1.0, 1.0e-14 );
        BOOST_CHECK_SMALL( firstDependentVariables.segment( 10, 5 ).norm( ), 1.0e-25 );
        checkConstantEquilibriumRelaxationSolution( propagationResults, moonRelaxationTime );
    }
}

BOOST_AUTO_TEST_CASE( testCircularSynchronousMoonGravityDeformation )
{
    spice_interface::loadStandardSpiceKernels( );
    const GravityDeformationPropagationResults propagationResults = propagateMoonGravityDeformation(
            circular_synchronous_moon_environment, moonRelaxationTime, moonElasticRelaxationTime, Eigen::VectorXd::Zero( 5 ) );
    checkSavedDependentVariables( propagationResults );

    const Eigen::VectorXd firstDependentVariables = propagationResults.dependentVariableHistory_.begin( )->second;
    const Eigen::VectorXd lastDependentVariables = propagationResults.dependentVariableHistory_.rbegin( )->second;
    checkVectorCloseFractionOrAbsolute( firstDependentVariables.segment( 5, 5 ), lastDependentVariables.segment( 5, 5 ), 5.0e-13, 1.0e-20 );
    BOOST_CHECK_SMALL( firstDependentVariables.segment( 10, 5 ).norm( ), 1.0e-18 );
    BOOST_CHECK_SMALL( lastDependentVariables.segment( 10, 5 ).norm( ), 1.0e-18 );

    checkConstantEquilibriumRelaxationSolution( propagationResults, moonRelaxationTime );
}

BOOST_AUTO_TEST_CASE( testSpiceMoonGravityDeformationParameterSweep )
{
    spice_interface::loadStandardSpiceKernels( );
    Eigen::VectorXd moonStaticCoefficients = Eigen::VectorXd::Zero( 5 );
    moonStaticCoefficients( 0 ) = -9.09e-5;
    moonStaticCoefficients( 2 ) = 3.47e-5;

    const std::vector< double > relaxationTimes = { moonRelaxationTime, 2.0 * moonRelaxationTime };
    const std::vector< double > elasticRelaxationTimes = { 0.5 * moonElasticRelaxationTime, moonElasticRelaxationTime };
    const std::vector< Eigen::VectorXd > staticCoefficientSets = { Eigen::VectorXd::Zero( 5 ), moonStaticCoefficients };

    std::vector< GravityDeformationPropagationResults > propagationResults;
    for( const double relaxationTime : relaxationTimes )
    {
        for( const double elasticRelaxationTime : elasticRelaxationTimes )
        {
            for( const Eigen::VectorXd& staticCoefficients : staticCoefficientSets )
            {
                propagationResults.push_back( propagateMoonGravityDeformation(
                        spice_moon_environment, relaxationTime, elasticRelaxationTime, staticCoefficients ) );
                checkSavedDependentVariables( propagationResults.back( ) );
            }
        }
    }

    const Eigen::VectorXd firstSpiceDerivative = propagationResults.front( ).dependentVariableHistory_.begin( )->second.segment( 0, 5 );
    const Eigen::VectorXd changedRelaxationTimeDerivative =
            propagationResults.at( 4 ).dependentVariableHistory_.begin( )->second.segment( 0, 5 );
    const Eigen::VectorXd changedElasticRelaxationTimeDerivative =
            propagationResults.at( 2 ).dependentVariableHistory_.begin( )->second.segment( 0, 5 );
    const Eigen::VectorXd changedStaticCoefficientsDerivative =
            propagationResults.at( 1 ).dependentVariableHistory_.begin( )->second.segment( 0, 5 );
    const Eigen::VectorXd equilibriumCoefficientDerivative =
            propagationResults.front( ).dependentVariableHistory_.begin( )->second.segment( 10, 5 );
    const Eigen::VectorXd normalizedStaticCoefficients =
            ( Eigen::VectorXd( 5 ) << moonStaticCoefficients( 0 ) * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
              moonStaticCoefficients( 1 ) * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ),
              moonStaticCoefficients( 2 ) * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ),
              moonStaticCoefficients( 3 ) * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ),
              moonStaticCoefficients( 4 ) * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) )
                    .finished( );

    BOOST_CHECK_GT( firstSpiceDerivative.norm( ), 0.0 );

    // The global relaxation time scales the entire derivative, the Maxwell relaxation time multiplies only the
    // equilibrium-coefficient derivative, and static coefficients enter with the opposite sign in the nominal field.
    const Eigen::VectorXd expectedChangedRelaxationTimeDerivative = 0.5 * firstSpiceDerivative;
    const Eigen::VectorXd expectedElasticRelaxationTimeDerivativeDifference =
            0.5 * moonElasticRelaxationTime / moonRelaxationTime * equilibriumCoefficientDerivative;
    const Eigen::VectorXd elasticRelaxationTimeDerivativeDifference = changedElasticRelaxationTimeDerivative - firstSpiceDerivative;
    const Eigen::VectorXd expectedStaticCoefficientsDerivativeDifference = normalizedStaticCoefficients / moonRelaxationTime;
    const Eigen::VectorXd staticCoefficientsDerivativeDifference = changedStaticCoefficientsDerivative - firstSpiceDerivative;

    checkVectorCloseCombined( changedRelaxationTimeDerivative, expectedChangedRelaxationTimeDerivative, 5.0e-14, 1.0e-26 );
    checkVectorCloseCombined(
            elasticRelaxationTimeDerivativeDifference, expectedElasticRelaxationTimeDerivativeDifference, 5.0e-14, 1.0e-26 );
    checkVectorCloseCombined( staticCoefficientsDerivativeDifference, expectedStaticCoefficientsDerivativeDifference, 5.0e-14, 1.0e-26 );
}

BOOST_AUTO_TEST_CASE( testCoupledRotationalGravityDirectAndIterativeAgreement )
{
    spice_interface::loadStandardSpiceKernels( );

    const auto propagateCoupledDynamics = []( const bool useDirectSolution ) {
        const double coupledInitialTime = 0.0;
        const double coupledFinalTime = 600.0;
        const double coupledStep = 60.0;
        BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Moon" }, coupledInitialTime, coupledFinalTime, "SSB", "J2000" );
        Eigen::Vector6d moonState = spice_interface::getBodyCartesianStateAtEpoch( "Moon", "Earth", "J2000", "NONE", coupledInitialTime );
        moonState.segment( 3, 3 ).setZero( );
        bodySettings.at( "Earth" )->ephemerisSettings =
                std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
        bodySettings.at( "Moon" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( moonState, "SSB", "J2000" );
        bodySettings.at( "Moon" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000", "MoonFixed", Eigen::Quaterniond::Identity( ), coupledInitialTime, 0.0 );
        const std::shared_ptr< SphericalHarmonicsGravityFieldSettings > moonGravitySettings =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Moon" )->gravityFieldSettings );
        BOOST_REQUIRE( moonGravitySettings != nullptr );
        moonGravitySettings->setScaledMeanMomentOfInertia( 0.4 );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        SelectedGravityDeformationModelMap deformationSettings;
        deformationSettings[ "Moon" ] = { maxwellDeformationSettings(
                moonElasticRelaxationTime, moonRelaxationTime, moonLoveNumber, 2, 2, "Earth", Eigen::VectorXd::Zero( 5 ), true, true ) };
        const basic_astrodynamics::GravityDeformationModelMap deformationModels =
                createGravityDeformationModelsMap( bodies, deformationSettings );

        SelectedTorqueMap torqueSettings;
        torqueSettings[ "Moon" ] = {};
        const basic_astrodynamics::TorqueModelMap torqueModels = createTorqueModelsMap( bodies, torqueSettings, { "Moon" } );

        Eigen::VectorXd initialRotation = Eigen::VectorXd::Zero( 7 );
        initialRotation.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond::Identity( ) );
        initialRotation.segment( 4, 3 ) = ( Eigen::Vector3d( ) << 1.0e-4, -2.0e-4, 3.0e-4 ).finished( );

        const std::shared_ptr< IntegratorSettings<> > integratorSettings =
                std::make_shared< IntegratorSettings<> >( rungeKutta4, coupledInitialTime, coupledStep );
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings =
                std::make_shared< PropagationTimeTerminationSettings >( coupledFinalTime );
        const std::shared_ptr< RotationalStatePropagatorSettings<> > rotationalSettings =
                std::make_shared< RotationalStatePropagatorSettings<> >( torqueModels,
                                                                         std::vector< std::string >( { "Moon" } ),
                                                                         initialRotation,
                                                                         coupledInitialTime,
                                                                         integratorSettings,
                                                                         terminationSettings );
        const std::shared_ptr< GravityDeformationPropagatorSettings<> > gravitySettings =
                std::make_shared< GravityDeformationPropagatorSettings<> >( std::vector< std::string >( { "Moon" } ),
                                                                            deformationModels,
                                                                            Eigen::VectorXd::Zero( 5 ),
                                                                            coupledInitialTime,
                                                                            integratorSettings,
                                                                            terminationSettings );
        const std::shared_ptr< MultiTypePropagatorSettings<> > settings = std::make_shared< MultiTypePropagatorSettings<> >(
                std::vector< std::shared_ptr< SingleArcPropagatorSettings<> > >( { rotationalSettings, gravitySettings } ),
                integratorSettings,
                coupledInitialTime,
                terminationSettings );
        settings->setCoupledStateDerivativeSolverSettings( std::make_shared< CoupledStateDerivativeSolverSettings >(
                useDirectSolution, 1.0e-12, 1.0e-14, 50, throw_exception_on_coupled_derivative_failure ) );

        SingleArcDynamicsSimulator<> simulator( bodies, settings );
        return simulator.getEquationsOfMotionNumericalSolution( );
    };

    const std::map< double, Eigen::VectorXd > directHistory = propagateCoupledDynamics( true );
    const std::map< double, Eigen::VectorXd > iterativeHistory = propagateCoupledDynamics( false );
    BOOST_REQUIRE_EQUAL( directHistory.size( ), 11 );
    BOOST_REQUIRE_EQUAL( iterativeHistory.size( ), directHistory.size( ) );
    BOOST_CHECK_GT( directHistory.rbegin( )->second.tail( 5 ).norm( ), 0.0 );

    for( const auto& directEntry : directHistory )
    {
        checkVectorCloseCombined( directEntry.second, iterativeHistory.at( directEntry.first ), 2.0e-10, 1.0e-16 );
    }
}

}  // namespace unit_tests
}  // namespace tudat
