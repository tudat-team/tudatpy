/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/equinoctialElementConversions.h"
#include "tudat/astro/ephemerides/tleFitting.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

//! Convert a Tudat aerodynamic area-to-mass ratio to the corresponding SGP4 B* convention.
double calculateBStarFromAerodynamicProperties( const double referenceArea, const double dragCoefficient, const double mass )
{
    // Vallado et al., Revisiting Spacetrack Report No. 3, AIAA 2006-6753, Appendix B gives
    // B* = 0.5 rho_0 C_D A / m, with rho_0 = 2.461e-5 kg/m^2 in the SGP4 Earth-radius convention.
    constexpr double sgp4ReferenceDensity = 2.461e-5;
    return 0.5 * sgp4ReferenceDensity * dragCoefficient * referenceArea / mass;
}

//! Propagate a TLE-seeded LEO with a force model intentionally richer than SGP4.
/*!
 * Earth gravity is evaluated to degree/order 16. Aerodynamic drag, cannonball solar-radiation pressure, and point-mass
 * Sun/Moon perturbations make the resulting history independent of the SGP4 forward model used by the fitter.
 */
std::map< double, Eigen::Vector6d > propagateTleWithExtendedDynamics( const std::shared_ptr< ephemerides::Tle >& initialTle,
                                                                      const double finalEpoch,
                                                                      const double aerodynamicReferenceArea,
                                                                      const double spacecraftMass,
                                                                      const double dragCoefficient,
                                                                      const double radiationPressureCoefficient )
{
    using namespace basic_astrodynamics;
    using namespace ephemerides;
    using namespace numerical_integrators;
    using namespace propagators;
    using namespace simulation_setup;

    const double initialEpoch = initialTle->getEpoch( );

    // Create Earth, Sun and Moon ephemerides over the complete integration interval.
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Sun", "Moon" }, initialEpoch - 300.0, finalEpoch + 300.0, "Earth", "J2000" );

    // Use a simple thermosphere model with approximately 4e-12 kg/m^3 density at 400 km.
    bodySettings.at( "Earth" )->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( 60.0e3, 1000.0, 3.1e-9 );

    // Apply identical reference area to drag and cannonball radiation pressure while retaining separate coefficients.
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->constantMass = spacecraftMass;
    bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings =
            constantAerodynamicCoefficientSettings( aerodynamicReferenceArea, dragCoefficient * Eigen::Vector3d::UnitX( ) );
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            cannonballRadiationPressureTargetModelSettings( aerodynamicReferenceArea, radiationPressureCoefficient );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Combine the extended force model requested for both independent numerical-propagation fitting tests.
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Vehicle" ][ "Earth" ] = { sphericalHarmonicAcceleration( 16, 16 ), aerodynamicAcceleration( ) };
    accelerationSettings[ "Vehicle" ][ "Sun" ] = { pointMassGravityAcceleration( ), radiationPressureAcceleration( ) };
    accelerationSettings[ "Vehicle" ][ "Moon" ] = { pointMassGravityAcceleration( ) };
    const std::vector< std::string > bodiesToPropagate = { "Vehicle" };
    const std::vector< std::string > centralBodies = { "Earth" };
    const AccelerationMap accelerationModels =
            createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    // Seed the numerical propagation from the exact Cartesian state represented by the input TLE at its own epoch.
    TleEphemeris initialTleEphemeris( "Earth", "J2000", initialTle );
    const Eigen::Vector6d initialState = initialTleEphemeris.getCartesianState( initialEpoch );

    // Integrate at 30-second fixed steps and retain a state every 30 minutes for the batch fit.
    const auto propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies, accelerationModels, bodiesToPropagate, initialState, finalEpoch );
    const auto integratorSettings = std::make_shared< IntegratorSettings<> >( rungeKutta4, initialEpoch, 30.0 );
    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, integratorSettings, propagatorSettings );

    std::map< double, Eigen::Vector6d > stateHistory;
    unsigned int outputIndex = 0;
    const auto& numericalStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( auto stateIterator = numericalStateHistory.begin( ); stateIterator != numericalStateHistory.end( );
         ++stateIterator, ++outputIndex )
    {
        if( outputIndex % 60 == 0 || std::next( stateIterator ) == numericalStateHistory.end( ) )
        {
            stateHistory[ stateIterator->first ] = stateIterator->second;
        }
    }
    return stateHistory;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_tle_fitting )

//! Verify that the extracted generic equinoctial conversions are mutual inverses in both inclination branches.
/*!
 * The fitter relies on the prograde branch for ordinary inclinations and the cotangent-based retrograde branch near
 * pi. Testing both here isolates the element conversion from SGP4 propagation and nonlinear least squares.
 */
BOOST_AUTO_TEST_CASE( testEquinoctialElementConversionRoundTrip )
{
    using namespace orbital_element_conversions;

    for( const bool useRetrogradeElements : { false, true } )
    {
        // Use the same nonsingular classical parameters in each branch, changing only inclination to select the
        // hemisphere for which that branch is intended.
        Eigen::Vector6d orbitalElements;
        orbitalElements << 0.0676, 0.01, useRetrogradeElements ? 2.8 : 0.9, 0.5, 0.7, 0.2;

        const Eigen::Vector6d equinoctialElements =
                convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements( orbitalElements, useRetrogradeElements );
        const Eigen::Vector6d recoveredOrbitalElements =
                convertEquinoctialToOrbitalElementsWithMeanMotionAndMeanAnomaly( equinoctialElements, useRetrogradeElements );

        // Check that mean motion, shape, inclination and all normalized angles survive the selected conversion round trip.
        BOOST_CHECK_SMALL( ( recoveredOrbitalElements - orbitalElements ).norm( ), 1.0e-14 );
    }
}

//! Verify an exact position-only round trip using the same full-precision SGP4 model for data generation and fitting.
/*!
 * The reference trajectory is generated from a numerical (unserialized) TLE. The initial solve-for TLE is deliberately
 * perturbed in all six orbital parameters while B* is held fixed. This isolates the batch least-squares algorithm and
 * confirms that neither formatted-TLE quantization nor velocity weighting is involved.
 */
BOOST_AUTO_TEST_CASE( testNumericalTlePositionOnlyRoundTrip )
{
    using namespace ephemerides;

    // CSPICE SGP4 requires the leap-second kernel even though all test epochs are already expressed as seconds since J2000.
    spice_interface::loadStandardSpiceKernels( );

    // Define an ordinary, mildly eccentric LEO using the full-precision numerical TLE constructor.
    const double tleEpoch = 1.0e6;
    const std::shared_ptr< Tle > referenceTle = std::make_shared< Tle >( tleEpoch, 2.0e-5, 0.9, 0.5, 0.01, 0.7, 0.2, 0.0676 );

    // Sample approximately one and a half orbital periods at ten-minute intervals.
    TleEphemeris referenceEphemeris( "Earth", "J2000", referenceTle );
    std::map< double, Eigen::Vector6d > stateHistory;
    for( unsigned int i = 0; i < 16; ++i )
    {
        const double epoch = tleEpoch + static_cast< double >( i ) * 600.0;
        stateHistory[ epoch ] = referenceEphemeris.getCartesianState( epoch );
    }

    // Perturb every orbital element enough to exercise nonlinear correction, while preserving the reference B*.
    TleFitSettings settings;
    settings.tleEpoch_ = tleEpoch;
    settings.initialTle_ = std::make_shared< Tle >( tleEpoch, 2.0e-5, 0.905, 0.505, 0.0105, 0.695, 0.205, 0.0676 * 1.0001 );
    settings.maximumNumberOfIterations_ = 50;
    settings.convergenceTolerance_ = 1.0e-10;

    const TleFitResult fitResult = fitTleToCartesianStateHistory( stateHistory, "J2000", settings );

    // Check that the fitter returned a usable full-precision TLE object.
    BOOST_REQUIRE( fitResult.fittedTle_ != nullptr );

    // Check that the least-squares iterations reduced the position RMS relative to the deliberately perturbed initial TLE.
    BOOST_CHECK( fitResult.positionRms_ < fitResult.initialPositionRms_ );

    // Check that the continuous fitted trajectory reproduces the synthetic positions to better than 2 cm RMS.
    BOOST_CHECK_SMALL( fitResult.positionRms_, 2.0e-2 );

    // Check that B* remains exactly fixed when estimateBStar_ is left disabled.
    BOOST_CHECK_SMALL( fitResult.fittedTle_->getBStar( ) - referenceTle->getBStar( ), 1.0e-15 );

    // Independently re-propagate the returned TLE and check every input epoch, not only the aggregate RMS diagnostic.
    TleEphemeris fittedEphemeris( "Earth", "J2000", fitResult.fittedTle_ );
    for( const auto& stateAtEpoch : stateHistory )
    {
        // Check that the fitted SGP4 position at this individual epoch differs from the reference by less than 5 cm.
        BOOST_CHECK_SMALL(
                ( fittedEphemeris.getCartesianState( stateAtEpoch.first ).segment( 0, 3 ) - stateAtEpoch.second.segment( 0, 3 ) ).norm( ),
                5.0e-2 );
    }
}

//! Verify the Cartesian/MEE initializer and both public input-frame orientations.
/*!
 * No initial TLE is supplied, so the fitter must construct one from the Cartesian state nearest the requested epoch by
 * way of Tudat's modified-equinoctial conversions. Repeating the same case in J2000 and ECLIPJ2000 verifies that model
 * positions and input positions are compared in the same selected frame.
 */
BOOST_AUTO_TEST_CASE( testCartesianInitializerAndSupportedFrames )
{
    using namespace ephemerides;

    // Load kernels once for the CSPICE SGP4 evaluations in both frame sub-cases.
    spice_interface::loadStandardSpiceKernels( );

    // Use a nonsingular inclined LEO as the common reference for both orientations.
    const double tleEpoch = 2.0e6;
    const std::shared_ptr< Tle > referenceTle = std::make_shared< Tle >( tleEpoch, 0.0, 1.1, 1.5, 0.02, 2.0, 0.8, 0.066 );

    for( const std::string& frameOrientation : { std::string( "J2000" ), std::string( "ECLIPJ2000" ) } )
    {
        // Generate each state history directly in the orientation passed to the fitter.
        TleEphemeris referenceEphemeris( "Earth", frameOrientation, referenceTle );
        std::map< double, Eigen::Vector6d > stateHistory;
        for( unsigned int i = 0; i < 20; ++i )
        {
            const double epoch = tleEpoch - 1800.0 + static_cast< double >( i ) * 300.0;
            stateHistory[ epoch ] = referenceEphemeris.getCartesianState( epoch );
        }

        // Deliberately leave initialTle_ null and request an epoch between samples. This exercises both the
        // Cartesian/MEE initializer and its requirement to attach the derived elements to the nearest state's epoch.
        TleFitSettings settings;
        settings.tleEpoch_ = tleEpoch + 121.0;
        settings.maximumNumberOfIterations_ = 60;
        settings.convergenceTolerance_ = 1.0e-10;
        const TleFitResult fitResult = fitTleToCartesianStateHistory( stateHistory, frameOrientation, settings );

        // Check that automatic initialization and fitting returned a TLE in this frame sub-case.
        BOOST_REQUIRE( fitResult.fittedTle_ != nullptr );

        // Check that the requested between-sample epoch was snapped to the nearest state at the reference TLE epoch.
        BOOST_CHECK_SMALL( fitResult.fittedTle_->getEpoch( ) - tleEpoch, 1.0e-15 );

        // Check that the batch fit improves on the single-state Cartesian/MEE initial approximation.
        BOOST_CHECK( fitResult.positionRms_ < fitResult.initialPositionRms_ );

        // Check that the automatically initialized fit reproduces the synthetic arc to sub-metre position RMS.
        BOOST_CHECK_SMALL( fitResult.positionRms_, 5.0e-1 );
    }
}

//! Verify that malformed inputs are rejected before any SGP4 or least-squares work is attempted.
BOOST_AUTO_TEST_CASE( testTleFitInputValidation )
{
    using namespace ephemerides;

    // A single position supplies only three scalar observations and cannot determine six orbital parameters.
    std::map< double, Eigen::Vector6d > stateHistory;
    stateHistory[ 0.0 ] = Eigen::Vector6d::Zero( );

    // Check that fewer than two Cartesian states produces the documented invalid-argument failure.
    BOOST_CHECK_THROW( fitTleToCartesianStateHistory( stateHistory, "J2000" ), std::invalid_argument );

    // Add a second finite state so that only the unsupported frame orientation remains invalid.
    stateHistory[ 1.0 ] = Eigen::Vector6d::Ones( );

    // Check that an unsupported TEME input request is rejected; the public interface accepts J2000 and ECLIPJ2000 only.
    BOOST_CHECK_THROW( fitTleToCartesianStateHistory( stateHistory, "TEME" ), std::invalid_argument );
}

//! Verify recovery of a physically consistent LEO TLE from an extended Tudat propagation.
/*!
 * A numerical TLE defines the initial Cartesian state near 400 km. Its B* is calculated from the same drag coefficient,
 * reference area and mass used by a one-day propagation with 16x16 Earth gravity, drag, cannonball radiation pressure,
 * and Sun/Moon perturbations. The history is fitted both with fixed B* and estimated B*. Both fits should remain close
 * to the seed mean elements, and adding B* as a solve-for should improve the trajectory fit.
 */
BOOST_AUTO_TEST_CASE( testFitRemainsCloseToTleForExtendedLeoPropagation )
{
    using namespace ephemerides;
    using namespace orbital_element_conversions;

    spice_interface::loadStandardSpiceKernels( );

    // Define one consistent spacecraft in both the SGP4 B* convention and the Tudat force model.
    const double initialEpoch = 2.5e6;
    const double aerodynamicReferenceArea = 5.0;
    const double spacecraftMass = 100.0;
    const double dragCoefficient = 2.2;
    const double radiationPressureCoefficient = 1.2;
    const double consistentBStar = calculateBStarFromAerodynamicProperties( aerodynamicReferenceArea, dragCoefficient, spacecraftMass );

    // Set mean motion for a roughly circular 400 km LEO using the WGS-72 constants employed by SGP4.
    const double initialSemiMajorAxis = 6378.135e3 + 400.0e3;
    const double initialMeanMotion = std::sqrt( 398600.8e9 / std::pow( initialSemiMajorAxis, 3.0 ) ) * 60.0;
    const std::shared_ptr< Tle > referenceTle =
            std::make_shared< Tle >( initialEpoch, consistentBStar, 0.9, 0.4, 0.001, 0.7, 0.2, initialMeanMotion );

    // Generate an independent one-day state history using all requested perturbations.
    const std::map< double, Eigen::Vector6d > stateHistory =
            propagateTleWithExtendedDynamics( referenceTle,
                                              initialEpoch + physical_constants::JULIAN_DAY,
                                              aerodynamicReferenceArea,
                                              spacecraftMass,
                                              dragCoefficient,
                                              radiationPressureCoefficient );

    // First fit the six mean orbital elements while holding the physically derived seed B* fixed.
    TleFitSettings fixedBStarSettings;
    fixedBStarSettings.initialTle_ = referenceTle;
    fixedBStarSettings.maximumNumberOfIterations_ = 50;
    fixedBStarSettings.convergenceTolerance_ = 1.0e-9;
    const TleFitResult fixedBStarResult = fitTleToCartesianStateHistory( stateHistory, "J2000", fixedBStarSettings );

    // Repeat the same fit with B* added as the seventh solve-for parameter.
    TleFitSettings estimatedBStarSettings = fixedBStarSettings;
    estimatedBStarSettings.estimateBStar_ = true;
    const TleFitResult estimatedBStarResult = fitTleToCartesianStateHistory( stateHistory, "J2000", estimatedBStarSettings );

    // Check that both versions of the extended-dynamics fit returned usable full-precision TLE objects.
    BOOST_REQUIRE( fixedBStarResult.fittedTle_ != nullptr && estimatedBStarResult.fittedTle_ != nullptr );

    Eigen::Vector6d referenceOrbitalElements;
    referenceOrbitalElements << referenceTle->getMeanMotion( ), referenceTle->getEccentricity( ), referenceTle->getInclination( ),
            referenceTle->getRightAscension( ), referenceTle->getArgOfPerigee( ), referenceTle->getMeanAnomaly( );
    Eigen::Vector6d fixedBStarOrbitalElements;
    fixedBStarOrbitalElements << fixedBStarResult.fittedTle_->getMeanMotion( ), fixedBStarResult.fittedTle_->getEccentricity( ),
            fixedBStarResult.fittedTle_->getInclination( ), fixedBStarResult.fittedTle_->getRightAscension( ),
            fixedBStarResult.fittedTle_->getArgOfPerigee( ), fixedBStarResult.fittedTle_->getMeanAnomaly( );
    Eigen::Vector6d estimatedBStarOrbitalElements;
    estimatedBStarOrbitalElements << estimatedBStarResult.fittedTle_->getMeanMotion( ), estimatedBStarResult.fittedTle_->getEccentricity( ),
            estimatedBStarResult.fittedTle_->getInclination( ), estimatedBStarResult.fittedTle_->getRightAscension( ),
            estimatedBStarResult.fittedTle_->getArgOfPerigee( ), estimatedBStarResult.fittedTle_->getMeanAnomaly( );
    const Eigen::Vector6d referenceEquinoctialElements =
            convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements( referenceOrbitalElements, false );
    const Eigen::Vector6d fixedBStarEquinoctialElements =
            convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements( fixedBStarOrbitalElements, false );
    const Eigen::Vector6d estimatedBStarEquinoctialElements =
            convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements( estimatedBStarOrbitalElements, false );

    BOOST_TEST_MESSAGE( "Extended-force-model TLE recovery: initial RMS = "
                        << fixedBStarResult.initialPositionRms_ << " m, fixed-B* RMS = " << fixedBStarResult.positionRms_
                        << " m, estimated-B* RMS = " << estimatedBStarResult.positionRms_ << " m, fixed-B* element change = "
                        << ( fixedBStarEquinoctialElements - referenceEquinoctialElements ).norm( ) << ", estimated-B* element change = "
                        << ( estimatedBStarEquinoctialElements - referenceEquinoctialElements ).norm( ) );

    // Check that the six-parameter fit reduces the seed-TLE position RMS by at least a factor of five.
    BOOST_CHECK( fixedBStarResult.positionRms_ < 2.0e-1 * fixedBStarResult.initialPositionRms_ );

    // Check that adding B* reduces position RMS by at least another factor of ten relative to the fixed-B* fit.
    BOOST_CHECK( estimatedBStarResult.positionRms_ < 1.0e-1 * fixedBStarResult.positionRms_ );

    // Check that the fixed-B* nonsingular mean-element vector changes by less than 5e-3 in absolute parameter norm.
    BOOST_CHECK_SMALL( ( fixedBStarEquinoctialElements - referenceEquinoctialElements ).norm( ), 5.0e-3 );

    // Check that estimating B* leaves the nonsingular mean-element vector within 1e-4 of the seed in absolute norm.
    BOOST_CHECK_SMALL( ( estimatedBStarEquinoctialElements - referenceEquinoctialElements ).norm( ), 1.0e-4 );
}

//! Verify that estimating B* captures drag in an independently numerically propagated Tudat LEO arc.
/*!
 * This is intentionally not an SGP4 self-fit. Tudat propagates a 400 km orbit for three days with 16x16 Earth gravity,
 * aerodynamic drag, cannonball radiation pressure and point-mass Sun/Moon gravity. The same state history is fitted
 * once with B* fixed at zero and once with B* estimated. Since B* is an empirical SGP4 drag parameter, the test checks
 * its sign and its effect on trajectory RMS, rather than claiming recovery of a physical ballistic coefficient.
 */
BOOST_AUTO_TEST_CASE( testBStarFitToNumericallyPropagatedDragArc )
{
    using namespace ephemerides;

    spice_interface::loadStandardSpiceKernels( );

    // Use a multi-day arc so the along-track signature of drag is distinguishable from a small mean-motion adjustment.
    const double initialEpoch = 3.0e6;
    const double finalEpoch = initialEpoch + 3.0 * physical_constants::JULIAN_DAY;
    const double aerodynamicReferenceArea = 5.0;
    const double spacecraftMass = 100.0;
    const double dragCoefficient = 2.2;
    const double radiationPressureCoefficient = 1.2;

    // Start the numerical propagation from an SGP4-consistent osculating state, with B* set to zero only in the seed TLE.
    const double initialSemiMajorAxis = 6378.135e3 + 400.0e3;
    const double initialMeanMotion = std::sqrt( 398600.8e9 / std::pow( initialSemiMajorAxis, 3.0 ) ) * 60.0;
    const std::shared_ptr< Tle > zeroDragInitialTle =
            std::make_shared< Tle >( initialEpoch, 0.0, 0.9, 0.4, 0.001, 0.7, 0.2, initialMeanMotion );
    const std::map< double, Eigen::Vector6d > dragStateHistory = propagateTleWithExtendedDynamics(
            zeroDragInitialTle, finalEpoch, aerodynamicReferenceArea, spacecraftMass, dragCoefficient, radiationPressureCoefficient );

    // Establish the baseline by fitting all six mean orbital parameters while B* is constrained to its zero seed value.
    TleFitSettings fixedBStarSettings;
    fixedBStarSettings.tleEpoch_ = initialEpoch;
    fixedBStarSettings.initialTle_ = zeroDragInitialTle;
    fixedBStarSettings.maximumNumberOfIterations_ = 50;
    fixedBStarSettings.convergenceTolerance_ = 1.0e-9;
    const TleFitResult fixedBStarResult = fitTleToCartesianStateHistory( dragStateHistory, "J2000", fixedBStarSettings );

    // Repeat from exactly the same seed and arc, adding only the scaled B* solve-for parameter.
    TleFitSettings estimatedBStarSettings = fixedBStarSettings;
    estimatedBStarSettings.estimateBStar_ = true;
    const TleFitResult estimatedBStarResult = fitTleToCartesianStateHistory( dragStateHistory, "J2000", estimatedBStarSettings );

    // Check that both baseline and seven-parameter fits returned usable full-precision TLE objects.
    BOOST_REQUIRE( fixedBStarResult.fittedTle_ != nullptr && estimatedBStarResult.fittedTle_ != nullptr );

    BOOST_TEST_MESSAGE( "Numerical drag-arc TLE fit: fixed-B* RMS = " << fixedBStarResult.positionRms_
                                                                      << " m, estimated-B* RMS = " << estimatedBStarResult.positionRms_
                                                                      << " m, fitted B* = " << estimatedBStarResult.fittedTle_->getBStar( )
                                                                      << ", converged = " << estimatedBStarResult.converged_
                                                                      << ", iterations = " << estimatedBStarResult.numberOfIterations_ );

    // Check that the fixed-B* fit preserves the supplied zero value exactly.
    BOOST_CHECK_SMALL( fixedBStarResult.fittedTle_->getBStar( ), 1.0e-15 );

    // Check that fitting a decaying drag arc produces the expected positive empirical B* value.
    BOOST_CHECK( estimatedBStarResult.fittedTle_->getBStar( ) > 0.0 );

    // Check that B* remains finite and below a deliberately broad limit, guarding against an unconstrained runaway solution.
    BOOST_CHECK( std::isfinite( estimatedBStarResult.fittedTle_->getBStar( ) ) && estimatedBStarResult.fittedTle_->getBStar( ) < 0.1 );

    // Check that the seven-parameter solver reached its parameter-step convergence criterion on the drag arc.
    BOOST_CHECK( estimatedBStarResult.converged_ );

    // Check that adding B* reduces position RMS by at least a factor of 100, strongly distinguishing drag recovery
    // from a marginal improvement caused by adding one more solve-for parameter.
    BOOST_CHECK( estimatedBStarResult.positionRms_ < 0.01 * fixedBStarResult.positionRms_ );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
