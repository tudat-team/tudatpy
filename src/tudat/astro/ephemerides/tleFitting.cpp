/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 */

#include "tudat/astro/ephemerides/tleFitting.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "tudat/astro/basic_astro/equinoctialElementConversions.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace ephemerides
{
namespace
{

// CSPICE evsgp4 uses the WGS-72 gravity model. Using the same gravitational parameter for Cartesian initialization
// avoids introducing an avoidable mean-motion offset before differential correction starts.
constexpr double wgs72EarthGravitationalParameter = 398600.8e9;

// Leave a small numerical margin below the elliptic-orbit boundary accepted by SGP4.
constexpr double maximumTleEccentricity = 0.999;

//! Map a numerical TLE into the generic equinoctial parameters and optional scaled B* used by the fitter.
Eigen::VectorXd getFitParametersFromTle( const Tle& tle,
                                         const bool estimateBStar,
                                         const double bStarScale,
                                         const bool useRetrogradeElements )
{
    Eigen::Vector6d orbitalElements;
    orbitalElements << tle.getMeanMotion( ), tle.getEccentricity( ), tle.getInclination( ), tle.getRightAscension( ),
            tle.getArgOfPerigee( ), tle.getMeanAnomaly( );

    Eigen::VectorXd parameters( estimateBStar ? 7 : 6 );
    parameters.segment( 0, 6 ) = orbital_element_conversions::convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements(
            orbitalElements, useRetrogradeElements );
    if( estimateBStar )
    {
        // Scale B* independently because its useful numerical magnitude differs strongly from the orbital parameters.
        parameters( 6 ) = tle.getBStar( ) / bStarScale;
    }
    return parameters;
}

//! Map generic equinoctial fit parameters and optional scaled B* into a full-precision numerical TLE.
/*!
 * Invalid trial vectors return nullptr so that the residual function can reject them without asking CSPICE to
 * propagate a nonphysical element set.
 */
std::shared_ptr< Tle > createTleFromFitParameters( const Eigen::VectorXd& parameters,
                                                   const double epoch,
                                                   const double fixedBStar,
                                                   const bool estimateBStar,
                                                   const double bStarScale,
                                                   const bool useRetrogradeElements )
{
    if( !parameters.allFinite( ) || parameters.rows( ) != ( estimateBStar ? 7 : 6 ) || parameters( 0 ) < -20.0 || parameters( 0 ) > 20.0 )
    {
        return nullptr;
    }

    const double inclinationParameterNorm = std::hypot( parameters( 3 ), parameters( 4 ) );
    if( inclinationParameterNorm > 1.0e6 )
    {
        return nullptr;
    }

    const Eigen::Vector6d orbitalElements = orbital_element_conversions::convertEquinoctialToOrbitalElementsWithMeanMotionAndMeanAnomaly(
            parameters.segment( 0, 6 ), useRetrogradeElements );
    const double bStar = estimateBStar ? parameters( 6 ) * bStarScale : fixedBStar;
    if( orbitalElements( 1 ) >= maximumTleEccentricity || !std::isfinite( bStar ) )
    {
        return nullptr;
    }
    return std::make_shared< Tle >( epoch,
                                    bStar,
                                    orbitalElements( 2 ),
                                    orbitalElements( 3 ),
                                    orbitalElements( 1 ),
                                    orbitalElements( 4 ),
                                    orbitalElements( 5 ),
                                    orbitalElements( 0 ) );
}

//! Construct an osculating-element first approximation when the caller did not supply an initial TLE.
/*!
 * The requested epoch is changed to the nearest available state epoch. This ensures that the Cartesian state used for
 * initialization and the epoch attached to its derived elements are identical.
 */
std::shared_ptr< Tle > createInitialTle( const std::map< double, Eigen::Vector6d >& stateHistory,
                                         const std::string& frameOrientation,
                                         double& tleEpoch,
                                         const TleFitSettings& settings )
{
    if( settings.initialTle_ )
    {
        // A supplied TLE is already tied to its own epoch; never relabel its mean elements at another epoch.
        tleEpoch = settings.initialTle_->getEpoch( );
        return settings.initialTle_;
    }

    // Snap the requested fitting epoch to a real observation before deriving elements from that observation.
    auto nearestState = std::min_element( stateHistory.begin( ), stateHistory.end( ), [ tleEpoch ]( const auto& left, const auto& right ) {
        return std::fabs( left.first - tleEpoch ) < std::fabs( right.first - tleEpoch );
    } );
    tleEpoch = nearestState->first;

    // Reuse the frame conversion provided by TleEphemeris so initialization and forward propagation cannot diverge.
    const Eigen::Matrix3d inputFrameToTeme = getRotationMatrixFromFrameToTeme( tleEpoch, frameOrientation );
    Eigen::Vector6d temeState;
    temeState.segment( 0, 3 ) = inputFrameToTeme * nearestState->second.segment( 0, 3 );
    temeState.segment( 3, 3 ) = inputFrameToTeme * nearestState->second.segment( 3, 3 );

    // Pass through Tudat's modified-equinoctial conversion to handle circular/equatorial Cartesian states robustly.
    // The result is osculating rather than SGP4 mean elements, but the batch fit corrects that initial discrepancy.
    const Eigen::Vector6d osculatingKeplerianElements =
            orbital_element_conversions::convertCartesianToKeplerianElements( temeState, wgs72EarthGravitationalParameter );
    const bool isRetrograde =
            osculatingKeplerianElements( orbital_element_conversions::inclinationIndex ) > mathematical_constants::PI / 2.0;
    const Eigen::Vector6d modifiedEquinoctialElements =
            orbital_element_conversions::convertKeplerianToModifiedEquinoctialElements( osculatingKeplerianElements, isRetrograde );
    Eigen::Vector6d keplerianElements =
            orbital_element_conversions::convertModifiedEquinoctialToKeplerianElements( modifiedEquinoctialElements, isRetrograde );
    const double eccentricity = keplerianElements( orbital_element_conversions::eccentricityIndex );
    if( !keplerianElements.allFinite( ) || eccentricity < 0.0 || eccentricity >= maximumTleEccentricity ||
        keplerianElements( orbital_element_conversions::semiMajorAxisIndex ) <= 0.0 )
    {
        throw std::runtime_error( "Cannot form a bound, elliptic initial TLE from the supplied Cartesian state." );
    }
    const double meanAnomaly = orbital_element_conversions::convertTrueAnomalyToMeanAnomaly(
            eccentricity, keplerianElements( orbital_element_conversions::trueAnomalyIndex ) );

    // Tle stores mean motion in radians per minute, whereas Tudat's Cartesian conversion uses SI units.
    const double meanMotion = std::sqrt( wgs72EarthGravitationalParameter /
                                         std::pow( keplerianElements( orbital_element_conversions::semiMajorAxisIndex ), 3.0 ) ) *
            60.0;
    return std::make_shared< Tle >( tleEpoch,
                                    settings.initialBStar_,
                                    keplerianElements( orbital_element_conversions::inclinationIndex ),
                                    keplerianElements( orbital_element_conversions::longitudeOfAscendingNodeIndex ),
                                    eccentricity,
                                    keplerianElements( orbital_element_conversions::argumentOfPeriapsisIndex ),
                                    meanAnomaly,
                                    meanMotion );
}

//! Propagate one trial numerical TLE and stack input-minus-model position residuals for all epochs.
Eigen::VectorXd getPositionResidualVector( const Eigen::VectorXd& parameters,
                                           const std::map< double, Eigen::Vector6d >& stateHistory,
                                           const std::string& frameOrientation,
                                           const double tleEpoch,
                                           const double fixedBStar,
                                           const bool estimateBStar,
                                           const double bStarScale,
                                           const bool useRetrogradeElements )
{
    const std::shared_ptr< Tle > tle =
            createTleFromFitParameters( parameters, tleEpoch, fixedBStar, estimateBStar, bStarScale, useRetrogradeElements );
    if( !tle )
    {
        // A large finite residual lets the damping logic reject invalid trials without contaminating the QR solve
        // with NaN or infinity values.
        return Eigen::VectorXd::Constant( 3 * stateHistory.size( ), 1.0e16 );
    }

    // TleEphemeris evaluates the existing CSPICE SGP4 implementation directly from the full-precision Tle fields.
    TleEphemeris ephemeris( "Earth", frameOrientation, tle, false );
    Eigen::VectorXd residuals( 3 * stateHistory.size( ) );
    unsigned int index = 0;
    for( const auto& stateAtEpoch : stateHistory )
    {
        // Propagation exceptions are intentionally not swallowed: a valid numerical TLE that CSPICE cannot evaluate
        // indicates a model/configuration failure that the caller must see, rather than an ordinary rejected trial.
        residuals.segment( 3 * index, 3 ) =
                stateAtEpoch.second.segment( 0, 3 ) - ephemeris.getCartesianState( stateAtEpoch.first ).segment( 0, 3 );
        ++index;
    }
    return residuals;
}

//! Assemble central-difference perturbations in the scaled solve-for parameter space.
Eigen::VectorXd getFiniteDifferenceSteps( const Eigen::VectorXd& parameters, const TleFitSettings& settings )
{
    Eigen::VectorXd steps = Eigen::VectorXd::Constant( parameters.rows( ), settings.equinoctialElementStep_ );
    steps( 0 ) = settings.logarithmicMeanMotionStep_;
    steps( 5 ) = settings.meanLongitudeStep_;
    if( parameters.rows( ) == 7 )
    {
        steps( 6 ) = settings.bStarStep_;
    }
    if( !( steps.array( ) > 0.0 ).all( ) || !steps.allFinite( ) )
    {
        throw std::invalid_argument( "All TLE-fit finite-difference steps must be finite and positive." );
    }
    return steps;
}

}  // namespace

//! Fit an SGP4 numerical TLE to the positions in a Cartesian state history.
TleFitResult fitTleToCartesianStateHistory( const std::map< double, Eigen::Vector6d >& cartesianStateHistory,
                                            const std::string& frameOrientation,
                                            const TleFitSettings& settings )
{
    // Validate all public inputs before constructing a TLE or entering the iterative solver.
    if( cartesianStateHistory.size( ) < 2 )
    {
        throw std::invalid_argument( "At least two Cartesian states are required to fit a TLE." );
    }
    if( frameOrientation != "J2000" && frameOrientation != "ECLIPJ2000" )
    {
        throw std::invalid_argument( "TLE fitting supports only J2000 and ECLIPJ2000 input states." );
    }
    for( const auto& stateAtEpoch : cartesianStateHistory )
    {
        if( !std::isfinite( stateAtEpoch.first ) || !stateAtEpoch.second.allFinite( ) )
        {
            throw std::invalid_argument( "TLE fitting requires finite epochs and Cartesian states." );
        }
    }
    if( settings.maximumNumberOfIterations_ == 0 || settings.convergenceTolerance_ <= 0.0 || settings.initialDamping_ <= 0.0 ||
        settings.bStarScale_ <= 0.0 || !std::isfinite( settings.bStarScale_ ) )
    {
        throw std::invalid_argument( "TLE fit iteration settings must be positive." );
    }

    // Default to a central observation so forward and backward propagation contribute comparably to the fit.
    const auto middleState = std::next( cartesianStateHistory.begin( ), cartesianStateHistory.size( ) / 2 );
    double tleEpoch = std::isfinite( settings.tleEpoch_ ) ? settings.tleEpoch_ : middleState->first;
    const std::shared_ptr< Tle > initialTle = createInitialTle( cartesianStateHistory, frameOrientation, tleEpoch, settings );

    // Select the nonsingular inclination branch once. Switching branches during finite differences would make the
    // residual function discontinuous.
    const bool useRetrogradeElements = initialTle->getInclination( ) > mathematical_constants::PI / 2.0;
    const Eigen::VectorXd initialParameters =
            getFitParametersFromTle( *initialTle, settings.estimateBStar_, settings.bStarScale_, useRetrogradeElements );
    const double fixedBStar = initialTle->getBStar( );

    // Capture one forward-model function so the nominal residual, finite differences and trial acceptance use exactly
    // the same full-precision SGP4 path.
    const auto residualFunction = [ & ]( const Eigen::VectorXd& parameters ) {
        return getPositionResidualVector( parameters,
                                          cartesianStateHistory,
                                          frameOrientation,
                                          tleEpoch,
                                          fixedBStar,
                                          settings.estimateBStar_,
                                          settings.bStarScale_,
                                          useRetrogradeElements );
    };
    const Eigen::VectorXd initialResiduals = residualFunction( initialParameters );

    // Form a central-difference Jacobian of the stacked position residuals. SGP4 is reinitialized for every perturbed
    // numerical TLE, with no text serialization or field quantization.
    const auto jacobianFunction = [ & ]( const Eigen::VectorXd& parameters ) {
        const Eigen::VectorXd residuals = residualFunction( parameters );
        const Eigen::VectorXd steps = getFiniteDifferenceSteps( parameters, settings );
        Eigen::MatrixXd jacobian( residuals.rows( ), parameters.rows( ) );
        for( int column = 0; column < parameters.rows( ); ++column )
        {
            Eigen::VectorXd positiveParameters = parameters;
            Eigen::VectorXd negativeParameters = parameters;
            positiveParameters( column ) += steps( column );
            negativeParameters( column ) -= steps( column );
            jacobian.col( column ) =
                    ( residualFunction( positiveParameters ) - residualFunction( negativeParameters ) ) / ( 2.0 * steps( column ) );
        }
        return jacobian;
    };

    // Reuse Tudat's Levenberg-Marquardt implementation rather than maintaining a TLE-specific nonlinear solver. Its
    // damping update follows Madsen, Nielsen and Tingleff, Methods for Non-Linear Least Squares Problems, 2nd ed. (2004).
    const auto residualAndJacobianFunction = [ & ]( const Eigen::VectorXd& parameters ) {
        return std::make_pair( residualFunction( parameters ), jacobianFunction( parameters ) );
    };
    unsigned int numberOfIterations = 0;
    bool converged = false;
    const Eigen::VectorXd fittedParameters = linear_algebra::nonLinearLeastSquaresFit( residualAndJacobianFunction,
                                                                                       initialParameters,
                                                                                       Eigen::VectorXd::Zero( initialResiduals.rows( ) ),
                                                                                       settings.initialDamping_,
                                                                                       settings.convergenceTolerance_,
                                                                                       settings.maximumNumberOfIterations_,
                                                                                       numberOfIterations,
                                                                                       converged );
    const Eigen::VectorXd fittedResiduals = residualFunction( fittedParameters );
    const std::shared_ptr< Tle > fittedTle = createTleFromFitParameters(
            fittedParameters, tleEpoch, fixedBStar, settings.estimateBStar_, settings.bStarScale_, useRetrogradeElements );
    if( !fittedTle )
    {
        throw std::runtime_error( "TLE fit ended at invalid TLE parameters." );
    }

    // Return unweighted physical diagnostics; positionRms_ is the RMS of per-epoch three-dimensional error norms.
    TleFitResult result;
    result.fittedTle_ = fittedTle;
    result.initialPositionRms_ = std::sqrt( initialResiduals.squaredNorm( ) / static_cast< double >( cartesianStateHistory.size( ) ) );
    result.positionRms_ = std::sqrt( fittedResiduals.squaredNorm( ) / static_cast< double >( cartesianStateHistory.size( ) ) );
    result.numberOfIterations_ = numberOfIterations;
    result.converged_ = converged;
    for( unsigned int i = 0; i < cartesianStateHistory.size( ); ++i )
    {
        result.positionResiduals_.push_back( fittedResiduals.segment( 3 * i, 3 ) );
    }
    return result;
}

}  // namespace ephemerides
}  // namespace tudat
