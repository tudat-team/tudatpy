/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms, with or without modification, are permitted exclusively under the
 *    terms of the Modified BSD license.
 */

#ifndef TUDAT_TLEFITTING_H
#define TUDAT_TLEFITTING_H

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/tleEphemeris.h"

namespace tudat
{
namespace ephemerides
{

//! Settings for fitting a numerical TLE to an Earth-centred Cartesian state history.
/*!
 * The fit minimizes Cartesian position residuals only. States must use SI units and be expressed in either J2000 or
 * ECLIPJ2000. The TLE epoch and all epochs are seconds since J2000 on the same time scale used by TleEphemeris.
 *
 * The optimizer works on a full-precision numerical Tle; it never serializes or reparses TLE text during the fit.
 */
struct TleFitSettings {
    //! Requested fitted-TLE epoch. The Cartesian initializer snaps it to the nearest state epoch; NaN selects the middle sample.
    double tleEpoch_ = std::numeric_limits< double >::quiet_NaN( );

    //! Optional full-precision initial TLE. A Cartesian/MEE-based initializer is used when this is null.
    std::shared_ptr< Tle > initialTle_ = nullptr;

    //! Estimate B* in addition to the six mean orbital parameters. It is false by default.
    bool estimateBStar_ = false;

    //! B* used by the Cartesian initializer. Ignored when initialTle_ is supplied.
    double initialBStar_ = 0.0;

    //! Maximum number of Levenberg-Marquardt iterations; matches Tudat's nonLinearLeastSquaresFit default.
    unsigned int maximumNumberOfIterations_ = 25;

    //! Parameter-step tolerance; matches the 1.0e-8 default in Tudat's nonLinearLeastSquaresFit utility.
    double convergenceTolerance_ = 1.0e-8;

    //! Initial damping; matches the 1.0e-6 initialScaling default in Tudat's nonLinearLeastSquaresFit utility.
    double initialDamping_ = 1.0e-6;

    //! Tudat-specific B* scale, selected for numerical conditioning and covered by the independent drag-arc regression test.
    /*!
     * This value is an implementation choice, not a physical constant or a value prescribed by Vallado and Crawford.
     * With the B* finite-difference step below, it gives a physical perturbation of 1.0e-7.
     */
    double bStarScale_ = 1.0e-4;

    //! Tudat-specific central-difference step for log mean motion, validated by the synthetic round-trip tests.
    /*!
     * Vallado and Crawford, AIAA 2008-6770, Section II.C motivates element-wise finite differences and reports a
     * relative perturbation of 1.0e-3 with a 1.0e-7 absolute fallback. The fixed steps here act on Tudat's scaled
     * equinoctial solve-for parameters, so their numerical values are local implementation defaults rather than values
     * copied directly from that paper.
     */
    double logarithmicMeanMotionStep_ = 1.0e-6;

    //! Tudat-specific central-difference step for the four nonsingular eccentricity/inclination parameters.
    double equinoctialElementStep_ = 1.0e-6;

    //! Tudat-specific central-difference step for mean longitude.
    double meanLongitudeStep_ = 1.0e-6;

    //! Scaled B* step; with bStarScale_=1.0e-4 this is 1.0e-7 in B*, Vallado and Crawford's absolute fallback size.
    double bStarStep_ = 1.0e-3;
};

//! Result of a Cartesian-state-history-to-TLE fit.
struct TleFitResult {
    //! Full-precision fitted TLE. Its raw line strings are empty because no TLE text serialization is performed.
    std::shared_ptr< Tle > fittedTle_;

    //! Position residuals (input state position minus fitted SGP4 position), in metres.
    std::vector< Eigen::Vector3d > positionResiduals_;

    //! RMS norm of position residuals, in metres.
    double positionRms_ = std::numeric_limits< double >::quiet_NaN( );

    //! RMS position residual of the supplied or generated initial TLE, in metres.
    double initialPositionRms_ = std::numeric_limits< double >::quiet_NaN( );

    //! Number of Levenberg-Marquardt iterations performed.
    unsigned int numberOfIterations_ = 0;

    //! Whether the final computed update met the configured convergence tolerance.
    bool converged_ = false;
};

//! Fit a full-precision TLE to an Earth-centred Cartesian state history.
/*!
 * \param cartesianStateHistory Map of epoch to Cartesian state, in metres and metres per second. Only the position is fitted.
 * \param frameOrientation Orientation of all input states: \c "J2000" or \c "ECLIPJ2000".
 * \param settings Fitter settings.
 * \return Fitted full-precision numerical TLE and residual diagnostics.
 *
 * This is a first, near-Earth SGP4 implementation. It uses the existing CSPICE evsgp4 forward model through
 * TleEphemeris and modified-equinoctial coordinates to avoid the circular/equatorial singularities of classical
 * elements.
 */
TleFitResult fitTleToCartesianStateHistory( const std::map< double, Eigen::Vector6d >& cartesianStateHistory,
                                            const std::string& frameOrientation,
                                            const TleFitSettings& settings = TleFitSettings( ) );

}  // namespace ephemerides
}  // namespace tudat

#endif  // TUDAT_TLEFITTING_H
