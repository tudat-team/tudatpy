/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_KEPLEREPHEMERIS_H
#define TUDAT_KEPLEREPHEMERIS_H

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"

namespace tudat
{

namespace ephemerides
{

//! Ephemeris derived class that calculates the Cartesian position as a function of time assuming
//! a Kepler orbit.
/*!
 *  Ephemeris derived class that calculates the Cartesian position as a function of time assuming
 *  a Kepler orbit. In this class, both the Kepler propagation and conversion to Cartesian
 *  coordinates is done internally using precomputed values of the initial Kepler elements,
 *  thereby reducing computation time compared to succesively using the propagateKeplerOrbit or
 *  convertKeplerianToCartesianElements functions.
 */
class KeplerEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianState;

    //! Class constructor.
    /*!
     *  Class constructor, sets the characteristics of the Kepler orbit and the root finder
     *  used for converting mean to eccentric anomalies.
     *  \param initialStateInKeplerianElements Kepler elements at time epochOfInitialState.
     *  \param epochOfInitialState Time at which initialStateInKeplerianElements represents
     *  the Keplerian state.
     *  \param centralBodyGravitationalParameter Gravitational parameter of the central body
     *  that is used in the computations.
     *  \param referenceFrameOrigin Origin of reference frame (string identifier) (default SSB).
     *  \param referenceFrameOrientation Orientation of reference frame (string identifier)
     *  (default ECLIPJ000
     *  \param rootFinderAbsoluteTolerance Convergence tolerance for root finder used to
     *  convert mean to eccentric anomaly on each call to getCartesianState
     *  (default 200*epsilon).
     *  \param rootFinderMaximumNumberOfIterations Maximum iteration for root finder used to
     *  convert mean to eccentric anomaly on each call to getCartesianState
     *  (default 1000).
     */
    KeplerEphemeris( const Eigen::Vector6d& initialStateInKeplerianElements,
                     const double epochOfInitialState,
                     const double centralBodyGravitationalParameter,
                     const std::string& referenceFrameOrigin = "SSB",
                     const std::string& referenceFrameOrientation = "ECLIPJ2000",
                     const double rootFinderAbsoluteTolerance =
                         200.0 * std::numeric_limits< double >::epsilon( ),
                     const double rootFinderMaximumNumberOfIterations = 1000.0 );

    //! Function to get state from ephemeris.
    /*!
     *  Returns state from ephemeris at given time, assuming a purely Keplerian orbit
     *  \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Keplerian orbit Cartesian state at given time.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch );

private:

    //! Kepler elements at time epochOfInitialState.
    Eigen::Vector6d initialStateInKeplerianElements_;

    //! Semi-latus rectum of orbit.
    double semiLatusRectum_;

    //! Mean anomaly at epochOfInitialState.
    double initialMeanAnomaly_;

    //! Semi-major axis of orbit.
    double semiMajorAxis_;

    //! Eccentricity of orbit.
    double eccentricity_;

    //! Rotation from orbital plane to frame in which the orbit is defined.
    Eigen::Quaterniond rotationFromOrbitalPlane_;

    //! Root finder used to convert mean to eccentric anomalies.
    std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder_;

    //! Initial epoch from which propagation of Kepler orbit is performed.
    double epochOfInitialState_;

    //! Gravitational parameter of central body about which the Kepler orbit is defined.
    double centralBodyGravitationalParameter_;

    //! Boolean denoting whether orbit is hyperbolic or elliptical (parabola not supported).
    bool isOrbitHyperbolic_;
};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_KEPLEREPHEMERIS_H
