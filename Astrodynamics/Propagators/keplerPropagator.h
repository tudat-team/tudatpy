/*! \file keplerPropagator.h
 *    Header file that defines the kepler propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 3 February, 2011
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The code at present does not work for near-parabolic orbits
 *      ( 0.8 < eccentricity < 1.2 ). In future, this neeeds to be included
 *      and perhaps a universal method to solve Kepler's equation needs to be
 *      employed. Presently, the code will output an error if the eccentricity
 *      of the orbit to be propagated lies within this range.
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110203    K. Kumar          File created.
 *      110207    E. Iorfida        Minor changes.
 */

#ifndef KEPLERPROPAGATOR_H
#define KEPLERPROPAGATOR_H

// Include statements.
#include <cmath>
#include "propagator.h"
#include "basicMathematicsFunctions.h"
#include "body.h"
#include "cartesianElements.h"
#include "keplerianElements.h"
#include "newtonRaphson.h"
#include "newtonRaphsonAdaptor.h"
#include "orbitalElementConversions.h"

//! Kepler propagator class.
/*!
 * Definition of Kepler propagator class that propagates Kepler orbits
 * analytically.
 */
class KeplerPropagator : public Propagator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    KeplerPropagator( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~KeplerPropagator( );

    //! Set central body.
    /*!
     * Sets the central body for given body to propagate.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToCentralBody Central body given as pointer to
     *          CelestialBody object.
     */
    void setCentralBody( Body* pointerToBody,
                         CelestialBody* pointerToCentralBody );

    //! Set Newton-Raphson method.
    /*!
     * Sets the Newton-Raphson method used.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphson( NewtonRaphson* pointerToNewtonRaphson );

    //! Propagate.
    /*!
     * This function executes propagation.
     */
    void propagate( );

protected:

private:

    //! Number of propagation steps.
    /*!
     * Number of propagation steps.
     */
    unsigned int numberOfPropagationSteps_;

    //! Eccentric anomaly.
    /*!
     * Eccentric anomaly.
     */
    double eccentricAnomaly_;

    //! Hyperbolic eccentric anomaly.
    /*!
     * Hyperbolic eccentric anomaly.
     */
    double hyperbolicEccentricAnomaly_;

    //! Hyperbolic mean anomaly.
    /*!
     * Hyperbolic mean anomaly.
     */
    double hyperbolicMeanAnomaly_;

    //! Mean anomaly.
    /*!
     * Mean anomaly.
     */
    double meanAnomaly_;

    //! True anomaly.
    /*!
     * True anomaly.
     */
    double trueAnomaly_;

    //! Pointer to Newton-Raphson.
    /*!
     * Pointer to Newton-Raphson method.
     */
    NewtonRaphson* pointerToNewtonRaphson_;

    //! Pointer to adaptor NewtonRaphsonAdaptor.
    /*!
     * Pointer to adaptor NewtonRaphsonAdaptor class.
     */
    NewtonRaphsonAdaptor< KeplerPropagator >
            newtonRaphsonAdaptorForKeplerPropagator_;

    //! Pointer to Keplerian elements.
    /*!
     * Pointer to Keplerian elements.
     */
    KeplerianElements* pointerToKeplerianElements_;

    //! Keplerian elements.
    /*!
     * Keplerian elements.
     */
    KeplerianElements keplerianElements_;

    //! Compute Kepler's equation for elliptical orbits.
    /*!
     * Computes Kepler's equation for elliptical orbits. Kepler's equation for
     * elliptical orbits is given as:
     * \f[
     *      f_{ E } = E - e * sin( E ) - M
     * \f]
     * where \f$ E \f$ is the eccentric anomaly, \f$ e \f$ is the eccentricity
     * and \f$ M \f$ is the mean anomaly.
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of Kepler's equation for elliptical orbits.
     */
    double computeKeplerEquationForEllipticalOrbits_( double&
                                                      eccentricAnomaly );

    //! Compute first-derivative of Kepler's equation for elliptical orbits.
    /*!
     * Computes the first-derivative of Kepler's equation with respect to
     * \f$ E \f$, the eccentric anomaly, for elliptical orbits. The
     * first-derivative of Kepler's for elliptical orbits equation is given as:
     * \f[
     *      \frac{ df_{ E } } { dE } = 1 - e * cos( E )
     * \f]
     * where \f$ e \f$ is the eccentricity.
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of first-derivative of Kepler's equation for elliptical
     *          orbits.
     */
    double computeFirstDerivativeKeplerEquationForEllipticalOrbits_(
            double& eccentricAnomaly );

    //! Compute Kepler's equation for hyperbolic orbits.
    /*!
     * Computes Kepler's equation for hyperbolic orbits. Kepler's equation for
     * hyperbolic orbits is given as:
     * \f[
     *      g_{ H } = e * sinh( H ) - H - Mh
     * \f]
     * where \f$ H \f$ is the hyperbolic eccentric anomaly, \f$ e \f$ is the
     * eccentricity and \f$ Mh \f$ is the hyperbolic mean anomaly.
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
     * \return Value of Kepler's equation for hyperbolic orbits.
     */
    double computeKeplerEquationForHyperbolicOrbits_(
            double& hyperbolicEccentricAnomaly );

    //! Compute first-derivative of Kepler's equation for hyperbolic orbits.
    /*!
     * Computes the first-derivative of Kepler's equation with respect to
     * \f$ E \f$, the eccentric anomaly, for hyperbolic orbits. The
     * first-derivative of Kepler's for hyperbolic orbits equation is given as:
     * \f[
     *      \frac{ dg_{ H } } { dH } = e * cos( H ) - 1
     * \f]
     * where \f$ e \f$ is the eccentricity.
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
     * \return Value of first-derivative of Kepler's equation for hyperbolic
     *          orbits.
     */
    double computeFirstDerivativeKeplerEquationForHyperbolicOrbits_(
            double& hyperbolicEccentricAnomaly );
};

#endif // KEPLERPROPAGATOR_H

// End of file.
