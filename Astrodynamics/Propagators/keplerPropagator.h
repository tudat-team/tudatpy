/*! \file keplerPropagator.h
 *    Header file that defines the kepler propagator class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 4
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
 *    Last modified     : 20 September, 2011
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
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110214    K. Kumar          Updated code to use orbital element conversion functions.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

#ifndef KEPLERPROPAGATOR_H
#define KEPLERPROPAGATOR_H

// Include statements.
#include <iostream>
#include <cmath>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Propagators/propagator.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/convertMeanAnomalyToEccentricAnomaly.h"
#include "Astrodynamics/States/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
#include "Astrodynamics/States/keplerianElements.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"

// Using declarations.
using orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly;
using orbital_element_conversions::ConvertMeanAnomalyToHyperbolicEccentricAnomaly;

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
    KeplerPropagator( ) : eccentricAnomaly_( -0.0 ), eccentricAnomalyChange_( -0.0 ),
        hyperbolicEccentricAnomaly_( -0.0 ), meanAnomaly_( -0.0 ), meanAnomalyChange_( -0.0 ),
        trueAnomaly_( -1.0 ), pointerToNewtonRaphson_( NULL ) { }

    //! Set central body.
    /*!
     * Sets the central body for given body to propagate.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToCentralBody Central body given as pointer to
     *          CelestialBody object.
     */
    void setCentralBody( Body* pointerToBody, CelestialBody* pointerToCentralBody )
    { bodiesToPropagate_[ pointerToBody ].pointerToCentralBody = pointerToCentralBody; }

    //! Set Newton-Raphson method.
    /*!
     * Sets the Newton-Raphson method used.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphson( NewtonRaphson* pointerToNewtonRaphson )
    { pointerToNewtonRaphson_ = pointerToNewtonRaphson; }

    //! Propagate.
    /*!
     * This function executes propagation.
     */
    void propagate( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param keplerPropagator Kepler propagator.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, KeplerPropagator& keplerPropagator );

protected:

private:

    //! Eccentric anomaly.
    /*!
     * Eccentric anomaly.
     */
    double eccentricAnomaly_;

    //! Change of eccentric anomaly.
    /*!
     * Change of eccentric anomaly.
     */
    double eccentricAnomalyChange_;

    //! Hyperbolic eccentric anomaly.
    /*!
     * Hyperbolic eccentric anomaly.
     */
    double hyperbolicEccentricAnomaly_;

    //! Mean anomaly.
    /*!
     * Mean anomaly.
     */
    double meanAnomaly_;

    //! Mean anomaly change.
    /*!
     * Change of mean anomaly between start time and end time.
     */
    double meanAnomalyChange_;

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

    //! Mean anomaly to eccentric anomaly conversion.
    /*!
     * Mean anomaly to eccentric anomaly conversion
     */
    ConvertMeanAnomalyToEccentricAnomaly convertMeanAnomalyToEccentricAnomaly_;

    //! Mean anomaly to hyperbolic eccentric anomaly conversion.
    /*!
     * Mean anomaly to hyperbolic eccentric anomaly conversion.
     */
    ConvertMeanAnomalyToHyperbolicEccentricAnomaly convertMeanAnomalyToHyperbolicEccentricAnomaly_;

    //! Cartesian elements object.
    /*!
     * Cartesian elements object.
     */
    CartesianElements cartesianElements_;

    //! Keplerian elements object.
    /*!
     * Keplerian elements object.
     */
    KeplerianElements keplerianElements_;
};

#endif // KEPLERPROPAGATOR_H

// End of file.
