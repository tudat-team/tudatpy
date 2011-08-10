/*! \file convertMeanAnomalyToEccentricAnomaly.h
 *    This header file contains a class to convert mean anomly to eccentric
 *    anomaly for elliptical orbits.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 3
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
 *    Date created      : 10 February, 2011
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series,
 *          VA, 2002.
 *
 *    Notes
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
 *      110210    K. Kumar          First creation of code.
 *      110215    E. Iorfida        Minor changes made.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef CONVERTMEANANOMALYTOECCENTRICANOMALY_H
#define CONVERTMEANANOMALYTOECCENTRICANOMALY_H

// Include statements.
#include "convertMeanAnomalyBase.h"

//! Orbital element conversions namespace.
/*!
 *  Orbital element conversions namespace.
 */
namespace orbital_element_conversions
{

//! Definition of mean anomaly to eccentric anomaly converter class.
/*!
 * Definition of mean anomaly to eccentric anomaly converter class.
 */
class ConvertMeanAnomalyToEccentricAnomaly : public ConvertMeanAnomalyBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConvertMeanAnomalyToEccentricAnomaly( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConvertMeanAnomalyToEccentricAnomaly( );

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Converts mean anomaly to eccentric anomaly for elliptical orbits.
     * Currently, the conversion does not work for near-parabolic
     * orbits ( 0.8 < eccentricity < 1.2 ).
     * \return Eccentric anomaly.
     */
    double convert( );

protected:

private:

    //! Eccentric anomaly.
    /*!
     * Eccentric anomaly.
     */
    double eccentricAnomaly_;

    //! Pointer to adaptor NewtonRaphsonAdaptor.
    /*!
     * Pointer to adaptor NewtonRaphsonAdaptor class.
     */
    NewtonRaphsonAdaptor< ConvertMeanAnomalyToEccentricAnomaly >
            newtonRaphsonAdaptor_;

    //! Compute Kepler's function for elliptical orbits.
    /*!
     * Computes Kepler's function, given as:
     * \f[
     *      f( E ) = E - e * sin( E ) - M
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly,
     * \f$ e \f$ is the eccentricity, \f$ M \f$ is the mean anomaly. Currently,
     * the case for near-parabolic orbits is not included
     * ( 0.8 < eccentricity < 1.2 ).
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of Kepler's function for elliptical orbits.
     */
    double computeKeplersFunctionForEllipticalOrbits_( double&
                                                       eccentricAnomaly );

    //! Compute first-derivative of Kepler's function for elliptical orbits.
    /*!
     * Computes the first-derivative of Kepler's function, given as:
     * \f[
     *      \frac{ df( E ) } { dE } = 1 - e * cos( E )
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, and
     * \f$ e \f$ is the eccentricity. Currently, the case for near-parabolic
     * orbits is not included ( 0.8 < eccentricity < 1.2 ).
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of first-derivative of Kepler's function for elliptical
     *          orbits.
     */
    double computeFirstDerivativeKeplersFunctionForEllipticalOrbits_(
            double& eccentricAnomaly );
};

}

#endif // CONVERTMEANANOMALYTOECCENTRICANOMALY_H

// End of file.
