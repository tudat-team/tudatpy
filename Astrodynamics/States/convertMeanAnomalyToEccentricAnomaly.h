/*! \file convertMeanAnomalyToEccentricAnomaly.h
 *    This header file contains a class to convert mean anomly to eccentric anomaly for elliptical
 *    orbits.
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
 *    Checker           : S. Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Date created      : 10 February, 2011
 *    Last modified     : 26 January, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 *    Notes
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
 *      110210    K. Kumar          First creation of code.
 *      110215    E. Iorfida        Minor changes made.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef TUDAT_CONVERTMEANANOMALYTOECCENTRICANOMALY_H
#define TUDAT_CONVERTMEANANOMALYTOECCENTRICANOMALY_H

// Include statements.
#include "Astrodynamics/States/convertMeanAnomalyBase.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Orbital element conversions namespace.
/*!
 * The orbital element conversions namespace that contains free functions and classes to
 * convert between orbital elements e.g., mean anomaly to eccentric anomaly, eccentric anomaly to
 * true anomaly etc.
 */
namespace orbital_element_conversions
{

//! Converter for Mean to Eccentric anomaly.
/*!
 * Conversion utility to convert the mean anomaly to eccentric anomaly for elliptical orbits.
 */
class ConvertMeanAnomalyToEccentricAnomaly : public ConvertMeanAnomalyBase
{
public:

    //! Construct converter with eccentricity and Eccentric Anomaly
    /*!
     * This constructor initializes the eccentricity and the Eccentric Anomaly for the conversion 
     * routine.
     * \param eccentricity The eccentricity of the orbit (default=0) [-].
     * \param eccentricAnomaly The eccentric anomaly to convert to Mean anomaly (default=0) [rad].
     */
    ConvertMeanAnomalyToEccentricAnomaly( double eccentricity = 0.0,
                                          double eccentricAnomaly = 0.0 )
        :  eccentricAnomaly_( eccentricAnomaly ) { setEccentricity( eccentricity ); }

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Converts mean anomaly to eccentric anomaly for elliptical orbits. Currently, the conversion
     * does not work for near-parabolic orbits ( only for eccentricity < 0.98 ).
     * If the conversion fails or the eccentricity falls outside the valid range, then double::NaN
     * is returned.
     * \return Eccentric anomaly [rad].
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
    NewtonRaphsonAdaptor< ConvertMeanAnomalyToEccentricAnomaly > newtonRaphsonAdaptor_;

    //! Compute Kepler's function for elliptical orbits.
    /*!
     * Computes Kepler's function, given as:
     * \f[
     *      f( E ) = E - e * sin( E ) - M
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, \f$ e \f$ is the
     * eccentricity, \f$ M \f$ is the mean anomaly. Currently, the case for near-parabolic orbits
     * is not included ( 0.8 < eccentricity < 1.2 ).
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of Kepler's function for elliptical orbits.
     */
    double computeKeplersFunctionForEllipticalOrbits_( double& eccentricAnomaly )
    {     return eccentricAnomaly - eccentricity_ * sin( eccentricAnomaly ) - meanAnomaly_; }

    //! Compute first-derivative of Kepler's function for elliptical orbits.
    /*!
     * Computes the first-derivative of Kepler's function, given as:
     * \f[
     *      \frac{ df( E ) } { dE } = 1 - e * cos( E )
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, and \f$ e \f$ is the
     * eccentricity. Currently, the case for near-parabolic orbits is not included
     * ( 0.8 < eccentricity < 1.2 ).
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of first-derivative of Kepler's function for elliptical orbits.
     */
    double computeFirstDerivativeKeplersFunctionForEllipticalOrbits_( double& eccentricAnomaly )
    { return 1.0 - eccentricity_ * cos( eccentricAnomaly ); }
};

} // Namespace orbital_element_conversions.

} // Namespace tudat.

#endif // TUDAT_CONVERTMEANANOMALYTOECCENTRICANOMALY_H

// End of file.
