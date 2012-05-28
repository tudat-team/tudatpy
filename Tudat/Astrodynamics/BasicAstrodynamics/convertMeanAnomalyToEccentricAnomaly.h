/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110210    K. Kumar          Creation of code.
 *      110215    E. Iorfida        Minor changes made.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120421    K. Kumar          Removed base class; updated to set values through constructor.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 */

#ifndef TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
#define TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H

#include <cmath>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

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
class ConvertMeanAnomalyToEccentricAnomaly
{
public:

    //! Construct converter with eccentricity and mean anomaly.
    /*!
     * This constructor initializes the eccentricity and the mean anomaly for the conversion
     * routine.
     * \param eccentricity The eccentricity of the orbit [-].
     * \param eccentricAnomaly The mean anomaly to convert to eccentric anomaly [rad].
     */
    ConvertMeanAnomalyToEccentricAnomaly( const double eccentricity, const double meanAnomaly,
                                          boost::shared_ptr< NewtonRaphson > newtonRaphson )
        : eccentricity_( eccentricity ),
          meanAnomaly_( meanAnomaly ),
          newtonRaphson_( newtonRaphson )
    { }

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

    //! Eccentricity.
    /*!
     * Eccentricity.
     */
    double eccentricity_;

    //! Mean anomaly.
    /*!
     * Mean anomaly.
     */
    double meanAnomaly_;

    //! Shared pointer to Newton-Raphson.
    /*!
     * Shared pointer to Newton-Raphson method.
     */
    boost::shared_ptr< NewtonRaphson > newtonRaphson_;

    //! NewtonRaphsonAdaptor.
    /*!
     * NewtonRaphsonAdaptor class.
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
    {
        return eccentricAnomaly - eccentricity_ * std::sin( eccentricAnomaly ) - meanAnomaly_;
    }

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
    {
        return 1.0 - eccentricity_ * std::cos( eccentricAnomaly );
    }
};

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
