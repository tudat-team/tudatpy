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
 *      110214    K. Kumar          Creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120421    K. Kumar          Removed base class; updated to set values through constructor.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 */

#ifndef TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H
#define TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H

#include <cmath>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace orbital_element_conversions
{

//! Definition of mean anomaly to hyperbolic eccentric anomaly converter class.
/*!
 * Definition of mean anomaly to hyperbolic eccentric anomaly converter class.
 */
class ConvertMeanAnomalyToHyperbolicEccentricAnomaly
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConvertMeanAnomalyToHyperbolicEccentricAnomaly(
            const double eccentricity, const double hyperbolicMeanAnomaly,
            boost::shared_ptr< NewtonRaphson > newtonRaphson )
        : eccentricity_( eccentricity ),
          hyperbolicMeanAnomaly_( hyperbolicMeanAnomaly ),
          newtonRaphson_( newtonRaphson )
    { }

    //! Convert mean anomaly to hyperbolic eccentric anomaly.
    /*!
     * Converts mean anomaly to hyperbolic eccentric anomaly for hyperbolic orbits. Currently, the
     * conversion does not work for near-parabolic orbits ( 0.8 < eccentricity < 1.2 ).
     * \return Hyperbolic eccentric anomaly.
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
    double hyperbolicMeanAnomaly_;

    //! Shared pointer to Newton-Raphson.
    /*!
     * Shared pointer to Newton-Raphson method.
     */
    boost::shared_ptr< NewtonRaphson > newtonRaphson_;

    //! NewtonRaphsonAdaptor.
    /*!
     *  NewtonRaphsonAdaptor class.
     */
    NewtonRaphsonAdaptor < ConvertMeanAnomalyToHyperbolicEccentricAnomaly > newtonRaphsonAdaptor_;

    //! Compute Kepler's function for hyperbolic orbits.
    /*!
     * Computes Kepler's function, given as:
     * \f[
     *      f( F ) = e * sinh( F ) - F - M
     * \f]
     * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, \f$ e \f$ is the
     * eccentricity, \f$ M \f$ is the mean anomaly. Currently, the case for near-parabolic orbits
     * is not included ( 0.8 < eccentricity < 1.2 ).
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
     * \return Value of Kepler's function for hyperbolic orbits.
     */
    double computeKeplersFunctionForHyperbolicOrbits_( double& hyperbolicEccentricAnomaly )
    {
        return eccentricity_ * std::sinh( hyperbolicEccentricAnomaly )
                - hyperbolicEccentricAnomaly - hyperbolicMeanAnomaly_;
    }

    //! Compute first-derivative of Kepler's function for hyperbolic orbits.
    /*!
     * Computes the first-derivative of Kepler's function, given as:
     * \f[
     *      \frac{ df( F ) } { dF } = e * cosh( F ) - 1
     * \f]
     * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, and \f$ e \f$ is
     * the eccentricity. Currently, the case for near-parabolic orbits is not included
     * ( 0.8 < eccentricity < 1.2 ).
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly
     * \return Value of first-derivative of Kepler's function for hyperbolic orbits.
     */
    double computeFirstDerivativeKeplersFunctionForHyperbolicOrbits_(
        double& hyperbolicEccentricAnomaly )
    {
        return eccentricity_ * std::cosh( hyperbolicEccentricAnomaly ) - 1.0;
    }
};

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H
