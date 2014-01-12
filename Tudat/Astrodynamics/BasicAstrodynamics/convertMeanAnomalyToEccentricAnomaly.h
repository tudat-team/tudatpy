/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110210    K. Kumar          File created.
 *      110215    E. Iorfida        Minor changes made.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120421    K. Kumar          Removed base class; updated to set values through constructor.
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      120822    P. Musegaas       Added functionality for near-parabolic cases. Added option for
 *                                  user to specifiy initial guess.
 *      120903    P. Musegaas       Added additional comments on the scope of the method and
 *                                  the 'conversion' of mean anomaly to 0 to 2*PI regime.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      130116    E. Heeren         Minor changes to comments.
 *      130120    K. Kumar          Added shared-pointer typedef.
 *      130123    K. Kumar          Added note about near-parabolic cases.
 *      140110    E. Brandon        Fixed Doxygen comments.
 *
 *    References
 *      Regarding method in general:
 *          Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *          Wakker, K. F. Astrodynamics I + II. Lecture Notes AE4-874, Delft University of
 *              Technology, Delft, Netherlands.
 *      Regarding the choice of initial guess:
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far]. Section available on tudat website (tudat.tudelft.nl)
 *              under issue #539.
 *
 *    Notes
 *      There are known to be some issues on some systems with near-parabolic orbits that are very
 *      close to eccentricity=1.0. For these orbits, the unit tests establish that they are able
 *      to generate consistent results by converting to and from eccentric anomaly to a precision
 *      of 1.0e-9. If this quality of solution is not adequate for specific applications, the user
 *      should investigate the code further.
 *
 */

#ifndef TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
#define TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H

#include <cmath>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace basic_astrodynamics
{
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
     * routine. Optionally also the rootfinder (with termination conditions specified inside) can
     * be given.
     * \param anEccentricity Eccentricity of the orbit [-].
     * \param aMeanAnomaly Mean anomaly to convert to eccentric anomaly [rad].
     * \param useDefaultInitialGuess_ Boolean specifying whether to use default initial guess [-].
     * \param userSpecifiedInitialGuess_ Initial guess for rootfinder [rad].
     * \param aRootFinder Shared-pointer to the rootfinder that is to be used. Default is
     *          Newton-Raphson using 1000 iterations as maximum and 1.0e-13 absolute X-tolerance.
     *          Higher precision may invoke machine precision problems for some values.
     */
    ConvertMeanAnomalyToEccentricAnomaly(
            const double anEccentricity, const double aMeanAnomaly,
            const bool useDefaultInitialGuess_ = true,
            const double userSpecifiedInitialGuess_ = TUDAT_NAN,
            root_finders::RootFinderPointer aRootFinder = root_finders::RootFinderPointer( ) );

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Converts mean anomaly to eccentric anomaly for elliptical orbits for all eccentricities >=
     * 0.0 and < 1.0. If the conversion fails or the eccentricity falls outside the valid range,
     * then double::NaN is returned. Calculated with an accuracy of 1.0e-13 for all cases, except
     * for some near-parabolic cases in which macine precision problems occur. These are tested
     * against an accuracy of 1.0e-9. Near-parabolic in this sense means e > 1.0-1.0e-11. Also
     * note that your mean anomaly is automatically transformed to fit within the 0 to 2.0*PI
     * spectrum.
     * \return Eccentric anomaly [rad].
     */
    double convert( );

protected:

private:

    //! Eccentricity.
    /*!
     * Eccentricity.
     */
    const double eccentricity;

    //! Mean anomaly.
    /*!
     * Mean anomaly.
     */
    double meanAnomaly;

    //! Boolean indicating whether default initial guess is used.
    /*!
     * A boolean that indicates whether the default initial guess is to be used. This may be useful,
     * because the current guess may not be ideal for certain regions in the eccentricity <-> mean
     * anomaly domain.
     */
    bool useDefaultInitialGuess;

    //! User-specified initial guess for root finder.
    /*!
     * The initial guess for the root finder can be passed as a double. This may be useful, because
     * the current guess may not be ideal for certain regions in the eccentricity<->mean anomaly
     * domain.
     */
    double userSpecifiedInitialGuess;

    //! Shared pointer to the rootfinder.
    /*!
     * Shared pointer to the rootfinder. The rootfinder contains termination conditions inside.
     */
    root_finders::RootFinderPointer rootFinder;

    //! Compute Kepler's function for elliptical orbits.
    /*!
     * Computes Kepler's function, given as:
     * \f[
     *      f( E ) = E - e * sin( E ) - M
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, \f$ e \f$ is the
     * eccentricity, \f$ M \f$ is the mean anomaly. All eccentricities >= 0.0 and < 1.0 are valid.
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of Kepler's function for elliptical orbits.
     */
    double computeKeplersFunctionForEllipticalOrbits( const double eccentricAnomaly )
    {
        return eccentricAnomaly - eccentricity * std::sin( eccentricAnomaly ) - meanAnomaly;
    }

    //! Compute first-derivative of Kepler's function for elliptical orbits.
    /*!
     * Computes the first-derivative of Kepler's function, given as:
     * \f[
     *      \frac{ df( E ) } { dE } = 1 - e * cos( E )
     * \f]
     * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, and \f$ e \f$ is the
     * eccentricity. All eccentricities >= 0.0 and < 1.0 are valid.
     * \param eccentricAnomaly Eccentric anomaly.
     * \return Value of first-derivative of Kepler's function for elliptical orbits.
     */
    double computeFirstDerivativeKeplersFunctionForEllipticalOrbits(
            const double eccentricAnomaly )
    {
        return 1.0 - eccentricity * std::cos( eccentricAnomaly );
    }
};

//! Typedef for shared-pointer to ConvertMeanAnomalyToEccentricAnomaly object.
typedef boost::shared_ptr< ConvertMeanAnomalyToEccentricAnomaly >
ConvertMeanAnomalyToEccentricAnomalyPointer;

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
