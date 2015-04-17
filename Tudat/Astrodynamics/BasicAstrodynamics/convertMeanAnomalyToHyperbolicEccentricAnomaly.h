/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110214    K. Kumar          File created.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120421    K. Kumar          Removed base class; updated to set values through constructor.
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      120822    P. Musegaas       Added functionality for near-parabolic cases. Added option for
 *                                  user to specifiy initial guess.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      130116    E. Heeren         Minor changes to comments.
 *      130120    K. Kumar          Added shared-ptr typedef.
 *      130123    K. Kumar          Added note about near-parabolic cases.
 *      140110    E. Brandon        Fixed Doxygen comments.
 *
 *    References
 *      Regarding method in general, including starter values used:
 *          Wakker, K. F. Astrodynamics I + II. Lecture Notes AE4-874, Delft University of
 *              Technology, Delft, Netherlands.
 *      Regarding the performance of different starter values and performance near-parabolic.
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far]. Section available on tudat website (tudat.tudelft.nl)
 *              under issue #539.
 *      Regarding older versions of the code (120813 and before):
 *          Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *          http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 *    Notes
 *      There are known to be some issues on some systems with near-parabolic orbits that are very
 *      close to eccentricity=1.0. For these orbits, the unit tests establish that they are able
 *      to generate consistent results by converting to and from hyperbolic eccentric anomaly to a
 *      precision of 1.0e-9. If this quality of solution is not adequate for specific applications,
 *      the user should investigate the code further.
 *
 */

#ifndef TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H
#define TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H

#include <cmath>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

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

    //! Construct converter with eccentricity and mean anomaly.
    /*!
     * This constructor initializes the eccentricity and the mean anomaly for the conversion
     * routine. Optionally also the rootfinder (with termination conditions specified inside) can
     * be given.
     * \param anEccentricity Eccentricity of the orbit [-].
     * \param aHyperbolicMeanAnomaly Hyperbolic mean anomaly to convert to eccentric anomaly [rad].
     * \param useDefaultInitialGuess_ Boolean specifying whether to use default initial guess [-].
     * \param userSpecifiedInitialGuess_ Initial guess for rootfinder [rad].
     * \param aRootFinder Shared-pointer to the rootfinder that is to be used. Default is
     *          Newton-Raphson using 1000 iterations as maximum and 5.0e-15 absolute X-tolerance.
     *          Higher precision may invoke machine precision problems for some values.
     */
    ConvertMeanAnomalyToHyperbolicEccentricAnomaly(
            const double anEccentricity, const double aHyperbolicMeanAnomaly,
            const bool useDefaultInitialGuess_ = true,
            const double userSpecifiedInitialGuess_ = TUDAT_NAN,
            root_finders::RootFinderPointer aRootFinder = root_finders::RootFinderPointer( ) );

    //! Convert mean anomaly to hyperbolic eccentric anomaly.
    /*!
     * Converts mean anomaly to hyperbolic eccentric anomaly for hyperbolic orbits for all
     * eccentricities > 1.0. If the conversion fails or the eccentricity falls outside the valid
     * range, then double::NaN is returned. Calculated with an accuracy of 1.0e-14 for all
     * reasonable cases (eccentricities up to 1.0e15, mean anomalies -1.2e12 to 1.2e12).
     * \return Hyperbolic eccentric anomaly.
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
    const double hyperbolicMeanAnomaly;

    //! Boolean indicating whether default initial guess is used.
    /*!
     * A boolean that indicates whether the default initial guess is to be used. This may be
     * useful, because the current guess may not be ideal for certain regions in the
     * eccentricity<->hyperbolic mean anomaly domain.
     */
    bool useDefaultInitialGuess;

    //! User-specified initial guess for root finder.
    /*!
     * The initial guess for the root finder can be passed as a double. This may be useful, because
     * the current guess may not be ideal for certain regions in the eccentricity<->hyperbolic mean
     * anomaly spectrum.
     */
    double userSpecifiedInitialGuess;

    //! Shared pointer to the rootfinder.
    /*!
     * Shared pointer to the rootfinder. The rootfinder contains termination conditions inside.
     */
    root_finders::RootFinderPointer rootFinder;

    //! Compute Kepler's function for hyperbolic orbits.
    /*!
     * Computes Kepler's function, given as:
     * \f[
     *      f( F ) = e * sinh( F ) - F - M
     * \f]
     * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, \f$ e \f$ is the
     * eccentricity, \f$ M \f$ is the mean anomaly. All eccentricities > 1.0 are valid.
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
     * \return Value of Kepler's function for hyperbolic orbits.
     */
    double computeKeplersFunctionForHyperbolicOrbits( const double hyperbolicEccentricAnomaly )
    {
        return eccentricity * std::sinh( hyperbolicEccentricAnomaly )
                - hyperbolicEccentricAnomaly - hyperbolicMeanAnomaly;
    }

    //! Compute first-derivative of Kepler's function for hyperbolic orbits.
    /*!
     * Computes the first-derivative of Kepler's function, given as:
     * \f[
     *      \frac{ df( F ) } { dF } = e * cosh( F ) - 1
     * \f]
     * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, and \f$ e \f$ is
     * the eccentricity. All eccentricities > 1.0 are valid.
     * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly
     * \return Value of first-derivative of Kepler's function for hyperbolic orbits.
     */
    double computeFirstDerivativeKeplersFunctionForHyperbolicOrbits(
        const double hyperbolicEccentricAnomaly )
    {
        return eccentricity * std::cosh( hyperbolicEccentricAnomaly ) - 1.0;
    }
};

//! Typedef for shared-pointer to ConvertMeanAnomalyToHyperbolicEccentricAnomaly object.
typedef boost::shared_ptr< ConvertMeanAnomalyToHyperbolicEccentricAnomaly >
ConvertMeanAnomalyToHyperbolicEccentricAnomalyPointer;

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_CONVERT_MEAN_ANOMALY_TO_HYPERBOLIC_ECCENTRIC_ANOMALY_H
