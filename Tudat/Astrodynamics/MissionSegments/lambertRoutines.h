/*    Copyright (c) 2010-2015, Del  ft University of Technology
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
 *      101111    E. Iorfida        File created.
 *      101126    E. Iorfida        Virtual void "solve" added.
 *      101206    E. Iorfida        setInitialLambertState added, protected variables added.
 *      101207    E. Iorfida        Set single variables, change variables names in more
 *                                  understandable ones.
 *      110113    E. Iorfida        Added necessary elements to build pointer-to-member-function to
 *                                  RootFinderAlgorithms and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to use the last
 *                                  version of Newton-Raphson code.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfVelocityVector' became
 *                                  'speed'. Also corrected 'tangential' to 'transverse'.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names
 *                                  (from heliocentric, to inertial).
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120619    T. Secretin       Converted to free functions. Added Izzo's approach.
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      140117    E. Brandon        Corrected doxygen documentation.
 *
 *    References
 *      Battin, R.H. An Introduction to the Mathematics and Methods of Astrodynamics,
 *          AIAA Education Series, 1999.
 *      Izzo, D. lambert_problem.h, keptoolbox.
 *      Gooding, R.H. A procedure for the solution of Lambert's orbital boundary-value problem,
 *          Celestial Mechanics and Dynamical Astronomy, 48:145-165, 1990.
 *
 *    Notes
 *
 */

#ifndef TUDAT_LAMBERT_ROUTINES_H
#define TUDAT_LAMBERT_ROUTINES_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace mission_segments
{

//! Solve Lambert Problem using Izzo's algorithm.
/*!
 * Solves the Lambert Problem using Izzo's algorithm. This code is an implementation of the method
 * developed by Dario Izzo from ESA/ACT and publicly available at:
 * http://keptoolbox.sourceforge.net/
 * After verification and validation, it was proven that this algorithm is faster and more robust
 * than the implemented Lancaster & Blanchard and Gooding method. Notably, this method does not
 * suffer from the near-pi singularity (pi-transfers are by nature singular). This method works in
 * adimensional units, meaning that position, time-of-flight and gravitational parameter can be
 * provided in any units, as long as they are coherent across all quantities. Results will be
 * returned in the same units as the input variables.
 * Note that while this implementation does not support multi-revolution transfers, the original
 * algorithm (see link above) does. The root-finder (Secant Method) is currently hard-coded.
 * \param cartesianPositionAtDeparture Cartesian position at departure. [Input]
 * \param cartesianPositionAtArrival Cartesian position at arrival. [Input]
 * \param timeOfFlight Time-of-flight between departure and arrival. [Input]
 * \param gravitationalParameter Gravitational parameter of the central body. [Input]
 * \param cartesianVelocityAtDeparture Velocity at departure. [Output]
 * \param cartesianVelocityAtArrival Velocity at arrival. [Output]
 * \param isRetrograde Boolean flag to indicate direction of motion. [Input, Optional]
 * \param convergenceTolerance Convergence tolerance for the root-finding process.
 *          [Input, Optional]
 * \param maximumNumberOfIterations Maximum number of iterations of the root-finding process.
 *          [Input, Optional]
 */
void solveLambertProblemIzzo( const Eigen::Vector3d& cartesianPositionAtDeparture,
                              const Eigen::Vector3d& cartesianPositionAtArrival,
                              const double timeOfFlight,
                              const double gravitationalParameter,
                              Eigen::Vector3d& cartesianVelocityAtDeparture,
                              Eigen::Vector3d& cartesianVelocityAtArrival,
                              const bool isRetrograde = false,
                              const double convergenceTolerance = 1e-9,
                              const unsigned int maximumNumberOfIterations = 50 );

//! Compute time-of-flight using Lagrange's equation.
/*!
 * Computes the time-of-flight according to Lagrange's equation as a function of the x-parameter.
 * \param xParameter x parameter in Izzo's algorithm.
 * \param semiPerimeter Semi-perimeter: \f$ s = \frac{ r_1 + r_2 + c }{ 2 } \f$.
 * \param chord Chord: \f$ c = \sqrt{ r_1^2 + r_2^2 - 2 * r_1 * r_2 \cos( \theta ) } \f$.
 * \param isLongway Boolean flag to indicate if the transfer is long-way (\f$ \theta > \pi \f$) or
 *          short-way (\f$ \theta < \pi \f$).
 * \param semiMajorAxisOfTheMinimumEnergyEllipse Semi-major axis of the minimum energy ellipse:
 *          \f$ a_m = s / 2 \f$.
 * \return timeOfFlight Computed time-of-flight.
 */
double computeTimeOfFlightIzzo( const double xParameter, const double semiPerimeter,
                                const double chord, const bool isLongway,
                                const double semiMajorAxisOfTheMinimumEnergyEllipse );

//! Solve Lambert Problem using Gooding's algorithm.
/*!
 * Solves the Lambert Problem using Lancaster and Blanchard's algorithm with further improvements
 * by Gooding.
 * The number of revolutions from departure to arrival body is zero by definition in this
 * routine. This can be made user-defined later on. The resulting trajectories are in
 * anti-clockwise (prograde) direction.
 * \param cartesianPositionAtDeparture Cartesian position at departure [m]. [Input]
 * \param cartesianPositionAtArrival Cartesian position at arrival [m]. [Input]
 * \param timeOfFlight Time-of-flight between departure and arrival [s]. [Input]
 * \param gravitationalParameter Gravitational parameter of the central body [m^3/s^2]. [Input]
 * \param cartesianVelocityAtDeparture Velocity at departure [m]. [Output]
 * \param cartesianVelocityAtArrival Velocity at arrival [m]. [Output]
 * \param rootFinder Shared-pointer to the rootfinder that is to be used. Default is Newton-Raphson
 *          using 1000 iterations as maximum and 1.0e-12 relative X-tolerance. [Input, optional]
 */
void solveLambertProblemGooding( const Eigen::Vector3d& cartesianPositionAtDeparture,
                                 const Eigen::Vector3d& cartesianPositionAtArrival,
                                 const double timeOfFlight,
                                 const double gravitationalParameter,
                                 Eigen::Vector3d& cartesianVelocityAtDeparture,
                                 Eigen::Vector3d& cartesianVelocityAtArrival,
                                 root_finders::RootFinderPointer rootFinder
                                    = root_finders::RootFinderPointer( ) );

//! Gooding Lambert functions class.
/*!
 * This class contains the auxiliary functions required by the root-finder in the Lambert routine
 * from Gooding.
 */
class LambertFunctionsGooding
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor that sets all the parameters in the Lambert functions for use in the
     * Newton-Raphson rootfinder.
     * \param aQParameter The value of the Lambert q-parameter.                                 [-]
     * \param aNormalizedTimeOfFlight The normalized time-of-flight for the Lambert problem.    [-]
     */
    LambertFunctionsGooding( const double aQParameter, const double aNormalizedTimeOfFlight )
        : qParameter( aQParameter ),
          normalizedTimeOfFlight( aNormalizedTimeOfFlight )
    { }

    //! Define general Lambert function.
    /*!
     * Defines the general Lambert function.
     * \param xParameter x-parameter.
     * \return Lambert function value.
     */
    double computeLambertFunctionGooding( const double xParameter );

    //! Define first derivative of general Lambert function.
    /*!
     * Defines the first derivative of general Lambert function.
     * \param xParameter x-parameter.
     * \return Derivative of the Lambert function.
     */
    double computeFirstDerivativeLambertFunctionGooding( const double xParameter );

    //! Define Lambert function for positive lambertEccentricAnomaly.
    /*!
     * Defines the Lambert function for a positive lambertEccentricAnomaly.
     * \param xParameter x-parameter.
     * \return Lambert function value for a positive lambertEccentricAnomaly.
     */
    double lambertFunctionPositiveGooding( const double xParameter );

    //! Define Lambert function for negative lambertEccentricAnomaly.
    /*!
     * Defines the Lambert function for a negative lambertEccentricAnomaly.
     * \param xParameter x-parameter.
     * \return Lambert function value for a negative lambertEccentricAnomaly.
     */
    double lambertFunctionNegativeGooding( const double xParameter );

    //! Define first derivative of Lambert function for positive lambertEccentricAnomaly.
    /*!
     * Defines the first derivative of Lambert function for a positive lambertEccentricAnomaly.
     * \param xParameter x-parameter.
     * \return Derivative of the Lambert function for a positive lambertEccentricAnomaly.
     */
    double lambertFirstDerivativeFunctionPositiveGooding( const double xParameter );

    //! Define first derivative of Lambert function for negative lambertEccentricAnomaly.
    /*!
     * Defines first derivative of Lambert function for negative lambertEccentricAnomaly.
     * \param xParameter x-parameter.
     * \return Derivative of the Lambert function for a negative lambertEccentricAnomaly.
     */
    double lambertFirstDerivativeFunctionNegativeGooding( const double xParameter );

protected:

private:

    //! Lambert q-parameter.
    /*!
     * Lambert q-parameter.
     */
    const double qParameter;

    //! Normalized time of flight for Lambert implementation.
    /*!
     * Normalized time of flight for Lambert implementation.
     */
    const double normalizedTimeOfFlight;
};

//! Typedef for shared-pointer to LambertFunctionsGooding object.
typedef boost::shared_ptr< LambertFunctionsGooding > LambertFunctionsGoodingPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_ROUTINES_H
