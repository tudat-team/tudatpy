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
 *      120210    T. Secretin       First creation of code.
 *
 *    References
 *      Battin, R.H. An Introduction to the Mathematics and Methods of Astrodynamics,
 *          AIAA Education Series, 1999.
 *      Izzo, D. lambert_problem.h, keptoolbox.
 *
 *    Notes
 *      This code is an implementation of the method developed by Dario Izzo from ESA/ACT and
 *      publicly available at: http://keptoolbox.sourceforge.net/.
 *      After verification and validation, it was proven that this algorithm is faster and more
 *      robust than the implemented Lancaster & Blanchard and Gooding method. Notably, this method
 *      does not suffer from the near-pi singularity (pi-transfers are by nature singular).
 *
 */

#ifndef TUDAT_LAMBERT_TARGETER_IZZO_H
#define TUDAT_LAMBERT_TARGETER_IZZO_H

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"

namespace tudat
{
namespace mission_segments
{

//! Izzo Lambert targeting algorithm class.
/*!
 * Implementation of the Izzo Lambert targeting algorithm in Tudat.
 */
class LambertTargeterIzzo : public LambertTargeter
{
public:

    //! Constructor with immediate definition of parameters and execution of the algorithm.
    /*!
     * Constructor with immediate definition of parameters and execution of the algorithm.
     */
    LambertTargeterIzzo( const Eigen::Vector3d cartesianPositionAtDeparture,
                         const Eigen::Vector3d cartesianPositionAtArrival,
                         const double timeOfFlight,
                         const double gravitationalParameter,
                         const bool isRetrograde = false,
                         const double convergenceTolerance = 1e-9,
                         const int maximumNumberOfIterations = 50.0 )
        : LambertTargeter( cartesianPositionAtDeparture,
                           cartesianPositionAtArrival,
                           timeOfFlight,
                           gravitationalParameter ),
          isRetrograde_( isRetrograde ),
          convergenceTolerance_( convergenceTolerance ),
          maximumNumberOfIterations_( maximumNumberOfIterations )
    {
        // Execute algorithm.
        execute( );
    }

    //! Get radial velocity at departure.
    /*!
     * Returns the radial velocity at departure.
     * \return Radial velocity at departure.
     */
    double getRadialVelocityAtDeparture( );

    //! Get radial velocity at arrival.
    /*!
     * Returns the radial velocity at arrival.
     * \return Radial velocity at arrival.
     */
    double getRadialVelocityAtArrival( );

    //! Get transverse velocity at departure.
    /*!
     * Returns the transverse velocity at departure.
     * \return Transverse velocity at departure.
     */
    double getTransverseVelocityAtDeparture( );

    //! Get transverse velocity at arrival.
    /*!
     * Returns the transverse velocity at arrival.
     * \return Transverse velocity at arrival.
     */
    double getTransverseVelocityAtArrival( );

    //! Get semi-major axis.
    /*!
     * Returns the semi-major axis of the computed conic.
     * \return Semi-major axis.
     */
    double getSemiMajorAxis( );

protected:

    //! Execute Lambert targeting algorithm.
    /*!
     * Executes the Lambert targeting algorithm.
     */
    void execute( );

private:

    //! Retrograde motion flag.
    /*!
     * Retrograde motion flag.
     */
    const bool isRetrograde_;

    //! Convergence tolerance.
    /*!
     * Convergence tolerance for the root-finding process.
     */
    const double convergenceTolerance_;

    //! Maximum number of iterations.
    /*!
     * Maximum number of iterations for the root-finding process.
     */
    const double maximumNumberOfIterations_;
};

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_IZZO_H
