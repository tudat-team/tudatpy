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
 *                                  'velocity'. Also corrected 'tangential' to 'transverse'.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names
 *                                  (from heliocentric, to inertial).
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120620    T. Secretin       Removed inheritance from TrajectoryDesignMethod class.
 *                                  Turned into base class for Lambert Targeters. Previous code
 *                                  adapted and moved to LambertTargeterGooding.h.
 *      120704    P. Musegaas       Moved getInertialVelocityVectors to header. Minor changes.
 *      130117    R.C.A. Boon       Added solved boolean flag to data members and public member
 *                                  functions (with help from S. Billemont). Updated copyright
 *                                  notice to 2013.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      140117    E. Brandon        Corrected doxygen documentation.
 *
 *    References
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 *    Notes
 *
 */

#ifndef TUDAT_LAMBERT_TARGETER_H
#define TUDAT_LAMBERT_TARGETER_H

#include <iostream>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace mission_segments
{

//! Lambert targeting algorithm class.
/*!
 * Implementation of Lambert targeting algorithm in Tudat.
 */
class LambertTargeter
{
public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Default constructor.
    /*!
     * Default constructor that only initializes parameters. Execution of solving routine happens
     * on member function call.
     * \param aCartesianPositionAtDeparture The position at departure in Cartesian coordinates. [m]
     * \param aCartesianPositionAtArrival The position at arrival in Cartesian coordinates.     [m]
     * \param aTimeOfFlight The time-of-flight between departure and arrival.                   [s]
     * \param aGravitationalParameter The gravitational parameter of the main body.      [m^3 s^-2]
     */
    LambertTargeter( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                     const Eigen::Vector3d& aCartesianPositionAtArrival,
                     const double aTimeOfFlight,
                     const double aGravitationalParameter )
        : cartesianPositionAtDeparture( aCartesianPositionAtDeparture ),
          cartesianPositionAtArrival( aCartesianPositionAtArrival ),
          timeOfFlight( aTimeOfFlight ),
          gravitationalParameter( aGravitationalParameter ),
          cartesianVelocityAtDeparture( Eigen::Vector3d::Zero( ) ),
          cartesianVelocityAtArrival( Eigen::Vector3d::Zero( ) ),
          solved( false )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~LambertTargeter( ){ }

    //! Get inertial velocity at departure.
    /*!
     * Returns the inertial velocity at departure ( heliocentric or planetocentric ).
     * \return Inertial velocity at departure.
     */
    Eigen::Vector3d getInertialVelocityAtDeparture( )
    {
        // If no solution has yet been computed, solve.
        if ( !solved ) execute( );

        // Return computed Cartesian velocity at departure.
        return cartesianVelocityAtDeparture;
    }

    //! Get inertial velocity at arrival.
    /*!
     * Returns the inertial velocity at arrival ( heliocentric or planetocentric ).
     * \return Inertial velocity at arrival.
     */
    Eigen::Vector3d getInertialVelocityAtArrival( )
    {
        // If no solution has yet been computed, solve.
        if ( !solved ) execute( );

        // Return computed Cartesian velocity at arrival.
        return cartesianVelocityAtArrival;
    }

    //! Get inertial velocity vectors.
    /*!
     * Returns a pair of vectors, the first corresponding to the velocity vector along the transfer
     * ellipse, the second to the the velocity vector at arrival.
     * \return Pair of velocity vectors.
     */
    std::pair< Eigen::Vector3d, Eigen::Vector3d > getInertialVelocityVectors( )
    {
        // If no solution has yet been computed, solve.
        if ( !solved ) execute( );

        // Return computed Cartesian velocities at departure and arrival.
        return std::pair< Eigen::Vector3d, Eigen::Vector3d >(
                    cartesianVelocityAtDeparture, cartesianVelocityAtArrival );
    }

protected:

    //! Execute Lambert targeting algorithm.
    /*!
     * Executes the Lambert targeting algorithm. Since the parameters of the Lambert routine are
     * set directly in the constructor, the same LambertTargeter object cannot be reused for a
     * different problem. This makes the class more robust, as all parameters are consistent with a
     * single problem at all times. Since each object corresponds to a unique problem, there is no
     * need to call this function more than once. Therefore, this function is protected and run at
     * construction.
     */
    virtual void execute( ) = 0;

    //! Cartesian position at departure.
    /*!
     * Cartesian position at departure.
     */
    const Eigen::Vector3d cartesianPositionAtDeparture;

    //! Cartesian position at arrival.
    /*!
     * Cartesian position at arrival.
     */
    const Eigen::Vector3d cartesianPositionAtArrival;

    //! Time-of-flight.
    /*!
     * Time-of-flight.
     */
    const double timeOfFlight;

    //! Gravitational parameter.
    /*!
     * Gravitational parameter.
     */
    const double gravitationalParameter;

    //! Cartesian velocity at departure.
    /*!
     * Cartesian velocity at departure.
     */
    Eigen::Vector3d cartesianVelocityAtDeparture;

    //! Cartesian velocity at arrival.
    /*!
     * Cartesian velocity at arrival.
     */
    Eigen::Vector3d cartesianVelocityAtArrival;

    //! Solved boolean flag.
    /*!
     * Boolean flag to indicate the problem has a computed outcome.
     */
    bool solved;

private:
};

//! Typedef for shared-pointer to LambertTargeter.
typedef boost::shared_ptr< LambertTargeter > LambertTargeterPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_H
