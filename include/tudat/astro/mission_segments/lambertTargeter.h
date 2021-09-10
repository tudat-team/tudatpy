/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 */

#ifndef TUDAT_LAMBERT_TARGETER_H
#define TUDAT_LAMBERT_TARGETER_H

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/mathematicalConstants.h"

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

    Eigen::Vector6d getDepartureState( )  const
    {
        return ( Eigen::Vector6d( ) << cartesianPositionAtDeparture, cartesianVelocityAtDeparture ).finished( );
    }

    double getCentralBodyGravitationalParameter( ) const
    {
        return gravitationalParameter;
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

Eigen::Vector6d getLambertTargeterInitialKeplerianState(
        const LambertTargeter& lambertTargeter );

Eigen::Vector6d getLambertTargeterKeplerianStateDuringTransfer(
        const LambertTargeter& lambertTargeter,
        const double timeAfterDeparture );

Eigen::Vector6d getLambertTargeterCartesianStateDuringTransfer(
        const LambertTargeter& lambertTargeter,
        const double timeAfterDeparture );

//! Typedef for shared-pointer to LambertTargeter.
typedef std::shared_ptr< LambertTargeter > LambertTargeterPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_H
