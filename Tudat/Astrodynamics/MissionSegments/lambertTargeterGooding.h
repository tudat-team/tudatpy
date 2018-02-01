/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_LAMBERT_TARGETER_GOODING_H
#define TUDAT_LAMBERT_TARGETER_GOODING_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace mission_segments
{

//! Gooding Lambert targeting algorithm class.
/*!
 * Implementation of the Gooding Lambert targeting algorithm in Tudat.
 */
class LambertTargeterGooding : public LambertTargeter
{
public:

    //! Constructor with immediate definition of parameters and execution of the algorithm.
    /*!
     * Constructor with immediate definition of parameters and execution of the algorithm.
     * \param aCartesianPositionAtDeparture The position at departure in Cartesian coordinates. [m]
     * \param aCartesianPositionAtArrival The position at arrival in Cartesian coordinates.     [m]
     * \param aTimeOfFlight The time-of-flight between departure and arrival.                   [s]
     * \param aGravitationalParameter The gravitational parameter of the main body.      [m^3 s^-2]
     * \param aRootFinder The shared-pointer to the rootfinder to be used to solve the problem. [-]
     * \sa LambertTargeter.
     */
    LambertTargeterGooding( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                            const Eigen::Vector3d& aCartesianPositionAtArrival,
                            const double aTimeOfFlight,
                            const double aGravitationalParameter,
                            root_finders::RootFinderPointer aRootFinder = 
                                root_finders::RootFinderPointer( ) );

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

    //! Shared pointer to the rootfinder.
    /*!
     * Shared pointer to the rootfinder. The rootfinder contains termination conditions inside.
     */
    root_finders::RootFinderPointer rootFinder;
};

//! Typedef for shared-pointer to LambertTargeterGooding object.
typedef boost::shared_ptr< LambertTargeterGooding > LambertTargeterGoodingPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_GOODING_H
