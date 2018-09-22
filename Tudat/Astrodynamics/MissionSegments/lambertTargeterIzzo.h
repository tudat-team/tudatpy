/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include <memory>

#include <Eigen/Core>

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
     * \param aCartesianPositionAtDeparture The position at departure in Cartesian coordinates. [m]
     * \param aCartesianPositionAtArrival The position at arrival in Cartesian coordinates.     [m]
     * \param aTimeOfFlight The time-of-flight between departure and arrival.                   [s]
     * \param aGravitationalParameter The gravitational parameter of the main body.      [m^3 s^-2]
     * \param isRetrograde Flag to indicate retrograde motion, (default is false).
     * \param convergenceTolerance Convergence tolerance for the root-finding process, (default
     *          is 1e-9).
     * \param maximumNumberOfIterations The maximum number of iterations for the root-finding
     *          process, (default is 50).
     * \sa LambertTargeter.
     */
    LambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                         const Eigen::Vector3d& aCartesianPositionAtArrival,
                         const double aTimeOfFlight,
                         const double aGravitationalParameter,
                         const bool isRetrograde = false,
                         const double convergenceTolerance = 1e-9,
                         const int maximumNumberOfIterations = 50 )
        : LambertTargeter( aCartesianPositionAtDeparture,
                           aCartesianPositionAtArrival,
                           aTimeOfFlight,
                           aGravitationalParameter ),
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

//! Typedef for shared-pointer to LambertTargeterIzzo object.
typedef std::shared_ptr< LambertTargeterIzzo > LambertTargeterIzzoPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_IZZO_H
