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
 *      lambertTargeterIzzo.h/.cpp source files, tudat revision 455.
 *      PyKEP toolbox, Dario Izzo, ESA Advanced Concepts Team.
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 *    Notes
 *      This is a new implementation of the lambertTargeterIzzo class, for better adaptability and
 *      extension towards subclasses and future improvements/additions. Therefore, it replaces the
 *      lambertTargeterIzzo class while still providing the same functionality.
 *
 */

#ifndef TUDAT_ZERO_REVOLUTION_LAMBERT_TARGETER_IZZO_H
#define TUDAT_ZERO_REVOLUTION_LAMBERT_TARGETER_IZZO_H

#include <Eigen/Core> // For Vector3d

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h" // base class

namespace tudat
{
namespace mission_segments
{

//! Izzo Lambert targeting algorithm class (zero revolution).
/*!
 * Solves the zero revolution Lambert problem using Izzo's algorithms, i.e. the two point boundary 
 * value problem of what kind of unperturbed Keplerian orbit connects two positions in a certain 
 * time of flight. 
 * 
 * This class uses lazy evaluation, solving the algorithm is postponed until the function call 
 * relying on the computations.
 * 
 * The algorithm used is based on the algorithm used in Izzo's PyKEP Keplerian toolbox, but is
 * implemented such that only one solution gets computed at a time. The approach from Izzo
 * transforms the problem into dimensionless parameters in terms of the initial position and
 * computes the solution in a 2D plane solving for Battin's x-variable. Contrasting the Izzo
 * approach, the non-dimensionalized parameters are stored internally in order to speed up future
 * recalculations.
 * 
 * \throws std::runtime_error For negative time of flight and gravitational parameter values.
 * \throws basic_mathematics::ConvergenceException When the internal rootfinding algorithm
 *              fails to converge within a reasonable amount of iterations.
 */
class ZeroRevolutionLambertTargeterIzzo : public LambertTargeter
{
public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Zero revolution Lambert problem constructor for problem definition.
    /*!
     * Constructs a lambert problem with all relevant parameters for a zero revolution problem (i.e.
     * without a full revolution in the solution). Constructor does not compute a solution in order
     * to save on computations; solving is done on first member function call by user.
     * \param aCartesianPositionAtDeparture The initial position in xyz-coordinates.            [m]
     * \param aCartesianPositionAtArrival The final position in xyz-coordinates.                [m]
     * \param aTimeOfFlight The time-of-flight between the two positions.                       [s]
     * \param aGravitationalParameter The gravitational parameter of the central body for the
     *                                computed solution orbit.                            [m^3/s^2]
     * \param aIsRetrograde A boolean flag to indicate retrograde motion with (default false).  [-]
     * \param aConvergenceTolerance The required tolerance for the root solving procedure.      [-]
     * \param aMaximumNumberOfIterations The required number of iterations the root solving
     *                                   procedure is allowed to go through (default 50).       [-]
     */
    ZeroRevolutionLambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                                       const Eigen::Vector3d& aCartesianPositionAtArrival,
                                       const double aTimeOfFlight,
                                       const double aGravitationalParameter,
                                       const bool aIsRetrograde = false,
                                       const double aConvergenceTolerance = 1.0e-9,
                                       const int aMaximumNumberOfIterations = 50 )
        : LambertTargeter( aCartesianPositionAtDeparture,
                           aCartesianPositionAtArrival,
                           aTimeOfFlight,
                           aGravitationalParameter ),
          isRetrograde( aIsRetrograde ),
          convergenceTolerance( aConvergenceTolerance ),
          maximumNumberOfIterations( aMaximumNumberOfIterations ),
          transformed( false )
    { }

    //! Get radial velocity at departure.
    /*!
     * Calculates the radial velocity at departure position in S.I. units.
     * \return radialVelocityAtDeparture Radial velocity at departure.                        [m/s]
     */
    double getRadialVelocityAtDeparture( );
    
    //! Get transverse velocity at departure.
    /*!
     * Calculates and returns transverse velocity at departure position in S.I. units.
     * \return transverseVelocityAtDeparture Transverse velocity at departure.                [m/s]
     */
    double getTransverseVelocityAtDeparture( );
    
    //! Get radial velocity at arrival.
    /*!
     * Calculates and returns radial velocity at arrival position in S.I. units.
     * \return radialVelocityAtArrival Radial velocity at arrival.                            [m/s]
     */
    double getRadialVelocityAtArrival( );
    
    //! Get transverse velocity at arrival.
    /*!
     * Calculates and returns transverse velocity at arrival position in S.I. units.
     * \return transverseVelocityAtDeparture Transverse velocity at departure.                [m/s]
     */
    double getTransverseVelocityAtArrival( );
    
    //! Get semi-major axis.
    /*!
     * Calculates and returns semi-major axis in S.I. units.
     * \return semiMajorAxis Semi-major axis of the problem's solution.                       [m/s]
     */
    double getSemiMajorAxis( );
    
    // what about retrograde flag, tolerance and max number of iterations? Might be useful to check.

protected:

    //! Execute the solving procedure.
    /*!
     * Solves the problem using dimensionless parameters by finding the root for the x-paramater and
     * recovering the dimensions of the corresponding orbit parameters.
     */
    virtual void execute( );

    //! Sanity check time of flight.
    /*!
     * Checks whether time of flight is positive: if negative, throws an error.
     */
    void sanityCheckTimeOfFlight( );

    //! Sanity check gravitational parameter.
    /*!
     * Checks whether gravitational parameter is positive: if negative, throws an error.
     */
    void sanityCheckGravitationalParameter( );

    //! Transform input to sub results in adimensional units.
    /*!
     * Transforms the input to the adimensional units required for the Izzo theory. It also
     * stores relevant (no longer changing) sub results for recovering velocities.
     */
    void transformDimensions( );
    
    //! Compute time-of-flight using Lagrange's equation.
    /*!
     * Computes the time-of-flight according to Lagrange's equation as a function of the
     * x-parameter as used by Battin, using the problem's sub results.
     * \param xParameter x-parameter in Izzo's theory.
     * \return timeOfFlight Computed dimensionless time-of-flight.
     */
    virtual double computeTimeOfFlight( const double xParameter );
    
    //! Solve the time of flight equation for x.
    /*!
     * Initializes and finds the root of the equation \f$ 0 = t_f - F( x ) \f$, where \f$ F( x ) \f$
     * is a function to compute the time-of-flight as a function of \f$ x \f$, for a  zero
     * revolution solution.
     * \return xParameter Computed x parameter.
     */
    virtual double computeRootTimeOfFlight( );
    
    //! Compute velocities at departure and arrival.
    /*!
     * Computes the inertial departure and arrival velocities according to Izzo's theory and
     * reassigns the dimensions of the input to the velocities.
     * \param xParameter x-parameter in Izzo's theory.
     */
    void computeVelocities( const double xParameter );
    
    //! Retrograde motion flag.
    /*!
     * Boolean indicator for retrograde motion.
     */
    const bool isRetrograde;

    //! Convergence tolerance.
    /*!
     * Convergence tolerance for the root-finding process.
     */
    const double convergenceTolerance;

    //! Maximum number of iterations.
    /*!
     * Maximum number of iterations allowed for the root-finding process.
     */
    const double maximumNumberOfIterations;

    //! Transformed flag.
    /*!
     * Boolean flag to indicate whether transformation of input dimensions has already been
     * performed or not.
     */
    bool transformed;

    //! Normalizing value for velocity.
    /*!
     * Calculated parameter to make velocities dimensionless.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double velocityNormalizingValue;

    //! Normalized time-of-flight.
    /*!
     * Calculated parameter expressing time-of-flight in dimensionless units.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double normalizedTimeOfFlight;

    //! Normalized radius at arrival.
    /*!
     * Calculated parameter expressing radius at arrival in dimensionless units.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double normalizedRadiusAtArrival;

    //! Normalized chord.
    /*!
     * Calculated parameter expressing chord (norm of difference between both position vectors) in
     * dimensionless units.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double normalizedChord;

    //! Normalized semi-perimeter.
    /*!
     * Calculated parameter expressing semi-perimeter of the triangle departure position, arrival
     * position and orbit focus in dimensionless units.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double normalizedSemiPerimeter;

    //! Longway flag.
    /*!
     * Calculated boolean flag to indicate long or short transfer angle between positions.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    bool isLongway;

    //! Normalized semi-major axis of minimum energy ellipse.
    /*!
     * Calculated parameter expressing the problem's minimum energy ellipse semi-major axis in
     * dimensionless units.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double normalizedMinimumEnergySemiMajorAxis;

    //! Transfer angle.
    /*!
     * Calculated angle between arrival and departure position.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double transferAngle;

    //! Lambda parameter.
    /*!
     * Calculated parameter necessary for recovering velocities.
     *
     * NB: This is calculated as a function of specified problem geometry, and
     * therefore not declared const.
     */
    double lambdaParameter;

private:
};

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_ZERO_REVOLUTION_LAMBERT_TARGETER_IZZO_H
