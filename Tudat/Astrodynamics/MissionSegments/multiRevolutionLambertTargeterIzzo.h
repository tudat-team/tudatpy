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
 *      PyKEP toolbox, Dario Izzo, ESA Advanced Concepts Team.
 *      Richard H. An Introduction to the Mathematics and Methods of Astrodynamics, Revised
 *          Edition.
 *      Battin, AIAA Education Series.
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 */

#ifndef TUDAT_MULTI_REVOLUTION_LAMBERT_TARGETER_IZZO_H
#define TUDAT_MULTI_REVOLUTION_LAMBERT_TARGETER_IZZO_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/MissionSegments/zeroRevolutionLambertTargeterIzzo.h"

namespace tudat
{
namespace mission_segments
{

//! Izzo Lambert targeting algorithm class (multiple revolution).
/*!
 * Solves the multi revolution (input parameter) Lambert problem using Izzo's algorithms, i.e., the
 * two point boundary value problem of what kind of unperturbed Keplerian orbit connects two 
 * positions in a certain time of flight. 
 *
 * This class uses lazy evaluation, solving the algorithm is postponed until the function call 
 * relying on the computations.
 *
 * The algorithm used is the same algorithms as the zero revolution base
 * class (Izzo's method): only the computation of the time of flight (based on Battin's x-variable)
 * now includes the number of revolutions. The branch is required for specifying the neighborhood
 * of the two roots of the problem. Furthermore, the non-dimensionalized parameters are stored
 * internally and reused for new calculations of new solutions for the same object. 
 *
 * \throws std::runtime_error For negative time of flight and gravitational parameter values.
 *              Also, the number of revolutions is checked for feasibility: no negative values 
 *              are allowed and is checked against an automatic computation of the maximum number
 *              of revolutions.
 * \throws basic_mathematics::ConvergenceException When the internal rootfinding algorithm
 *              fails to converge within a reasonable amount of iterations.
 */
class MultiRevolutionLambertTargeterIzzo : public ZeroRevolutionLambertTargeterIzzo
{
public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Multiple revolution Lambert problem constructor for problem definition.
    /*!
     * Constructs a lambert problem with all relevant parameters for a multiple revolution problem
     * (i.e. allowing multiple full revolutions in the solution). Constructor does not compute a
     * solution in order to save on computations; solving is done on first member function call by
     * user.
     * \param aCartesianPositionAtDeparture The initial position in xyz-coordinates.            [m]
     * \param aCartesianPositionAtArrival The final position in xyz-coordinates.                [m]
     * \param aTimeOfFlight The time-of-flight between the two positions.                       [s]
     * \param aGravitationalParameter The gravitational parameter of the central body for the
     *                                computed solution orbit.                            [m^3/s^2]
     * \param aNumberOfRevolutions The required number of revolutions in the problem solution
     *                             (default 0).                                                 [-]
     * \param aIsRightBranch A boolean flag to indicate whether the right or left branch
                             (corresponding to a low or high energy transfer arc) should be used in
                             the solution (default false).                                      [-]
     * \param aIsRetrograde A boolean flag to indicate retrograde motion with (default false).  [-]
     * \param aConvergenceTolerance The required tolerance for the root solving procedure.      [-]
     * \param aMaximumNumberOfIterations The required number of iterations the root solving
     *                                   procedure is allowed to go through (default 50).       [-]
     */
    MultiRevolutionLambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                                        const Eigen::Vector3d& aCartesianPositionAtArrival,
                                        const double aTimeOfFlight,
                                        const double aGravitationalParameter,
                                        // NB: numberOfRevolutions and isRightBranch placed here
                                        // to prevent retyping last three parameters when only
                                        // change in revolutions or branch is required.
                                        const int aNumberOfRevolutions = 0,
                                        const bool aIsRightBranch = false,
                                        const bool aIsRetrograde = false,
                                        const double aConvergenceTolerance = 1e-9,
                                        const int aMaximumNumberOfIterations = 50 )
        : ZeroRevolutionLambertTargeterIzzo( aCartesianPositionAtDeparture,
                                             aCartesianPositionAtArrival,
                                             aTimeOfFlight,
                                             aGravitationalParameter,
                                             aIsRetrograde,
                                             aConvergenceTolerance,
                                             aMaximumNumberOfIterations ),
          numberOfRevolutions( aNumberOfRevolutions ),
          isRightBranch( aIsRightBranch ),
          maximumNumberOfRevolutions( NO_MAXIMUM_REVOLUTIONS ) // Signifies it's not calculated yet
    { }

    //! Flag indicating that the maximum number of revolutions has not yet been calculated!
    /*!
     * When maximumNumberOfRevolutions is this amount of revolutions, it has not yet been
     * initialized. Initialization is done automatically for the user.
     */
    static const int NO_MAXIMUM_REVOLUTIONS = -1;

    //! Compute solution for N revolutions and branch.
    /*!
     * Using the constants from the problem specified in the constructor, this calculates an
     * alternative solution specified by the number of revolutions and branch selected.
     * \param aNumberOfRevolutions Number of revolutions for the alternative solution.          [-]
     * \param aIsRightBranch A boolean flag to indicate retrograde motion.                      [-]
     */
    void computeForRevolutionsAndBranch( const int aNumberOfRevolutions,
                                         const bool aIsRightBranch );

    //! Get maximum number of revolutions calculated.
    /*!
     * Returns the maximum number of revolutions possible in a solution to the problem specified.
     * \return maximumNumberOfRevolutions Maximum number of revolutions possible for the problem.
     */
    int getMaximumNumberOfRevolutions( );

protected:

    //! Sanity check number of revolutions.
    /*!
     * Computes (if not done before) the maximum number of revolutions possible and checks whether
     * the number of revolutions specified is within this solution space. This only yields a proper
     * outcome after dimensions have been transformed.
     */
    void sanityCheckNumberOfRevolutions( );

    //! Execute the solving procedure (for multiple revolutions).
    /*!
     * Executes the de-dimensionalization of the input parameter, solves the root for the time of
     * flight equation, and subsequently reconstructs the inertial velocities to obtain a solution
     * to the Lambert problem specified by the constructor of the class. If the number of
     * revolutions specified is zero, it automatically calls to the base class function to calculate
     * for zero revolutions.
     */
    void execute( );

    //! Compute time-of-flight using Lagrange's equation (for multiple revolutions).
    /*!
     * Computes the time-of-flight according to Lagrange's equation as a function of the
     * x-parameter as used by Battin, using the class's sub results, including multiple revolutions.
     * \param xParameter x parameter in Izzo's algorithm.
     * \return timeOfFlight Computed time-of-flight.
     */
    virtual double computeTimeOfFlight( const double xParameter );

    //! Solve the time of flight equation for x (for multiple revolutions).
    /*!
     * Initializes and finds the root of the equation \f$ 0 = t_f - F( x ) \f$, where \f$ F( x ) \f$
     * is a function to compute the time-of-flight as a function of \f$ x \f$, for a multiple
     * revolution solution.
     * \return xParameter Computed x parameter.
     */
    double computeRootTimeOfFlight( );

    //! Number of revolutions in problem.
    /*!
     * Number of full revolutions needed in the solution to the problem specified.
     *
     * NB: the following parameters can be readjusted by the appropriate public function and
     * therefore not declared const.
     */
    int numberOfRevolutions;

    //! Left/right branch flag.
    /*!
     * Boolean indicator for left or right branch solution, corresponding to a low or high energy
     * transfer of the multiple solutions.
     */
    bool isRightBranch;

    //! Maximum number of revolutions possible in problem.
    /*!
     * Number of full revolutions possible as a solution in the problem specified.
     *
     * NB: The following parameter is calculated on construction of the object and therefore not
     * declared const.
     */
    int maximumNumberOfRevolutions;

private:
};

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_MULTI_REVOLUTION_LAMBERT_TARGETER_IZZO_H
