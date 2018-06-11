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
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128147,
 *          2011.
 *      Wall, B.J. and D. Novak. A 3D Shape-Based Approximation Method for Low-Thrust Trajectory
 *          Design, Advances in the Astronautical Sciences 142(2), pp. 1163-1176, 2012.
 *      Wall, B.J., Pols, B. and B. Lanktree. Shape-Based Approximation Method for Low-Thrust
 *          Interception and Rendezvous Trajectory Design, Advances in the Astronautical Sciences
 *          136(2), pp. 1447-1458, 2010.
 *
 *    Notes
 *      This file contains the improved inverse polynomial function, which is described by Wall et
 *      al. [2010] (See also Wall and Novak [2012]). This function is a mathematical representation
 *      of the radial position, r (in spherical coordinates), of a thrusting spacecraft. The
 *      improved inverse polynomial function can be used to approximate a continuous-thrust
 *      trajectory (also referred to as a low-thrust trajectory), when combined with a function
 *      which represents the (out-of-plane) elevation angle, phi (in spherical coordinates), of a
 *      thrusting spacecraft. The spherical coordinate system that is used for the calculations and
 *      the descriptions is taken from Novak and Vasile [2011].
 */

#ifndef TUDAT_IMPROVED_INVERSE_POLYNOMIAL_WALL_H
#define TUDAT_IMPROVED_INVERSE_POLYNOMIAL_WALL_H

#include <functional>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/function.h"

namespace tudat
{
namespace mission_segments
{

//! Improved inverse polynomial function approximating the radial position of a thrusting
//! spacecraft.
/*!
 * This class contains an improved inverse polynomial function, and its exact first and second
 * derivative w.r.t. the azimuthal angle \f$ \theta \f$. It is a mathematical function, and can be
 * used to approximate the radial position of a thrusting spacecraft. The function is documented
 * by Wall et al. [2010] as the "improved inverse polynomial method". The radial position of the
 * spacecraft is parameterized in terms of an independent variable, the (in-plane) azimuthal angle
 * \f$ \theta \f$. The first and second derivative are taken with respect to this independent
 * variable.
 *
 * The function is completely described by the independent variable \f$ \theta \f$ and a set of
 * seven parameters. These parameters are related to the boundary conditions of the problem: three
 * initial conditions, three final conditions and the total time-of-flight. One of the parameters
 * is used to solve for the time-of-flight and is therefore a time-dependent parameter. These
 * seven parameters are passed to this class as std::functions to facilitate the flexible
 * external manipulation of their values.
 *
 */
class ImprovedInversePolynomialWall : public basic_mathematics::Function< >
{
public:

    //! Default constructor with immediate definition of parameters.
    /*!
     * Default constructor with immediate definition of parameters through std::functions.
     * This setup allows for a flexible external manipulation of the values of the parameters.
     *
     * \param aTimeDependentParameter The parameter that is used to solve for the time-of-flight.
     *      This parameter is equivalent to parameter d from Wall et al. [2010].
     * \param aSetOfBoundaryParameters A set of six parameters, related to the position boundary
     *      conditions.
     *      These parameters are equivalent to the parameters a, b, c, e, f, and g from Wall et al.
     *      [2010]. The order is important!
     *      aSetOfBoundaryParameters.first( 0 ) = a
     *      aSetOfBoundaryParameters.first( 1 ) = b
     *      aSetOfBoundaryParameters.first( 2 ) = c
     *      aSetOfBoundaryParameters.second( 0 ) = e
     *      aSetOfBoundaryParameters.second( 1 ) = f
     *      aSetOfBoundaryParameters.second( 2 ) = g
     */
    ImprovedInversePolynomialWall( const std::function< double(  ) > aTimeDependentParameter,
                                   const std::function< std::pair< Eigen::Vector3d ,
                                   Eigen::Vector3d >(  ) > aSetOfBoundaryParameters ) :
        timeDependentParameter_( aTimeDependentParameter ),
        boundaryParameters_( aSetOfBoundaryParameters ){  }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ImprovedInversePolynomialWall(  ){  }

    //! Evaluate the function value for a given (in-plane) azimuthal angle.
    /*!
     * Evaluates the inverse polynomial function for a given (in-plane) azimuthal angle
     * \f$ \theta \f$.
     *
     * The (improved version of the) inverse polynomial function is given by Wall et al. [2010] as:
     *
     * \f[
     *      r ( \theta ) = \frac{ 1 }
     *      { a + b \cos( \theta + c ) + d \theta^3 + e \theta^4 + f \theta^5 + g \theta^6 }
     * \f]
     *
     * in which \f$ r \f$ is the radial position of the spacecraft, \f$ \theta \f$ is the
     * (in-plane) azimuthal angle and \f$ ( a, b, c, d, e, f, g ) \f$ are parameters that define
     * the shape of the trajectory. The parameter \f$ d \f$ is the time-dependent parameter which
     * is used to solve for a required time-of-flight.
     *
     * \param anAzimuthalAngle Value of the (in-plane) azimuthal angle [rad].
     * \return The radial position at \f$ \theta \f$ [m].
     */
    double evaluate( const double anAzimuthalAngle );

    //! Compute the derivative of the function.
    /*!
     * Computes the derivative of a given order of the function, for a given (in-plane) azimuthal
     * angle \f$ \theta \f$. It only evaluates the first, second or third derivative, therefore the
     * order should be 1, 2 or 3.
     *
     * The first derivative of the function is:
     *
     * \f[
     *      \frac{ dr }{ d\theta } ( \theta ) =
     *      -r^2 \left( -b \sin \left( \theta + c \right) + 3 d \theta^2 + 4 e \theta^3 + 5 f
     *      \theta^4 + 6 g \theta^5 \right)
     * \f]
     *
     * The second derivative of the function is:
     *
     * \f[
     *      \frac{ d^2 r }{ d \theta^2 } ( \theta ) =
     *      2 r^3 \left( -b \sin \left( \theta + c \right) + 3 d \theta^2 + 4 e \theta^3 + 5 f
     *      \theta^4 + 6 g \theta^5 \right)^2 -
     *      r^2 \left( -b \cos \left( \theta + c \right) + 6 d \theta + 12 e \theta^2 + 20 f
     *      \theta^3 + 30 g \theta^4 \right)
     * \f]
     *
     * When substituting for the first derivative, the equation becomes:
     *
     * \f[
     *      \frac{ d^2 r }{ d \theta^2 } ( \theta ) =
     *      \frac{ 2 }{ r } r^{ \prime 2 } - r^2 \left( -b \cos \left( \theta + c \right)
     *      + 6 d \theta + 12 e \theta^2 + 20 f \theta^3 + 30 g \theta^4 \right)
     * \f]
     *
     * The third derivative of the function is:
     *
     * \f[
     *      \frac{ d^3 r }{ d \theta^3 } ( \theta ) =
     *      -6 r^4 \left( -b \sin \left( \theta + c \right) + 3 d \theta^2 + 4 e \theta^3
     *          + 5 f \theta^4 + 6 g \theta^5 \right)^3
     *      + 6 r^3 \left( -b \sin \left( \theta + c \right) + 3 d \theta^2 + 4 e \theta^3
     *          + 5 f \theta^4 + 6 g \theta^5 \right) \left( -b \cos \left( \theta + c \right)
     *          + d \theta + 12 e \theta^2 + 20 f \theta^3 + 30 g \theta^4 \right)
     *      - r^2 \left( b \sin \left( \theta + c \right) + 6 d + 24 e \theta + 60 f \theta^2
     *          + 120 g \theta^3 \right)
     * \f]
     *
     * When substituting for the first and second derivative, the equation becomes:
     *
     * \f[
     *      \frac{ d^3 r }{ d \theta^3 } ( \theta ) =
     *      \frac{ 6 r^{ \prime 3 } }{ r^2 } - \frac{ 6 r^{ \prime } }{ r }
     *          \left( \frac{ 2 r^{ \prime 2 } }{ r } - r^{ \prime \prime } \right)
     *      - r^2 \left( b \sin \left( \theta + c \right) + 6 d + 24 e \theta + 60 f \theta^2
     *          + 120 g \theta^3 \right)
     * \f]
     *
     * In these formulas \f$ r \f$ is the radial position of the spacecraft, \f$ \theta \f$ is the
     * (in-plane) azimuthal angle, the prime indicates the derivative with respect to
     * \f$ \theta \f$, and \f$ ( a, b, c, d, e, f, g ) \f$ are parameters that define the shape of
     * the trajectory. The parameter \f$ d \f$ is the time-dependent parameter which is used to
     * solve for a required time-of-flight.
     *
     * \param order             Order of the derivative to evaluate. Can be either 1, 2 or 3.
     * \param anAzimuthalAngle  Value of the (in-plane) azimuthal angle [rad].
     * \return Value of the derivative of order 1, 2 or 3 of the function [-, rad^-1, rad^-2].
     */
    double computeDerivative( const unsigned int order, const double anAzimuthalAngle );

    //! Compute the definite integral of the function.
    /*!
     * This function throws a runtime error when it is called. The integral of the inverse
     * polynomial function is not part of the shape-based method and has not been implemented.
     *
     * \param order Order of the integral to evaluate.
     * \param lowerBound Integration lower bound (integrate from this point).
     * \param upperBound Integration upper bound (integrate to this point).
     * \return Value of the integral.
     *
     * \throws std::runtime_error when this function is used.
     */
    double computeDefiniteIntegral( const unsigned int order,
                                    const double lowerBound,
                                    const double upperBound );

protected:

private:

    //! The time-dependent parameter.
    /*!
     * The time-dependent parameter, which is used to satisfy a required time-of-flight.
     */
    std::function< double(  ) > timeDependentParameter_;

    //! The parameters of the function, related to the boundary conditions.
    /*!
     * The six time-independent parameters of the function. These parameters determine the shape of
     * the trajectory and are dependent on the initial and final boundary conditions.
     *
     * The first three parameters are dependent on the initial boundary conditions and the last
     * three parameters are dependent on the final boundary conditions plus the time-dependent
     * parameter 'aTimeDependentParameter'. The parameters correspond to Wall et al. [2010] in the
     * following manner: the first three parameters are a, b and c; the last three parameters are
     * e, f, and g.
     *
     * The order is important!
     * boundaryParameters_.first( 0 ) = a
     * boundaryParameters_.first( 1 ) = b
     * boundaryParameters_.first( 2 ) = c
     * boundaryParameters_.second( 0 ) = e
     * boundaryParameters_.second( 1 ) = f
     * boundaryParameters_.second( 2 ) = g
     */
    std::function< std::pair< Eigen::Vector3d , Eigen::Vector3d >(  ) > boundaryParameters_;
};

//! Typedef for shared-pointer to ImprovedInversePolynomialWall object.
typedef std::shared_ptr< ImprovedInversePolynomialWall > ImprovedInversePolynomialWallPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_IMPROVED_INVERSE_POLYNOMIAL_WALL_H
