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
 *      Novak, D.M. Methods and Tools for Preliminary Low Thrust Mission Analysis. PhD Thesis,
 *          University of Glasgow, UK, 2012.
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128147,
 *          2011.
 *      Wall, B.J. and D. Novak. A 3D Shape-Based Approximation Method for Low-Thrust Trajectory
 *          Design, Advances in the Astronautical Sciences 142(2), pp. 1163-1176, 2012.
 *
 *    Notes
 *      This file contains a specific oscillating function, which is described by Novak and Vasile
 *      [2011] (See also Novak [2012] and Wall and Novak [2012]). This function is a mathematical
 *      representation of the (out-of-plane) elevation angle, phi (in spherical coordinates), of a
 *      thrusting spacecraft. This oscillating function can be used to approximate a
 *      continuous-thrust trajectory (also referred to as a low-thrust trajectory), when combined
 *      with a function which represents the radial position, r (in spherical coordinates), of a
 *      thrusting spacecraft. The spherical coordinate system that is used for the calculations and
 *      the descriptions is taken from Novak and Vasile [2011].
 */

#ifndef TUDAT_OSCILLATING_FUNCTION_NOVAK_H
#define TUDAT_OSCILLATING_FUNCTION_NOVAK_H

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/function.h"

namespace tudat
{
namespace mission_segments
{

//! Oscillating function approximating the (out-of-plane) elevation angle of a thrusting
//! spacecraft.
/*!
 * This class contains an oscillating function, and its exact first and second derivative w.r.t.
 * the azimuthal angle \f$ \theta \f$. It is a mathematical function, and can be used to
 * approximate the out-of-plane motion of a thrusting spacecraft. The function is documented by
 * Novak and Vasile [2011]. It approximates the (out-of-plane) elevation angle of the spacecraft.
 * The (out-of-plane) elevation angle is parameterized in terms of an independent variable, the
 * (in-plane) azimuthal angle \f$ \theta \f$. The first and second derivative are taken with
 * respect to this independent variable.
 *
 * The function is completely described by the independent variable \f$ \theta \f$ and a set of
 * four parameters. These parameters are related to the boundary conditions of the problem: two
 * initial conditions and two final conditions. These parameters are passed to this class as
 * boost::functions to facilitate the flexible external manipulation of their values.
 *
 */
class OscillatingFunctionNovak : public basic_mathematics::Function< >
{
public:

    //! Default constructor with immediate definition of parameters.
    /*!
     * Default constructor with immediate definition of parameters through boost::functions.
     * This setup allows for a flexible external manipulation of the values of the parameters.
     *
     * \param aSetOfBoundaryParameters A set of four parameters, related to the boundary conditions.
     *      These parameters are equivalent to the parameters b0, b1, b2, and b3 from Novak and
     *      Vasile [2011]. The order is important!
     *      aSetOfBoundaryParameters.first( 0 ) = b0
     *      aSetOfBoundaryParameters.first( 1 ) = b1
     *      aSetOfBoundaryParameters.second( 0 ) = b2
     *      aSetOfBoundaryParameters.second( 1 ) = b3
     */
    OscillatingFunctionNovak( const boost::function< std::pair< Eigen::Vector2d,
                              const Eigen::Vector2d >(  ) > aSetOfBoundaryParameters ) :
        boundaryParameters_( aSetOfBoundaryParameters ){  }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~OscillatingFunctionNovak(  ){  }

    //! Evaluate the function value for a given (in-plane) azimuthal angle.
    /*!
     * Evaluates the oscillating function for a given (in-plane) azimuthal angle \f$ \theta \f$.
     *
     * The function is given by Novak and Vasile [2011] as:
     *
     * \f[
     *      \phi ( \theta ) = \left( b_0 + b_1 \theta \right) \cos \theta +
     *      \left( b_2 + b_3 \theta \right) sin \theta
     * \f]
     *
     * in which \f$ \phi \f$ is the (out-of-plane) elevation angle of the spacecraft,
     * \f$ \theta \f$ is the (in-plane) azimuthal angle and \f$ ( b_0 ,b_1, b_2, b_3) \f$ are
     * parameters that define the shape of the trajectory.
     *
     * \param anAzimuthalAngle Value of the (in-plane) azimuthal angle [rad].
     * \return The (out-of-plane) elevation angle at \f$ \theta \f$ [m].
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
     *      \frac{ d \phi }{ d \theta } ( \theta ) =
     *      \left( b_1 + b_2 + b_3 \theta \right) \cos \theta -
     *      \left( b_0 + b_1 \theta - b_3 \right) \sin \theta
     * \f]
     *
     * The second derivative of the function is:
     *
     * \f[
     *      \frac{ d^2 \phi }{ d \theta^2 } ( \theta ) =
     *      \left( - b_0 - b_1 \theta + 2 b_3 \right) \cos \theta -
     *      \left( 2 b_1 + b_2 + b_3 \theta \right) \sin \theta
     * \f]
     *
     * The third derivative of the function is:
     *
     * \f[
     *      \frac{ d^3 \phi }{ d \theta^3 } ( \theta ) =
     *      - \left( 3 b_1 + b_2 + b_3 \theta \right) \cos \theta -
     *      \left( b_0 + b_1 \theta - 3 b_3 \right) \sin \theta
     * \f]
     *
     * In these formulas \f$ \phi \f$ is the (out-of-plane) elevation angle of the spacecraft,
     * \f$ \theta \f$ is the (in-plane) azimuthal angle, the prime indicates the derivative with
     * respect to \f$ \theta \f$, and \f$ ( b_0 ,b_1, b_2, b_3 ) \f$ are parameters that define the
     * shape of the trajectory.
     *
     * \param order              Order of the derivative to evaluate. Can be either 1, 2 or 3.
     * \param anAzimuthalAngle   Value of the (in-plane) azimuthal angle [rad].
     * \return Value of the derivative of order 1, 2 or 3 of the function [-, rad^-1, rad^-2].
     */
    double computeDerivative( const unsigned int order, const double anAzimuthalAngle );

    //! Compute the definite integral of the function.
    /*!
     * This function throws a runtime error when it is called. The integral of the inverse
     * polynomial function is not part of the shape-based method.
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

    //! The parameters of the function, related to the boundary conditions.
    /*!
     * The four parameters of the function. These parameters determine the shape of the trajectory
     * and are dependent on the initial and final boundary conditions.
     *
     * The first two parameters are dependent on the initial boundary conditions and the last two
     * parameters are dependent on the final boundary conditions. The parameters correspond to
     * Novak and Vasile [2011] in the following manner: the first two parameters are \f$ b_0 \f$
     * and \f$ b_1 \f$; the last two parameters are \f$ b_2 \f$ and \f$ b_3 \f$.
     *
     * The order is important!
     * boundaryParameters_.first( 0 ) = b_0
     * boundaryParameters_.first( 1 ) = b_1
     * boundaryParameters_.second( 0 ) = b_2
     * boundaryParameters_.second( 1 ) = b_3
     */
    boost::function< std::pair< Eigen::Vector2d, Eigen::Vector2d >(  ) > boundaryParameters_;
};

//! Typedef for shared-pointer to oscillatingFunctionNovak object.
typedef boost::shared_ptr< OscillatingFunctionNovak > OscillatingFunctionNovakPointer;

} // namespace mission_segments
} // namespace tudat

#endif // TUDAT_OSCILLATING_FUNCTION_NOVAK_H
