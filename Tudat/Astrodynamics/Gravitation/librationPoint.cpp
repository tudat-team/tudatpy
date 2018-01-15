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
 *      van der Ham, L. Interplanetary trajectory design using dynamical systems theory,
 *          MSc thesis, Delft University of Technology, Delft, The Netherlands, 2012.
 *      Mireles James, J.D. Celestial Mechanics Notes Set 4: The Circular Restricted Three Body
 *          Problem, 2006, http://www.math.utexas.edu/users/jjames/hw4Notes.pdf,
 *          last accessed: 18th May, 2012.
 *
 *    Notes
 *      WARNING: There seems to be a bug in the computation of the L3 location!
 *
 */

#include <cmath>
#include <stdexcept>

#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

using namespace root_finders;
using namespace basic_mathematics;

//! Compute location of Lagrange libration point.
void LibrationPoint::computeLocationOfLibrationPoint(
        LagrangeLibrationPoints lagrangeLibrationPoint )
{
    using std::pow;
    using std::sqrt;

    // Set functions for Newton-Raphson based on collinear libration point passed as input
    // parameter, or computed locations directly of equilateral libration points.
    switch( lagrangeLibrationPoint )
    {
    case l1:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = boost::make_shared< UnivariateProxy >(
                    boost::bind( &LibrationPoint::computeL1LocationFunction, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, boost::bind( &LibrationPoint::
                computeL1FirstDerivativeLocationFunction, this, _1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, 1.0 ), 0.0, 0.0;
    }
        break;

    case l2:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = boost::make_shared< UnivariateProxy >(
                    boost::bind( &LibrationPoint::computeL2LocationFunction, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, boost::bind( &LibrationPoint::
                computeL2FirstDerivativeLocationFunction, this, _1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, 1.0 ), 0.0, 0.0;
    }
        break;

    case l3:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = boost::make_shared< UnivariateProxy >(
                    boost::bind( &LibrationPoint::computeL3LocationFunction, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, boost::bind( &LibrationPoint::
                computeL3FirstDerivativeLocationFunction, this, _1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, -1.0 ), 0.0, 0.0;
    }
        break;

    case l4:

        // Set position vector of L4 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter;
        positionOfLibrationPoint_.y( ) = 0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case l5:

        // Set position vector of L5 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter;
        positionOfLibrationPoint_.y( ) = -0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    default:

        throw std::runtime_error(
                            "The Lagrange libration point requested does not exist." );
    };
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat
