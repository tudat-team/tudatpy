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
 *      110527    L. van der Ham    File created.
 *      110602    L. van der Ham    Made LibrationPoints class and added comments.
 *      110625    K. Kumar          Minor modifications to layout and comments.
 *      110629    L. van der Ham    Modifications according to comments first code check.
 *      110710    K. Kumar          Removed duplicated code; added missing include statement;
 *                                  modified libration point compute-functions; changed filename
 *                                  and class; added enum.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110927    L. van der Ham    Reverted to full equations of motion for determination location
 *                                  of colinear libration points.
 *      111027    K. Kumar          Moved 1-line functions from source file.
 *      120307    K. Kumar          Moved file.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120813    P. Musegaas       Changed code to new root finding structure. Added option to
 *                                  specify which rootfinder and termination conditions to use.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References:
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

#ifndef TUDAT_LIBRATION_POINT_H
#define TUDAT_LIBRATION_POINT_H

#include <cmath>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute mass parameter.
/*!
 * Computes dimensionless mass parameter for the CRTBP.
 * \param primaryGravitationalParameter Gravitational parameter of primary body.
 * \param secondaryGravitationalParameter Gravitational parameter of secondary body.
 * \param Dimensionless mass parameter.
 */
inline double computeMassParameter( const double primaryGravitationalParameter,
                                    const double secondaryGravitationalParameter )
{
    return secondaryGravitationalParameter
            / ( primaryGravitationalParameter + secondaryGravitationalParameter );
}

//! Libration point class.
/*!
 * This class includes functions to compute the location of Lagrange libration points in the CRTBP.
 */
class LibrationPoint
{

public:

    //! Lagrange libration points.
    enum LagrangeLibrationPoints { L1, L2, L3, L4, L5 };

    //! Default constructor.
    LibrationPoint( const double aPrimaryGravitationalParameter,
                    const double aSecondaryGravitationalParameter,
                    const root_finders::RootFinderPointer aRootFinder )
        : massParameter( computeMassParameter( aPrimaryGravitationalParameter,
                                               aSecondaryGravitationalParameter ) ),
          rootFinder( aRootFinder )
    { }

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LibrationPoint( const double aMassParameter,
                    const root_finders::RootFinderPointer aRootFinder )
        : massParameter( aMassParameter ),
          rootFinder( aRootFinder )
    { }

    //! Get dimensionless mass parameter.
    /*!
     * Returns the dimensionless mass parameter based on the gravitational
     * parameters of the primary and secondary bodies.
     */
    double getMassParameter( ) { return massParameter; }

    //! Get location of Lagrange libration point.
    /*!
     * Returns the position vector in Cartesian elements of a Lagrange libration point. The
     * libration point is selected when the computeLocationOfLibrationPoint( ) function is
     * called.
     * \return Cartesian position elements of Lagrange libration point.
     */
    Eigen::Vector3d getLocationOfLagrangeLibrationPoint( )
    {
        return positionOfLibrationPoint_;
    }

    //! Compute location of Lagrange libration point.
    /*!
     * Computes the location as a position vector in Cartesian elements of a given Lagrange
     * libration point.
     * \param lagrangeLibrationPoint Lagrange libration point.
     */
    void computeLocationOfLibrationPoint( LagrangeLibrationPoints lagrangeLibrationPoint );

protected:

private:

    //! Dimensionless mass parameter.
    /*!
    * Dimensionless mass parameter of the smaller of the massive bodies in the CRTBP.
    */
    const double massParameter;

    //! Shared pointer to the rootfinder.
    /*!
     * Shared pointer to the rootfinder. The rootfinder contains termination conditions inside.
     */
    const root_finders::RootFinderPointer rootFinder;

    //! Cartesian position elements of a Lagrange libration point.
    /*!
     * Cartesian position elements of a Lagrange libration point.
     */
    Eigen::Vector3d positionOfLibrationPoint_;

    //! Compute equation of motion for location of L1 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of location function.
     */
    double computeL1LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate -
                ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                + massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L1 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL1FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                + 2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L2 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of location function.
     */
    double computeL2LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate
                - ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                - massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L2 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL2FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                - 2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L3 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of location function.
     */
    double computeL3LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate
                + ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                - massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L3 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL3FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 - 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                -2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }
};

// Typedef for shared-pointer to LibrationPoint object.
typedef boost::shared_ptr< LibrationPoint > LibrationPointPointer;

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat

#endif // TUDAT_LIBRATION_POINT_H
