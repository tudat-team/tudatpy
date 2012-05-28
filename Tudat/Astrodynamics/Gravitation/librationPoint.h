/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110527    L. van der Ham    Creation of code.
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
 *
 *    References:
 *      Mireles James, J.D. Celestial Mechanics Notes Set 4: The Circular Restricted Three Body
 *          Problem, 2006, http://www.math.utexas.edu/users/jjames/hw4Notes.pdf,
 *          last accessed: 26 May, 2012.
 *      van der Ham, L., TBD
 */

#ifndef TUDAT_LIBRATION_POINT_H
#define TUDAT_LIBRATION_POINT_H

#include <cmath>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Libration point class.
/*!
 * This class includes functions to compute the location of a Lagrange
 * libration point in the CRTBP.
 */
class LibrationPoint
{
public:

    //! Typedef for shared pointer to celestial body.
    /*!
     * Typedef for shared pointer to celestial body.
     */
    typedef boost::shared_ptr< bodies::CelestialBody > CelestialBodyPointer;

    //! Typedef for shared pointer to newton-raphson method.
    /*!
     * Typedef for shared pointer to newton-raphson method.
     */
    typedef boost::shared_ptr< NewtonRaphson > NewtonRaphsonPointer;

    //! Lagrange libration points.
    /*!
     * Lagrange libration points.
     */
    enum LagrangeLibrationPoints { L1, L2, L3, L4, L5 };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LibrationPoint( )
        : massParameter_( TUDAT_NAN ),
          massParameterSquared_( TUDAT_NAN ),
          oneMinusMassParameterSquared_( TUDAT_NAN )
    { }

    //! Set primary celestial body.
    /*!
     * Sets primary celestial body in the CRTBP.
     * \param primaryCelestialBody Pointer to primary celestial body.
     */
    void setPrimaryCelestialBody( CelestialBodyPointer primaryCelestialBody )
    {
        primaryCelestialBody_ = primaryCelestialBody;
    }

    //! Set secondary celestial body.
    /*!
     * Sets secondary celestial body in the CRTBP.
     * \param secondaryCelestialBody Pointer to secondary celestial body.
     */
    void setSecondaryCelestialBody( CelestialBodyPointer secondaryCelestialBody )
    {
        secondaryCelestialBody_ = secondaryCelestialBody;
    }

    //! Set Newton-Raphson method for Lagrange libration points algorithm.
    /*!
     * Sets Newton-Raphson method used to compute the location of the Lagrange
     * libration points.
     * \param newtonRaphson Pointer to Newton-Raphson method.
     */
    void setNewtonRaphsonMethod( NewtonRaphsonPointer newtonRaphson )
    {
        newtonRaphson_ = newtonRaphson;
    }

    //! Set mass parameter.
    /*!
     * Sets mass parameter for the CRTBP.
     * \param massParameter Mass parameter.
     */
    void setMassParameter( const double massParameter ) { massParameter_ = massParameter; }

    //! Get dimensionless mass parameter.
    /*!
     * Returns the dimensionless mass parameter based on the gravitational
     * parameters of the primary and secondary bodies.
     */
    double getMassParameter( ) { return massParameter_; }

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

    //! Compute mass parameter.
    /*!
     * Computes dimensionless mass parameter. This requires that the setPrimaryCelestialBody( ) and
     * setSecondaryCelestialBody( ) functions have been called.
     */
    void computeMassParameter( )
    {
        massParameter_ = secondaryCelestialBody_->getGravityFieldModel( )->
                getGravitationalParameter( ) / ( primaryCelestialBody_->
                getGravityFieldModel( )->getGravitationalParameter( )
                    + secondaryCelestialBody_->getGravityFieldModel( )->
                                                 getGravitationalParameter( ) );
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
    double massParameter_;

    //! Dimensionless mass parameter squared.
    /*!
     * Dimensionless mass parameter squared. This is used locally in the get-functions for the
     * Lagrange libration points.
     */
    double massParameterSquared_;

    //! Square of one minus dimensionless mass parameter.
    /*!
     * Square of one minus dimensionless mass parameter \f$( 1 -  \mu )^{2}\f$. This is used
     * locally in the get-functions for the Lagrange libration points.
     */
    double oneMinusMassParameterSquared_;

    //! Pointer to object of NewtonRaphson class.
    /*!
     * Pointer to object of NewtonRaphson class.
     */
    NewtonRaphsonPointer newtonRaphson_;

    //! Pointer to adaptor object of NewtonRaphsonAdaptor class.
    /*!
     * Pointer to adaptor object of NewtonRaphsonAdaptor class. The template
     * parameter passed is this class.
     */
    NewtonRaphsonAdaptor< LibrationPoint > newtonRaphsonAdaptorForLibrationPoint_;

    //! Cartesian position elements of a Lagrange libration point.
    /*!
     * Cartesian position elements of a Lagrange libration point.
     */
    Eigen::Vector3d positionOfLibrationPoint_;

    //! Pointer to CelestialBody class for primary body.
    /*!
     * Pointer to CelestialBody class for primary body in CRTBP.
     */
    CelestialBodyPointer primaryCelestialBody_;

    //! Pointer to CelestialBody class for secondary body.
    /*!
     * Pointer to CelestialBody class for secondary body in CRTBP.
     */
    CelestialBodyPointer secondaryCelestialBody_;

    //! Compute equation of motion for location of L1 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of location function.
     */
    double computeL1LocationFunction_( double& xLocationEstimate )
    {
        return xLocationEstimate -
                ( 1.0 - massParameter_ ) / pow( massParameter_ + xLocationEstimate, 2.0 )
                + massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L1 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL1FirstDerivativeLocationFunction_( double& xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter_ )
                / pow( massParameter_ + xLocationEstimate, 3.0 )
                + 2.0 * massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L2 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of location function.
     */
    double computeL2LocationFunction_( double& xLocationEstimate )
    {
        return xLocationEstimate
                - ( 1.0 - massParameter_ ) / pow( massParameter_ + xLocationEstimate, 2.0 )
                - massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L2 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL2FirstDerivativeLocationFunction_( double& xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter_ )
                / pow( massParameter_ + xLocationEstimate, 3.0 )
                - 2.0 * massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L3 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of location function.
     */
    double computeL3LocationFunction_( double& xLocationEstimate )
    {
        return xLocationEstimate
                + ( 1.0 - massParameter_ ) / pow( massParameter_ + xLocationEstimate, 2.0 )
                - massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L3 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, TBD).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL3FirstDerivativeLocationFunction_( double& xLocationEstimate )
    {
        return 1.0 - 2.0 * ( 1.0 - massParameter_ )
                / pow( massParameter_ + xLocationEstimate, 3.0 )
                -2.0 * massParameter_ / pow( 1.0 - massParameter_ - xLocationEstimate, 3.0 );
    }
};

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_LIBRATION_POINT_H
