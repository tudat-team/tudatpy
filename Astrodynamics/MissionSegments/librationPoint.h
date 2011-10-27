/*! \file librationPoint.h
 *    Header file that defines the computation of the location of a Lagrange libration point in
 *    the Circular Restricted Three-Body Problem (CRTBP).
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 27 May, 2011
 *    Last modified     : 27 October, 2011
 *
 *    References:
 *      van der Ham, L., TBD
 *      Mireles James, J.D., "Celestial Mechanics Notes Set 4: The Circular Restricted Three Body
 *          Problem", http://www.math.utexas.edu/users/jjames/celestMech, 2006.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110527    L. van der Ham    First creation of code.
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
 */

#ifndef LIBRATIONPOINT_H
#define LIBRATIONPOINT_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/States/cartesianPositionElements.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Libration point class.
/*!
 * This class includes functions to compute the location of a Lagrange
 * libration point in the CRTBP.
 */
class LibrationPoint
{
public:

    //! Lagrange libration points.
    /*!
     * Lagrange libration points.
     */
    enum LagrangeLibrationPoints { L1, L2, L3, L4, L5 };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LibrationPoint( ) : massParameter_( -0.0 ), massParameterSquared_( -0.0 ),
        oneMinusMassParameterSquared_( -0.0 ), pointerToNewtonRaphson_( NULL ),
        pointerToPrimaryCelestialBody_( NULL ), pointerToSecondaryCelestialBody_( NULL ) { }

    //! Set primary celestial body.
    /*!
     * Sets primary celestial body in the CRTBP.
     * \param pointerToPrimaryCelestialBody Pointer to primary celestial body.
     */
    void setPrimaryCelestialBody( CelestialBody* pointerToPrimaryCelestialBody )
    { pointerToPrimaryCelestialBody_ = pointerToPrimaryCelestialBody; }

    //! Set secondary celestial body.
    /*!
     * Sets secondary celestial body in the CRTBP.
     * \param pointerToSecondaryCelestialBody Pointer to secondary celestial
     *          body.
     */
    void setSecondaryCelestialBody( CelestialBody* pointerToSecondaryCelestialBody )
    { pointerToSecondaryCelestialBody_ = pointerToSecondaryCelestialBody; }

    //! Set Newton-Raphson method for Lagrange libration points algorithm.
    /*!
     * Sets Newton-Raphson method used to compute the location of the Lagrange
     * libration points.
     * \param pointerToNewtonRaphson Pointer to Newton-Raphson method.
     */
    void setNewtonRaphsonMethod( NewtonRaphson *pointerToNewtonRaphson )
    { pointerToNewtonRaphson_ = pointerToNewtonRaphson; }

    //! Set mass parameter.
    /*!
     * Sets mass parameter for the CRTBP.
     * \param massParameter Mass parameter.
     */
    void setMassParameter( const double& massParameter ) { massParameter_ = massParameter; }

    //! Get dimensionless mass parameter.
    /*!
     * Returns the dimensionless mass parameter based on the gravitational
     * parameters of the primary and secondary bodies.
     */
    const double& getMassParameter( ) { return massParameter_; }

    //! Get location of Lagrange libration point.
    /*!
     * Returns the position vector in Cartesian elements of a Lagrange libration point. The
     * libration point is selected when the computeLocationOfLibrationPoint() function is
     * called.
     * \return Cartesian position elements of Lagrange libration point.
     */
    CartesianPositionElements& getLocationOfLagrangeLibrationPoint( )
    { return positionOfLibrationPoint_; }

    //! Compute mass parameter.
    /*!
     * Computes dimensionless mass parameter. This requires that the setPrimaryCelestialBody() and
     * setSecondaryCelestialBody() functions have been called.
     */
    void computeMassParameter( )
    {
        massParameter_ = pointerToSecondaryCelestialBody_->getGravitationalParameter( )
                / ( pointerToPrimaryCelestialBody_->getGravitationalParameter( )
                    + pointerToSecondaryCelestialBody_->getGravitationalParameter( ) );
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
    NewtonRaphson* pointerToNewtonRaphson_;

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
    CartesianPositionElements positionOfLibrationPoint_;

    //! Pointer to CelestialBody class for primary body.
    /*!
     * Pointer to CelestialBody class for primary body in CRTBP.
     */
    CelestialBody* pointerToPrimaryCelestialBody_;

    //! Pointer to CelestialBody class for secondary body.
    /*!
     * Pointer to CelestialBody class for secondary body in CRTBP.
     */
    CelestialBody* pointerToSecondaryCelestialBody_;

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

}

#endif // LIBRATIONPOINT_H

// End of file.
