/*! \file cartesianStateNumericalPropagator.h
 *    Header file that defines the Cartesian state numerical propagator class
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 14 May, 2011
 *    Last modified     : 14 May, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD    Author              Comment
 *      110514    K. Kumar            File created.
 */

#ifndef CARTESIANSTATENUMERICALPROPAGATOR_H
#define CARTESIANSTATENUMERICALPROPAGATOR_H

// Include statements.
#include "numericalPropagator.h"

//! Cartesian state numerical propagator.
/*!
 * Definition of numerical propagator that propagates Cartesian state
 * consisting of three position components (x, y, z) and three velocity
 * components (xdot, ydot, zdot).
 */
class CartesianStateNumericalPropagator : public NumericalPropagator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CartesianStateNumericalPropagator( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CartesianStateNumericalPropagator( );

    //! Propagate.
    /*!
     * Executes numerical propagation.
     */
    void propagate( );

    //! Set initial state of body.
    /*!
     * Sets the initial state of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToInitialState Initial state given as pointer to a
     *          CartesianElements object.
     */
    void setInitialState( Body* pointerToBody,
                          State* pointerToInitialState );

    //! Compute state derivative.
    /*!
     * Computes the state derivative given to the integrator. This is an
     * implementation of the computeStateDerivative() virtual function in the
     * StateDerivativeBase class.
     * \param indepedendentVariable Value of independent variable.
     * \param pointerToAssembledState Assembled state given as a pointer to
     *          State object.
     * \param pointerToAssembledStateDerivative Computed assembled state
     *          derivative given as a pointer to a State object.
     */
    void computeStateDerivative( double& independentVariable,
                                 State* pointerToAssembledState,
                                 State* pointerToAssembledStateDerivative );

protected:

private:

    //! Size of assembled state.
    /*!
     * Size of assembled state.
     */
    unsigned int sizeOfAssembledState_;

    //! Assembled state in Cartesian elements.
    /*!
     * Assembled state in Cartesian elements.
     */
    CartesianElements assembledState_;
};

#endif // CARTESIANSTATENUMERICALPROPAGATOR_H

// End of file.
