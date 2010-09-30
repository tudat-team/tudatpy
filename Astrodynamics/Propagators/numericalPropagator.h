/*! \file numericalPropagator.h
 *    Header file that defines the numerical propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagator/
 *    Version           : 5
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
 *    Date created      : 14 September, 2010
 *    Last modified     : 29 September, 2010
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
 *      YYMMDD    author              comment
 *      100915    K. Kumar            File created.
 *      100926    K. Kumar            Doxygen comments added.
 *      100928    K. Kumar            Completed missing comments.
 *      100929    J. Melman           Added a dot behind "Include statements" :)
 *      100929    K. Kumar            Minor comment modifications.
 */

#ifndef NUMERICALPROPAGATOR_H
#define NUMERICALPROPAGATOR_H

// Include statements.
#include "propagator.h"
#include "integrator.h"
#include "bodyContainer.h"
#include "body.h"
#include "basicMathematicsFunctions.h"

// Forward declarations.
class Integrator;
class BodyContainer;

//! Numerical propagator class.
/*!
 * Numerical propagator class.
 */
class NumericalPropagator : public Propagator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NumericalPropagator( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~NumericalPropagator( );

    //! Set integrator for propagation.
    /*!
     * This function sets an integrator for propagation.
     * \param pointerToIntegrator Pointer to Integrator object.
     */
    void setIntegrator( Integrator* pointerToIntegrator );

    //! Propagate.
    /*!
     * This function executes propagation.
     */
    void propagate( );

protected:

private:

    //! Pointer to Integrator object.
    /*!
     * Pointer to Integrator object.
     */
    Integrator* pointerToIntegrator_;
};

#endif // NUMERICALPROPAGATOR_H

// End of file.
