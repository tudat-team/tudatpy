/*! \file numericalPropagator.h
 *    Header file that defines the numerical propagator class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 8
 *    Check status      : Checked
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
 *    Last modified     : 20 September, 2011
 *
 *    References
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
 *      100915    K. Kumar          File created.
 *      100926    K. Kumar          Doxygen comments added.
 *      100928    K. Kumar          Completed missing comments.
 *      100929    J. Melman         Added a dot behind "Include statements".
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to use Integrator adaptor instead of
 *                                  pointers-to-member functions; updated code to use State class.
 *      110512    K. Kumar          Updated code not to use dynamics memory allocation; split into
 *                                  base and derived class.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

#ifndef NUMERICALPROPAGATOR_H
#define NUMERICALPROPAGATOR_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Propagators/propagator.h"
#include "Mathematics/NumericalIntegrators/integrator.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Numerical propagator class.
/*!
 * Numerical propagator class.
 */
class NumericalPropagator : public Propagator, public StateDerivativeBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NumericalPropagator( ) : sizeOfAssembledState_( 0 ), pointerToIntegrator_( NULL ) { }

    //! Set integrator for propagation.
    /*!
     * Sets an integrator for numerical propagation.
     * \param pointerToIntegrator Pointer to Integrator object.
     */
    void setIntegrator( Integrator* pointerToIntegrator )
    { pointerToIntegrator_ = pointerToIntegrator; }

    //! Add force model for propagation of a body.
    /*!
     * Adds a force model for propagation of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToForceModel Pointer to ForceModel object.
     */
    void addForceModel( Body* pointerToBody, ForceModel* pointerToForceModel )
    { bodiesToPropagate_[ pointerToBody ].vectorContainerOfPointersToForceModels
                .push_back( pointerToForceModel ); }

    //! Set initial state of body.
    /*!
     * Sets the initial state of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToInitialState Initial state given as pointer to a
     *          CartesianElements object.
     */
    void setInitialState( Body* pointerToBody, State* pointerToInitialState );

    //! Propagate.
    /*!
     * Executes numerical propagation.
     */
    void propagate( );

protected:

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

    //! Pointer to Integrator object.
    /*!
     * Pointer to Integrator object.
     */
    Integrator* pointerToIntegrator_;

private:
};

}

#endif // NUMERICALPROPAGATOR_H

// End of file.
