/*! \file numericalPropagator.h
 *    Header file that defines the numerical propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 7
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
 *    Last modified     : 7 February, 2011
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
 *      100915    K. Kumar            File created.
 *      100926    K. Kumar            Doxygen comments added.
 *      100928    K. Kumar            Completed missing comments.
 *      100929    J. Melman           Added a dot behind "Include statements"
 *      100929    K. Kumar            Minor comment modifications.
 *      110201    K. Kumar            Updated code to use Integrator adaptor
 *                                    instead of pointers-to-member functions;
 *                                    Update code to use State class.
 */

#ifndef NUMERICALPROPAGATOR_H
#define NUMERICALPROPAGATOR_H

// Include statements.
#include "basicMathematicsFunctions.h"
#include "body.h"
#include "propagator.h"
#include "integrator.h"
#include "integratorAdaptor.h"

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
     * Sets an integrator for numerical propagation.
     * \param pointerToIntegrator Pointer to Integrator object.
     */
    void setIntegrator( Integrator* pointerToIntegrator );

    //! Add force model for propagation of a body.
    /*!
     * Adds a force model for propagation of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToForceModel Pointer to ForceModel object.
     */
    void addForceModel( Body* pointerToBody, ForceModel* pointerToForceModel );

    //! Propagate.
    /*!
     * Executes numerical propagation.
     */
    void propagate( );

protected:

private:

    //! Size of assembled state.
    /*!
     * Size of assmebled state.
     */
    unsigned int sizeOfAssembledState_;

    //! Pointer to Integrator object.
    /*!
     * Pointer to Integrator object.
     */
    Integrator* pointerToIntegrator_;

    //! Assembled state.
    /*!
     * Assembled state given as a State object.
     */
    State assembledState_;

    //! Assembled state derivative.
    /*!
     * Assembled state derivative given as a State object.
     */
    State assembledStateDerivative_;

    //! Pointer to adaptor object of IntegratorAdaptor class.
    /*!
     * Pointer to adaptor object of IntegratorAdaptor class. The template
     * parameter passed is this class.
     */
    IntegratorAdaptor < NumericalPropagator >
            integratorAdaptorForNumericalPropagator_;

    //! Iterator for vector container of pointers to force models.
    /*!
     * Iterator for vector container of pointers to force models.
     */
    std::vector < ForceModel* >::
            iterator iteratorContainerOfPointersToForceModels_;

    //! Map of integration history.
    /*!
     * Map of integration history.
     */
    std::map < double, State* > integrationHistory_;

    //! Compute sum of state derivatives.
    /*!
     * Computes the sum of state derivatives.
     * \param pointerToAssembledState Assembled state given as a pointer to
     *          State object.
     * \return Assembled state derivative given as a pointer to a State object.
     */
    State* computeSumOfStateDerivatives_( State* pointerToAssembledState );
};

#endif // NUMERICALPROPAGATOR_H

// End of file.
