/*! \file propagator.h
 *    Header file that defines the base class for all propagators included in
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
 *      100914    K. Kumar            File created.
 *      100926    K. Kumar            Added Doxygen comments.
 *      100928    K. Kumar            Completed missing comments.
 *      100929    J. Melman           Reference to Integrator::
 *                                    setFixedOutputInterval deleted,
 *                                    'class' in comments changed into 'object'.
 *      100929    K. Kumar            EigenRoutines.h replaced in include
 *                                    statements by linearAlgebra.h.
 */

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

// Include statements.
#include <map>
#include <vector>
#include "forceModel.h"
#include "bodyContainer.h"
#include "body.h"
#include "linearAlgebra.h"

// Forward declarations.
class BodyContainer;

//! Propagator class.
/*!
 * Base class for all propagators included in Tudat.
 */
class Propagator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    Propagator( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Propagator( );

    //! Set start of propagation interval.
    /*!
     * This function sets the start of the propagation interval.
     * \param propagationIntervalStart Start of propagation interval.
     */
    void setPropagationIntervalStart( const double& propagationIntervalStart );

    //! Set end of propagation interval.
    /*!
     * This function sets the end of the propagation interval.
     * \param propagationIntervalEnd End of propagation interval.
     */
    void setPropagationIntervalEnd( const double& propagationIntervalEnd );

    //! Add body to be propagated.
    /*!
     * This function adds a body to be propagated.
     * \param pointerToBody Pointer to Body object.
     */
    void addBody( Body* pointerToBody );

    //! Add force model for propagation of a body.
    /*!
     * This function adds a force model for propagation of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToForceModel Pointer to ForceModel object.
     */
    void addForceModel( Body* pointerToBody, ForceModel* pointerToForceModel );

    //! Set propagator for propagation of a body.
    /*!
     * This function sets a propagator for propagation of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToPropagator Pointer to Propagator object.
     */
    void setPropagator( Body* pointerToBody, Propagator* pointerToPropagator );

    //! Set initial state of body.
    /*!
     * This function sets the initial state of given body.
     * \param pointerToBody Pointer to Body object.
     * \param initialStateVector Initial state vector.
     */
    void setInitialState( Body* pointerToBody, VectorXd& initialStateVector );

    //! Set fixed output interval.
    /*!
     * This function sets the fixed output interval at which propagation output
     * should be generated and stored in propagationHistory_. Calls to this
     * function are optional.
     * \param fixedOutputInterval Fixed output interval.
     */
    void setFixedOutputInterval( const double& fixedOutputInterval );

    //! Get final state of body.
    /*!
     * This function gets the final state of given body.
     * \param pointerToBody Pointer to Body object.
     * \return Final state vector.
     */
    VectorXd& getFinalState( Body* pointerToBody );

    //! Get propagation history of body at fixed output intervals.
    /*!
     * This function gets the propagation history of given body at specified
     * fixed output intervals.
     * \param pointerToBody Pointer to Body object.
     * \return Map of propagation history.
     */
    std::map < double, VectorXd >
            getPropagationHistoryAtFixedOutputIntervals( Body* pointerToBody );

    //! Propagate.
    /*!
     * This function executes propagation.
     */
    virtual void propagate( ) = 0;

protected:

    //! Size of assembled state vector.
    /*!
     * Size of assmebled state vector.
     */
    int sizeOfAssembledStateVector_;

    //! Set start of propagation interval.
    /*!
     * Set start of propagation interval.
     */
    double propagationIntervalStart_;

    //! Set end of propagation interval.
    /*!
     * Set end of propagation interval.
     */
    double propagationIntervalEnd_;

    //! Fixed output interval.
    /*!
     * Fixed interval for output of state.
     */
    double fixedOutputInterval_;

    //! Assembled state vector.
    /*!
     * Assembled state vector.
     */
    VectorXd assembledStateVector_;

    //! Assembled state derivative vector.
    /*!
     * Assembled state derivative vector.
     */
    VectorXd assembledStateDerivativeVector_;

    //! Map of bodies to be propagated.
    /*!
     * Map of bodies to be propagated.
     */
    std::map < Body*, BodyContainer* > bodiesToBePropagated_;

    //! Iterator for map of bodies to be propagated.
    /*!
     * Iterator for map of bodies to be propagated.
     */
    std::map < Body*, BodyContainer* >::iterator
            iteratorBodiesToBePropagated_;

    //! Iterator for vector container of pointers to force models.
    /*!
     * Iterator for vector container of pointers to force models.
     */
    std::vector < ForceModel* >
            ::iterator iteratorContainerOfPointersToForceModels_;

    //! Map of propagation history.
    /*!
     * Map of propagation history.
     */
    std::map < double, VectorXd > propagationHistory_;

    //! Iterator for map of propagation history.
    /*!
     * Iterator for map of propagation history.
     */
    std::map < double, VectorXd >::iterator iteratorPropagationHistory_;

    //! Map of integration history.
    /*!
     * Map of integration history. // JM: Niet nodig voor analytical propagator
     */
    std::map < double, VectorXd > integrationHistory_;

    //! Pointer to Body object.
    /*!
     * Pointer to Body object.
     */
    Body* pointerToBody_;

    //! Compute sum of state derivatives.
    /*!
     * This function computes the sum of state derivatives.
     * \param stateVector State vector.
     * \param stateDerivativeVector State derivative vector.
     */
    void computeSumOfStateDerivatives_
            ( VectorXd& stateVector, VectorXd& stateDerivativeVector );

private:
};

#endif // PROPAGATOR_H

// End of file.
