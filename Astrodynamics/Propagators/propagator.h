/*! \file propagator.h
 *    Header file that defines the base class for all propagators included in Tudat.
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
 *      It may be possible that there memory leaks in the library due to the
 *      Propagator class. This must be tested in future using a tool such as
 *      Valgrind ( www.valgrind.org ).
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
 *      100914    K. Kumar          File created.
 *      100926    K. Kumar          Added Doxygen comments.
 *      100928    K. Kumar          Completed missing comments.
 *      100929    J. Melman         Reference to Integrator::setFixedOutputInterval deleted,
 *                                  'class' in comments changed into 'object'.
 *      100929    K. Kumar          EigenRoutines.h replaced in include statements by
 *                                  linearAlgebra.h.
 *      110201    K. Kumar          Updated code to make use of State class; moved
 *                                  computeSumOfStateDerivatives_() to NumericalPropagator class.
 *      110207    K. Kumar          Added ostream overload.
 *      110216    K. Kumar          Added get function for fixedOutputInterval_.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

// Include statements.
#include <map>
#include <vector>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/ForceModels/forceModel.h"
#include "Astrodynamics/Propagators/propagatorDataContainer.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/state.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

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
    Propagator( ) : propagationIntervalStart_( -0.0 ), propagationIntervalEnd_( -0.0 ),
        bodiesToPropagate_( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Propagator( ) { }

    //! Set start of propagation interval.
    /*!
     * Sets the start of the propagation interval.
     * \param propagationIntervalStart Start of propagation interval.
     */
    void setPropagationIntervalStart( double propagationIntervalStart )
    { propagationIntervalStart_ = propagationIntervalStart; }

    //! Set end of propagation interval.
    /*!
     * Sets the end of the propagation interval.
     * \param propagationIntervalEnd End of propagation interval.
     */
    void setPropagationIntervalEnd( double propagationIntervalEnd )
    { propagationIntervalEnd_ = propagationIntervalEnd; }

    //! Add body to be propagated.
    /*!
     * Adds a body to be propagated.
     * \param pointerToBody Pointer to Body object.
     */
    void addBody( Body* pointerToBody )
    { bodiesToPropagate_[ pointerToBody ].pointerToBody = pointerToBody; }

    //! Set initial state of body.
    /*!
     * Sets the initial state of given body.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToInitialState Initial state given as pointer to a
     *          CartesianElements object.
     */
    virtual void setInitialState( Body* pointerToBody, State* pointerToInitialState )
    {
        bodiesToPropagate_[ pointerToBody ].pointerToInitialState = pointerToInitialState;
        bodiesToPropagate_[ pointerToBody ].sizeOfState = pointerToInitialState->state.size( );
    }

    //! Set bodies to propagate.
    /*!
     * Set map of bodies to propagate.
     * \param bodiesToProgate Map of bodies to propagate.
     */
    void setBodiesToPropagate( std::map< Body*, PropagatorDataContainer > bodiesToPropagate )
    { bodiesToPropagate_ = bodiesToPropagate; }

    //! Get start of propagation interval.
    /*!
     * Gets the start of the propagation interval.
     * \return Start of propagation interval.
     */
    double& getPropagationIntervalStart( ) { return propagationIntervalStart_; }

    //! Get end of propagation interval.
    /*!
     * Gets the end of the propagation interval.
     * \return End of propagation interval.
     */
    double& getPropagationIntervalEnd( ) { return propagationIntervalEnd_; }

    //! Get initial state of body.
    /*!
     * Returns the initial state of given body.
     * \param pointerToBody Pointer to Body object.
     * \return Initial state.
     */
    State* getInitialState( Body* pointerToBody )
    { return bodiesToPropagate_[ pointerToBody ].pointerToInitialState; }

    //! Get final state of body.
    /*!
     * Returns the final state of given body.
     * \param pointerToBody Pointer to Body object.
     * \return Final state.
     */
    State* getFinalState( Body* pointerToBody )
    { return &bodiesToPropagate_[ pointerToBody ].finalState; }

    //! Get map of bodies to propagate.
    /*!
     * Returns map of bodies to propagate.
     * \return Map of bodies to propagate.
     */
    std::map< Body*, PropagatorDataContainer>& getBodiesToPropagate( )
    { return bodiesToPropagate_; }

    //! Propagate.
    /*!
     * This function executes propagation.
     */
    virtual void propagate( ) = 0;

protected:

    //! Type definition of map of body propagator data.
    /*!
     * Type definition of map of body propagator data.
     */
    typedef std::map< Body*, PropagatorDataContainer > BodyPropagatorDataMap;

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

    //! Map of bodies to propagate and associated data.
    /*!
     * Map of bodies to be propagated and associated data.
     */
    BodyPropagatorDataMap bodiesToPropagate_;

private:
};

}

#endif // PROPAGATOR_H

// End of file.
