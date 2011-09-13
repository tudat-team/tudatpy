/*! \file seriesPropagator.h
 *    Header file that defines a class that executes a series of propagation
 *    steps.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 1
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
 *    Date created      : 6 May, 2011
 *    Last modified     : 6 May, 2011
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
 *      110506    K. Kumar          File created.
 */

#ifndef SERIESPROPAGATOR_H
#define SERIESPROPAGATOR_H

// Include statements.
#include <cmath>
#include <map>
#include "body.h"
#include "cartesianElements.h"
#include "propagator.h"
#include "propagatorDataContainer.h"

// Using declarations.
using std::map;

// Forward declarations.
class Propagator;
class PropagatorDataContainer;

//! Definition of series propagator.
/*!
 * Definition of series propagator, which loops calls to a propagator to
 * generate a series of propagation steps. For example, if the independent
 * variable is time, this class loops calls to a given propagator to generate
 * timeseries data.
 */
class SeriesPropagator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SeriesPropagator( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SeriesPropagator( );

    //! Set start of series propagation.
    /*!
     * Sets the start of the series propagation.
     * \param seriesPropagationStart Start of series propagation.
     */
    void setSeriesPropagationStart( const double& seriesPropagationStart );

    //! Set end of series propagation.
    /*!
     * Sets the end of the series propagation.
     * \param seriesPropagationEnd End of series propagation.
     */
    void setSeriesPropagationEnd( const double& seriesPropagationEnd );

    //! Set fixed output interval.
    /*!
     * Sets the fixed output interval at which propagation output should be
     * generated and stored for each body being propagated. This function
     * can only be called after both setSeriesPropagationStart() and
     * setSeriesPropagationEnd() are called.
     * \param fixedOutputInterval Fixed output interval.
     */
    void setFixedOutputInterval( const double& fixedOutputInterval );

    //! Set initial state of body for series propagation.
    /*!
     * Sets the initial state of given body at start of series propagation.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToInitialState Initial state at start of series
     *          propagation.
     */
    void setInitialState( Body* pointerToBody,
                          State* pointerToInitialState );

    //! Set propagator.
    /*!
     * Sets propagator used for series propagation.
     * \param pointerToPropagator Pointer to propagator.
     */
    void setPropagator( Propagator* pointerToPropagator );

    //! Get fixed output interval.
    /*!
     * Gets the fixed output interval at which propagation output should be
     * generated and stored in propagationHistory_.
     * \return Fixed output interval.
     */
    double& getFixedOutputInterval( );

    //! Get start of series propagation.
    /*!
     * Returns the start of the series propagation.
     * \return Start of series propagation.
     */
    double& getSeriesPropagationStart( );

    //! Get end of series propagation.
    /*!
     * Returns the end of the series propagation.
     * \return End of series propagation.
     */
    double& getSeriesPropagationEnd( );

    //! Get propagation history of body at fixed output intervals.
    /*!
     * Returns the propagation history of given body at specified fixed output
     * intervals.
     * \param pointerToBody Pointer to Body object.
     * \return Map of propagation history.
     */
    map< double, State >&
            getPropagationHistoryAtFixedOutputIntervals( Body* pointerToBody );

    //! Execute.
    /*!
     * Executes series propagation.
     */
    void execute( );

protected:

private:

    //! Number of propagation steps.
    /*!
     * Number of propagation steps. This value is computed internally
     * using the start and end of the series propagation and the fixed output
     * interval.
     */
    unsigned int numberOfPropagationSteps_;

    //! Start of series propagation.
    /*!
     * Start of series propagation.
     */
    double seriesPropagationStart_;

    //! End of series propagation.
    /*!
     * End of series propagation.
     */
    double seriesPropagationEnd_;

    //! Fixed output interval.
    /*!
     * Fixed interval for output of state.
     */
    double fixedOutputInterval_;

    //! Pointer to propagator.
    /*!
     * Pointer to propagator used for series propagation.
     */
    Propagator* pointerToPropagator_;

    //! Map of propagated bodies and associated data.
    /*!
     * Map of propagated bodies and associated data.
     */
    map< Body*, PropagatorDataContainer > propagatedBodies_;

    //! Iterator for map of propagated bodies and associated data.
    /*!
     * Iterator for map of propagated bodies and associated data.
     */
    map< Body*, PropagatorDataContainer >::iterator iteratorPropagatedBodies_;
};

#endif // SERIESPROPAGATOR_H

// End of file.
