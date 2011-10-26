/*! \file seriesPropagator.h
 *    Header file that defines a class that executes a series of propagation steps.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 2
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
 *      110506    K. Kumar          File created.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

#ifndef SERIESPROPAGATOR_H
#define SERIESPROPAGATOR_H

// Include statements.
#include <cmath>
#include <iostream>
#include <map>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Propagators/propagator.h"
#include "Astrodynamics/Propagators/propagatorDataContainer.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Mathematics/basicMathematicsFunctions.h"

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
    SeriesPropagator( ) : numberOfPropagationSteps_( 0 ), seriesPropagationStart_( -0.0 ),
        seriesPropagationEnd_( -0.0 ), fixedOutputInterval_( -0.0 ), pointerToPropagator_( NULL ),
        bodiesToPropagate_( ) { }

    //! Set start of series propagation.
    /*!
     * Sets the start of the series propagation.
     * \param seriesPropagationStart Start of series propagation.
     */
    void setSeriesPropagationStart( const double& seriesPropagationStart )
    { seriesPropagationStart_ = seriesPropagationStart; }

    //! Set end of series propagation.
    /*!
     * Sets the end of the series propagation.
     * \param seriesPropagationEnd End of series propagation.
     */
    void setSeriesPropagationEnd( const double& seriesPropagationEnd )
    { seriesPropagationEnd_ = seriesPropagationEnd; }

    //! Set fixed output interval.
    /*!
     * Sets the fixed output interval at which propagation output should be
     * generated and stored for each body being propagated. This function
     * can only be called after both setSeriesPropagationStart() and
     * setSeriesPropagationEnd() are called. Also computes number of propagation
     * steps minus one. This will only lead to a sensible result if setSeriesPropagationStart()
     * and setSeriesPropagationEnd() have been called. To prevent numerical instabilities from
     * occuring (e.g., ceil( 4 / 2 ) != ceil( ( 4 + 1e-16 ) / 2 )), the square root of the machine
     * precision is subtracted before applying the ceil function.
     * \param fixedOutputInterval Fixed output interval.
     */
    void setFixedOutputInterval( const double& fixedOutputInterval )
    {
        fixedOutputInterval_ = fixedOutputInterval;
        numberOfPropagationSteps_ = static_cast< unsigned int >(
                    std::ceil( ( seriesPropagationEnd_ - seriesPropagationStart_ )
                               / fixedOutputInterval_
                               - std::sqrt( mathematics::MACHINE_PRECISION_DOUBLES ) ) );
    }

    //! Set initial state of body for series propagation.
    /*!
     * Sets the initial state of given body at start of series propagation.
     * \param pointerToBody Pointer to Body object.
     * \param pointerToInitialState Initial state at start of series
     *          propagation.
     */
    void setInitialState( Body* pointerToBody, State* pointerToInitialState )
    { pointerToPropagator_->setInitialState( pointerToBody, pointerToInitialState ); }

    //! Set propagator.
    /*!
     * Sets propagator used for series propagation.
     * \param pointerToPropagator Pointer to propagator.
     */
    void setPropagator( Propagator* pointerToPropagator )
    { pointerToPropagator_ = pointerToPropagator; }

    //! Get fixed output interval.
    /*!
     * Gets the fixed output interval at which propagation output should be
     * generated and stored in propagationHistory_.
     * \return Fixed output interval.
     */
    double& getFixedOutputInterval( ) { return fixedOutputInterval_; }

    //! Get start of series propagation.
    /*!
     * Returns the start of the series propagation.
     * \return Start of series propagation.
     */
    double& getSeriesPropagationStart( ) { return seriesPropagationStart_; }

    //! Get end of series propagation.
    /*!
     * Returns the end of the series propagation.
     * \return End of series propagation.
     */
    double& getSeriesPropagationEnd( ) { return seriesPropagationEnd_; }

    //! Get propagation history of body at fixed output intervals.
    /*!
     * Returns the propagation history of given body at specified fixed output
     * intervals.
     * \param pointerToBody Pointer to Body object.
     * \return Map of propagation history.
     */
    std::map< double, State >& getPropagationHistoryAtFixedOutputIntervals( Body* pointerToBody )
    { return pointerToPropagator_->getBodiesToPropagate( )[ pointerToBody ].propagationHistory; }

    //! Execute.
    /*!
     * Executes series propagation.
     */
    void execute( );

protected:

private:

    //! Type definition of map of body propagator data.
    /*!
     * Type definition of map of body propagator data.
     */
    typedef std::map< Body*, PropagatorDataContainer > BodyPropagatorDataMap;

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

    //! Map of bodies to propagate and associated data.
    /*!
     * Map of bodies to be propagated and associated data.
     */
    BodyPropagatorDataMap bodiesToPropagate_;

};

#endif // SERIESPROPAGATOR_H

// End of file.
