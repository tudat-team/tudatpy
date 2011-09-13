/*! \file deepSpaceManeuver.h
 *    This header file contains a base class for deep space maneuver elements.
 *
 *    Path              : /Astrodynamics/MissionSegments/DeepSpaceManeuver/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 6 April, 2011
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
 *      110224    E. Iorfida        First creation of code.
 *      110406    K. Kumar          Minor modifications.
 */

#ifndef DEEPSPACEMANEUVER_H
#define DEEPSPACEMANEUVER_H

// Include statements.
#include "state.h"

//! Deep space maneuver base class.
/*!
 * Deep space maneuver class.
 */
class DeepSpaceManeuver
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    DeepSpaceManeuver( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~DeepSpaceManeuver( );

    //! Set time of deep space maneuver event.
    /*!
     * Sets time of deep space maneuver event.
     * \param timeOfDeepSpaceManeuver Time of deep space maneuver event.
     */
    void setTime( const double& timeOfDeepSpaceManeuver );

    //! Set state at deep space maneuver event.
    /*!
     * Sets pointer to state at deep space maneuver event.
     * \param pointerToState Pointer to state at deep space maneuver event.
     */
    void setState( State* pointerToState );

    //! Set delta-V of deep space maneuver event.
    /*!
     * Sets delta-V of deep space maneuver event.
     * \param deltaV Delta-V of deep space maneuver event.
     */
    void setDeltaV( const double& deltaV );

    //! Get time of deep space maneuver event.
    /*!
     * Returns the the time of deep space maneuver event.
     * \return Time of deep space maneuver event.
     */
    double& getTime( );

    //! Get state at deep space maneuver event.
    /*!
     * Returns a pointer to state at deep space maneuver event.
     * \return Pointer to state at deep space maneuver event.
     */
    State* getState( );

    //! Get delta-V of deep space maneuver event.
    /*!
     * Returns delta-V of deep space maneuver event.
     * \return Delta-V of deep space maneuver event.
     */
    double& getDeltaV( );

protected:

private:

    //! Delta-V of deep space maneuver event.
    /*!
     * Delta-V of deep space maneuver event.
     */
    double deltaV_;

    //! Time of deep space maneuver event.
    /*!
     * Time of deep space maneuver event.
     */
    double timeOfDeepSpaceManeuver_;

    //! Pointer to state at deep space maneuver event.
    /*!
     * Pointer to state at deep space maneuver event.
     */
    State* pointerToState_;
};

#endif // DEEPSPACEMANEUVER_H

// End of file.
