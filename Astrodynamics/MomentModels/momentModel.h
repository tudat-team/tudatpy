/*! \file momentModel.h
 *    This header file contains the moment model base class included in Tudat.
 *
 *    Path              : /Astrodynamics/MomentModels/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Checker           : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 10 May, 2011
 *    Last modified     : 24 August, 2011
 *
 *    References
 *
 *    Notes
 *      This baseclass already has the architecture to calculate the extra
 *      moment due to a force, in other words:
 *     total moment = arm x force + moment.
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
 *      110510    F.M. Engelen      First creation of code.
 *      110810    K. Kumar          Minor corrections; changed function names and removed redundant
 *                                  functions.
 *      110824    J. Leloux         Corrected doxygen documentation.
 */

#ifndef MOMENTMODEL_H
#define MOMENTMODEL_H

// Include statements.
#include <Eigen/Core>
#include <iostream>
#include "Astrodynamics/ForceModels/forceModel.h"

//! Tudat library namespace.
/*!
 *  The Tudat library namespace.
 */
namespace tudat
{

//! Base class for moment models.
/*!
 * Base class for the moment models used in Tudat.
 */
class MomentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    MomentModel( ) : pointerToForceModel_( NULL ) { forceApplicationArm_.setZero( ); }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~MomentModel( ) { }

    //! Set force application point.
    /*!
     * Sets force application arm, i.e., the vector from the origin of the reference frame in
     * which the moment is calculated, to the application point of the input force.
     * \param forceApplicationArm Vector arm to application point of force.
     */
    void setForceApplicationArm( Eigen::Vector3d& forceApplicationArm )
    { forceApplicationArm_ = forceApplicationArm; }

    //! Get force application arm.
    /*!
     * Returns force application arm, i.e., returns the vector from the origin of the reference
     * frame in which the moment is calculated, to the application point of the input force.
     * \return Vector arm to application point of force.
     */
    Eigen::Vector3d& getForceApplicationArm( ) { return forceApplicationArm_; }

    //! Set force model.
    /*!
     * Sets the force model. This force should be calculated in a reference frame that has the same
     * orientation as the frame the moment is calculated in (normally the body frame).
     * \param pointerToForceModel Pointer to force model.
     */
    void setForceModel( ForceModel* pointerToForceModel )
    { pointerToForceModel_ = pointerToForceModel; }

    //! Get force model.
    /*!
     * Returns force model. This force is calculated in a reference frame that has the same
     * orientation as the frame the moment is calculated in (normally the body frame).
     * \return Pointer to force model.
     */
    ForceModel* getForceModel( ) { return pointerToForceModel_; }

    //! Get moment.
    /*!
     * Returns the moment.
     * \return Moment.
     */
    Eigen::Vector3d& getMoment( ) { return moment_; }

    //! Set moment.
    /*!
     * Sets the moment.
     * \param moment Moment.
     */
    void setMoment( Eigen::Vector3d& moment ) { moment_ = moment; }

    //! Compute moment.
     /*!
      * Computes the moment.
      */
    virtual void computeMoment( State* pointerToState ) = 0;

protected:

    //! Pointer to force model.
    /*!
     * Pointer to force model.
     */
    ForceModel* pointerToForceModel_;

    //! Moment.
    /*!
     * Moment.
     */
    Eigen::Vector3d moment_;

    //! Force application point.
    /*!
     * Force application point, which can also be considered the force arm.
     */
    Eigen::Vector3d forceApplicationArm_;

private:
};

}

#endif // MOMENTMODEL_H

// End of file.
