/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120204    D. Dirkx          Split MomentModel class into PureMomentModel and  this class.
 *
 *    References
 *
 */

#ifndef TUDAT_MOMENT_DUE_TO_FORCEMODEL_H
#define TUDAT_MOMENT_DUE_TO_FORCEMODEL_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/generalizedForceModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"

namespace tudat
{

//! Base class for all generalized forces.
/*!
 *  Class for moment due to an resultant force. Interface for computing moment takes arm and force,
 *  so that only the cross product is calculated. An additional interface with a ForceModel for
 *  instance can be added in the future.
 */
class MomentDueToForceModel : public GeneralizedForceModel< Eigen::Vector3d, 3 >
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    MomentDueToForceModel( ):momentDueToForce_( Eigen::Vector3d::Zero( ) ) {}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~MomentDueToForceModel( ) { }

    //! Get moment.
    /*!
     * Returns the moment.
     * \return Moment.
     */
    Eigen::Vector3d getMomentDueToForce( ) { return momentDueToForce_; }

    //! Get moment.
    /*!
     * Returns the moment.
     * \return Moment.
     */
    virtual Eigen::Vector3d getGeneralizedForce(  ) { return getMomentDueToForce( ); }

     //! Compute moment.
     /*!
      * Computes the moment.
      */
    virtual void computeMoment( Eigen::Vector3d force, Eigen::Vector3d arm )
    {
            momentDueToForce_ = arm.cross( force );
    }

protected:

private:

    //! Moment.
    /*!
     * Moment.
     */
    Eigen::Vector3d momentDueToForce_;
};

} // namespace tudat

#endif // TUDAT_MOMENT_DUE_TO_FORCEMODEL_H
