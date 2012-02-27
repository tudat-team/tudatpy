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
 *      120204    D. Dirkx          Split MomentModel class into this class and
 *                                  MomentDueToForceModel.
 *
 *    References
 *
 */

#ifndef TUDAT_MOMENT_MODEL_H
#define TUDAT_MOMENT_MODEL_H

#include <Eigen/Core>
#include <iostream>
#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/generalizedForceModel.h"

namespace tudat
{

//! Base class for moment models.
/*!
 * Base class for the moment models used in Tudat.
 */
class PureMomentModel: public GeneralizedForceModel< Eigen::Vector3d, 3 >
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    PureMomentModel( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~PureMomentModel( ) { }

    //! Get moment.
    /*!
     * Returns the moment.
     * \return Moment.
     */
    Eigen::Vector3d getMoment( ) { return moment_; }

    virtual Eigen::Vector3d  getGeneralizedForce(  ) { return getMoment( ); }

    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    virtual void computeMoment( State* pointerToState, double time ) = 0;

protected:

    //! Moment.
    /*!
     * Moment.
     */
    Eigen::Vector3d moment_;

private:
};

} // namespace tudat

#endif // TUDAT_MOMENT_MODEL_H
