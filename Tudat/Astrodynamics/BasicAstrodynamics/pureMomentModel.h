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
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#ifndef TUDAT_MOMENT_MODEL_H
#define TUDAT_MOMENT_MODEL_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/generalizedForceModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace moment_models
{

//! Base class for moment models.
/*!
 * Base class for the moment models used in Tudat.
 */
class PureMomentModel: public GeneralizedForceModel< Eigen::Vector3d, 3 >
{
public:

    //! Typedef for shared pointer to state.
    /*!
     * Typedef for shared pointer to state.
     */
    typedef boost::shared_ptr< states::State > StatePointer;

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

    //! Get generalized force, which is a moment.
    /*!
     * Returns the generalized force, which in this case is a moment.
     * \return Moment.
     */
    virtual Eigen::Vector3d  getGeneralizedForce(  ) { return getMoment( ); }

    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param state Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    virtual void computeMoment( StatePointer state, const double time ) = 0;

protected:

    //! Moment.
    /*!
     * Moment.
     */
    Eigen::Vector3d moment_;

private:
};

} // namespace moment_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_MOMENT_MODEL_H
