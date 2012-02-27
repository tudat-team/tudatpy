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
 *      100914    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor corrections to include statements and comments.
 *      110113    K. Kumar          Modified CelestialBody object to pointer; minor comment
 *                                  changes.
 *      110119    K. Kumar          Changed computeStateDerivatives( ) to computeForce( ).
 *      110202    K. Kumar          Updated code to make use of State class.
 *      110707    F.M. Engelen      Replaced code with new code.
 *      110809    K. Kumar          Split code into base class and derived class
 *                                  (SixDegreeOfFreedomForceModel).
 *
 *    References
 *
 */

#ifndef TUDAT_FORCE_MODEL_H
#define TUDAT_FORCE_MODEL_H

#include <Eigen/Core>
#include "Tudat/Astrodynamics/States/state.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/generalizedForceModel.h"

namespace tudat
{

//! Force model class.
/*!
 * Base class for all force models.
 */
class ForceModel : public GeneralizedForceModel< Eigen::Vector3d, 3 >
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ForceModel( ) : force_( Eigen::Vector3d::Zero( ) ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~ForceModel( ) { }

   //! Get force.
    /*!
     * Returns the force.
     * \return Force.
     */
    Eigen::Vector3d getGeneralizedForce( ) { return getGeneralizedForce( ); }

    //! Get force.
     /*!
      * Returns the force.
      * \return Force.
      */
     Eigen::Vector3d getForce( ) { return force_; }

    //! Compute force.
    /*!
     * Compute the force.
     * \param pointerToState Pointer to an object of the State class.
     * \param time Time (or other independent variable).
     */
    virtual void computeForce( State* pointerToState, double time ) = 0;

protected:

    //! Force.
    /*!
     * Force given in [N].
     */
    Eigen::Vector3d force_;

private:

};

} // namespace tudat

#endif // TUDAT_FORCE_MODEL_H
