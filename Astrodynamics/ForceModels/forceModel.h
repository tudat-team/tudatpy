/*! \file forceModel.h
 *    Header file that defines the base class for all force models included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 14 September, 2010
 *    Last modified     : 9 August, 2011
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
 *      100914    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor corrections to include statements and comments.
 *      110113    K. Kumar          Modified CelestialBody object to pointer; minor comment
 *                                  changes.
 *      110119    K. Kumar          Changed computeStateDerivatives() to computeForce().
 *      110202    K. Kumar          Updated code to make use of State class.
 *      110707    F.M. Engelen      Replaced code with new code.
 *      110809    K. Kumar          Split code into base class and derived class
 *                                  (SixDegreeOfFreedomForceModel).
 */

#ifndef FORCEMODEL_H
#define FORCEMODEL_H

// Include statements.
#include <Eigen/Core>
#include "Astrodynamics/States/state.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Force model class.
/*!
 * Base class for all force models.
 */
class ForceModel
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

    //! Set force.
    /*!
     * Sets the force.
     * \param force Force.
     */
    void setForce( const Eigen::Vector3d& force ) { force_ = force; }

   //! Get force.
    /*!
     * Returns the force.
     * \return Force.
     */
    Eigen::Vector3d& getForce( ) { return force_; }

    //! Compute force.
    /*!
     * Compute the force.
     * \param pointerToState Pointer to an object of the State class.
     */
    virtual void computeForce( State* pointerToState ) = 0;

protected:

    //! Force.
    /*!
     * Force given in [N].
     */
    Eigen::Vector3d force_;

private:
};

}

#endif // FORCEMODEL_H

// End of file.
