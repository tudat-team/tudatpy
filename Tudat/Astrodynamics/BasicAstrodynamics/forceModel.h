/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#ifndef TUDAT_FORCE_MODEL_H
#define TUDAT_FORCE_MODEL_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/States/state.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/generalizedForceModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace force_models
{

//! Force model class.
/*!
 * Base class for all force models.
 */
class ForceModel : public GeneralizedForceModel< Eigen::Vector3d, 3 >
{
public:

    //! Typedef for shared pointer to state.
    /*!
     * Typedef for shared pointer to state.
     */
    typedef boost::shared_ptr< states::State > StatePointer;

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
     * \param state Pointer to an object of the State class.
     * \param time Time (or other independent variable).
     */
    virtual void computeForce( StatePointer state, const double time ) = 0;

protected:

    //! Force.
    /*!
     * Force given in [N].
     */
    Eigen::Vector3d force_;

private:
};

} // namespace force_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_FORCE_MODEL_H
