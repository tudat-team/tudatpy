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
 *      110510    F.M. Engelen      Creation of code.
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
