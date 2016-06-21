/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120710    D. Dirkx          File created.
 *      120716    A. Ronse          Corrected spelling errors. Addition of file header,
 *                                  renamed virtual function. Added template functionality.
 *      120723    K. Kumar          Added template specialization for 1D cases.
 *      120818    A. Ronse          Made template specializations of updateAndGetAcceleration
 *                                  inline.
 *      120821    K. Kumar          Changed template parameters for AccelerationModel class to
 *                                  DataType only; removed obsolete template specializations of
 *                                  class and updateAndGetAcceleration() function.
 *      130225    K. Kumar          Modified updateMembers() function to pure virtual with void
 *                                  return-type.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_ACCELERATION_MODEL_H
#define TUDAT_ACCELERATION_MODEL_H

#include <map>
#include <unordered_map>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{  
namespace basic_astrodynamics
{

//! Base class for (translational) acceleration models.
/*!
 * Base class for (translational) acceleration models. Derived classes should contain
 * implementations to perform calculations of accelerations. Therefore, the getAcceleration()
 * function has no arguments.
 * \tparam AccelerationDataType Data type used to represent accelerations
 *          (default=Eigen::Vector3d).
 */
template < typename AccelerationDataType = Eigen::Vector3d >
class AccelerationModel
{
public:

    //! Constructor.
    AccelerationModel( ):
        currentTime_( TUDAT_NAN ){ }

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~AccelerationModel( ) { }

    //! Get acceleration.
    /*!
     * Returns the acceleration. No arguments are passed to this function for generality.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc which are to be set in a derived class and evaluated by the
     * updateMembers() function below.
     * \return Acceleration.
     * \sa updateMembers().
     */
    virtual AccelerationDataType getAcceleration( ) = 0;

    //! Update member variables used by the acceleration model.
    /*!
     * Updates member variables used by the acceleration model. In the case of acceleration models
     * containing varying parameters, function-pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function-pointers and updates member variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used by the getAcceleration() function.
     *
     * N.B.: This pure virtual function must be overridden by derived classes!
     * \param currentTime Time at which acceleration model is to be updated.
     */
    virtual void updateMembers( const double currentTime = TUDAT_NAN ) = 0;

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

protected:

    //! Previous time to which acceleration model was updated.
    double currentTime_;

protected:

private:
};

//! Typedef to a 3D acceleration model.
typedef AccelerationModel< > AccelerationModel3d;

//! Typedef for shared-pointer to a 3D acceleration model.
typedef boost::shared_ptr< AccelerationModel3d > AccelerationModel3dPointer;

//! Typedef to a 2D acceleration model.
typedef AccelerationModel< Eigen::Vector2d > AccelerationModel2d;

//! Typedef for shared-pointer to a 2D acceleration model.
typedef boost::shared_ptr< AccelerationModel2d > AccelerationModel2dPointer;

//! Update the members of an acceleration model and evaluate the acceleration.
/*!
 * Updates the member variables of an acceleration model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration.
 * \tparam AccelerationDataType Data type used to represent accelerations
 *          (default=Eigen::Vector3d).
 * \param accelerationModel Acceleration model that is to be evaluated.
 * \param currentTime Time at which acceleration model is to be updated.
 * \return Acceleration that is obtained following the member update.
 */
template < typename AccelerationDataType >
AccelerationDataType updateAndGetAcceleration(
        const boost::shared_ptr< AccelerationModel< AccelerationDataType > > accelerationModel,
        const double currentTime = TUDAT_NAN )
{
    // Update members.
    accelerationModel->updateMembers( currentTime );

    // Evaluate and return acceleration.
    return accelerationModel->getAcceleration( );
}

//! Typedef defining a list of accelerations acting on a single body, key is the name of each
//! body exerting a acceletation, value is a list of accelerations exerted by that body.
typedef std::unordered_map< std::string, std::vector<
boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >
SingleBodyAccelerationMap;

//! Typedef defining a list of accelerations acting on a set of bodies, key is the name of each
//! body undergoing a acceletation, value is SingleBodyAccelerationMap, defining all accelerations
//! acting on it.
typedef std::unordered_map< std::string, SingleBodyAccelerationMap > AccelerationMap;

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_ACCELERATION_MODEL_H
