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
 *      120710    D. Dirkx          Creation of code.
 *      120716    A. Ronse          Corrected spelling errors. Addition of file header,
 *                                  renamed virtual function. Added template functionality.
 *      120723    K. Kumar          Added template specialization for 1D cases.
 *      120818    A. Ronse          Made template specializations of updateAndGetAcceleration
 *                                  inline.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_ACCELERATION_MODEL_H
#define TUDAT_ACCELERATION_MODEL_H

#include <boost/shared_ptr.hpp>
#include <boost/exception/all.hpp>

#include <TudatCore/Mathematics/BasicMathematics/linearAlgebra.h>

namespace tudat
{  
namespace basic_astrodynamics
{
namespace acceleration_models
{

//! Base class for (translational) acceleration models.
/*!
 * Base class for (translational) acceleration models. Derived classes should contain pointers to
 * data strutures from where calculations of accelerations are to be performed. Therefore, the
 * getAcceleration() function has no arguments.
 * \tparam SpatialDimensions Size of VectorXd that is used for position, velocity and acceleration
 * i.e. number of spatial dimensions in simulation.
 * \tparam DataType Type used for the entries of the position, velocity and acceleration (i.e.
 * double, float, etc.)
 */
template < int SpatialDimensions = 3, typename DataType = double >
class AccelerationModel
{
public:

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
     * \return Acceleration as Eigen::Matrix.
     * \sa updateMembers().
     */
    virtual Eigen::Matrix< DataType, SpatialDimensions, 1 > getAcceleration( ) = 0;

    //! Update member variables used by the acceleration model.
    /*!
     * Updates member variables used by the acceleration model. In the case of acceleration models
     * containing varying parameters, function pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function pointers and sets member variables of the 'current'
     * values of these parameters. Only these current values, not the function pointers are then
     * used by the getAcceleration function.
     * \return True if the update was successful. NOTE: This could be modified to throw
     *          an exception in the future.
     */
    virtual bool updateMembers( ) { return true; }

protected:

private:
};

//! Base class for 1D (translational) acceleration models.
/*!
 * Base class for 1-dimensional (translational) acceleration models. This is a template
 * specialization of the general AccelerationModel base class.
 * \tparam DataType Type used for the entries of the position, velocity and acceleration (i.e.
 * double, float, etc.)
 */
template < typename DataType >
class AccelerationModel< 1, DataType >
{
public:

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~AccelerationModel( ) { }

    //! Get acceleration.
    /*!
     * Returns the 1D acceleration. No arguments are passed to this function for generality.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc which are to be set in a derived class and evaluated by the
     * updateMembers() function below.
     * \return Acceleration as scalar-type (double, float, int, etc.).
     * \sa updateMembers().
     */
    virtual DataType getAcceleration( ) = 0;

    //! Update member variables used by the acceleration model.
    /*!
     * Updates member variables used by the acceleration model. In the case of acceleration models
     * containing varying parameters, function pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function pointers and sets member variables of the 'current'
     * values of these parameters. Only these current values, not the function pointers are then
     * used by the getAcceleration function.
     * \return True if the update was successful. NOTE: This could be modified to throw
     *          an exception in the future.
     */
    virtual bool updateMembers( ) { return true; }

protected:

private:
};

//! Update the members of an acceleration model and evaluate the acceleration.
/*!
 * Updates the member variables of an acceleration model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration.
 * In case the updateMembers function fails (returns false), this function throws a runtime_error.
 * \tparam SpatialDimensions Size of VectorXd that is used for position, velocity and acceleration
 *          i.e. number of spatial dimensions in simulation.
 * \tparam DataType Type used for the entries of the position, velocity and acceleration (i.e.
 *          double, float, etc.)
 * \param accelerationModel Acceleration model that is to be evaluated.
 * \return Acceleration that is obtained following the member update.
 */
template < int SpatialDimensions, typename DataType >
inline Eigen::Matrix< DataType, SpatialDimensions, 1 > updateAndGetAcceleration(
        boost::shared_ptr< AccelerationModel< SpatialDimensions, DataType > > accelerationModel )
{
    // Update members and throw exception if it fails.
    if( !accelerationModel->updateMembers( ) )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "Unable to update acceleration model." ) ) );
    }

    // Evaluate and return acceleration.
    return accelerationModel->getAcceleration( );
}

//! Update the members of an acceleration model and evaluate the acceleration.
/*!
 * Updates the member variables of an acceleration model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration. This is an overload for 1-D int-types of the general
 * updateAndGetAcceleration() function (specialization of function templates is not permitted).
 * In case the updateMembers function fails (returns false), this function throws a runtime_error.
 * \param accelerationModel Acceleration model that is to be evaluated.
 * \return Acceleration that is obtained following the member update.
 */
inline int updateAndGetAcceleration(
        boost::shared_ptr< AccelerationModel< 1, int > > accelerationModel )
{
    // Update members and throw exception if it fails.
    if( !accelerationModel->updateMembers( ) )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "Unable to update acceleration model." ) ) );
    }

    // Evaluate and return acceleration.
    return accelerationModel->getAcceleration( );
}

//! Update the members of an acceleration model and evaluate the acceleration.
/*!
 * Updates the member variables of an acceleration model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration. This is an overload for 1-D double-types of the general
 * updateAndGetAcceleration() function (specialization of function templates is not permitted).
 * In case the updateMembers function fails (returns false), this function throws a runtime_error.
 * \param accelerationModel Acceleration model that is to be evaluated.
 * \return Acceleration that is obtained following the member update.
 */
inline double updateAndGetAcceleration(
        boost::shared_ptr< AccelerationModel< 1, double > > accelerationModel )
{
    // Update members and throw exception if it fails.
    if( !accelerationModel->updateMembers( ) )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "Unable to update acceleration model." ) ) );
    }

    // Evaluate and return acceleration.
    return accelerationModel->getAcceleration( );
}

//! Update the members of an acceleration model and evaluate the acceleration.
/*!
 * Updates the member variables of an acceleration model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration. This is an overload for 1-D float-types of the general
 * updateAndGetAcceleration() function (specialization of function templates is not permitted).
 * In case the updateMembers function fails (returns false), this function throws a runtime_error.
 * \param accelerationModel Acceleration model that is to be evaluated.
 * \return Acceleration that is obtained following the member update.
 */
inline float updateAndGetAcceleration(
        boost::shared_ptr< AccelerationModel< 1, float > > accelerationModel )
{
    // Update members and throw exception if it fails.
    if( !accelerationModel->updateMembers( ) )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "Unable to update acceleration model." ) ) );
    }

    // Evaluate and return acceleration.
    return accelerationModel->getAcceleration( );
}

} // namespace acceleration_models
} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_ACCELERATION_MODEL_H
