/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ACCELERATION_MODEL_H
#define TUDAT_ACCELERATION_MODEL_H

#include <vector>
#include <map>
#include <unordered_map>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"

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
template< typename AccelerationDataType = Eigen::Vector3d >
class AccelerationModel
{
public:

    //! Constructor.
    AccelerationModel( ):
        currentTime_( TUDAT_NAN ),
        currentAcceleration_( Eigen::Vector3d::Constant( TUDAT_NAN ) ){ }

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~AccelerationModel( ) { }

//    //! Get acceleration.
//    /*!
//     * Returns the acceleration. No arguments are passed to this function for generality.
//     * Instead, all data required for computation is to be obtained from pointers to functions/
//     * classes/structs, etc which are to be set in a derived class and evaluated by the
//     * updateMembers() function below.
//     * \return Acceleration.
//     * \sa updateMembers().
//     */
//    virtual AccelerationDataType getAcceleration( ) = 0;

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

    AccelerationDataType& getAccelerationReference( )
    {
        return currentAcceleration_;
    }

    AccelerationDataType getAcceleration( )
    {
        return currentAcceleration_;
    }

    void getAccelerationByReference( AccelerationDataType& acceleration ) const
    {
        acceleration = currentAcceleration_;
    }

    void addCurrentAcceleration( AccelerationDataType& acceleration ) const
    {
        acceleration += currentAcceleration_;
    }

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

protected:

    //! Previous time to which acceleration model was updated.
    double currentTime_;

    AccelerationDataType currentAcceleration_;

protected:

private:

};

//! Typedef to a 3D acceleration model.
typedef AccelerationModel< > AccelerationModel3d;

//! Typedef for shared-pointer to a 3D acceleration model.
typedef std::shared_ptr< AccelerationModel3d > AccelerationModel3dPointer;


//! Typedef to a 2D acceleration model.
typedef AccelerationModel< Eigen::Vector2d > AccelerationModel2d;

//! Typedef for shared-pointer to a 2D acceleration model.
typedef std::shared_ptr< AccelerationModel2d > AccelerationModel2dPointer;

extern template class AccelerationModel< Eigen::Vector3d >;


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
        const std::shared_ptr< AccelerationModel< AccelerationDataType > > accelerationModel,
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
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >
SingleBodyAccelerationMap;


//! Typedef defining a list of accelerations acting on a set of bodies, key is the name of each
//! body undergoing an acceletation, value is SingleBodyAccelerationMap, defining all accelerations
//! acting on it.
typedef std::unordered_map< std::string, SingleBodyAccelerationMap > AccelerationMap;

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_ACCELERATION_MODEL_H
