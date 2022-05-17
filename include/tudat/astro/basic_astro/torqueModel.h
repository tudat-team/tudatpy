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

#ifndef TUDAT_TORQUEMODEL_H
#define TUDAT_TORQUEMODEL_H

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Pure virtual base class for a model that computes the torque acting on a body.
/*!
 *  Pure virtual base class for a model that computes the torque acting on a body. These models are used in the propagation
 *  of the rotational equations of motion. Every specific torque model requires its own derived class.
 */
class TorqueModel
{
public:

    //! Constructor
    TorqueModel( ):
        currentTime_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~TorqueModel( ) { }

    //! Function to retrieve the current value of the torque
    /*!
     * Returns the torque. No arguments are passed to this function for generality.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc which are to be set in a derived class and evaluated by the
     * updateMembers() function below.
     * \return Current torque, as computed by last call to updateMembers function.
     */
    virtual Eigen::Vector3d getTorque( ) = 0;

    //! Update member variables used by the torque model.
    /*!
     * Updates member variables used by the torque model. In the case of torque models
     * containing varying parameters, function-pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function-pointers and updates member variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used by the getTorque function.
     * This pure virtual function must be overridden by derived classes.
     * \param currentTime Time at which torque model is to be updated.
     */
    virtual void updateMembers( const double currentTime ) = 0;

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the torque model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

protected:

    //! Last time at which the updateMembers was evaluated (or NaN if next call is to recompute torque).
    double currentTime_;

private:

};

//! Function to compute the inertial torque term
/*!
 *  Function to compute the inertial torque term, which is included as the angular velocity rate is computed in a body-fixed frame.
 */
class InertialTorqueModel: public TorqueModel
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param angularVelocityFunction Function that returns body's body-fixed angular velocity vector
     * \param inertiaTensorFunction Function that returns body's inertia tensor
     */
    InertialTorqueModel(
            const std::function< Eigen::Vector3d( ) > angularVelocityFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction ):TorqueModel( ),
    angularVelocityFunction_( angularVelocityFunction ),
    inertiaTensorFunction_( inertiaTensorFunction ){ }

    //! Destructor
    ~InertialTorqueModel( ) { }


    //! Get inertial torque.
    /*!
     * Returns the inertial torque.
     * \return Inertial torque.
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the torque model.
    /*!
     * Updates member variables used by the torque model.
     * Function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    virtual void updateMembers( const double currentTime )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentTorque_ = -angularVelocityFunction_( ).cross( inertiaTensorFunction_( ) * angularVelocityFunction_( ) );
            currentTime_ = currentTime;
        }
    }

protected:

    //! Function that returns body's body-fixed angular velocity vector
    std::function< Eigen::Vector3d( ) > angularVelocityFunction_;

    //! Function that returns body's inertia tensor
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    //! Current torque, as computed by last call to updateMembers function
    Eigen::Vector3d currentTorque_;

private:

};

//! Typedef for list of torques acting on a body (map key is body exerting torque).
typedef std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > > > SingleBodyTorqueModelMap;

//! Typedef for list of torques acting on a set of bodies (map key is body undergoing torque).
typedef std::map< std::string, SingleBodyTorqueModelMap > TorqueModelMap;

//! Update the members of a torque model and evaluate the acceleration.
/*!
 * Updates the member variables of a torque model and subsequently evaluates the
 * acceleration. This allows the user to suffice with a single function call to both update the
 * members and evaluate the acceleration.
 * \param torqueModel Torque model that is to be evaluated.
 * \param currentTime Time at which torque model is to be updated.
 * \return Torque that is obtained following the member update.
 */
Eigen::Vector3d updateAndGetTorque(
        const std::shared_ptr< TorqueModel > torqueModel,
        const double currentTime = TUDAT_NAN );
}

}

#endif // TUDAT_TORQUEMODEL_H
