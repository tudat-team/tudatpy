/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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

    //! Function to retrieve the current value of the toruqe
    /*!
     * Returns the toruqe. No arguments are passed to this function for generality.
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
    virtual void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

protected:

    //! Last time at which the updateMembers was evaluated (or NaN if next call is to recompute torque).
    double currentTime_;

private:

};

//! Typedef for list of torques acting on a body (map key is body exerting torque).
typedef std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::TorqueModel > > > SingleBodyTorqueModelMap;

//! Typedef for list of torques acting on a set of bodies (map key is body undergoing torque).
typedef std::map< std::string, SingleBodyTorqueModelMap > TorqueModelMap;

}

}

#endif // TUDAT_TORQUEMODEL_H
