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

#ifndef TUDAT_CUSTOM_TORQUE_H
#define TUDAT_CUSTOM_TORQUE_H

#include <boost/function.hpp>
#include "tudat/astro/basic_astro/torqueModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Class to link the a custom torque model to a body.
/*!
 *  Class to link the a custom torque model to a body.
 */
class CustomTorqueModel : public TorqueModel
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param customTorqueFunction Function to retrieve the value of the custom torque.
     *  \param customUpdateFunction Function to update the value of the custom torque based on the current time (which is the
     *      only input).
     */
    CustomTorqueModel( const std::function< Eigen::Vector3d( const double ) >& customTorqueFunction ) :
        customTorqueFunction_( customTorqueFunction )
    { }

    //! Destructor
    ~CustomTorqueModel( ) { }

    //! Function to retrieve the current value of the torque.
    /*!
     *  Function to retrieve the current value of the torque, which is computed by the
     *  custom torque function, after having its members updated via the updateMembers function.
     *  \return Current torque, as computed by last call to updateMembers function.
     */
    Eigen::Vector3d getTorque( )
    {
        // Retrieve and return the torque based on the custom torque model
        return currentTorque_;
    }

    //! Update member variables used by the torque model.
    /*!
     *  Update member variables used by the torque model. This function simply calls the custom
     *  update function with the current time as input.
     */
    void updateMembers( const double currentTime )
    {
        // Update the custom torque model
        if ( currentTime_ != currentTime )
        {
            currentTorque_ = customTorqueFunction_( currentTime );
            currentTime_ = currentTime;
        }
    }

protected:

private:

    //! Function to be used to retrieve the torque.
    const std::function< Eigen::Vector3d( const double ) > customTorqueFunction_;

    Eigen::Vector3d currentTorque_;

};

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_CUSTOM_TORQUE_H
