/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ENGINEMODEL_H
#define TUDAT_ENGINEMODEL_H

#include <Eigen/Core>
#include <functional>

#include "tudat/astro/propulsion/thrustFunctions.h"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace system_models
{

//! Base class for a model of an engine, the mass rate/thrust computations are not implemented here.
/*!
 *  Base class for a model of an engine, the mass rate/thrust computations are not implemented here, and should be
 *  performed in a derived class.
 */
class EngineModel
{
public:

    //! Constructor
    /*!
     *  Constructor
     *  \param bodyFixedThrustDirection Vector denoting the direction of the thrust delivered by the engine in the body-fixed
     *  frame (for an ideal engine, this is the opposite direction of the nozzle direction). By default, this vector is set
     *  to the positive x-direction, along the longitudinal axis.
     */
    EngineModel(
            const std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper,
            const std::string engineName,
            const std::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection =
                []( const double ){ return Eigen::Vector3d::UnitX( ); } ):
       thrustMagnitudeWrapper_( thrustMagnitudeWrapper ),
       engineName_( engineName ),
       bodyFixedThrustDirection_( bodyFixedThrustDirection )
    { }

    //! Destructor.
    virtual ~EngineModel( ){ }


    //! Pure virtual function to update the engine model to the current time
    /*!
     *  Pure virtual function to update the engine model to the current time. This function must be implemented in the derived
     *  class, wehre it typically resets the currentThrust_ variable (and any associated variables).
     *  \param currentTime Current itme in simulation.
     */
    virtual void updateEngineModel( const double currentTime )
    {
        thrustMagnitudeWrapper_->update( currentTime );
        currentBodyFixedThrustDirection_ = bodyFixedThrustDirection_( currentTime ).normalized( );
    }


    //! Function to retrive the magnitude of the current engine thrust
    /*!
     *  Function to retrive the magnitude of the current engine thrust, which was set by the last call to the
     *  updateEngineModel function (to be defined in the derived class).
     *  \return Current engine thrust.
     */
    double getCurrentThrust( const double currentMass = TUDAT_NAN )
    {
        return thrustMagnitudeWrapper_->getCurrentThrustForceMagnitude( currentMass );
    }


    double getCurrentThrustAcceleration( const double currentMass = TUDAT_NAN )
    {
        return thrustMagnitudeWrapper_->getCurrentThrustAccelerationMagnitude( currentMass );
    }

    //! Pure virtual function to retrieve the propellant mass rate.
    /*!
     *  Pure virtual function to retrieve the propellant mass rate.
     *  \return Propellant mass rate.
     */
    double getCurrentMassRate( const double currentMass = TUDAT_NAN )
    {
        return thrustMagnitudeWrapper_->getCurrentMassRate( currentMass );
    }

    //! Function to retrieve the vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
    /*!
     * Function to retrieve the vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
     * \return Vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
     */
    Eigen::Vector3d getBodyFixedThrustDirection( )
    {
        return currentBodyFixedThrustDirection_;
    }

    void resetCurrentTime( )
    {
        thrustMagnitudeWrapper_->resetCurrentTime( );
        bodyFixedThrustDirection_( TUDAT_NAN );
    }

    const std::string getEngineName( )
    {
        return engineName_;
    }

    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > getThrustMagnitudeWrapper( )
    {
        return thrustMagnitudeWrapper_;
    }


protected:

    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper_;

    const std::string engineName_;

    std::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection_;

    Eigen::Vector3d currentBodyFixedThrustDirection_;


};


} // namespace system_models

} // namespace tudat

#endif // TUDAT_ENGINEMODEL_H
