/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Propulsion/thrustFunctions.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
       currentThrust_( TUDAT_NAN ),
       bodyFixedThrustDirection_( bodyFixedThrustDirection.normalized( ) )
    { }

    //! Destructor.
    virtual ~EngineModel( ){ }

    //! Function to retrive the magnitude of the current engine thrust
    /*!
     *  Function to retrive the magnitude of the current engine thrust, which was set by the last call to the
     *  updateEngineModel function (to be defined in the derived class).
     *  \return Current engine thrust.
     */
    double getCurrentThrust( )
    {
        return currentThrust_;
    }

    //! Pure virtual function to update the engine model to the current time
    /*!
     *  Pure virtual function to update the engine model to the current time. This function must be implemented in teh derived
     *  class, wehre it typically resets the currentThrust_ variable (and any associated variables).
     *  \param currentTime Current itme in simulation.
     */
    virtual void updateEngineModel( const double currentTime ) = 0;

    //! Pure virtual function to retrieve the propellant mass rate.
    /*!
     *  Pure virtual function to retrieve the propellant mass rate.
     *  \return Propellant mass rate.
     */
    virtual double getCurrentMassRate( ) = 0;

    //! Function to retrieve the vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
    /*!
     * Function to retrieve the vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
     * \return Vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
     */
    Eigen::Vector3d getBodyFixedThrustDirection( )
    {
        return bodyFixedThrustDirection_;
    }

protected:

    //! Current magnitude of the thrust delivered by the engine.
    double currentThrust_;

    //! Vector denoting the direction of the thrust delivered by the engine in the body-fixed frame
    /*!
     * Vector denoting the direction of the thrust delivered by the engine in the body-fixed
     * frame (for an ideal engine, this is the opposite direction of the nozzle direction). By default, this vector is set
     * to the positive x-direction, along the longitudinal axis
     */
    Eigen::Vector3d bodyFixedThrustDirection_;

};

//! Engine model derived class in which the thrust is computed directly from teh propellant mass flow and specific impulse
/*!
 *  Engine model derived class in which the thrust is computed directly from teh propellant mass flow and specific impulse
 *  These two variables may be either constant or variable (magnitudes controlled by associated class defining e.g. GNC
 *  system.
 */
class DirectEngineModel: public EngineModel
{
public:

    //! Constructor
    /*!
     *  Constructor
     *  \param specificImpulseFunction Variable specific impulse of engine. (no input arguments provided; must be updated by
     *  associated guidance law).
     *  \param massFlowFunction Variable mass flow function. (no input arguments provided; must be updated by associated
     *  guidance law).
     *  \param bodyFixedThrustDirection Vector denoting the direction of the thrust delivered by the engine in the body-fixed
     *  frame (for an ideal engine, this is the opposite direction of the nozzle direction). By default, this vector is set
     *  to the positive x-direction, along the longitudinal axis.
     */
    DirectEngineModel(
            const boost::function< double( ) > specificImpulseFunction,
            const boost::function< double( ) > massFlowFunction,
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        EngineModel( bodyFixedThrustDirection ),
        specificImpulseFunction_( specificImpulseFunction ),
        massFlowFunction_( massFlowFunction )
    { }

    //! Function to update the engine model to the current time
    /*!
     *  Function to update the engine model to the current time. Retrieves the current mass flow and specific impulse;
     *  consequently computes the current thrust.
     *  \param currentTime Current itme in simulation.
     */
    void updateEngineModel( const double currentTime )
    {
        currentThrust_ = propulsion::computeThrustFromSpecificImpulse( massFlowFunction_( ), specificImpulseFunction_( ) );
    }

    //! Function to retrieve the propellant mass rate.
    /*!
     *  Function to retrieve the propellant mass rate.
     *  \return Propellant mass rate.
     */
    double getCurrentMassRate( )
    {
        return massFlowFunction_( );
    }


protected:

    //! Variable specific impulse of engine (no input arguments provided; must be updated by associated guidance law).
    boost::function< double( ) > specificImpulseFunction_;

    //! Variable mass flow function (no input arguments provided; must be updated by associated guidance law).
    boost::function< double( ) > massFlowFunction_;
};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_ENGINEMODEL_H
