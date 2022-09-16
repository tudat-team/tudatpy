#ifndef TUDAT_THRUSTDIRECTIONGUIDANCE_H
#define TUDAT_THRUSTDIRECTIONGUIDANCE_H

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *
 *    Kluever (2010), Low-Thrust Trajectory Optimization Using Orbital Averaging and Control Parameterization, In: Conway,
 *    (editor) Spacecraft trajectory optimization. Cambridge University Press, 2010.
 *    Boudestijn (2014), Development of a Low-Thrust Earth-Centered Transfer Optimizer for the Preliminary Mission Design Phase,
 *    M.Sc. Thesis, Delft University of Technology
 */

#include <cmath>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/missionGeometry.h"
#include "tudat/astro/propulsion/thrustGuidance.h"

namespace tudat
{

namespace propulsion
{

////! Derived class which calculates the thrust direction based on costates for modified equinoctial elements
///*!
// *  Derived class which calculates the thrust direction based on costates for modified equinoctial elements. The values of the
// *  costates are provided as input to the class (as a function of time), and are used to determine the thrust direction,
// *  based on the algorithm by Kluever (2010), extended to Modified Equinoctial Elements by Boudestijn (2014)
// */
//class MeeCostateBasedThrustGuidance: public BodyFixedForceDirectionGuidance
//{
//public:


//    //! Constructor
//    /*!
//     * Constructor
//     * \param thrustingBodyStateFunction Function returning the state of the body under thrust as a function of time
//     * \param centralBodyStateFunction Function returning the state of the central body as a function of time
//     * \param centralBodyGravitationalParameterFunction Function returning the gravitational parameter of the central body as a
//     * function of time
//     * \param costateFunction Function which returns the costates as a function of time.
//     * \param bodyFixedForceDirection Function returning the direction of the force in the body-fixed frame.
//     */
//    MeeCostateBasedThrustGuidance(
//            const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
//            const std::function< Eigen::Vector6d( ) > centralBodyStateFunction,
//            const std::function< double( ) > centralBodyGravitationalParameterFunction,
//            const std::function< Eigen::VectorXd( const double ) > costateFunction,
//            const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
//            [ ]( ){ return  Eigen::Vector3d::UnitX( ); } );

//    //! Destructor
//    ~MeeCostateBasedThrustGuidance( ){ }

//    //! Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
//    /*!
//     *  Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection,
//     *  computed based on the MEE costates provided by the costateFunction_ member variable
//     *  \return Current force direction, expressed in propagation frame.
//     */
//    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
//    {
//        return currentForceDirection_;
//    }

//    //! Function to get the rotation from body-fixed to inertial frame.
//    /*!
//     *  Function to get the rotation from body-fixed to inertial frame. NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
//     *  \return NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
//     */
//    Eigen::Quaterniond getRotationToGlobalFrame( )
//    {
//        throw std::runtime_error( "Error, body-fixed frame to propagation frame not yet implemented for MeeCostateBasedThrustGuidance." );
//    }

//private:

//    //! Function to update the force direction to the current time.
//    /*!
//     *  Function to update the force direction to the current time.
//     *  \param time Time to which object is to be updated.
//     */
//    void updateForceDirection( const double time );


//    //!  Function returning the state of the body under thrust as a function of time
//    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction_;

//    //!  Function returning the state of the central body as a function of time
//    std::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

//    //!  Function returning the gravitational parameter of the central body as a function of time
//    std::function< double( ) > centralBodyGravitationalParameterFunction_;

//    //! General function which gives the costates vector as a function of time.
//    std::function< Eigen::VectorXd( const double ) > costateFunction_;

//    //! Direction of thrust force in propagation frame, as computed by last call to updateForceDirection function.
//    Eigen::Vector3d currentForceDirection_;

//};

} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTDIRECTIONGUIDANCE_H
