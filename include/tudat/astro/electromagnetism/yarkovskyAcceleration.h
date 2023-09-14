/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *          Pérez-Hernández, J. A., & Benet, L. (2022). Non-zero Yarkovsky acceleration for near-Earth
 *          asteroid (99942) Apophis. Communications Earth & Environment, 3(1), Article 1.
 *          DOI: https://doi.org/10.1038/s43247-021-00337-x
 */






#ifndef TUDAT_YARKOVSKYACCELERATION_H
#define TUDAT_YARKOVSKYACCELERATION_H
#include <functional>
#include <boost/lambda/lambda.hpp>
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/basicTypedefs.h"


namespace tudat
{
namespace electromagnetism
{

//! Compute Yarkovsky Acceleration using a simplified tangential model.
/*!
 * \param yarkovskyParameter Yarkovsky Parameter N2                                          [m/s^2]
 * \param stateVector is the state vector pointing from the source to the body
 *          undergoing the acceleration                                                         [m]
 * \return Acceleration due to Yarkovsky effect.                                            [m/s^2]
 */
Eigen::Vector3d computeYarkovskyAcceleration( double yarkovskyParameter, const Eigen::Vector6d& stateVector );

//! Class for calculating an Yarkovsky acceleration, based on (Pérez-Hernández & Benet, 2022).
/*!
 * Class for calculating an Yarkovsky acceleration, based on (Pérez-Hernández & Benet, 2022).
 * The acceleration is only considered in the tangential direction and is proportional to
 * a = A2 * (r0/rS)^2, where A2 is the Yarkovsky parameter, r0 = 1AU and rS is the heliocentric
 * distance in AU.
 */
class YarkovskyAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param yarkovskyParameter Yarkovsky parameter
     * \param bodyGravitationalParameterFunction Function that returns the state of the body.
     * \param centralBodyStateFunction Functon that returns the state of central body.
     */
    YarkovskyAcceleration( const double yarkovskyParameter,
                           const std::function< Eigen::Vector6d( ) >& bodyStateFunction,
                           const std::function< Eigen::Vector6d( ) >& centralBodyStateFunction = []( ) { return Eigen::Vector6d::Zero( ); } )
            : yarkovskyParameter_( yarkovskyParameter ), bodyStateFunction_( bodyStateFunction ),
              centralBodyStateFunction_( centralBodyStateFunction )
    {
    }

    //! Destructor
    ~YarkovskyAcceleration( ) override = default;

    //! Function to update constituent elements of Yarkovsky acceleration to current time
    /*!
     * Function to update constituent elements of Yarkovsky acceleration to current time
     * \param currentTime Time to which Yarkovsky acceleration elements are to be updated.
     */
    void updateMembers( const double currentTime ) override
    {
        if ( this->currentTime_ != currentTime ) {
            // Calculate current relative state of accelerated body
            currentState_ = bodyStateFunction_( ) - centralBodyStateFunction_( );

            // Update
            this->currentAcceleration_ = computeYarkovskyAcceleration( yarkovskyParameter_, currentState_ );
            this->currentTime_ = currentTime;
        }
    }

    //! Function to retrieve current state of the body that is undergoing the Yarkovsky acceleration, relative to central body
    /*!
     *  Function to retrieve Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body,
     *  in propagation frame.
     *  \return Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body, in global frame.
     */
    Eigen::Vector6d getCurrentState( )
    {
        return currentState_;
    }

private:
    //! Yarkovsky Parameter
    double yarkovskyParameter_;

    //! State function of the body that is undergoing the Yarkovsky acceleration.
    std::function< Eigen::Vector6d( ) > bodyStateFunction_;

    //! State function of the body that is being orbited
    std::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //! Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body, in global frame.
    Eigen::Vector6d currentState_;
};

}
}
#endif // TUDAT_YARKOVSKYACCELERATION_H
