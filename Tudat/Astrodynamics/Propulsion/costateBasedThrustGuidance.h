#ifndef TUDAT_THRUSTDIRECTIONGUIDANCE_H
#define TUDAT_THRUSTDIRECTIONGUIDANCE_H

/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <cmath>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"

namespace tudat
{

namespace propulsion
{

//! Derived class which calculates the thrust direction in the body frame
//! based on modified equinoctial elements.
class MeeCostateBasedThrustGuidance: public BodyFixedForceDirectionGuidance
{
public:


    //! Constructor with a function as input
    /*!
     * \brief MeeCostateBasedThrustGuidance control parameterisation in terms of modifed equinoctial
     * elements, costates are inputted as a boost function.
     * \param bodyMap Obtained in the simulation setup, contains information about the satellite and the celestial
     * bodies.
     * \param vehicleName Name of satellite.
     * \param costateFunction boost function which returns the costates.
     * \param thrustMagnitude Magnitude of Thrust force.
     */
    MeeCostateBasedThrustGuidance(
            const boost::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
            const boost::function< Eigen::Vector6d( ) > centralBodyStateFunction,
            const boost::function< double( ) > centralBodyGravitationalParameterFunction,
            const boost::function< Eigen::VectorXd( const double ) > costateFunction,
            const boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) );

    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
    {
        return currentForceDirection_;
    }

    //! Function to get the rotation from body-fixed to inertial frame.
    /*!
     *  Function to get the rotation from body-fixed to inertial frame. NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
     *  \return NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
     */
    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        throw std::runtime_error( "Error, body-fixed frame to propagation frame not yet implemented for DirectionBasedForceGuidance." );
    }

private:

    void updateForceDirection( const double time );

    //! interpolator for the costates, independent parameter is time
    const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > costatesInterpolator_;

    //! General function which gives the costates vector.
    boost::function< Eigen::VectorXd( const double ) > costateFunction_;

    boost::function< Eigen::Vector6d( ) > thrustingBodyStateFunction_;

    boost::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    boost::function< double( ) > centralBodyGravitationalParameterFunction_;


    Eigen::Vector3d currentForceDirection_;

};

} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTDIRECTIONGUIDANCE_H
