/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CENTRALGRAVITYACCELERATIONPARTIALS_H
#define TUDAT_CENTRALGRAVITYACCELERATIONPARTIALS_H

#include <iostream>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

using namespace gravitation;

//! Calculates partial derivative of point mass gravitational acceleration wrt the position of body undergoing acceleration.
/*!
 *  Calculates partial derivative of point mass gravitational acceleration wrt the position of body undergoing acceleration.
 *  \param acceleratedBodyPosition Cartesian state of body being accelerated.
 *  \param acceleratingBodyPosition Cartesian state of body exerting acceleration.
 *  \param gravitationalParameter Gravitational parameter of gravitating body.
 *  \return Matrix with the Jacobian of the acceleration vector w.r.t. the position vector.
 */
Eigen::Matrix3d calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
        const Eigen::Vector3d& acceleratedBodyPosition,
        const Eigen::Vector3d& acceleratingBodyPosition,
        const double gravitationalParameter );

//! Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
/*!
 *  Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
 *  \param acceleratedBodyPosition Cartesian state of body being accelerated.
 *  \param acceleratingBodyPosition Cartesian state of body exerting acceleration.
 *  \return Vector with the partial of the acceleration vector w.r.t. ational parameter of the central body.
 */
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d& acceleratedBodyPosition,
                                                                         const Eigen::Vector3d& acceleratingBodyPosition);


//! Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
/*!
 *  Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
 *  \param gravitationalAcceleration Gravitational acceleration vector for which partial is to be computed.
 *  \param gravitationalParameter Gravitational parameter of gravitating body.
 *  \return Vector with the partial of the acceleration vector w.r.t. ational parameter of the central body.
 */
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d& gravitationalAcceleration,
                                                                         const double gravitationalParameter );


//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class CentralGravitationPartial: public AccelerationPartial
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    CentralGravitationPartial(
            const boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > gravitationalAcceleration,
            const std::string acceleratedBody,
            const std::string acceleratingBody );

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtPositionOfAcceleratedBody( )
    {
        return currentPositionPartial;
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtVelocityOfAcceleratedBody( )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtPositionOfAcceleratingBody( )
    {
        return -currentPositionPartial;
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    Eigen::Matrix3d wrtVelocityOfAcceleratingBody( )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    virtual std::pair< boost::function< Eigen::MatrixXd( ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< Eigen::MatrixXd( ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for updating partial to current state.
    /*!
     *  Function for updating partial to current state. For the central gravitational acceleration, position partial is
     *  computed and set.
     *  \param currentTime Time at which partials are to be calcu
     */
    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentPositionPartial = calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        acceleratedBodyState_( ),
                        centralBodyState_( ),
                        gravitationalParameterFunction_( ) );
            currentTime_ = currentTime;
        }
    }

protected:

    //! Function to create a function returning the current partial w.r.t. a gravitational parameter.
    /*!
     * Function to create a function returning the current partial w.r.t. a gravitational parameter.
     * \param parameterId Identified of parameter for which the partial is to be created.
     * \return Pair with partial function and paramater partial size. The partial function is non-empty only
     * if the parameterId input represents the gravitational parameter of acceleratingBody_ (or acceleratedBody_ if
     * accelerationUsesMutualAttraction_ is true).
     */
    std::pair< boost::function< Eigen::MatrixXd( ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function to calculate central gravity partial w.r.t. central body gravitational parameter.
    Eigen::Vector3d wrtGravitationalParameterOfCentralBody( );

    //! Function to retrieve current gravitational parameter of central body.
    boost::function< double( ) > gravitationalParameterFunction_;

    //! Function to retrieve current state of body exerting acceleration.
    boost::function< Eigen::Vector3d( ) > centralBodyState_;

    //! Function to retrieve current state of body undergoing acceleration.
    boost::function< Eigen::Vector3d( ) > acceleratedBodyState_;

    //! Boolean denoting whether the gravitational attraction of the central body on the accelerated body is included.
    bool accelerationUsesMutualAttraction_;


    //! Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPositionPartial;

};

}

}

}
#endif // TUDAT_CENTRALGRAVITYACCELERATIONPARTIALS_H
