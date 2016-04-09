/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ACCELERATIONPARTIALS_H
#define TUDAT_ACCELERATIONPARTIALS_H

#include <string>
#include <map>
#include <Eigen/Core>

#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


//! Base class for objects calculating partial derivatives of accelerations w.r.t. states, model parameters.
/*!
 *  Base class for objects calculating partial derivatives of accelerations  w.r.t. states, model parameters. Such
 *  calculations are used in orbit determination, for the computation of the state transition; sensitivity matrices.
 *  Derived classes implement derivative-calculating models for specific acceleration models, so that the calculation
 *  of all partials of a single type acceleration model is encompassed in a single derived class.
 */
class AccelerationPartial: public StateDerivativePartial
{

public:
    //! Base class constructor.
    /*!
     *  Constructor of base class, sets the base class member variables identifying the body undergoing and exerting the
     *  acceleration.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     *  \param accelerationType Type of acceleration w.r.t. which partial is taken.
     */
    AccelerationPartial( const std::string& acceleratedBody, const std::string& acceleratingBody,
                         const basic_astrodynamics::AvailableAcceleration accelerationType ):
        StateDerivativePartial( propagators::transational_state, std::make_pair( acceleratedBody, "" ) ),
        acceleratedBody_( acceleratedBody ), acceleratingBody_( acceleratingBody ),accelerationType_( accelerationType ) { }

    //! Virtual destructor.
    virtual ~AccelerationPartial( ) { }

    //! Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
    /*!
     * Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for translational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return Pair with function, returning partial derivative, and number of columns in partial vector,
     */
    std::pair< boost::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >  getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        // Initialize to empty function; 0 parameter size.
        std::pair< boost::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
                partialFunction = std::make_pair( boost::function< void( Eigen::Block< Eigen::MatrixXd > ) >( ), 0 );

        // Check if state is translational.
        if( integratedStateType == propagators::transational_state )
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for dynamics" );
            }
            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == acceleratedBody_ )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAcceleratedBody, this, _1 ), 3 );
            }
            else if( stateReferencePoint.first == acceleratingBody_ )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAcceleratingBody, this, _1 ), 3 );
            }
            else if( isAccelerationPartialWrtAdditionalBodyNonNull( stateReferencePoint.first ) )
            {
                partialFunction = std::make_pair( boost::bind( &AccelerationPartial::wrtStateOfAdditionalBody,
                                                               this, _1, stateReferencePoint.first ), 3 );
            }
        }


        return partialFunction;
    }

    //! Function to check whether a partial w.r.t. some integrated state is non-zero.
    /*!
     * Function to check whether a partial w.r.t. some integrated state is non-zero.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for translational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return True if dependency exists, false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        bool isDependent = 0;

        // Check if state is translational.
        if( integratedStateType == propagators::transational_state )
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when checking state derivative partial dependency of acceleration model, cannot have reference point on body for dynamics" );
            }

            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == acceleratedBody_ || stateReferencePoint.first == acceleratingBody_ ||
                     isAccelerationPartialWrtAdditionalBodyNonNull( stateReferencePoint.first ) )
            {
                isDependent = 1;
            }
        }
        return isDependent;
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    virtual Eigen::Matrix3d wrtPositionOfAcceleratedBody( ) = 0;

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body.
     *  \return Partial derivative of acceleration w.r.t. velocity of body undergoing acceleration.
     */
    virtual Eigen::Matrix3d wrtVelocityOfAcceleratedBody( ) = 0;

    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body undergoing acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body
     *  undergoing acceleration.
     *  \return Partial derivative of acceleration w.r.t. state of body exerting acceleration.
     */
    void wrtStateOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        partialMatrix.block( 0, 0, 3, 3 ) += wrtPositionOfAcceleratedBody( );
        partialMatrix.block( 3, 0, 3, 3 ) += wrtVelocityOfAcceleratedBody( );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the body exerting
     *  acceleration.
     *  \return Partial derivative of acceleration w.r.t. position of body exerting acceleration.
     */
    virtual Eigen::Matrix3d wrtPositionOfAcceleratingBody( ) = 0;

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting a
     *  acceleration.
     *  \return Partial derivative of acceleration w.r.t. velocity of body exerting acceleration.
     */
    virtual Eigen::Matrix3d wrtVelocityOfAcceleratingBody( ) = 0;

    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body exerting
     *  acceleration.
     *  \return Partial derivative of acceleration w.r.t. state of body exerting acceleration.
     */
    void wrtStateOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        partialMatrix.block( 0, 0, 3, 3 ) += wrtPositionOfAcceleratingBody( );
        partialMatrix.block( 3, 0, 3, 3 ) += wrtVelocityOfAcceleratingBody( );
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the third body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the third body.
     *  \param bodyName Name of third body.
     *  \return Partial derivative of acceleration w.r.t. position of third body.
     */
    virtual Eigen::Matrix3d wrtPositionOfAdditionalBody(
            const std::string& bodyName )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the third body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the velocity of the third body.
     *  \param bodyName Name of third body.
     *  \return Partial derivative of acceleration w.r.t. velocity of third body.
     */
    virtual Eigen::Matrix3d wrtVelocityOfAdditionalBody( const std::string& bodyName )
    {
        return Eigen::Matrix3d::Zero( );
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the Cartesian state of the third body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the Cartesian state of the third body.
     *  \param bodyName Name of third body.
     *  \return Partial derivative of acceleration w.r.t. Cartesian state of third body.
     */
    void wrtStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix, const std::string& bodyName )
    {
        partialMatrix.block( 0, 0, 3, 3 ) += wrtPositionOfAdditionalBody( bodyName );
        partialMatrix.block( 3, 0, 3, 3 ) += wrtVelocityOfAdditionalBody( bodyName );
    }

    //! Function to check whether the partial derivative w.r.t. the translational state of a third body is non-zero.
    /*!
     * Function to check whether the partial derivative w.r.t. the translational state of a third body is non-zero. This
     * function returns false by default, should be redefined in derived class if any third-bodyd dependencies exist.
     * \param bodyName Name of third body.
     * \return True if third body dependency exists, false otherwise.
     */
    virtual bool isAccelerationPartialWrtAdditionalBodyNonNull( const std::string& bodyName )
    {
        return 0;
    }

    //! Function to retrieve the name of the body undergoing acceleration.
    /*!
     *  Function to retrieve the name of the body undergoing acceleration.
     *  \return Name of the body undergoing acceleration.
     */
    std::string getAcceleratedBody( ) { return acceleratedBody_; }

    //! Function to retrieve the name of the body exerting acceleration.
    /*!
     *  Function to retrieve the name of the body exerting acceleration.
     *  \return Name of the body exerting acceleration.
     */
    std::string getAcceleratingBody( ) { return acceleratingBody_; }

    //! Function to retrieve the type of acceleration w.r.t. which partial is taken..
    /*!
     *  Function to retrieve the type of acceleration w.r.t. which partial is taken..
     *  \return Type of acceleration w.r.t. which partial is taken..
     */
    basic_astrodynamics::AvailableAcceleration getAccelerationType( )
    {
        return accelerationType_;
    }

protected:    
    //! Name of the body undergoing acceleration.
    std::string acceleratedBody_;

    //! Name of the body exerting acceleration.
    std::string acceleratingBody_;

    //! Type of acceleration w.r.t. which partial is taken..
    basic_astrodynamics::AvailableAcceleration accelerationType_;
};


}

}

}
#endif // TUDAT_ACCELERATIONPARTIALS_H
