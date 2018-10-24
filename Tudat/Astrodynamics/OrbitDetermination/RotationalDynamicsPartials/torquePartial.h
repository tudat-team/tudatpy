/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TORQUE_PARTIAL_H
#define TUDAT_TORQUE_PARTIAL_H

#include <string>
#include <map>
#include <Eigen/Core>

#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/stateDerivativePartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Base class for objects calculating partial derivatives of torques w.r.t. states, model parameters.
/*!
 *  Base class for objects calculating partial derivatives of torques  w.r.t. states, model parameters. Such
 *  calculations are used in orbit determination, for the computation of the state transition; sensitivity matrices.
 *  Derived classes implement derivative-calculating models for specific torque models, so that the calculation
 *  of all partials of a single type torque model is encompassed in a single derived class.
 */
class TorquePartial: public orbit_determination::StateDerivativePartial
{

public:
    //! Base class constructor.
    /*!
     *  Constructor of base class, sets the base class member variables identifying the body undergoing and exerting the
     *  torque.
     *  \param bodyUndergoingTorque Body undergoing torque.
     *  \param bodyExertingTorque Body exerting torque.
     *  \param torqueType Type of torque w.r.t. which partial is taken.
     */
    TorquePartial( const std::string& bodyUndergoingTorque, const std::string& bodyExertingTorque,
                   const basic_astrodynamics::AvailableTorque torqueType ):
        StateDerivativePartial( propagators::translational_state, std::make_pair( bodyUndergoingTorque, "" ) ),
        bodyUndergoingTorque_( bodyUndergoingTorque ), bodyExertingTorque_( bodyExertingTorque ),torqueType_( torqueType ) { }

    //! Virtual destructor.
    virtual ~TorquePartial( ) { }

    //! Function for determining if the torque is dependent on a non-rotational integrated state.
    /*!
     *  Function for determining if the torque is dependent on a non-rotational integrated state. Default none, may be
     *  overridden by derived class
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    virtual bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return false;
    }

    //! Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
    /*!
     * Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for rotational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return Pair with function, returning partial derivative, and number of columns in partial vector,
     */
    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
    getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        // Initialize to empty function; 0 parameter size.
        std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
                partialFunction = std::make_pair( std::function< void( Eigen::Block< Eigen::MatrixXd > ) >( ), 0 );

        // Check if state dependency exists
        switch( integratedStateType )
        {
        case propagators::rotational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial torque model, cannot have reference point on body for dynamics" );
            }
            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == bodyUndergoingTorque_ )
            {
                partialFunction = std::make_pair( std::bind( &TorquePartial::wrtRotationalStateOfAcceleratedBody, this, std::placeholders::_1 ), 7 );
            }
            else if( stateReferencePoint.first == bodyExertingTorque_ )
            {
                partialFunction = std::make_pair( std::bind( &TorquePartial::wrtRotationalStateOfAcceleratingBody, this, std::placeholders::_1 ), 7 );
            }
            else if( isTorquePartialWrtAdditionalBodyNonNull( stateReferencePoint.first ) )
            {
                partialFunction = std::make_pair( std::bind( &TorquePartial::wrtRotationalStateOfAdditionalBody,
                                                               this, std::placeholders::_1, stateReferencePoint.first ), 3 );
            }
            break;
        }
        case propagators::translational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting torque partial, cannot have reference point on body for body mass" );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &TorquePartial::wrtNonRotationalStateOfAdditionalBody,
                                                               this, std::placeholders::_1, stateReferencePoint, integratedStateType ), 1 );
            }
        }
        case propagators::custom_state:
        {
            break;
        }
        default:
            std::string errorMessage =
                    "Error when getting state derivative partial torque model, dynamics type " +
                    std::to_string( integratedStateType ) + "not recognized" ;
            throw std::runtime_error( errorMessage );
            break;
        }


        return partialFunction;
    }

    //! Pure virtual function for calculating the partial of the torque w.r.t. the orientation of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the torque w.r.t. the orientation of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtOrientationOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ) = 0;

    //! Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body and
     *  adding it to the existing partial block. Function may be overridden in derived class, default dependency is none.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. angular velocity of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtRotationalVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    { }

    //! Function for calculating the partial of the torque w.r.t. the rotational state of the body undergoing torque.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the rotational state of the body
     *  undergoing torque  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Cartesian state of body
     *  undergoing torque where current partial is to be added.
     */
    void wrtRotationalStateOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        wrtOrientationOfAcceleratedBody( partialMatrix, true, 0, 0 );
        wrtRotationalVelocityOfAcceleratedBody( partialMatrix, true, 0, 4 );
    }

    //! Pure virtual function for calculating the partial of the torque w.r.t. the orientation of the body exerting torque.
    /*!
     *  Pure virtual function for calculating the partial of the torque w.r.t. the orientation of the body exerting
     *  torque and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of body
     *  exerting torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtOrientationOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    { }

    //! Function for calculating the partial of the torque w.r.t. the angular velocity of the body exerting torque.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the angular velocity of the body exerting a
     *  torque. . Function may be overridden in derived class, default dependency is none.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. angular velocity of body
     *  exerting torque where current partial is to be added.
     */
    virtual void wrtRotationalVelocityOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 ){ }

    //! Function for calculating the partial of the torque w.r.t. the Cartesian state of the body exerting torque.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the Cartesian state of the body exerting
     *  torque and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Cartesian state of body
     *  exerting torque where current partial is to be added.
     */
    void wrtRotationalStateOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        wrtOrientationOfAcceleratingBody( partialMatrix, true, 0, 0 );
        wrtRotationalVelocityOfAcceleratingBody( partialMatrix, true, 0, 3 );
    }

    //! Function for calculating the partial of the torque w.r.t. the orientation of the third body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the orientation of the third body
     *  and adding it to the existing partial block. Function may be overridden in derived class, default dependency is none.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtOrientationOfAdditionalBody(
            const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

    //! Function for calculating the partial of the torque w.r.t. the angular velocity of the third body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the angular velocity of the third body
     *  and adding it to the existing partial block. . Function may be overridden in derived class, default dependency is
     *  none.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. angular velocity of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtRotationalVelocityOfAdditionalBody(
            const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 ){ }

    //! Function for calculating the partial of the torque w.r.t. the Cartesian state of the third body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the Cartesian state of the third body
     *  and adding it to the existing partial block.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Cartesian state of third body where
     *  current partial is to be added.
     */
    void wrtRotationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix, const std::string& bodyName )
    {
        wrtOrientationOfAdditionalBody( bodyName, partialMatrix, true, 0, 0 );
        wrtRotationalVelocityOfAdditionalBody( bodyName, partialMatrix, true, 0, 3 );
    }

    //! Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
    /*!
     *  Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
     *  and adding it to the existing partial block. Function may be overridden in derived class, default dependency is
     *  none.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     */
    virtual void wrtNonRotationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ){ }

    //! Function to check whether the partial derivative w.r.t. the rotational state of a third body is non-zero.
    /*!
     * Function to check whether the partial derivative w.r.t. the rotational state of a third body is non-zero. This
     * function returns false by default, should be redefined in derived class if any third-bodyd dependencies exist.
     * \param bodyName Name of third body.
     * \return True if third body dependency exists, false otherwise.
     */
    virtual bool isTorquePartialWrtAdditionalBodyNonNull( const std::string& bodyName )
    {
        return 0;
    }

    //! Function to retrieve the name of the body undergoing torque.
    /*!
     *  Function to retrieve the name of the body undergoing torque.
     *  \return Name of the body undergoing torque.
     */
    std::string getBodyUndergoingTorque( ) { return bodyUndergoingTorque_; }

    //! Function to retrieve the name of the body exerting torque.
    /*!
     *  Function to retrieve the name of the body exerting torque.
     *  \return Name of the body exerting torque.
     */
    std::string getBodyExertingTorque( ) { return bodyExertingTorque_; }

    //! Function to retrieve the type of torque w.r.t. which partial is taken..
    /*!
     *  Function to retrieve the type of torque w.r.t. which partial is taken..
     *  \return Type of torque w.r.t. which partial is taken..
     */
    basic_astrodynamics::AvailableTorque getTorqueType( )
    {
        return torqueType_;
    }

protected:    
    //! Name of the body undergoing torque.
    std::string bodyUndergoingTorque_;

    //! Name of the body exerting torque.
    std::string bodyExertingTorque_;

    //! Type of torque w.r.t. which partial is taken..
    basic_astrodynamics::AvailableTorque torqueType_;
};


} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_TORQUE_PARTIAL_H
