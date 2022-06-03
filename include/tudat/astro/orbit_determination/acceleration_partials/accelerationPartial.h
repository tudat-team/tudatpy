/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/stateDerivativePartial.h"

namespace tudat
{

namespace acceleration_partials
{

inline Eigen::Vector3d computeDerivativeOfForceBasedAccelerationWrtMass(
        const Eigen::Vector3d& acceleration,
        const double mass )
{
    if( mass == 0.0 )
    {
        throw std::runtime_error( "Error calculating derivative of force-based acceleration w.r.t. mass; cannot calculate for zero mass." );
    }
    return -acceleration / mass;
}

//! Base class for objects calculating partial derivatives of accelerations w.r.t. states, model parameters.
/*!
 *  Base class for objects calculating partial derivatives of accelerations  w.r.t. states, model parameters. Such
 *  calculations are used in orbit determination, for the computation of the state transition; sensitivity matrices.
 *  Derived classes implement derivative-calculating models for specific acceleration models, so that the calculation
 *  of all partials of a single type acceleration model is encompassed in a single derived class.
 */
class AccelerationPartial: public orbit_determination::StateDerivativePartial
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
        StateDerivativePartial( propagators::translational_state, std::make_pair( acceleratedBody, "" ) ),
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
        case propagators::translational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for dynamics" );
            }
            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == acceleratedBody_ )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAcceleratedBody, this, std::placeholders::_1 ), 6 );
            }
            else if( stateReferencePoint.first == acceleratingBody_ )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAcceleratingBody, this, std::placeholders::_1 ), 6 );
            }
            else if( isAccelerationPartialWrtAdditionalBodyNonnullptr( stateReferencePoint.first ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAdditionalBody,
                                                               this, std::placeholders::_1, stateReferencePoint.first ), 6 );
            }
            break;
        }
        case propagators::rotational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for body mass" );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtNonTranslationalStateOfAdditionalBody,
                                                               this, std::placeholders::_1, stateReferencePoint, integratedStateType, true ), 7 );
            }
            break;
        }
        case propagators::body_mass_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for body mass" );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtNonTranslationalStateOfAdditionalBody,
                                                               this, std::placeholders::_1, stateReferencePoint, integratedStateType, true ), 1 );
            }
            break;
        }
        case propagators::custom_state:
        {
            break;
        }
        default:
            std::string errorMessage =
                    "Error when getting state derivative partial acceleration model, dynamics type " +
                    std::to_string( integratedStateType ) + "not recognized" ;
            throw std::runtime_error( errorMessage );
            break;
        }


        return partialFunction;
    }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state (default none).
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented is returned in this base class function, but may be overriden by derived class.
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
        switch( integratedStateType )
        {
        case propagators::translational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when checking state derivative partial dependency of acceleration model, cannot have reference point on body for dynamics" );
            }

            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == acceleratedBody_ || stateReferencePoint.first == acceleratingBody_ ||
                     isAccelerationPartialWrtAdditionalBodyNonnullptr( stateReferencePoint.first ) )
            {
                isDependent = 1;
            }
            break;
        }
        default:
            isDependent = isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType );
            break;
        }
        return isDependent;
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ) = 0;

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of the accelerated body and
     *  adding it to the existing partial block. Function may be overridden in derived class, default dependency is none.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    { }

    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body
     *  undergoing acceleration  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian state of body
     *  undergoing acceleration where current partial is to be added.
     */
    void wrtStateOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        wrtPositionOfAcceleratedBody( partialMatrix, true, 0, 0 );
        wrtVelocityOfAcceleratedBody( partialMatrix, true, 0, 3 );
    }

    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the body exerting acceleration.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the body exerting
     *  acceleration and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtPositionOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ) = 0;

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of the body exerting a
     *  acceleration. . Function may be overridden in derived class, default dependency is none.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  exerting acceleration where current partial is to be added.
     */
    virtual void wrtVelocityOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body exerting
     *  acceleration and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian state of body
     *  exerting acceleration where current partial is to be added.
     */
    void wrtStateOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        wrtPositionOfAcceleratingBody( partialMatrix, true, 0, 0 );
        wrtVelocityOfAcceleratingBody( partialMatrix, true, 0, 3 );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of the third body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of the third body
     *  and adding it to the existing partial block. Function may be overridden in derived class, default dependency is none.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtPositionOfAdditionalBody(
            const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of the third body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of the third body
     *  and adding it to the existing partial block. . Function may be overridden in derived class, default dependency is
     *  none.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtVelocityOfAdditionalBody(
            const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the third body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the third body
     *  and adding it to the existing partial block.
     *  \param bodyName Name of third body.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian state of third body where
     *  current partial is to be added.
     */
    void wrtStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix, const std::string& bodyName )
    {
        wrtPositionOfAdditionalBody( bodyName, partialMatrix, true, 0, 0 );
        wrtVelocityOfAdditionalBody( bodyName, partialMatrix, true, 0, 3 );
    }

    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block. Function may be overridden in derived class, default dependency is
     *  none.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    virtual void wrtNonTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType,
            const bool addContribution = true ){ }

    //! Function to check whether the partial derivative w.r.t. the translational state of a third body is non-zero.
    /*!
     * Function to check whether the partial derivative w.r.t. the translational state of a third body is non-zero. This
     * function returns false by default, should be redefined in derived class if any third-bodyd dependencies exist.
     * \param bodyName Name of third body.
     * \return True if third body dependency exists, false otherwise.
     */
    virtual bool isAccelerationPartialWrtAdditionalBodyNonnullptr( const std::string& bodyName )
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


} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_ACCELERATIONPARTIALS_H
