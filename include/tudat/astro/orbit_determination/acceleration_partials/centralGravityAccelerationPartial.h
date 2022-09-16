/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/gravitation/centralGravityModel.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

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
            const std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > gravitationalAcceleration,
            const std::string acceleratedBody,
            const std::string acceleratingBody );

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                const std::pair< std::string, std::string >& stateReferencePoint,
                const propagators::IntegratedStateType integratedStateType )
    {
        if( ( ( stateReferencePoint.first == acceleratingBody_ ||
              ( stateReferencePoint.first == acceleratedBody_  && accelerationUsesMutualAttraction_ ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            std::cout<<stateReferencePoint.first<<" "<<acceleratingBody_<<" "<<acceleratedBody_<<" "<<
                       accelerationUsesMutualAttraction_<<std::endl;
            std::cerr<<"Warning, dependency of central gravity on body masses (only on gravitational parameter) not yet implemented"<<std::endl;
        }
        return 0;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the central gravitational acceleration,
     *  position partial is computed and set.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN )
    {
        accelerationUpdateFunction_( currentTime );

        if( !( currentTime_ == currentTime ) )
        {
             acceleratedBodyState_( currentAcceleratedBodyState_ );
             centralBodyState_( currentCentralBodyState_ );
            currentGravitationalParameter_ = gravitationalParameterFunction_( );

            currentPartialWrtPosition_ = calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        currentAcceleratedBodyState_,
                        currentCentralBodyState_,
                        currentGravitationalParameter_ );

            currentTime_ = currentTime;
        }
    }

protected:

    //! Function to create a function returning the current partial w.r.t. a gravitational parameter.
    /*!
     * Function to create a function returning the current partial w.r.t. a gravitational parameter.
     * \param parameterId Identifier of parameter for which the partial is to be created.
     * \return Pair with partial function and paramater partial size. The partial function is non-empty only
     * if the parameterId input represents the gravitational parameter of acceleratingBody_ (or acceleratedBody_ if
     * accelerationUsesMutualAttraction_ is true).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function to calculate central gravity partial w.r.t. central body gravitational parameter.
    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& gravitationalParameterPartial );

    //! Function to retrieve current gravitational parameter of central body.
    std::function< double( ) > gravitationalParameterFunction_;

    //! Function to retrieve current state of body exerting acceleration.
    std::function< void( Eigen::Vector3d& ) > centralBodyState_;

    //! Function to retrieve current state of body undergoing acceleration.
    std::function< void( Eigen::Vector3d& ) > acceleratedBodyState_;

    //! Boolean denoting whether the gravitational attraction of the central body on the accelerated body is included.
    bool accelerationUsesMutualAttraction_;

    //! Current state of the body undergoing the acceleration (as set by update function).
    Eigen::Vector3d currentAcceleratedBodyState_;

    //! Current state of the body exerting the acceleration (as set by update function).
    Eigen::Vector3d currentCentralBodyState_;

    //! Current gravitational parameetr of the body exerting the acceleration (as set by update function).
    double currentGravitationalParameter_;

    //! Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! Function to update the gravitational acceleration model.
    std::function< void( const double ) > accelerationUpdateFunction_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_CENTRALGRAVITYACCELERATIONPARTIALS_H
