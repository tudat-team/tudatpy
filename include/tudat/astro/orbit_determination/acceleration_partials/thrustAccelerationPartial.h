/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MOMENTUMWHEELDESATURATIONPARTIALS_H
#define TUDAT_MOMENTUMWHEELDESATURATIONPARTIALS_H

#include "tudat/astro/propulsion/thrustAccelerationModel.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"

namespace tudat
{

namespace acceleration_partials
{

class ThrustMagnitudePartial
{
    ThrustMagnitudePartial( ){ }

    virtual ~ThrustMagnitudePartial( ){ }

    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getStatePartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        return std::make_pair( nullptr, 0 );
    }

    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        return std::make_pair( nullptr, 0 );
    }

    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        return std::make_pair( nullptr, 0 );
    }
};

class ConstantThrustMagnitudePartial
{
    ConstantThrustMagnitudePartial(
            const std::shared_ptr< propulsion::ConstantThrustMagnitudeWrapper > constantThrustWrapper ):
    constantThrustWrapper_( constantThrustWrapper ){ }

    virtual ~ConstantThrustMagnitudePartial( ){ }

    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        return std::make_pair( nullptr, 0 );
    }

    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        return std::make_pair( nullptr, 0 );
    }
protected:

    std::shared_ptr< propulsion::ConstantThrustMagnitudeWrapper > constantThrustWrapper_;
};

class ThrustAccelerationPartial: public AccelerationPartial
{
public:

    ThrustAccelerationPartial(
            const std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration,
            const std::string acceleratedBody,
            const std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > >& rotationMatrixPartials  );

    ~ThrustAccelerationPartial( ){ }

    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

    void wrtPositionOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType );

    //! Function to calculate an acceleration partial wrt a rotational parameter.
    void wrtRotationModelParameter(
            Eigen::MatrixXd& accelerationPartial,
            const estimatable_parameters::EstimatebleParametersEnum parameterType,
            const std::string& secondaryIdentifier );

    void wrtBodyMass( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                      const bool addContribution = true );

    void wrtThrustMagnitude( Eigen::MatrixXd& partialMatrix,
                             const int engineIndex );


    void wrtNonTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType,
            const bool addContribution = true );

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
                std::make_pair( nullptr, 0 );

        // Check if any partials w.r.t. rotation matrix parameters are found matching the request
        if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                           parameter->getSecondaryIdentifier( ) ) ) != 0 )
        {
            partialFunction = std::make_pair(
                        std::bind( &ThrustAccelerationPartial::wrtRotationModelParameter, this,
                                         std::placeholders::_1,
                                         parameter->getParameterName( ).first,
                                         parameter->getSecondaryIdentifier( ) ), 1 );
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::constant_thrust_magnitude_parameter &&
                 parameter->getParameterName( ).second.first == acceleratedBody_ &&
                 getEngineModelIndex( parameter->getParameterName( ).second.second ) >= 0 )
        {
            int engineIndex = getEngineModelIndex( parameter->getParameterName( ).second.second );
            partialFunction = std::make_pair(
                        std::bind( &ThrustAccelerationPartial::wrtThrustMagnitude, this,
                                         std::placeholders::_1,
                                         engineIndex ), 1 );

        }
        return partialFunction;
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
                std::make_pair( nullptr, 0 );

        // Check if any partials w.r.t. rotation matrix parameters are found matching the request
        if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first,
                                                           parameter->getSecondaryIdentifier( ) ) ) != 0 )
        {
            partialFunction = std::make_pair(
                        std::bind( &ThrustAccelerationPartial::wrtRotationModelParameter, this,
                                         std::placeholders::_1,
                                         parameter->getParameterName( ).first,
                                         parameter->getSecondaryIdentifier( ) ),
                        parameter->getParameterSize( ) );
        }
        return partialFunction;
    }


    void update( const double currentTime = TUDAT_NAN )
    {

        if( !( currentTime_ == currentTime ) )
        {
            thrustAcceleration_->updateMembers( currentTime );
            currentTime_ = currentTime;
        }
    }

protected:

    int getEngineModelIndex(
            const std::string& engineId )
    {
        int engineModelIndex = -1;
        for( unsigned int i = 0; i < thrustSources_.size( ); i++ )
        {
            if( thrustSources_.at( i )->getEngineName( ) == engineId )
            {
                engineModelIndex = i;
            }
        }
        return engineModelIndex;
    }
    std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration_;

    std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > > rotationMatrixPartials_;

    std::vector< std::shared_ptr< system_models::EngineModel > > thrustSources_;

    std::vector< unsigned int > massDependentThrustSources_;

    bool isAccelerationDependentOnMass_;

    bool isAccelerationDependentOnTranslationalState_;

    bool isRotationDirectionBased_;
};


//! Class to calculate the partials of the momentum wheel desaturation acceleration w.r.t. parameters and states.
class MomentumWheelDesaturationPartial: public AccelerationPartial
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param thrustAcceleration Momentum wheel desaturation thrust acceleration model.
     * \param acceleratedBody Name of the body undergoing acceleration.
     */
    MomentumWheelDesaturationPartial(
            const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > thrustAcceleration,
            const std::string acceleratedBody );

    //! Destructor.
    ~MomentumWheelDesaturationPartial( ){ }

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
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

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
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

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
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the momentum wheel
     *  desaturation acceleration, only the thrust acceleration model is updated.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN )
    {

        if( !( currentTime_ == currentTime ) )
        {
            thrustAcceleration_->updateMembers( currentTime );
            currentTime_ = currentTime;
        }
    }

protected:

    //! Function to compute the partial derivative w.r.t. the deltaV values of the momentum desaturation maneuvers
    /*!
     * Function to compute the partial derivative w.r.t. the deltaV values of the momentum desaturation maneuvers
     * \param accelerationPartial Partial derivative w.r.t. deltaV values of the momentum desaturation maneuvers
     */
    void wrtDesaturationDeltaVValues( Eigen::MatrixXd& accelerationPartial );

    //! Momentum wheel desaturation thrust acceleration.
    std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > thrustAcceleration_;


};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_MOMENTUMWHEELDESATURATIONPARTIALS_H
