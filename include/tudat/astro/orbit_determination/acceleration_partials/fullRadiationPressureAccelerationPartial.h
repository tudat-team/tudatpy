/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLRADIATIONPRESSUREACCELERATIONPARTIALS_H
#define TUDAT_FULLRADIATIONPRESSUREACCELERATIONPARTIALS_H

#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/radiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{



//! Class to calculate the partials of the custom acceleration w.r.t. parameters and states.
class RadiationPressureAccelerationPartial: public AccelerationPartial
{
public:

    RadiationPressureAccelerationPartial(
            const std::string acceleratedBody,
            const std::string acceleratingBody,
            const std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > radiationPressureAcceleration,
            const std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > customAccelerationPartialSet = nullptr );

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
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtUndergoingState_.block( 0, 0, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtUndergoingState_.block( 0, 0, 3, 3 );
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
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtExertingState_.block( 0, 0, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtExertingState_.block( 0, 0, 3, 3 );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtUndergoingState_.block( 0, 3, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtUndergoingState_.block( 0, 3, 3, 3 );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtExertingState_.block( 0, 3, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtExertingState_.block( 0, 3, 3, 3 );
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
        return 0;
    }

    void createCustomParameterPartialFunction(
        Eigen::MatrixXd& partialBlock,
        const std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > customPartialCalculator )
    {
        partialBlock = customPartialCalculator->computePartial( currentTime_, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
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
        int parameterSize = 0;
        if( customAccelerationPartialSet_->customDoubleParameterPartials_.count( parameter->getParameterName() )!= 0 )
        {
            partialFunction = std::bind( &RadiationPressureAccelerationPartial::createCustomParameterPartialFunction, this,
                                         std::placeholders::_1,
                                         customAccelerationPartialSet_->customDoubleParameterPartials_.at( parameter->getParameterName() ) );
            parameterSize = 1;
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::radiation_pressure_coefficient &&
            parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            if ( std::dynamic_pointer_cast<electromagnetism::CannonballRadiationPressureTargetModel>(
                radiationPressureAcceleration_->getTargetModel( ) ) != nullptr )
            {
                partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtRadiationPressureCoefficient,
                                             this, std::placeholders::_1, std::dynamic_pointer_cast<electromagnetism::CannonballRadiationPressureTargetModel>(
                        radiationPressureAcceleration_->getTargetModel( ) ) );
                parameterSize = 1;
            }
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::specular_reflectivity &&
                 parameter->getParameterName( ).second.first == acceleratedBody_)
        {
            if(std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                    radiationPressureAcceleration_->getTargetModel( ) ) != nullptr){
                throw std::runtime_error( "Error when creating specular reflectivity partial, PaneledRadiationPressureTargetModel not specified" );
            }
            if(parameter->getParameterName( ).second.second == ""){
                throw std::runtime_error( "Error when creating specular reflectivity partial, panel group name not specified" );
            }
            else{
                partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtSpecularReflectivity,
                                             this,
                                             std::placeholders::_1,
                                             std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                                                     radiationPressureAcceleration_->getTargetModel( ) ),
                                             parameter->getParameterName( ).second.second );
                parameterSize = 1;
            };
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::diffuse_reflectivity &&
                 parameter->getParameterName( ).second.first == acceleratedBody_)
        {
            if(std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                    radiationPressureAcceleration_->getTargetModel( ) ) != nullptr){
                throw std::runtime_error( "Error when creating diffuse reflectivity partial, PaneledRadiationPressureTargetModel not specified" );
            }
            if(parameter->getParameterName( ).second.second == ""){
                throw std::runtime_error( "Error when creating diffuse reflectivity partial, panel group name not specified" );
            }
            else{
                partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtDiffuseReflectivity,
                                             this,
                                             std::placeholders::_1,
                                             std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                                                     radiationPressureAcceleration_->getTargetModel( ) ),
                                             parameter->getParameterName( ).second.second );
                parameterSize = 1;
            };
        }
        // Check if parameter dependency exists.
        else if( parameter->getParameterName( ).second.first == acceleratedBody_ && parameter->getParameterName( ).second.second == acceleratingBody_ )
        {
            switch( parameter->getParameterName( ).first )
            {
            case estimatable_parameters::source_direction_radiation_pressure_scaling_factor:

                partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourceDirectionScaling,
                                             radiationPressureAcceleration_,
                                             std::placeholders::_1 );
                parameterSize = 1;

                break;
            case estimatable_parameters::source_perpendicular_direction_radiation_pressure_scaling_factor:

                partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourcePerpendicularDirectionScaling,
                                             radiationPressureAcceleration_,
                                             std::placeholders::_1 );
                parameterSize = 1;

                break;
            default:
                break;
            }
        }

        return std::make_pair( partialFunction, parameterSize );
    }

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
        int parameterSize = 0;
        if( customAccelerationPartialSet_->customVectorParameterPartials_.count( parameter->getParameterName() )!= 0 )
        {
            partialFunction = std::bind( &RadiationPressureAccelerationPartial::createCustomParameterPartialFunction, this,
                                         std::placeholders::_1,
                                         customAccelerationPartialSet_->customVectorParameterPartials_.at( parameter->getParameterName() ) );
            parameterSize = parameter->getParameterSize( );
        }
        return std::make_pair( partialFunction, parameterSize );
    }

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the central gravitational acceleration,
     *  position partial is computed and set.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );

protected:

    void wrtRadiationPressureCoefficient(
        Eigen::MatrixXd& partial, std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > targetModel );

    void wrtSpecularReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId);

    void wrtDiffuseReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId);

    std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > radiationPressureAcceleration_;

    std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > customAccelerationPartialSet_;

    std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > bodyUndergoingPositionPartial_;

    std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > bodyExertingPositionPartial_;

    Eigen::Matrix< double, 3, 6 > currentPartialWrtUndergoingState_;

    Eigen::Matrix< double, 3, 6 > currentPartialWrtExertingState_;


};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_FULLRADIATIONPRESSUREACCELERATIONPARTIALS_H
