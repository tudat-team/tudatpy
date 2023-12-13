/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H
#define TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H

#include <memory>

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"

namespace tudat
{

namespace acceleration_partials
{


//! Class to calculate the partials of the cannnonball radiation pressure acceleration w.r.t. parameters and states.
class CannonBallRadiationPressurePartial: public AccelerationPartial
{
public:

    CannonBallRadiationPressurePartial(
        const std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > cannonballTargetModel,
        const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > accelerationModel,
        const std::string& acceleratedBody, const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody,
                             basic_astrodynamics::cannon_ball_radiation_pressure ),
        sourceBodyState_( accelerationModel->getSourcePositionFunction() ),
        acceleratedBodyState_( accelerationModel->getTargetPositionFunction( ) ),
        areaFunction_( std::bind( &electromagnetism::CannonballRadiationPressureTargetModel::getArea, cannonballTargetModel ) ),
        radiationPressureCoefficientFunction_(  std::bind( &electromagnetism::CannonballRadiationPressureTargetModel::getCoefficient,
                                                           cannonballTargetModel ) ),
        radiationPressureFunction_( std::bind( &electromagnetism::RadiationPressureAcceleration::getCurrentRadiationPressure,
                                               accelerationModel ) ),
        acceleratedBodyMassFunction_( accelerationModel->getTargetMassFunction() ),
        accelerationUpdateFunction_( std::bind( &basic_astrodynamics::AccelerationModel3d::updateMembers, accelerationModel, std::placeholders::_1 ) ),
        cannonballTargetModel_( cannonballTargetModel ),
        accelerationModel_( accelerationModel )
    {

    }

    //! Destructor.
    ~CannonBallRadiationPressurePartial( ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to existing partial block.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
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

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration and
     *  adding it to exting partial block.
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

    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    void wrtNonTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType,
            const bool addContribution = true )
    {
        if( stateReferencePoint.first == acceleratedBody_ && integratedStateType == propagators::body_mass_state )
        {
            partialMatrix.block( 0, 0, 3, 1 ) +=
                   ( addContribution ? 1.0 : -1.0 ) * ( radiationPressureFunction_( ) * areaFunction_( ) * radiationPressureCoefficientFunction_( ) *
                    ( sourceBodyState_( ) - acceleratedBodyState_( ) ).normalized( ) /
                    ( acceleratedBodyMassFunction_( ) * acceleratedBodyMassFunction_( ) ) );
        }
    }

    //! Function for determining if the azcceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        bool isDependent = 0;

        // Acceleration is dependent on mass of body undergoing acceleration.
        if( stateReferencePoint.first == acceleratedBody_ && integratedStateType == propagators::body_mass_state )
        {
            isDependent = 1;
        }
        return isDependent;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function to compute the partial derivative w.r.t. a constant radiation pressure coefficient
    /*!
     * Function to compute the partial derivative w.r.t. a constant radiation pressure coefficient
     * \param partial Partial derivative w.r.t. a constant radiation pressure coefficient (returned by reference)
     */
    void wrtRadiationPressureCoefficient( Eigen::MatrixXd& partial );

    //! Function to compute the partial derivative w.r.t. an arcwise radiation pressure coefficient
    /*!
     * Function to compute the partial derivative w.r.t. an arcwise radiation pressure coefficient
     * \param parameter Parameter of arcwise radiation pressure coefficient w.r.t. which partial is to be taken
     * \param partial Partial derivative w.r.t. an arcwise radiation pressure coefficient (returned by reference)
     */
    void wrtArcWiseRadiationPressureCoefficient(
            Eigen::MatrixXd& partial,
            const std::shared_ptr< estimatable_parameters::ArcWiseRadiationPressureCoefficient > parameter )
    {
        // Get partial w.r.t. radiation pressure coefficient
        Eigen::MatrixXd partialWrtSingleParameter = Eigen::Vector3d::Zero( );
        this->wrtRadiationPressureCoefficient( partialWrtSingleParameter );

        // Retrieve current arc
        std::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp =
                parameter->getArcTimeLookupScheme( );
        partial.setZero( );
        if( currentArcIndexLookUp->getMinimumValue( ) <= currentTime_ )
        {
            int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

            if( currentArc >= partial.cols( ) )
            {
                throw std::runtime_error( "Error when getting arc-wise radiation pressure coefficient partials, data not consistent" );
            }

            // Set partial
            partial.block( 0, currentArc, 3, 1 ) = partialWrtSingleParameter;
        }

    }

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the radiation pressure acceleration,
     *  position partial is computed and set.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = 0.0 )
    {
//        std::cout<<"Computing pre "<<currentTime<<" "<<currentTime_<<std::endl;

        if( !( currentTime_ == currentTime ) )
        {
//            std::cout<<"Computing post "<<currentTime<<" "<<currentTime_<<std::endl;

            accelerationUpdateFunction_( currentTime );

            // Compute helper quantities.
            Eigen::Vector3d rangeVector = ( acceleratedBodyState_( ) - sourceBodyState_( ) );
            double range = rangeVector.norm( );
            double rangeInverse = 1.0 / ( range );

            // Compute position partial.
            currentPartialWrtPosition_ =
                    ( radiationPressureCoefficientFunction_( ) * areaFunction_( ) * radiationPressureFunction_( ) /
                      acceleratedBodyMassFunction_( ) ) * ( Eigen::Matrix3d::Identity( ) * rangeInverse - 3.0 *
                                                            rangeVector * rangeVector.transpose( ) * rangeInverse / (
                                                                range * range ) );
            currentTime_ = currentTime;
        }
    }

private:

    //! Function returning position of radiation source.
    std::function< Eigen::Vector3d( ) > sourceBodyState_;

    //! Function returning position of body undergoing acceleration.
    std::function< Eigen::Vector3d( )> acceleratedBodyState_;

    //! Function returning reflecting (or reference) area of radiation pressure on acceleratedBody_
    std::function< double( ) > areaFunction_;

    //! Function returning current radiation pressure coefficient (usually denoted C_{r}).
    std::function< double( ) > radiationPressureCoefficientFunction_;

    //! Function returning current radiation pressure (in N/m^{2})
    std::function< double( ) > radiationPressureFunction_;

    //! Function returning the mass of the body undergoing the acceleration.
    std::function< double( ) > acceleratedBodyMassFunction_;

    std::function< void( const double ) > accelerationUpdateFunction_;

    //! Current partial of acceleration w.r.t. position of body undergoing acceleration (equal to minus partial w.r.t.
    //! position of body exerting acceleration).
    Eigen::Matrix3d currentPartialWrtPosition_;

    std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > cannonballTargetModel_;
    std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > accelerationModel_;
};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_RADIATIONPRESSUREACCELERATIONPARTIAL_H
