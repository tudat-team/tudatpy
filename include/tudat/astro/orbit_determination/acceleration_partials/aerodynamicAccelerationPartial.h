/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AERODYNAMICACCELERATIONPARTIALS_H
#define TUDAT_AERODYNAMICACCELERATIONPARTIALS_H

#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate the partials of the aerodynamic acceleration w.r.t. parameters and states.
/*!
 * Class to calculate the partials of the aerodynamic acceleration w.r.t. parameters and states. Note that the state partials
 * are computed numerically by 2nd-order central difference with perturbations hard-coded in the constructor
 */
class AerodynamicAccelerationPartial: public AccelerationPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param aerodynamicAcceleration Object that computes the aerodynamic acceleration
     * \param flightConditions Object that computes the current atmospheric and flight conditions, as well as associated angles,
     * for the body undergoing acceleration
     * \param vehicleStateGetFunction Function to retrieve the state of the body undergoing the acceleration.
     * \param vehicleStateSetFunction Function to set the state of the body undergoing the acceleration.
     * \param acceleratedBody Body undergoing acceleration.
     * \param acceleratingBody Body exerting acceleration.
     */
    AerodynamicAccelerationPartial(
            const std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration,
            const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
            const std::function< Eigen::Vector6d( ) > vehicleStateGetFunction,
            const std::function< void( const Eigen::Vector6d& ) > vehicleStateSetFunction,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::aerodynamic ),
        aerodynamicAcceleration_( aerodynamicAcceleration ), flightConditions_( flightConditions ),
        vehicleStateGetFunction_( vehicleStateGetFunction ), vehicleStateSetFunction_( vehicleStateSetFunction )
    {
        bodyStatePerturbations_ << 10.0, 10.0, 10.0, 1.0E-2, 1.0E-2, 1.0E-2;
    }

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
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentAccelerationStatePartials_.block( 0, 0, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentAccelerationStatePartials_.block( 0, 0, 3, 3 );
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
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentAccelerationStatePartials_.block( 0, 3, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentAccelerationStatePartials_.block( 0, 3, 3, 3 );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration and
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
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentAccelerationStatePartials_.block( 0, 0, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentAccelerationStatePartials_.block( 0, 0, 3, 3 );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
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
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentAccelerationStatePartials_.block( 0, 3, 3, 3 );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentAccelerationStatePartials_.block( 0, 3, 3, 3 );
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
                ( stateReferencePoint.first == acceleratedBody_  ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            throw std::runtime_error( "Warning, dependency of aerodynamic acceleration on body masses not yet implemented" );
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
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfColumns = 0;

        // Check if parameter is gravitational parameter.
        if( parameter->getParameterName( ).first ==  estimatable_parameters::constant_drag_coefficient )
        {
            // Check if parameter body is accelerated body,
            if( parameter->getParameterName( ).second.first == acceleratedBody_ )
            {
                partialFunction = std::bind(
                            &AerodynamicAccelerationPartial::computeAccelerationPartialWrtCurrentDragCoefficient,
                            this, std::placeholders::_1 );
                numberOfColumns = 1;
            }
        }

        return std::make_pair( partialFunction, numberOfColumns );
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
        int numberOfColumns = 0;

        // Check if parameter is arc-wise constant drag coefficient
        if( parameter->getParameterName( ).first ==  estimatable_parameters::arc_wise_constant_drag_coefficient )
        {
            // Check if parameter body is accelerated body,
            if( parameter->getParameterName( ).second.first == acceleratedBody_ )
            {
                if( std::dynamic_pointer_cast< estimatable_parameters::ArcWiseConstantDragCoefficient >( parameter ) != nullptr )
                {
                    partialFunction = std::bind(
                                &AerodynamicAccelerationPartial::computeAccelerationPartialWrtArcwiseDragCoefficient,
                                this, std::placeholders::_1,
                                std::dynamic_pointer_cast< estimatable_parameters::ArcWiseConstantDragCoefficient >( parameter ) );
                    numberOfColumns = parameter->getParameterSize( );

                }
                else
                {
                    throw std::runtime_error( "Error when making radiation drag, arcwise drag parameter not consistent" );
                }


            }
        }

        return std::make_pair( partialFunction, numberOfColumns );
    }

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. The partial of the acceleration w.r.t. the current
     *  state (in inertial frame) is computed numerically.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );

protected:

    //! Function to compute the partial derivative of the acceleration w.r.t. the drag coefficient
    /*!
     * Function to compute the partial derivative of the acceleration w.r.t. the drag coefficient
     * \param accelerationPartial Derivative of acceleration w.r.t. drag coefficient (returned by reference).
     */
    void computeAccelerationPartialWrtCurrentDragCoefficient( Eigen::MatrixXd& accelerationPartial )
    {
        Eigen::Quaterniond rotationToInertialFrame =
                flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                    reference_frames::aerodynamic_frame, reference_frames::inertial_frame );

        double currentAirspeed = flightConditions_->getCurrentAirspeed( );
        accelerationPartial =
                rotationToInertialFrame * Eigen::Vector3d::UnitX( ) * (
                    -0.5 * flightConditions_->getCurrentDensity( ) * currentAirspeed * currentAirspeed *
                    flightConditions_->getAerodynamicCoefficientInterface( )->getReferenceArea( ) ) /
                aerodynamicAcceleration_->getCurrentMass( );

    }

    //! Function to compute the partial derivative of the acceleration w.r.t. the arc-wise constant drag coefficient
    /*!
     * Function to compute the partial derivative of the acceleration w.r.t. the arc-wise constant drag coefficient
     * \param accelerationPartial Derivative of acceleration w.r.t. arc-wise constant drag coefficient (returned by reference).
     * \param parameter Parameter object containing information on arcwise drag coefficient that is to be estimated
     */
    void computeAccelerationPartialWrtArcwiseDragCoefficient(
            Eigen::MatrixXd& accelerationPartial,
            const std::shared_ptr< estimatable_parameters::ArcWiseConstantDragCoefficient > parameter )
    {
        // Get partial w.r.t. rdrag coefficient
        Eigen::MatrixXd partialWrtSingleParameter = Eigen::Vector3d::Zero( );
        this->computeAccelerationPartialWrtCurrentDragCoefficient( partialWrtSingleParameter );

        // Retrieve current arc
        std::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp =
                parameter->getArcTimeLookupScheme( );
        accelerationPartial.setZero( 3, parameter->getNumberOfArcs( ) );
        if( currentArcIndexLookUp->getMinimumValue( ) <= currentTime_ )
        {
            int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

            if( currentArc >= accelerationPartial.cols( ) )
            {
                throw std::runtime_error( "Error when getting arc-wise radiation pressure coefficient partials, data not consistent" );
            }

            // Set partial
            accelerationPartial.block( 0, currentArc, 3, 1 ) = partialWrtSingleParameter;
        }
    }

    //! Perturbations of Cartesian state used in the numerical (central difference) computation of
    //! currentAccelerationStatePartials_
    Eigen::Vector6d bodyStatePerturbations_;

    //! Partial derivative of aerodynamic acceleration w.r.t. current state, numerically computed by update function
    Eigen::Matrix< double, 3, 6 > currentAccelerationStatePartials_;

    //! Object that computes the aerodynamic acceleration
    std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration_;

    //! Object that computes the current atmospheric and flight conditions, as well as associated angles, for the body undergoing
    //! acceleration
    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions_;

    //! Function to retrieve the state of the body undergoing the acceleration.
    std::function< Eigen::Vector6d( ) > vehicleStateGetFunction_;

    //! Function to set the state of the body undergoing the acceleration
    std::function< void( const Eigen::Vector6d& ) > vehicleStateSetFunction_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_AERODYNAMICACCELERATIONPARTIALS_H
