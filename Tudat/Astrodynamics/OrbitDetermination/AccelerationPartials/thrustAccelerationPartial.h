/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

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
     * \param partial Partial derivative w.r.t. deltaV values of the momentum desaturation maneuvers
     */
    void wrtDesaturationDeltaVValues( Eigen::MatrixXd& accelerationPartial );

    //! Momentum wheel desaturation thrust acceleration.
    std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > thrustAcceleration_;


};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_MOMENTUMWHEELDESATURATIONPARTIALS_H
