/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RTGACCELERATIONPARTIAL_H
#define TUDAT_RTGACCELERATIONPARTIAL_H

#include <functional>
#include <memory>

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/rtgForceVector.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/system_models/rtgAccelerationModel.h"

namespace tudat
{

namespace acceleration_partials
{

class RTGAccelerationPartial : public AccelerationPartial
{
public:
    RTGAccelerationPartial( std::shared_ptr< system_models::RTGAccelerationModel > rtgAcceleration,
                            std::string acceleratedBody,
                            std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, rtgAcceleration, basic_astrodynamics::rtg_acceleration ),
        rtgAcceleration_( rtgAcceleration )  // pos/vel partials declared to be const
    {
        if( acceleratedBody != acceleratingBody )
        {
            throw std::runtime_error(
                    "Error when setting up parameter partial of rtg acceleration - body undergoing and excerting are not the same (but are "
                    "required to be the same for given acceleration model)" );
        }
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
    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
    {}

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
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 )
    {}

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
    {}

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1,
                                        const int startRow = 0,
                                        const int startColumn = 0 )
    {}

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function to compute the partial w.r.t. time-independent empirical acceleration components
    /*!
     * Function to compute the partial w.r.t. time-independent empirical acceleration components from list of components and
     * functional shapes.
     * \param numberOfAccelerationComponents Total number of empirical acceleration components w.r.t. which partials are to
     * be computed.
     * \param accelerationIndices Map denoting list of components of accelerations that are to be computed. Key: functional
     * shape of empirical accelerations. Value: list of acceleration vaector entries that are to be used (0: radial (R),
     * 1: along-track (S), 2: cross-track (W)).
     * \param partialDerivativeMatrix Matrix of partial derivatives of accelerations w.r.t. empirical accelerations (returned
     * by reference)
     */
    // void wrtEmpiricalAccelerationCoefficientFromIndices(
    //        const int numberOfAccelerationComponents,
    //        const std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >& accelerationIndices,
    //        Eigen::MatrixXd& partialDerivativeMatrix );

    //! Function for updating common blocks of partial to current state.
    /*!
     *  Function for updating common blocks of partial to current state. Position and velocity partials are computed and set.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );

    //! Function to compute the partial w.r.t. arcwise empirical acceleration components
    /*!
     * Function to compute the partial w.r.t. arcwise empirical acceleration components
     * \param parameter Object defining the properties of the arcwise components that are to be estimated.
     * \param partialDerivativeMatrix Matrix of partial derivatives of accelerations w.r.t. empirical accelerations (returned
     * by reference)
     */
    void wrtRTGForceVector( Eigen::MatrixXd& partialDerivativeMatrix );

    void wrtRTGForceVectorMagnitude( Eigen::MatrixXd& partialDerivativeMatrix );

private:
    //! Acceleration w.r.t. which partials are to be computed.
    std::shared_ptr< system_models::RTGAccelerationModel > rtgAcceleration_;
};

}  // namespace acceleration_partials
}  // namespace tudat

#endif  // TUDAT_RTGACCELERATIONPARTIAL_H
