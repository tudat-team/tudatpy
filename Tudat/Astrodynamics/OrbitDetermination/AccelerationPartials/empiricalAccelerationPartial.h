/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EMPIRICALACCELERATIONPARTIAL_H
#define TUDAT_EMPIRICALACCELERATIONPARTIAL_H

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <memory>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/empiricalAcceleration.h"

namespace tudat
{

namespace acceleration_partials
{
//! Function determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state
/*!
 *  unction determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state,
 *  from the Cartesian state as input.
 *  A first-order central difference with a user-defined Cartesian state perturbation vector is used.
 *  \param cartesianElements Nominal Cartesian elements at which the partials are to be computed
 *  \param gravitationalParameter Gravitational parameter of central body around which Keplerian orbit is given
 *  \param cartesianStateElementPerturbations Numerical perturbations of Cartesian state that are to be used
 *  \return Partial of Cartesian state wrt true anomaly of orbit.
 */
Eigen::Matrix< double, 1, 6 > calculateNumericalPartialOfTrueAnomalyWrtState(
        const Eigen::Vector6d& cartesianElements, const double gravitationalParameter,
        const Eigen::Vector6d& cartesianStateElementPerturbations );

class EmpiricalAccelerationPartial: public AccelerationPartial
{
public:
    using AccelerationPartial::getParameterPartialFunction;

    EmpiricalAccelerationPartial(
            std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration,
            std::string acceleratedBody,
            std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::empirical_acceleration ),
        empiricalAcceleration_( empiricalAcceleration ){ cartesianStateElementPerturbations << 0.1, 0.1, 0.1, 0.001, 0.001, 0.001; }

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
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPositionPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPositionPartial_;
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
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPositionPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPositionPartial_;
        }
    }

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
    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentVelocityPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentVelocityPartial_;
        }
    }

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
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentVelocityPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentVelocityPartial_;
        }
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

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
    void wrtArcWiseEmpiricalAccelerationCoefficient(
            std::shared_ptr< estimatable_parameters::ArcWiseEmpiricalAccelerationCoefficientsParameter > parameter,
            Eigen::MatrixXd& partialDerivativeMatrix );

    //! Function to compute the partial w.r.t. time-independent empirical acceleration components
    /*!
     * Function to compute the partial w.r.t. time-independent empirical acceleration components
     * \param parameter Object defining the properties of the components that are to be estimated.
     * \param partialDerivativeMatrix Matrix of partial derivatives of accelerations w.r.t. empirical accelerations (returned
     * by reference)
     */
    void wrtEmpiricalAccelerationCoefficient(
            std::shared_ptr< estimatable_parameters::EmpiricalAccelerationCoefficientsParameter > parameter,
            Eigen::MatrixXd& partialDerivativeMatrix )
    {
        return wrtEmpiricalAccelerationCoefficientFromIndices(
                    parameter->getParameterSize( ), parameter->getIndices( ), partialDerivativeMatrix );
    }

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
    void wrtEmpiricalAccelerationCoefficientFromIndices(
            const int numberOfAccelerationComponents,
            const std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >& accelerationIndices,
            Eigen::MatrixXd& partialDerivativeMatrix );


private:

    //! Acceleration w.r.t. which partials are to be computed.
    std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration_;

    //! Current partial of empirical acceleration w.r.t. position of body undergoing acceleration.
    Eigen::Matrix3d currentPositionPartial_;

    //! Current partial of empirical acceleration w.r.t. velocity of body undergoing acceleration.
    Eigen::Matrix3d currentVelocityPartial_;

    //! Perturbations to use on Cartesian state elements when computing partial of true anomaly w.r.t. state.
    Eigen::Matrix< double, 1, 6 > cartesianStateElementPerturbations;

};

}
}

#endif // TUDAT_EMPIRICALACCELERATIONPARTIAL_H
