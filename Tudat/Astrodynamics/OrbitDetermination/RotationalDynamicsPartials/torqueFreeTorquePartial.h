/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TORQUEFREETORQUEPARTIALS_H
#define TUDAT_TORQUEFREETORQUEPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

class TorqueFreeTorquePartial: public TorquePartial
{
public:

    TorqueFreeTorquePartial(
            boost::function< Eigen::Vector3d( ) > angularVelocityFunction,
            boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const std::string acceleratedBody ):
        TorquePartial( acceleratedBody, acceleratedBody, basic_astrodynamics::torque_free ),
        angularVelocityFunction_( angularVelocityFunction ),
        inertiaTensorFunction_( inertiaTensorFunction ){ }

    ~TorqueFreeTorquePartial( ){ }

    //! Function for determining if the acceleration is dependent on a non-rotational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-rotational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedNonRotationalState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return 0;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        boost::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    void wrtOrientationOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 ){ }

    void wrtRotationalVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialDerivativeWrtAngularVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialDerivativeWrtAngularVelocity_;
        }
    }

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return false;
    }

    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentAngularVelocityVector_ = angularVelocityFunction_( );
            currentInertiaTensor_ = inertiaTensorFunction_( );
            currentAngularMomentumVector_ = currentInertiaTensor_ * currentAngularVelocityVector_;

            currentPartialDerivativeWrtAngularVelocity_ =
                    currentInertiaTensor_.inverse( ) * (
                        linear_algebra::getCrossProductMatrix( currentAngularMomentumVector_ ) -
                        linear_algebra::getCrossProductMatrix( currentAngularVelocityVector_ ) * currentInertiaTensor_ );
        }
    }

protected:

    boost::function< Eigen::Vector3d( ) > angularVelocityFunction_;

    boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    Eigen::Vector3d currentAngularVelocityVector_;

    Eigen::Matrix3d currentInertiaTensor_;

    Eigen::Vector3d currentAngularMomentumVector_;

    Eigen::Matrix3d currentPartialDerivativeWrtAngularVelocity_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_TORQUEFREETORQUEPARTIALS_H
