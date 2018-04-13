/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
#define TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
        const double premultiplier,
        const Eigen::Matrix3d& inertiaTensor,
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const Eigen::Vector3d& intertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions );

//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class SecondDegreeGravitationalTorquePartial: public TorquePartial
{
public:

    SecondDegreeGravitationalTorquePartial(
            const boost::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > torqueModel,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::second_order_gravitational_torque ),
        torqueModel_( torqueModel )
    {
        currentRotationMatrixDerivativesWrtQuaternion_.resize( 4 );
    }

    ~SecondDegreeGravitationalTorquePartial( ){ }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
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

    virtual void wrtOrientationOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialDerivativeWrtQuaternion_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialDerivativeWrtQuaternion_;
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
            torqueModel_->updateMembers( currentTime );
            currentQuaternionVector_ = linear_algebra::convertQuaternionToVectorFormat(
                        ( torqueModel_->getCurrentRotationToBodyFixedFrame( ) ).inverse( ) );
            linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
                        currentQuaternionVector_,  currentRotationMatrixDerivativesWrtQuaternion_ );
            currentPartialDerivativeWrtQuaternion_ = getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
                        torqueModel_->getCurrentTorqueMagnitudePremultiplier( ),
                        torqueModel_->getCurrentInertiaTensorOfRotatingBody( ),
                        torqueModel_->getCurrentRelativeBodyFixedPositionOfBodySubjectToTorque(),
                        torqueModel_->getCurrentRelativePositionOfBodySubjectToTorque( ),
                        currentRotationMatrixDerivativesWrtQuaternion_ );
        }
    }

protected:

    boost::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > torqueModel_;

    Eigen::Vector4d currentQuaternionVector_;

    std::vector< Eigen::Matrix3d > currentRotationMatrixDerivativesWrtQuaternion_;

    Eigen::Matrix< double, 3, 4 > currentPartialDerivativeWrtQuaternion_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
