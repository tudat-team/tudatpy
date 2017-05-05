#ifndef RELATIVISTICACCELERATIONPARTIAL_H
#define RELATIVISTICACCELERATIONPARTIAL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/relativisticAccelerationCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
namespace tudat
{

namespace acceleration_partials
{


Eigen::Matrix3d computePartialOfSchwardschildAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d relativeState, Eigen::Vector3d currentAcceleration,
        const double gravitationalParameter, const double ppnParameterGamma = 1.0, const double ppnParameterBeta = 1.0 );

Eigen::Matrix3d computePartialOfSchwardschildAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter, const double ppnParameterGamma = 1.0 );

void computePartialOfSchwardschildAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partial, const double ppnParameterGamma = 1.0, const double ppnParameterBeta = 1.0 );

Eigen::Vector3d computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter );

Eigen::Vector3d computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterBeta(
        const Eigen::Vector6d relativeState,
        const double gravitationalParameter );

class RelativisticAccelerationPartial: public AccelerationPartial
{
public:
    using AccelerationPartial::getParameterPartialFunction;

    RelativisticAccelerationPartial(
            const boost::shared_ptr< relativity::RelativisticAccelerationCorrection > accelerationModel,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::relativistic_correction_acceleration )
    {
        if( accelerationModel->getCalculateDeSitterCorrection( ) || accelerationModel->getCalculateLenseThirringCorrection( ) )
        {
            throw std::runtime_error( "Error when creating relativistic acceleration correction partial, only Schwarzschild term implemented" );
        }
        centralBodyState_  = accelerationModel->getStateFunctionOfCentralBody( );
        acceleratedBodyState_ = accelerationModel->getStateFunctionOfAcceleratedBody( );

        ppnGammaParameterFunction_ = accelerationModel->getPpnParameterGammaFunction_( );
        ppnBetaParameterFunction_ = accelerationModel->getPpnParameterBetaFunction_( );
        centralBodyGravitationalParameterFunction_ = accelerationModel->getGravitationalParameterFunctionOfCentralBody( );
        currentAccelerationFunction_ = boost::bind( &relativity::RelativisticAccelerationCorrection::getAcceleration,
                                                    accelerationModel );
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

    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
    }

    void wrtVelocityOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
    }


    bool isStateDerivativeDependentOnIntegratedNonTranslationalState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        if( ( ( stateReferencePoint.first == acceleratingBody_ ||
                ( stateReferencePoint.first == acceleratedBody_  ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            throw std::runtime_error( "Warning, dependency of relativistic acceleration on body masses not yet implemented" );
        }
        return 0;
    }

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix )
    {
        computePartialOfSchwardschildAccelerationCorrectionWrtGravitationalParameter(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ), partialMatrix,
                    ppnGammaParameterFunction_( ), ppnBetaParameterFunction_( ) );
    }


    Eigen::MatrixXd wrtPpnParameterGamma( )
    {
        return computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterGamma(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ) );
    }

    Eigen::MatrixXd wrtPpnParameterBeta( )
    {
        return computePartialOfSchwardschildAccelerationCorrectionWrtPpnParameterBeta(
                    currentRelativeState_, centralBodyGravitationalParameterFunction_( ) );
    }

    void update( const double currentTime = TUDAT_NAN );

private:

    boost::function< Eigen::Vector6d( ) > centralBodyState_;

    boost::function< Eigen::Vector6d( ) > acceleratedBodyState_;

    boost::function< double( ) > ppnGammaParameterFunction_;

    boost::function< double( ) > ppnBetaParameterFunction_;

    boost::function< double( ) > centralBodyGravitationalParameterFunction_;

    boost::function< Eigen::Vector3d( ) > currentAccelerationFunction_;

    Eigen::Matrix3d currentPartialWrtPosition_;

    Eigen::Matrix3d currentPartialWrtVelocity_;

    Eigen::Vector6d currentRelativeState_;

    Eigen::Vector3d currentAcceleration_;

};

}

}

#endif // RELATIVISTICACCELERATIONPARTIAL_H
