#include "tudat/astro/relativity/solarSystemMetric.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace relativity
{

void SolarSystemMetric::update( const Eigen::Matrix< double, 6, 1 >& state,
                                 const double time,
                                 const bool updateCurrentMetric,
                                 const bool updateCurrentChristoffelSymbols )
{
    currentEvaluationState_ = state;
    currentTime_ = time;

    updateCurrentState( true );

    if ( updateCurrentMetric )
    {
        updateMetric();
    }

    if ( updateCurrentChristoffelSymbols )
    {
        updateChristoffelSymbols();
    }
}

void SolarSystemMetric::updateMetric()
{
    currentCovariantMetricContribution_.setZero();

    currentExternalScalarPotentials_ = calculateScalarPotentials();

    currentTotalExternalScalarPotential_ = 0.0;
    currentTotalExternalScalarPotentialSquaredTerms_ = 0.0;
    currentSecondOrderPotentialSquaredTerms_ = 0.0;

    for ( unsigned int i = 0; i < currentExternalScalarPotentials_.size(); ++i )
    {
        const double potential = currentExternalScalarPotentials_.at( i );
        currentTotalExternalScalarPotential_ += potential;
        currentTotalExternalScalarPotentialSquaredTerms_ += potential * potential;

        if ( std::find( secondOrderBodyList_.begin(), secondOrderBodyList_.end(), i ) != secondOrderBodyList_.end() )
        {
            currentSecondOrderPotentialSquaredTerms_ += potential * potential;
        }
    }

    currentVectorPotential_ = calculateVectorPotential( currentExternalScalarPotentials_ );
    currentAnisotropicMatrixPotential_ = calculateAnisotropicMatrixPotential();

    currentCovariantMetricContribution_( 0, 0 ) =
        2.0 * ( currentTotalExternalScalarPotential_ * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT -
                ppnParameterSet_->getParameterBeta() * currentTotalExternalScalarPotentialSquaredTerms_ *
                physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT );

    Eigen::Vector3d spaceTimeCoupling =
        -( ppnParameterSet_->getParameterGamma() + 1.0 ) * 2.0 * currentVectorPotential_ *
        physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT;

    currentCovariantMetricContribution_.block( 1, 0, 3, 1 ) = spaceTimeCoupling;
    currentCovariantMetricContribution_.block( 0, 1, 1, 3 ) = spaceTimeCoupling.transpose();

    double secondOrderPpnParameter = ppnParameterSet_->getParameterGamma() * ppnParameterSet_->getParameterGamma() +
                                     ppnParameterSet_->getParameterBeta() - 1.0;

    double spaceSpaceDiagonalContribution =
        2.0 * ppnParameterSet_->getParameterGamma() * currentTotalExternalScalarPotential_ *
        physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT +
        2.0 * secondOrderPpnParameter * currentSecondOrderPotentialSquaredTerms_ *
        physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT;

    for ( int i = 1; i < 4; ++i )
    {
        currentCovariantMetricContribution_( i, i ) = spaceSpaceDiagonalContribution;
    }

    currentCovariantMetricContribution_.block( 1, 1, 3, 3 ) +=
        2.0 * physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT *
        ( ppnParameterSet_->getParameterGamma() + 1.0 ) * currentAnisotropicMatrixPotential_;
}

void SolarSystemMetric::updateCurrentState( const bool updateEvaluationPointIndependentVariables )
{
    for( const auto& updater : rotationUpdateFunctions_ )
    {
        updater.second( currentTime_ );
    }

    for ( size_t i = 0; i < bodyGravitationalParameterFunctions_.size(); ++i )
    {
        if ( updateEvaluationPointIndependentVariables )
        {
            currentBodyGravitationalParameters_[ i ] = bodyGravitationalParameterFunctions_[ i ]();
            currentBodyStates_[ i ] = bodyStateFunctions_[ i ]();
        }

        currentBodyRelativePositions_[ i ] = ( currentEvaluationState_ - currentBodyStates_[ i ] ).segment( 0, 3 );
        currentBodyDistances_[ i ] = currentBodyRelativePositions_[ i ].norm();
    }

    for ( const auto& [index, wrapper] : bodySphericalHarmonicGravityWrappers_ )
    {
        currentSphericalHarmonicPotentials_[ index ] =
            wrapper->getCurrentPotential( currentEvaluationState_.segment( 0, 3 ) );
    }
}

void SolarSystemMetric::updatePotentialPartials()
{
    currentScalarPotentialTimePartial_ = 0.0;
    currentScalarPotentialPositionPartial_.setZero();

    currentVectorPotentialTimePartial_.setZero();
    currentVectorPotentialPositionPartial_.setZero();

    currentInherentScalarPotentialTimePartial_ = 0.0;
    currentInherentVectorPotentialTimePartial_.setZero();

    for ( size_t i = 0; i < bodyList_.size(); ++i )
    {
        double inherentTimeDeriv = 0.0;

        if ( bodySphericalHarmonicGravityWrappers_.count( static_cast<int>(i) ) > 0 )
        {
            auto wrapper = bodySphericalHarmonicGravityWrappers_.at( i );
            currentScalarPotentialGradients_[ i ] =
                wrapper->getSphericalHarmonicPotentialGradient( currentEvaluationState_.segment( 0, 3 ), currentTime_ );

            Eigen::Matrix3d rotationInverse = wrapper->getRotationToBodyFixedFrame().inverse().toRotationMatrix();
            Eigen::Matrix3d rotationDerivative = wrapper->getRotationDerivativeToBodyFixedFrame();
            Eigen::Vector3d relativePosition = currentBodyRelativePositions_[ i ];

            Eigen::Matrix3d temp = rotationInverse * rotationDerivative;
            Eigen::Vector3d intermediate = temp * relativePosition;

            Eigen::Vector3d centralAcceleration = Eigen::Vector3d::Zero( );
            const double r = currentBodyDistances_[ i ];
            if( r > 0.0 )
            {
                centralAcceleration =
                    -currentBodyGravitationalParameters_[ i ] / ( r * r * r ) * currentBodyRelativePositions_[ i ];
            }
            const Eigen::Vector3d nonCentralAcceleration =
                currentScalarPotentialGradients_[ i ] - centralAcceleration;

            inherentTimeDeriv = nonCentralAcceleration.dot( intermediate );
            // size error (?)
            //inherentTimeDeriv =
            //    ( currentScalarPotentialGradients_[ i ].transpose() *
            //      wrapper->getRotationToBodyFixedFrame().inverse() *
            //      wrapper->getRotationDerivativeToBodyFixedFrame() *
            //      currentBodyRelativePositions_[ i ] )( 0 );
        }
        else
        {
            double r = currentBodyDistances_[ i ];
            currentScalarPotentialGradients_[ i ] =
                -currentBodyGravitationalParameters_[ i ] / ( r * r * r ) * currentBodyRelativePositions_[ i ];
        }

        currentInherentScalarPotentialTimePartial_ += inherentTimeDeriv;

        currentBodyAccelerations_[ i ] = bodyAccelerationFunctions_.at( i )( currentTime_ );

        double scalarTimeDeriv =
            currentScalarPotentialGradients_[ i ].dot( -currentBodyStates_[ i ].segment( 3, 3 ) ) + inherentTimeDeriv;

        currentScalarPotentialTimePartial_ += scalarTimeDeriv;
        currentScalarPotentialPositionPartial_ += currentScalarPotentialGradients_[ i ];

        currentVectorPotentialTimePartial_ +=
            currentExternalScalarPotentials_[ i ] * currentBodyAccelerations_[ i ] +
            currentBodyStates_[ i ].segment( 3, 3 ) * scalarTimeDeriv;

        currentInherentVectorPotentialTimePartial_ +=
            currentBodyStates_[ i ].segment( 3, 3 ) * inherentTimeDeriv;

        currentVectorPotentialPositionPartial_ +=
            currentBodyStates_[ i ].segment( 3, 3 ) * currentScalarPotentialGradients_[ i ].transpose();
    }
}

void SolarSystemMetric::updateChristoffelSymbols()
{
    updatePotentialPartials();

    for ( int i = 0; i < 4; ++i )
    {
        currentChristoffelSymbols_[ i ].setZero();
    }

    currentChristoffelSymbols_[ 0 ]( 0, 0 ) =
        -physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * currentScalarPotentialTimePartial_;

    currentChristoffelSymbols_[ 0 ].block( 0, 1, 1, 3 ) =
        -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * currentScalarPotentialPositionPartial_.transpose();

    currentChristoffelSymbols_[ 0 ].block( 1, 0, 3, 1 ) =
        -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * currentScalarPotentialPositionPartial_;

    currentChristoffelSymbols_[ 0 ].block( 1, 1, 3, 3 ) =
        physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT *
        ( ( ppnParameterSet_->getParameterGamma() + 1.0 ) *
          ( currentVectorPotentialPositionPartial_ + currentVectorPotentialPositionPartial_.transpose() ) +
          ppnParameterSet_->getParameterGamma() * Eigen::Matrix3d::Identity() * currentScalarPotentialTimePartial_ );

    Eigen::Vector3d timeTimeTerms =
        -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * currentScalarPotentialPositionPartial_ +
        physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT *
        ( -2.0 * ( ppnParameterSet_->getParameterBeta() + 1.0 ) * currentVectorPotentialTimePartial_ +
          2.0 * ( ppnParameterSet_->getParameterBeta() + ppnParameterSet_->getParameterGamma() ) *
          currentTotalExternalScalarPotential_ * currentScalarPotentialPositionPartial_ );

    Eigen::Matrix3d spaceTimeTerms =
        physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT *
        ( -( ppnParameterSet_->getParameterGamma() + 1.0 ) *
            ( currentVectorPotentialPositionPartial_ - currentVectorPotentialPositionPartial_.transpose() ) +
          ppnParameterSet_->getParameterGamma() *
            Eigen::Matrix3d::Identity() * currentScalarPotentialTimePartial_ );

    for ( int i = 1; i < 4; ++i )
    {
        currentChristoffelSymbols_[ i ]( 0, 0 ) = timeTimeTerms( i - 1 );
        currentChristoffelSymbols_[ i ].block( 0, 1, 1, 3 ) = spaceTimeTerms.block( i - 1, 0, 1, 3 );
        currentChristoffelSymbols_[ i ].block( 1, 0, 3, 1 ) = spaceTimeTerms.block( i - 1, 0, 1, 3 ).transpose();

        currentChristoffelSymbols_[ i ].block( 1, 1, 3, 3 ) -=
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            ppnParameterSet_->getParameterGamma() *
            ( Eigen::Matrix3d::Identity() * currentScalarPotentialPositionPartial_( i - 1 ) );

        currentChristoffelSymbols_[ i ].block( i, 1, 1, 3 ) +=
            ppnParameterSet_->getParameterGamma() *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            currentScalarPotentialPositionPartial_.transpose();

        currentChristoffelSymbols_[ i ].block( 1, i, 3, 1 ) +=
            ppnParameterSet_->getParameterGamma() *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            currentScalarPotentialPositionPartial_;
    }
}

} // namespace relativity

} // namespace tudat
