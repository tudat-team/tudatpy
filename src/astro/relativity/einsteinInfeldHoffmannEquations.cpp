
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity//einsteinInfeldHoffmannEquations.h"

namespace tudat
{

namespace gravitation
{


EinsteinInfeldHoffmannEquations::EinsteinInfeldHoffmannEquations(
        const std::vector< std::string > acceleratedBodies,
        const std::vector< std::function< double( ) > > gravitationalParameterFunction,
        const std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions,
        const std::function< double( ) > ppnGammaFunction,
        const std::function< double( ) > ppnBetaFunction ):
    acceleratedBodies_( acceleratedBodies ), gravitationalParameterFunction_( gravitationalParameterFunction ),
    bodyStateFunctions_( bodyStateFunctions ), ppnGammaFunction_( ppnGammaFunction ), ppnBetaFunction_( ppnBetaFunction )
{
    for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
    {
        currentGravitationalParameters_.resize( acceleratedBodies_.size( ) );
        currentPositions_.resize( acceleratedBodies_.size( ) );
        currentVelocities_.resize( acceleratedBodies_.size( ) );
        currentSquareSpeeds_.resize( acceleratedBodies_.size( ) );

        currentRelativePositions_.resize( acceleratedBodies_.size( ) );
        currentRelativeDistances_.resize( acceleratedBodies_.size( ) );
        secondaryCentralBodyTerms_.resize( acceleratedBodies_.size( ) );
        relativePositionVelocityProduct_.resize( acceleratedBodies_.size( ) );
        velocityInnerProducts_.resize( acceleratedBodies_.size( ) );
        totalPointMassAccelerations_.resize( acceleratedBodies_.size( ) );
        currentAccelerations_.resize( acceleratedBodies.size( ) );
        for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
        {
            currentRelativePositions_[ i ].resize( acceleratedBodies_.size( ) );
            currentRelativeDistances_[ i ].resize( acceleratedBodies_.size( ) );
            secondaryCentralBodyTerms_[ i ].resize( acceleratedBodies_.size( ) );
            relativePositionVelocityProduct_[ i ].resize( acceleratedBodies_.size( ) );
            velocityInnerProducts_.resize( acceleratedBodies_.size( ) );
        }
    }
}

void EinsteinInfeldHoffmannEquations::update( const double currentTime )
{
    if( currentTime_ != currentTime )
    {
        currentPpnGamma_ = ppnGammaFunction_( );
        currentPpnBeta_ = ppnBetaFunction_( );

        Eigen::Matrix< double, 6, 1 > currentBodyState;
        for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
        {
            totalPointMassAccelerations_[ i ].setZero( );
            currentBodyState = bodyStateFunctions_[ i ]( );

            currentGravitationalParameters_[ i ] = gravitationalParameterFunction_[ i ]( );
            currentPositions_[ i ] = currentBodyState.segment( 0, 3 );
            currentVelocities_[ i ] = currentBodyState.segment( 3, 3 );
            currentSquareSpeeds_[ i ] = currentVelocities_[ i ].dot( currentVelocities_[ i ] );

            for( unsigned int j = 0; j < acceleratedBodies_.size( ); j++ )
            {
                if( i != j )
                {
                    currentRelativePositions_[ i ][ j ] = currentPositions_[ j ] - currentPositions_[ i ];
                    currentRelativeDistances_[ i ][ j ] = currentRelativePositions_[ i ][ j ].norm( );
                    secondaryCentralBodyTerms_[ i ][ j ] = currentGravitationalParameters_[ j ] /
                            currentRelativeDistances_[ i ][ j ];
                    relativePositionVelocityProduct_[ i ][ j ] = ( currentPositions_[ i ] - currentPositions_[ j ] ).dot( currentVelocities_[ j ] );
                    velocityInnerProducts_[ i ][ j ] = currentVelocities_[ i ].dot( currentVelocities_[ j ] );
                    singlePointMassAccelerations_[ i ][ j ] = currentGravitationalParameters_[ j ] * currentRelativePositions_[ i ][ j ] /
                                                              ( currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] );
                    totalPointMassAccelerations_[ i ] += singlePointMassAccelerations_[ i ][ j ];
                }
                else
                {
                    currentRelativePositions_[ i ][ j ].setZero( );
                    currentRelativeDistances_[ i ][ j ] = 0.0;
                    secondaryCentralBodyTerms_[ i ][ j ] = 0.0;
                    relativePositionVelocityProduct_[ i ][ j ] = 0.0;
                    velocityInnerProducts_[ i ][ j ] = 0.0;
                }
            }

        }
        calculateAccelerations( );
    }
    currentTime_ = currentTime;
}

void EinsteinInfeldHoffmannEquations::calculateAccelerations( )
{
    using namespace tudat::physical_constants;

    double summedTerm1, summedTerm2;
    Eigen::Vector3d singleTerm;

    for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
    {
        currentAccelerations_[ i ].setZero( );

        summedTerm1 = 0;
        for( unsigned int l = 0; l < acceleratedBodies_.size( ); l++ )
        {
            if( l!= i )
            {
                summedTerm1 += secondaryCentralBodyTerms_[ i ][ l ];
            }
        }
        singleTerm.setZero( );

        for( unsigned int j = 0; j < acceleratedBodies_.size( ); j++ )
        {
            currentSingleAccelerations_[ i ][ j ].setZero( );
            if( i != j )
            {
                summedTerm2 = 0;
                for( unsigned int k = 0; k < acceleratedBodies_.size( ); k++ )
                {
                    if( k!= j )
                    {
                        summedTerm2 += secondaryCentralBodyTerms_[ j ][ k ];
                    }
                }

                currentSingleAccelerations_[ i ][ j ] += singlePointMassAccelerations_[ i ][ j ]  * ( 1.0 - INVERSE_SQUARE_SPEED_OF_LIGHT *
                        ( 2.0 * ( currentPpnBeta_ + currentPpnGamma_ ) * summedTerm1 - ( 2.0 * currentPpnBeta_ - 1.0 ) * summedTerm2 +
                          currentPpnGamma_ * currentSquareSpeeds_[ i ] + ( 1.0 + currentPpnGamma_ ) * currentSquareSpeeds_[ j ] -
                          2.0 * ( 1.0 + currentPpnGamma_ ) * velocityInnerProducts_[ i ][ j ] -
                          1.5 * relativePositionVelocityProduct_[ i ][ j ] / currentRelativeDistances_[ i ][ j ] +
                          0.5 * currentRelativePositions_[ i ][ j ].dot( totalPointMassAccelerations_[ j ] ) ) );

                currentSingleAccelerations_[ i ][ j ] += INVERSE_SQUARE_SPEED_OF_LIGHT * (
                            currentGravitationalParameters_[ j ] / ( currentRelativeDistances_[ i ][ j ]  * currentRelativeDistances_[ i ][ j ]  *
                                                                     currentRelativeDistances_[ i ][ j ]  ) *
                            ( relativePositionVelocityProduct_[ j ][ i ] * ( 2.0 + 2.0 * currentPpnGamma_ ) -
                              relativePositionVelocityProduct_[ i ][ j ] * ( 1.0 + 2.0 * currentPpnGamma_ ) ) ) *
                        ( currentVelocities_[ i ] - currentVelocities_[ j ] );
                currentSingleAccelerations_[ i ][ j ] += INVERSE_SQUARE_SPEED_OF_LIGHT * ( 3.0 + 4.0 * currentPpnGamma_ ) / 2.0 * currentGravitationalParameters_[ j ] * totalPointMassAccelerations_[ j ] /
                    currentRelativeDistances_[ i ][ j ];
            }
            currentAccelerations_[ i ] += currentSingleAccelerations_[ i ][ j ];
        }
    }
}

