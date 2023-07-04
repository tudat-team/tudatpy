#include <iostream>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"

namespace tudat
{

namespace relativity
{


EinsteinInfeldHoffmannEquations::EinsteinInfeldHoffmannEquations(
        const std::vector< std::string > acceleratedBodies,
        const std::vector< std::string > acceleratingBodies,
        const std::vector< std::function< double( ) > > gravitationalParameterFunction,
        const std::vector< std::function< Eigen::Matrix< double, 6, 1 >( ) > > bodyStateFunctions,
        const std::function< double( ) > ppnGammaFunction,
        const std::function< double( ) > ppnBetaFunction ):
    acceleratedBodies_( acceleratedBodies ), acceleratingBodies_( acceleratingBodies ),
    gravitationalParameterFunction_( gravitationalParameterFunction ),
    bodyStateFunctions_( bodyStateFunctions ),
    ppnGammaFunction_( ppnGammaFunction ),
    ppnBetaFunction_( ppnBetaFunction )
{
    for( unsigned int i = 0; i < acceleratingBodies_.size(); i++ )
    {
        currentGravitationalParameters_.resize( acceleratingBodies_.size() );
        currentPositions_.resize( acceleratingBodies_.size() );
        currentVelocities_.resize( acceleratingBodies_.size() );
        currentSquareSpeeds_.resize( acceleratingBodies_.size() );
        currentLocalPotentials_.resize( acceleratingBodies_.size() );
        totalPointMassAccelerations_.resize( acceleratingBodies_.size() );

        currentRelativePositions_.resize( acceleratedBodies_.size() );
        currentRelativeDistances_.resize( acceleratedBodies_.size() );
        relativePositionVelocityProduct_.resize( acceleratedBodies_.size() );
        velocityInnerProducts_.resize( acceleratedBodies_.size() );
        currentSingleSourceLocalPotential_.resize( acceleratedBodies_.size() );
        singlePointMassAccelerations_.resize( acceleratedBodies_.size() );
        currentSingleAccelerations_.resize( acceleratedBodies_.size() );

        for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
        {
            currentRelativePositions_[ i ].resize( acceleratingBodies_.size( ) );
            currentRelativeDistances_[ i ].resize( acceleratingBodies_.size( ) );
            relativePositionVelocityProduct_[ i ].resize( acceleratingBodies_.size( ) );
            velocityInnerProducts_[ i ].resize( acceleratingBodies_.size( ) );
            currentSingleSourceLocalPotential_[ i ].resize( acceleratingBodies_.size( ) );
            singlePointMassAccelerations_[ i ].resize( acceleratingBodies_.size( ) );
            currentSingleAccelerations_[ i ].resize( acceleratingBodies_.size( ) );

        }

        currentAccelerations_.resize( acceleratedBodies_.size( ) );

    }
    for( unsigned int i = 0; i < acceleratedBodies.size(); i++ )
    {
        acceleratedBodyMap_[ acceleratedBodies.at( i ) ] = i;
    }
    expansionMultipliers_.resize( 10 );
    recomputeExpansionMultipliers( );
}

void EinsteinInfeldHoffmannEquations::recomputeExpansionMultipliers( )
{
    currentPpnGamma_ = ppnGammaFunction_( );
    currentPpnBeta_ = ppnBetaFunction_( );

    expansionMultipliers_[ 0 ] = -2.0 * ( currentPpnGamma_ + currentPpnBeta_ );
    expansionMultipliers_[ 1 ] = -( 2.0 * currentPpnBeta_ + 1.0 );
    expansionMultipliers_[ 2 ] = currentPpnGamma_;
    expansionMultipliers_[ 3 ] = 1.0 + currentPpnGamma_;
    expansionMultipliers_[ 4 ] = -2.0 * ( 1.0 + currentPpnGamma_ );
    expansionMultipliers_[ 5 ] = -1.5;
    expansionMultipliers_[ 6 ] = 0.5;
    expansionMultipliers_[ 7 ] = 2.0 * ( 1.0 + currentPpnGamma_ );
    expansionMultipliers_[ 8 ] = -1.0 + 2.0 * currentPpnGamma_;
    expansionMultipliers_[ 9 ] = ( 3.0 + 4.0 * currentPpnGamma_ ) / 2.0;
}

void EinsteinInfeldHoffmannEquations::update( const double currentTime )
{
    if( currentTime_ != currentTime )
    {

        Eigen::Matrix< double, 6, 1 > currentBodyState;
        for( unsigned int i = 0; i < acceleratingBodies_.size( ); i++ )
        {
            // Extract data from environment
            currentBodyState = bodyStateFunctions_[ i ]( );
            currentGravitationalParameters_[ i ] = gravitationalParameterFunction_[ i ]( );

            // Set local variables
            currentPositions_[ i ] = currentBodyState.segment( 0, 3 );
            currentVelocities_[ i ] = currentBodyState.segment( 3, 3 );
            currentSquareSpeeds_[ i ] = currentVelocities_[ i ].dot( currentVelocities_[ i ] );

        }

        for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
        {
            // Reset value
            totalPointMassAccelerations_[ i ].setZero( );
            currentLocalPotentials_[ i ] = 0.0;

            for( unsigned int j = 0; j < acceleratingBodies_.size( ); j++ )
            {
                if( i != j )
                {

                    currentRelativePositions_[ i ][ j ] = currentPositions_[ j ] - currentPositions_[ i ];
                    currentRelativeDistances_[ i ][ j ] = currentRelativePositions_[ i ][ j ].norm( );
                    relativePositionVelocityProduct_[ i ][ j ] = currentRelativePositions_[ i ][ j ].dot( currentVelocities_[ j ] );
                    velocityInnerProducts_[ i ][ j ] = currentVelocities_[ i ].dot( currentVelocities_[ j ] );

                    currentSingleSourceLocalPotential_[ i ][ j ] = currentGravitationalParameters_[ j ] / currentRelativeDistances_[ i ][ j ];
                    singlePointMassAccelerations_[ i ][ j ] = currentGravitationalParameters_[ j ] * currentRelativePositions_[ i ][ j ] /
                                                              ( currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] );

                    currentLocalPotentials_[ i ] += currentSingleSourceLocalPotential_[ i ][ j ];
                    totalPointMassAccelerations_[ i ] += singlePointMassAccelerations_[ i ][ j ];
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

        for( unsigned int j = 0; j < acceleratingBodies_.size( ); j++ )
        {
            currentSingleAccelerations_[ i ][ j ].setZero( );
            if( i != j )
            {
                currentSingleAccelerations_[ i ][ j ] +=
                    singlePointMassAccelerations_[ i ][ j ]  * ( 1.0 + INVERSE_SQUARE_SPEED_OF_LIGHT *
                        ( expansionMultipliers_[ 0 ] * currentLocalPotentials_[ i ] + expansionMultipliers_[ 1 ] * currentLocalPotentials_[ j ] +
                          expansionMultipliers_[ 2 ] * currentSquareSpeeds_[ i ] + expansionMultipliers_[ 3 ] * currentSquareSpeeds_[ j ] +
                          expansionMultipliers_[ 4 ] * velocityInnerProducts_[ i ][ j ] + expansionMultipliers_[ 5 ] *
                              ( relativePositionVelocityProduct_[ i ][ j ] * relativePositionVelocityProduct_[ i ][ j ] ) /
                              ( currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] ) +
                          expansionMultipliers_[ 6 ] * currentRelativePositions_[ i ][ j ].dot( totalPointMassAccelerations_[ j ] ) ) );

                currentSingleAccelerations_[ i ][ j ] += INVERSE_SQUARE_SPEED_OF_LIGHT * (
                           singlePointMassAccelerations_[ i ][ j ].dot(
                               expansionMultipliers_[ 7 ] * currentVelocities_[ i ] + expansionMultipliers_[ 8 ] * currentVelocities_[ j ] ) *
                               ( currentVelocities_[ i ] - currentVelocities_[ j ] ) +
                               expansionMultipliers_[ 9 ] * currentSingleSourceLocalPotential_[ i ][ j ] * totalPointMassAccelerations_[ j ] );
            }
            currentAccelerations_[ i ] += currentSingleAccelerations_[ i ][ j ];
        }
    }
}

}

}