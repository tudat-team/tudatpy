#include <iostream>
#include <iomanip>

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
    ppnBetaFunction_( ppnBetaFunction ),
    omitMainTerm_( false ),
    currentTime_( TUDAT_NAN )
{
//    if( acceleratedBodies.size( ) > acceleratingBodies.size( ) )
//    {
//        throw std::runtime_error( "Error when creating EIH object, number of undergoing bodies exceed exerting bodies" );
//    }
//
//    for( unsigned int i = 0; i < acceleratingBodies_.size( ); i++ )
//    {
//        if( acceleratedBodies_.at( i ) != acceleratingBodies_.at( i ) )
//        {
//            throw std::runtime_error( "Error when creating EIH object, ordering of bodies is inconsistent" );
//        }
//    }

    scalarEihCorrections_.resize( 7 );
    vectorEihCorrections_.resize( 3 );

    currentGravitationalParameters_.resize( acceleratingBodies_.size() );
    currentPositions_.resize( acceleratingBodies_.size() );
    currentVelocities_.resize( acceleratingBodies_.size() );
    currentSquareSpeeds_.resize( acceleratingBodies_.size() );
    currentLocalPotentials_.resize( acceleratingBodies_.size() );
    totalPointMassAccelerations_.resize( acceleratingBodies_.size() );

    currentRelativePositions_.resize( acceleratingBodies_.size() );
    currentInverseSquareDistances_.resize( acceleratingBodies_.size() );
    currentRelativeDistances_.resize( acceleratingBodies_.size() );
    currentRelativeVelocities_.resize( acceleratingBodies_.size() );

    lineOfSightSpeed_.resize( acceleratingBodies_.size() );
    velocityInnerProducts_.resize( acceleratingBodies_.size() );
    currentSingleSourceLocalPotential_.resize( acceleratingBodies_.size() );
    singlePointMassAccelerations_.resize( acceleratingBodies_.size() );
    currentSingleAccelerations_.resize( acceleratingBodies_.size() );

    totalScalarTermCorrection.resize( acceleratingBodies_.size() );
    totalVectorTermCorrection_.resize( acceleratingBodies_.size() );

    for( int k = 0; k < 7; k++ )
    {
        scalarEihCorrections_[ k ].resize( acceleratingBodies_.size( ) );
    }

    for( int k = 0; k < 3; k++ )
    {
        vectorEihCorrections_[ k ].resize( acceleratingBodies_.size( ) );
    }

    for( unsigned int i = 0; i < acceleratingBodies_.size( ); i++ )
    {
        currentRelativePositions_[ i ].resize( acceleratingBodies_.size( ) );
        currentInverseSquareDistances_[ i ].resize( acceleratingBodies_.size( ) );
        currentRelativeDistances_[ i ].resize( acceleratingBodies_.size( ) );
        currentRelativeVelocities_[ i ].resize( acceleratingBodies_.size( ) );

        lineOfSightSpeed_[ i ].resize( acceleratingBodies_.size( ) );
        velocityInnerProducts_[ i ].resize( acceleratingBodies_.size( ) );
        currentSingleSourceLocalPotential_[ i ].resize( acceleratingBodies_.size( ) );
        singlePointMassAccelerations_[ i ].resize( acceleratingBodies_.size( ) );
        currentSingleAccelerations_[ i ].resize( acceleratingBodies_.size( ) );
        totalScalarTermCorrection[ i ].resize( acceleratingBodies_.size() );
        totalVectorTermCorrection_[ i ].resize( acceleratingBodies_.size() );

        for( int k = 0; k < 7; k++ )
        {
            scalarEihCorrections_[ k ][ i ].resize( acceleratingBodies_.size( ) );
        }

        for( int k = 0; k < 3; k++ )
        {
            vectorEihCorrections_[ k ][ i ].resize( acceleratingBodies_.size( ) );
        }
    }

    currentAccelerations_.resize( acceleratedBodies_.size( ) );


    for( unsigned int i = 0; i < acceleratedBodies.size(); i++ )
    {
        acceleratedBodyMap_[ acceleratedBodies.at( i ) ] = i;
    }

    for( unsigned int i = 0; i < acceleratingBodies_.size(); i++ )
    {
        acceleratingBodyMap_[ acceleratingBodies_.at( i ) ] = i;
    }
    scalarTermMultipliers_.resize( 7 );
    vectorTermMultipliers_.resize( 3 );

    recomputeExpansionMultipliers( );
}

void EinsteinInfeldHoffmannEquations::recomputeExpansionMultipliers( )
{
    currentPpnGamma_ = ppnGammaFunction_( );
    currentPpnBeta_ = ppnBetaFunction_( );

    scalarTermMultipliers_[ 0 ] = -2.0 * ( currentPpnGamma_ + currentPpnBeta_ );
    scalarTermMultipliers_[ 1 ] = - 2.0 * currentPpnBeta_ - 1.0 ;
    scalarTermMultipliers_[ 2 ] = currentPpnGamma_;
    scalarTermMultipliers_[ 3 ] = 1.0 + currentPpnGamma_;
    scalarTermMultipliers_[ 4 ] = -2.0 * ( 1.0 + currentPpnGamma_ );
    scalarTermMultipliers_[ 5 ] = -1.5;
    scalarTermMultipliers_[ 6 ] = 0.5;

    vectorTermMultipliers_[ 0 ] = 2.0 * ( 1.0 + currentPpnGamma_ );
    vectorTermMultipliers_[ 1 ] = -1.0 + 2.0 * currentPpnGamma_;
    vectorTermMultipliers_[ 2 ] = ( 3.0 + 4.0 * currentPpnGamma_ ) / 2.0;
}

void EinsteinInfeldHoffmannEquations::update( const double currentTime )
{
    if( currentTime_ != currentTime )
    {
        currentTime_ = currentTime;
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

        for( unsigned int i = 0; i < acceleratingBodies_.size( ); i++ )
        {
            // Reset value
            currentLocalPotentials_[ i ] = 0.0;
            totalPointMassAccelerations_[ i ].setZero( );

            for( unsigned int j = 0; j < acceleratingBodies_.size( ); j++ )
            {

                // v_{i} * v_{j}
                velocityInnerProducts_[ i ][ j ] = currentVelocities_[ i ].dot( currentVelocities_[ j ] );

                if( i != j )
                {
                    // r_{ij} = r_{j} - r_{i}
                    currentRelativePositions_[ i ][ j ] = currentPositions_[ j ] - currentPositions_[ i ];
                    currentRelativeDistances_[ i ][ j ] = currentRelativePositions_[ i ][ j ].norm( );
                    currentInverseSquareDistances_[ i ][ j ] = 1.0 / ( currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ]  );

                    currentRelativeVelocities_[ i ][ j ] = currentVelocities_[ j ] - currentVelocities_[ i ];

                    // r_{ij} * v_{j}
                    lineOfSightSpeed_[ i ][ j ] = currentRelativePositions_[ i ][ j ].dot( currentVelocities_[ j ] );

                    // mu_j / ||r_{ij}||
                    currentSingleSourceLocalPotential_[ i ][ j ] = currentGravitationalParameters_[ j ] / currentRelativeDistances_[ i ][ j ];

                    // mu_{j} * r_{ij} / ||r_{ij}||^3
                    singlePointMassAccelerations_[ i ][ j ] = currentGravitationalParameters_[ j ] * currentRelativePositions_[ i ][ j ] *
                        currentInverseSquareDistances_[ i ][ j ] / currentRelativeDistances_[ i ][ j ];


                    currentLocalPotentials_[ i ] += currentSingleSourceLocalPotential_[ i ][ j ];
                    totalPointMassAccelerations_[ i ] += singlePointMassAccelerations_[ i ][ j ];
                }
            }
        }

        for( unsigned int i = 0; i < acceleratedBodies_.size( ); i++ )
        {
            for( unsigned int j = 0; j < acceleratingBodies_.size( ); j++ )
            {
                if( i != j )
                {
                    scalarEihCorrections_[ 0 ][ i ][ j ] = scalarTermMultipliers_[ 0 ] * currentLocalPotentials_[ i ];
                    scalarEihCorrections_[ 1 ][ i ][ j ] = scalarTermMultipliers_[ 1 ] * currentLocalPotentials_[ j ];
                    scalarEihCorrections_[ 2 ][ i ][ j ] = scalarTermMultipliers_[ 2 ] * velocityInnerProducts_[ i ][ i ];

                    scalarEihCorrections_[ 3 ][ i ][ j ] = scalarTermMultipliers_[ 3 ] * currentVelocities_[ j ].dot( currentVelocities_[ j ] );
                    scalarEihCorrections_[ 4 ][ i ][ j ] = scalarTermMultipliers_[ 4 ] * velocityInnerProducts_[ i ][ j ];
                    scalarEihCorrections_[ 5 ][ i ][ j ] =
                        scalarTermMultipliers_[ 5 ] * lineOfSightSpeed_[ i ][ j ] * lineOfSightSpeed_[ i ][ j ] /
                         ( currentRelativeDistances_[ i ][ j ] * currentRelativeDistances_[ i ][ j ] );
                    scalarEihCorrections_[ 6 ][ i ][ j ] = scalarTermMultipliers_[ 6 ] * currentRelativePositions_[ i ][ j ].dot( totalPointMassAccelerations_[ j ] );

                    vectorEihCorrections_[ 0 ][ i ][ j ] = vectorTermMultipliers_[ 0 ] * currentRelativePositions_[ i ][ j ].dot( currentVelocities_[ i ] ) * currentInverseSquareDistances_[ i ][ j ] * currentRelativeVelocities_[ i ][ j ];
                    vectorEihCorrections_[ 1 ][ i ][ j ] = vectorTermMultipliers_[ 1 ] * lineOfSightSpeed_[ i ][ j ] * currentInverseSquareDistances_[ i ][ j ] * currentRelativeVelocities_[ i ][ j ];
                    vectorEihCorrections_[ 2 ][ i ][ j ] = vectorTermMultipliers_[ 2 ] * totalPointMassAccelerations_[ j ];
                }
            }
        }

        calculateAccelerations( );
    }
    else
    {
        currentTime_ = currentTime;
    }
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
                totalScalarTermCorrection[ i ][ j ] = 0;
                totalVectorTermCorrection_[ i ][ j ].setZero( );
                for(  int k = 0; k < 7; k++ )
                {
                    totalScalarTermCorrection[ i ][ j ] += scalarEihCorrections_[ k ][ i ][ j ];
                }
                currentSingleAccelerations_[ i ][ j ] = singlePointMassAccelerations_[ i ][ j ] *
                    ( 1.0 + totalScalarTermCorrection[ i ][ j ] * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT );

                for(  int k = 0; k < 3; k++ )
                {
                    totalVectorTermCorrection_[ i ][ j ] += vectorEihCorrections_[ k ][ i ][ j ];
                }
                currentSingleAccelerations_[ i ][ j ] += currentSingleSourceLocalPotential_[ i ][ j ] *
                    totalVectorTermCorrection_[ i ][ j ] * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;

            }
            currentAccelerations_[ i ] += currentSingleAccelerations_[ i ][ j ];
        }
    }
}

}

}