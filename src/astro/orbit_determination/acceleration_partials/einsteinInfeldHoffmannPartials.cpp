/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/einsteinInfeldHoffmannPartials.h"

namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix< double, 1, 3 > calculatePartialOfPointMassPotentialWrtBodyPosition(
    const Eigen::Vector3d& relativePosition,
    const double gravitationalParameter )
{
    double relativePositionNorm = relativePosition.norm( );
    return -gravitationalParameter * relativePosition.transpose( ) /
           ( relativePositionNorm * relativePositionNorm * relativePositionNorm );
}

EihEquationsPartials::EihEquationsPartials(
    const std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations ):
    eihEquations_( eihEquations ),
    numberOfExertingBodies_( eihEquations_->getBodiesExertingAcceleration( ).size( ) ),
    numberOfUndergoingBodies_( eihEquations_->getBodiesUndergoingAcceleration( ).size( ) ),
    currentTime_( TUDAT_NAN )
{
    currentTotalAccelerationsWrtPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    currentTotalAccelerationsWrtVelocity_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );

    currentSinglePointMassAccelerationsWrtExertingPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    currentTotalPointMassAccelerationsWrtPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    currentLocalPotentialsWrtPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    currentTotalPotentialWrtPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );

    currentTotalScalarTermWrtExertingPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    currentTotalScalarTermWrtUndergoingPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    currentTotalVectorTermWrtExertingPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    currentTotalVectorTermWrtUndergoingPosition_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );

    currentTotalScalarTermWrtExertingVelocity_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    currentTotalScalarTermWrtUndergoingVelocity_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 1, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 1, 3 >::Zero( ) );
    currentTotalVectorTermWrtExertingVelocity_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
    currentTotalVectorTermWrtUndergoingVelocity_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );

    currentTotalAccelerationsWrtPositionCrossTerms_ = utilities::getTwoDimensionalVector< Eigen::Matrix< double, 3, 3 > >(
        numberOfExertingBodies_, numberOfExertingBodies_, Eigen::Matrix< double, 3, 3 >::Zero( ) );
}

void EihEquationsPartials::update( const double currentTime )
{
    if( currentTime_ != currentTime )
    {
        eihEquations_->update( currentTime );
        for( int i = 0; i < numberOfExertingBodies_; i++ )
        {
            currentTotalPointMassAccelerationsWrtPosition_[ i ][ i ].setZero( );
            currentTotalPotentialWrtPosition_[ i ][ i ].setZero( );
            for ( int j = 0; j < numberOfExertingBodies_; j++ )
            {
                if( i != j )
                {
                    currentSinglePointMassAccelerationsWrtExertingPosition_[ i ][ j ] = -calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        eihEquations_->getRelativePositions( i, j ), eihEquations_->getGravitationalParameter( j ) );
                    currentTotalPointMassAccelerationsWrtPosition_[ i ][ j ] = currentSinglePointMassAccelerationsWrtExertingPosition_[ i ][ j ];
                    currentTotalPointMassAccelerationsWrtPosition_[ i ][ i ] -= currentTotalPointMassAccelerationsWrtPosition_[ i ][ j ];

                    currentLocalPotentialsWrtPosition_[ i ][ j ] = calculatePartialOfPointMassPotentialWrtBodyPosition(
                        eihEquations_->getRelativePositions( i, j ), eihEquations_->getGravitationalParameter( j ) );
                    currentTotalPotentialWrtPosition_[ i ][ j ] = currentLocalPotentialsWrtPosition_[ i ][ j ];
                    currentTotalPotentialWrtPosition_[ i ][ i ] -= currentLocalPotentialsWrtPosition_[ i ][ j ];
                }

                currentTotalScalarTermWrtExertingPosition_[ i ][ j ].setZero( );
                currentTotalScalarTermWrtUndergoingPosition_[ i ][ j ].setZero( );

                currentTotalVectorTermWrtExertingPosition_[ i ][ j ].setZero( );
                currentTotalVectorTermWrtUndergoingPosition_[ i ][ j ].setZero( );

                currentTotalScalarTermWrtExertingVelocity_[ i ][ j ].setZero( );
                currentTotalScalarTermWrtUndergoingVelocity_[ i ][ j ].setZero( );

                currentTotalVectorTermWrtExertingVelocity_[ i ][ j ].setZero( );
                currentTotalVectorTermWrtUndergoingVelocity_[ i ][ j ].setZero( );


                for( int k = 0; k < 7; k++ )
                {
                    addSingleScalarTermWrtPositionPartial(
                        currentTotalScalarTermWrtExertingPosition_[ i ][ j ], i, j, true, k );
                    addSingleScalarTermWrtPositionPartial(
                        currentTotalScalarTermWrtUndergoingPosition_[ i ][ j ], i, j, false, k );

                    addSingleScalarTermWrtVelocityPartial(
                        currentTotalScalarTermWrtExertingVelocity_[ i ][ j ], i, j, true, k );
                    addSingleScalarTermWrtVelocityPartial(
                        currentTotalScalarTermWrtUndergoingVelocity_[ i ][ j ], i, j, false, k );

                    if( k < 3 )
                    {
                        addSingleVectorTermWrtPositionPartial(
                            currentTotalVectorTermWrtExertingPosition_[ i ][ j ], i, j, true, k );
                        addSingleVectorTermWrtPositionPartial(
                            currentTotalVectorTermWrtUndergoingPosition_[ i ][ j ], i, j, false, k );

                        addSingleVectorTermWrtVelocityPartial(
                            currentTotalVectorTermWrtExertingVelocity_[ i ][ j ], i, j, true, k );
                        addSingleVectorTermWrtVelocityPartial(
                            currentTotalVectorTermWrtUndergoingVelocity_[ i ][ j ], i, j, false, k );

                    }
                }
            }
        }

        for( int i = 0; i < numberOfExertingBodies_; i++ )
        {
            currentTotalAccelerationsWrtPosition_[ i ][ i ].setZero( );
            currentTotalAccelerationsWrtVelocity_[ i ][ i ].setZero( );

            for( int m = 0; m < numberOfExertingBodies_; m++ )
            {
                currentTotalAccelerationsWrtPositionCrossTerms_[ i ][ m ].setZero( );
                if( i != m )
                {
                    for( int j = 0; j < numberOfExertingBodies_; j++ )
                    {
                        if( ( j != i ) && ( j != m ) )
                        {
                            currentTotalAccelerationsWrtPositionCrossTerms_[ i ][ m ] +=
                                physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                ( eihEquations_->getSinglePointMassAccelerations( i, j ) * (
                                    getSingleScalarCrossTermWrtPositionPartial( i, j, m, 0 ) +
                                    getSingleScalarCrossTermWrtPositionPartial( i, j, m, 1 ) +
                                    getSingleScalarCrossTermWrtPositionPartial( i, j, m, 6 ) ) +
                                  eihEquations_->getSingleSourceLocalPotential( i, j ) *
                                  getSingleVectorCrossTermWrtPositionPartial( i, j, m, 2 ) );
                        }
                    }


                    currentTotalAccelerationsWrtPosition_[ i ][ m ] =
                        currentSinglePointMassAccelerationsWrtExertingPosition_[ i ][ m ] *
                        ( 1.0 + physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getTotalScalarTermCorrection( i, m ) ) +
                        physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                        eihEquations_->getTotalVectorTermCorrection( i, m ) * currentLocalPotentialsWrtPosition_.at( i ).at( m );

                    currentTotalAccelerationsWrtPosition_[ i ][ i ] -= currentTotalAccelerationsWrtPosition_[ i ][ m ];

                    currentTotalAccelerationsWrtPosition_[ i ][ m ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                       ( eihEquations_->getSinglePointMassAccelerations( i, m ) * currentTotalScalarTermWrtExertingPosition_.at( i ).at( m ) +
                                                                         eihEquations_->getSingleSourceLocalPotential( i, m ) * currentTotalVectorTermWrtExertingPosition_.at( i ).at( m ) );

                    currentTotalAccelerationsWrtPosition_[ i ][ i ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                       ( eihEquations_->getSinglePointMassAccelerations( i, m ) * currentTotalScalarTermWrtUndergoingPosition_.at( i ).at( m ) +
                                                                         eihEquations_->getSingleSourceLocalPotential( i, m ) * currentTotalVectorTermWrtUndergoingPosition_.at( i ).at( m ) );

                    currentTotalAccelerationsWrtPosition_[ i ][ m ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                       currentTotalAccelerationsWrtPositionCrossTerms_[ i ][ m ];

                    currentTotalAccelerationsWrtVelocity_[ i ][ m ] = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                      ( eihEquations_->getSinglePointMassAccelerations( i, m ) * currentTotalScalarTermWrtExertingVelocity_.at( i ).at( m ) +
                                                                        eihEquations_->getSingleSourceLocalPotential( i, m ) * currentTotalVectorTermWrtExertingVelocity_.at( i ).at( m ) );

                    currentTotalAccelerationsWrtVelocity_[ i ][ i ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                       ( eihEquations_->getSinglePointMassAccelerations( i, m ) * currentTotalScalarTermWrtUndergoingVelocity_.at( i ).at( m ) +
                                                                         eihEquations_->getSingleSourceLocalPotential( i, m ) * currentTotalVectorTermWrtUndergoingVelocity_.at( i ).at( m ) );


                }
            }
        }
    }

    currentTime_ = currentTime;
}

void EihEquationsPartials::addSingleScalarTermWrtPositionPartial(
    Eigen::Matrix< double, 1, 3 >& scalarTermWrtPosition, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex )
{
    switch( termIndex )
    {
    case 0:
        if( wrtExerting )
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyUndergoing ][ bodyExerting ];
        }
        else
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyUndergoing ][ bodyUndergoing ];
        }
        break;
    case 1:
        if( wrtExerting )
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyExerting ][ bodyExerting ];
        }
        else
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyExerting ][ bodyUndergoing ];
        }
        break;
    case 2:
        break;
    case 3:
        break;
    case 4:
        break;
    case 5:
    {
        Eigen::Vector3d normalizedRelativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ) /
                                                     eihEquations_->getRelativeDistance( bodyUndergoing, bodyExerting );
        Eigen::Vector3d exertingVelocity = eihEquations_->getVelocity( bodyExerting );
        if( wrtExerting )
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * 2.0 * normalizedRelativePosition.dot( exertingVelocity ) * exertingVelocity.transpose( ) *
                                     ( Eigen::Matrix3d::Identity( ) - normalizedRelativePosition * normalizedRelativePosition.transpose( ) ) /
                                     eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).norm( );
        }
        else
        {
            scalarTermWrtPosition -= eihEquations_->getScalarTermMultiplier( termIndex ) * 2.0 * normalizedRelativePosition.dot( exertingVelocity ) * exertingVelocity.transpose( ) *
                                     ( Eigen::Matrix3d::Identity( ) - normalizedRelativePosition * normalizedRelativePosition.transpose( ) ) /
                                     eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).norm( );
        }
        break;
    }
    case 6:
        if( wrtExerting )
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * (
                eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) *
                currentTotalPointMassAccelerationsWrtPosition_[ bodyExerting ][ bodyExerting ] +
                eihEquations_->getTotalPointMassAcceleration( bodyExerting ).transpose( ) );
        }
        else
        {
            scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * (
                eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) *
                currentTotalPointMassAccelerationsWrtPosition_[ bodyExerting ][ bodyUndergoing ] -
                eihEquations_->getTotalPointMassAcceleration( bodyExerting ).transpose( ) );
        }
        break;
    default:
        throw std::runtime_error( "Error when getting EIH scalar term partial w.r.t. positon, index " + std::to_string( termIndex ) + " not allowed." );

    }
}

Eigen::Matrix< double, 1, 3 > EihEquationsPartials::getSingleScalarCrossTermWrtPositionPartial(
    const int bodyUndergoing, const int bodyExerting, const int bodyPartial, const int termIndex )
{
    Eigen::Matrix< double, 1, 3 > scalarTermWrtPosition;
    switch( termIndex )
    {
    case 0:
        scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyUndergoing ][ bodyPartial ];
        break;
    case 1:
        scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * currentTotalPotentialWrtPosition_[ bodyExerting ][ bodyPartial ];
        break;

    case 6:
        scalarTermWrtPosition += eihEquations_->getScalarTermMultiplier( termIndex ) * eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) *
                                 currentTotalPointMassAccelerationsWrtPosition_[ bodyExerting ][ bodyPartial ];
        break;
    default:
        throw std::runtime_error( "Error when getting EIH scalar cross term partial w.r.t. positon, index " + std::to_string( termIndex ) + " not allowed." );

    }
    return scalarTermWrtPosition;
}

void EihEquationsPartials::addSingleScalarTermWrtVelocityPartial(
    Eigen::Matrix< double, 1, 3 >& scalarTermWrtVelocity, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex )
{
    switch( termIndex )
    {
    case 0:
        break;
    case 1:
        break;
    case 2:
        if( !wrtExerting )
        {
            scalarTermWrtVelocity += 2.0 * eihEquations_->getScalarTermMultiplier( termIndex ) * eihEquations_->getVelocity( bodyUndergoing ).transpose( );
        }
        break;
    case 3:
        if( wrtExerting )
        {
            scalarTermWrtVelocity += eihEquations_->getScalarTermMultiplier( termIndex ) * 2.0 * eihEquations_->getVelocity( bodyExerting ).transpose( );
        }
        break;
    case 4:
        if( wrtExerting )
        {
            scalarTermWrtVelocity += eihEquations_->getScalarTermMultiplier( termIndex ) * eihEquations_->getVelocity( bodyUndergoing ).transpose( );
        }
        else
        {
            scalarTermWrtVelocity += eihEquations_->getScalarTermMultiplier( termIndex ) * eihEquations_->getVelocity( bodyExerting ).transpose( );
        }
        break;
    case 5:
    {
        if( wrtExerting )
        {
            Eigen::Vector3d normalizedRelativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ) /
                                                         eihEquations_->getRelativeDistance( bodyUndergoing, bodyExerting );
            Eigen::Vector3d exertingVelocity = eihEquations_->getVelocity( bodyExerting );
            scalarTermWrtVelocity += eihEquations_->getScalarTermMultiplier( termIndex ) *
                                     2.0 * normalizedRelativePosition.dot( exertingVelocity ) * normalizedRelativePosition.transpose( );
        }
        break;
    }
    case 6:
        break;
    default:
        throw std::runtime_error( "Error when getting EIH scalar term partial w.r.t. velocity, index " + std::to_string( termIndex ) + " not allowed." );

    }
}

void EihEquationsPartials::addSingleVectorTermWrtPositionPartial(
    Eigen::Matrix3d& vectorTermWrtPosition, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex )
{

    switch( termIndex )
    {
    case 0:
    {
        double inverseSquareDistance = eihEquations_->getInverseSquareDistance( bodyUndergoing, bodyExerting );
        Eigen::Vector3d relativePositionNormalized = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ) /
                                                     eihEquations_->getRelativeDistance( bodyUndergoing, bodyExerting );


        double sign = ( wrtExerting ? 1.0 : -1.0 );
        vectorTermWrtPosition +=
            eihEquations_->getVectorTermMultiplier( termIndex ) * sign * eihEquations_->getRelativeVelocity( bodyUndergoing, bodyExerting ) *
            ( eihEquations_->getVelocity( bodyUndergoing ).transpose( ) * inverseSquareDistance *
              ( Eigen::Matrix3d::Identity( ) - 2.0 * relativePositionNormalized * relativePositionNormalized.transpose( ) ) );
        break;
    }
    case 1:
    {
        double inverseSquareDistance = eihEquations_->getInverseSquareDistance( bodyUndergoing, bodyExerting );
        Eigen::Vector3d relativePositionNormalized = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ) /
                                                     eihEquations_->getRelativeDistance( bodyUndergoing, bodyExerting );


        double sign = ( wrtExerting ? 1.0 : -1.0 );
        vectorTermWrtPosition +=
            eihEquations_->getVectorTermMultiplier( termIndex ) * sign * eihEquations_->getRelativeVelocity( bodyUndergoing, bodyExerting ) *
            ( eihEquations_->getVelocity( bodyExerting ) .transpose( ) * inverseSquareDistance *
              ( Eigen::Matrix3d::Identity( ) - 2.0 * relativePositionNormalized * relativePositionNormalized.transpose( ) ) );
        break;
    }
    case 2:
        if( wrtExerting )
        {
            vectorTermWrtPosition += eihEquations_->getVectorTermMultiplier( termIndex ) * currentTotalPointMassAccelerationsWrtPosition_.at( bodyExerting ).at( bodyExerting );
        }
        else
        {
            vectorTermWrtPosition += eihEquations_->getVectorTermMultiplier( termIndex ) * currentTotalPointMassAccelerationsWrtPosition_.at( bodyExerting ).at( bodyUndergoing );
        }
        break;
    default:
        throw std::runtime_error( "Error when getting EIH vector term partial w.r.t. position, index " + std::to_string( termIndex ) + " not allowed." );

    }
}

Eigen::Matrix< double, 3, 3 > EihEquationsPartials::getSingleVectorCrossTermWrtPositionPartial(
    const int bodyUndergoing, const int bodyExerting, const int bodyPartial, const int termIndex )
{
    Eigen::Matrix< double, 3, 3 > vectorTermWrtPosition;
    switch( termIndex )
    {

    case 2:
        vectorTermWrtPosition += eihEquations_->getVectorTermMultiplier( termIndex ) * currentTotalPointMassAccelerationsWrtPosition_.at( bodyUndergoing ).at( bodyPartial );
        break;
    default:
        throw std::runtime_error( "Error when getting EIH vector cross term partial w.r.t. positon, index " + std::to_string( termIndex ) + " not allowed." );

    }
    return vectorTermWrtPosition;
}

void EihEquationsPartials::addSingleVectorTermWrtVelocityPartial(
    Eigen::Matrix3d& vectorTermWrtVelocity, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex )
{
    switch( termIndex )
    {
    case 0:
    {
        double inverseSquareDistance = eihEquations_->getInverseSquareDistance(
            bodyUndergoing, bodyExerting );
        double baseMultiplier = eihEquations_->getLineOfSighSpeed( bodyExerting, bodyUndergoing ) * inverseSquareDistance;
        {
            if( wrtExerting )
            {
                vectorTermWrtVelocity -= eihEquations_->getVectorTermMultiplier( termIndex ) *baseMultiplier * Eigen::Matrix3d::Identity( );

            }
            else
            {
                vectorTermWrtVelocity += eihEquations_->getVectorTermMultiplier( termIndex ) * ( baseMultiplier * Eigen::Matrix3d::Identity( ) + eihEquations_->getRelativeVelocity(
                    bodyUndergoing, bodyExerting ) * eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) * inverseSquareDistance );
            }
        }
        break;
    }
    case 1:
    {
        double inverseSquareDistance = eihEquations_->getInverseSquareDistance(
            bodyUndergoing, bodyExerting );
        double baseMultiplier = eihEquations_->getLineOfSighSpeed( bodyUndergoing, bodyExerting ) * inverseSquareDistance;
        {
            if( wrtExerting )
            {
                vectorTermWrtVelocity += eihEquations_->getVectorTermMultiplier( termIndex ) * ( baseMultiplier * Eigen::Matrix3d::Identity( ) + eihEquations_->getRelativeVelocity(
                    bodyUndergoing, bodyExerting ) * eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) * inverseSquareDistance );

            }
            else
            {
                vectorTermWrtVelocity -= eihEquations_->getVectorTermMultiplier( termIndex ) * baseMultiplier * Eigen::Matrix3d::Identity( );

            }
        }
        break;
    }
    case 2:
        break;
    default:
        throw std::runtime_error( "Error when getting EIH vector term partial w.r.t. position, index " + std::to_string( termIndex ) + " not allowed." );

    }
}

void EihEquationsPartials::getAccelerationWrtGravitationalParameter( Eigen::MatrixXd& muPartial, const int bodyIndex, const int muIndex )
{
    muPartial.setZero( 3, 1 );
    double currentMu = eihEquations_->getGravitationalParameter( muIndex     );
    if( currentMu == 0.0 )
    {
        throw std::runtime_error( "Error when computing EIH partial w.r.t. mu, value of mu is 0" );
    }

    if( bodyIndex != muIndex )
    {
        muPartial += eihEquations_->getSinglePointMassAccelerations( bodyIndex, muIndex )  *
            ( 1.0 + physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getTotalScalarTermCorrection( bodyIndex, muIndex ) );
        muPartial += eihEquations_->getSingleSourceLocalPotential( bodyIndex, muIndex ) *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getTotalVectorTermCorrection( bodyIndex, muIndex );
    }

    double scalarMultiplier0 = eihEquations_->getScalarTermMultiplier( 0 );
    double scalarMultiplier1 = eihEquations_->getScalarTermMultiplier( 1 );
    double scalarMultiplier6 = eihEquations_->getScalarTermMultiplier( 6 );
    double vectorMultiplier2 = eihEquations_->getVectorTermMultiplier( 2 );

    double scalarTerm0Partial = 0.0;
    if( bodyIndex != muIndex )
    {
        scalarTerm0Partial = eihEquations_->getSingleSourceLocalPotential( bodyIndex, muIndex );
    }

    for( int j = 0; j < numberOfExertingBodies_; j++ )
    {
        if ( j != bodyIndex && j != muIndex )
        {
            muPartial += (
                eihEquations_->getSinglePointMassAccelerations( bodyIndex, j ) * (
                scalarMultiplier0 * scalarTerm0Partial +
                scalarMultiplier1 * eihEquations_->getSingleSourceLocalPotential( j, muIndex ) +
                scalarMultiplier6 * eihEquations_->getRelativePositions( bodyIndex, j ).dot( eihEquations_->getSinglePointMassAccelerations( j, muIndex ) )
                ) +
                vectorMultiplier2 * eihEquations_->getSingleSourceLocalPotential( bodyIndex, j ) * eihEquations_->getSinglePointMassAccelerations( j, muIndex )
                                                                                           ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
        }
    }
    muPartial /= currentMu;
}



void EihEquationsPartials::getAccelerationWrtGamma( Eigen::MatrixXd& gammaPartial, const int bodyIndex )
{
    gammaPartial.setZero( 3, 1 );
    for( int j = 0; j < numberOfExertingBodies_; j++ )
    {
        if( j != bodyIndex )
        {
            double accelerationMultiplier = 0.0;
            accelerationMultiplier += -2.0 * eihEquations_->getScalarEihCorrection( 0, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 0 ) ;
            accelerationMultiplier += eihEquations_->getScalarEihCorrection( 2, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 2 ) ;
            accelerationMultiplier += eihEquations_->getScalarEihCorrection( 3, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 3 ) ;
            accelerationMultiplier += -2.0 * eihEquations_->getScalarEihCorrection( 4, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 4 ) ;

            Eigen::Vector3d potentialMultiplier = Eigen::Vector3d::Zero( );
            potentialMultiplier += 2.0 * eihEquations_->getVectorEihCorrection( 0, bodyIndex, j ) / eihEquations_->getVectorTermMultiplier( 0 );
            potentialMultiplier += 2.0 * eihEquations_->getVectorEihCorrection( 1, bodyIndex, j ) / eihEquations_->getVectorTermMultiplier( 1 );
            potentialMultiplier += 2.0 * eihEquations_->getVectorEihCorrection( 2, bodyIndex, j ) / eihEquations_->getVectorTermMultiplier( 2 );

            gammaPartial += accelerationMultiplier * eihEquations_->getSinglePointMassAccelerations( bodyIndex, j ) +
                            potentialMultiplier * eihEquations_->getSingleSourceLocalPotential( bodyIndex, j );

        }
    }
    gammaPartial *= physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

void EihEquationsPartials::getAccelerationWrtBeta( Eigen::MatrixXd& betaPartial, const int bodyIndex )
{
    betaPartial.setZero( 3, 1 );
    for( int j = 0; j < numberOfExertingBodies_; j++ )
    {
        if( j != bodyIndex )
        {
            double accelerationMultiplier = 0;
            accelerationMultiplier += -2.0 * eihEquations_->getScalarEihCorrection( 0, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 0 );
            accelerationMultiplier += -2.0 * eihEquations_->getScalarEihCorrection( 1, bodyIndex, j ) / eihEquations_->getScalarTermMultiplier( 1 );

            betaPartial += accelerationMultiplier * eihEquations_->getSinglePointMassAccelerations( bodyIndex, j );

        }
    }
    betaPartial *= physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}


} // namespace acceleration_partials

} // namespace tudat

