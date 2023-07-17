/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EIHPARTIALS_H
#define TUDAT_EIHPARTIALS_H

#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"

#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

inline Eigen::Matrix< double, 1, 3 > calculatePartialOfPointMassPotentialWrtBodyPosition(
    const Eigen::Vector3d& relativePosition,
    const double gravitationalParameter )
{
    double relativePositionNorm = relativePosition.norm( );
    return -gravitationalParameter * relativePosition.transpose( ) /
        ( relativePositionNorm * relativePositionNorm * relativePositionNorm );
}

//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class EihEquationsPartials
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    EihEquationsPartials(
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
//        currentTotalAccelerationsWrtGravitationalParameter_.resize( numberOfUndergoingBodies_ );

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

    void update( const double currentTime )
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

    void addSingleScalarTermWrtPositionPartial(
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

    Eigen::Matrix< double, 1, 3 > getSingleScalarCrossTermWrtPositionPartial(
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

    void addSingleScalarTermWrtVelocityPartial(
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

    void addSingleVectorTermWrtPositionPartial(
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

    Eigen::Matrix< double, 3, 3 > getSingleVectorCrossTermWrtPositionPartial(
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

    void addSingleVectorTermWrtVelocityPartial(
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

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalAccelerationsWrtPosition( )
    {
        return currentTotalAccelerationsWrtPosition_;
    }

    Eigen::Matrix3d getCurrentTotalAccelerationWrtPosition( const int bodyUndergoing, const int bodyExerting )
    {
        return currentTotalAccelerationsWrtPosition_.at( bodyUndergoing ).at( bodyExerting );
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalAccelerationsWrtVelocity( )
    {
        return currentTotalAccelerationsWrtVelocity_;
    }

    Eigen::Matrix3d getCurrentTotalAccelerationWrtVelocity( const int bodyUndergoing, const int bodyExerting )
    {
        return currentTotalAccelerationsWrtVelocity_.at( bodyUndergoing ).at( bodyExerting );
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > getCurrentTotalPotentialWrtPosition( )
    {
        return currentTotalPotentialWrtPosition_;
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > getCurrentLocalPotentialWrtPosition( )
    {
        return currentLocalPotentialsWrtPosition_;
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentSinglePointMassAccelerationWrtPosition( )
    {
        return currentSinglePointMassAccelerationsWrtExertingPosition_;
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalPointMassAccelerationWrtPosition( )
    {
        return currentTotalPointMassAccelerationsWrtPosition_;
    }

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > getEihEquations( )
    {
        return eihEquations_;
    }


    void getAccelerationWrtGamma( Eigen::MatrixXd& gammaPartial, const int bodyIndex )
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

    void getAccelerationWrtBeta( Eigen::MatrixXd& betaPartial, const int bodyIndex )
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

    void getAccelerationWrtGravitationalParameter( Eigen::MatrixXd& muPartial, const int bodyIndex, const int muIndex );

protected:

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations_;

    int numberOfExertingBodies_;

    int numberOfUndergoingBodies_;

    double currentTime_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtPosition_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtPositionCrossTerms_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtVelocity_;


    std::vector< std::vector< double > > currentTotalAccelerationsWrtGravitationalParameter_;




    std::vector< std::vector< Eigen::Matrix3d > > currentSinglePointMassAccelerationsWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalPointMassAccelerationsWrtPosition_;

//    std::vector< std::vector< Eigen::Vector3d > > currentSingleAccelerationsWrtGravitationalParameter_;



    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentTotalPotentialWrtPosition_;


//    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtGravitationalParameter_;


    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtUndergoingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtUndergoingPosition_;


    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtExertingVelocity_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtUndergoingVelocity_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtExertingVelocity_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtUndergoingVelocity_;



};


//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class EihAccelerationPartial: public AccelerationPartial
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    EihAccelerationPartial(
        const std::shared_ptr< acceleration_partials::EihEquationsPartials > fulEihPartials,
        const std::string acceleratedBody ):
        AccelerationPartial( acceleratedBody, "", basic_astrodynamics::einstein_infeld_hoffmann_acceleration ),
        fullEihPartials_( fulEihPartials )
    {
        std::vector< std::string > bodyList = fullEihPartials_->getEihEquations( )->getBodiesExertingAcceleration( );
        for( unsigned int i = 0; i < bodyList.size( ); i++ )
        {
            bodyIndices_[ bodyList.at( i ) ] = i;
        }
        acceleratedBodyIndex_ = bodyIndices_.at( acceleratedBody );
    }


    void wrtPositionOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
    }


    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        throw std::runtime_error( "Error when calculating EIH partial w.r.t. body exerting acceleration, no such single body exists" );
    }

    virtual void wrtPositionOfAdditionalBody(
        const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        int additionalBodyIndex = bodyIndices_.at( bodyName );
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
    }


    void wrtVelocityOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
    }


    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        throw std::runtime_error( "Error when calculating EIH partial w.r.t. body exerting acceleration, no such single body exists" );
    }

    virtual void wrtVelocityOfAdditionalBody(
        const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        int additionalBodyIndex = bodyIndices_.at( bodyName );
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
    }



    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
    {
        return 0;
    }


    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

        // Check dependencies.
        if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter )
        {
//            std::cout<<"Checking "<<parameter->getParameterName( ).second.first<<" "
//            <<fullEihPartials_->getEihEquations( )->getAcceleratingBodyMap( ).count( parameter->getParameterName( ).second.first )<<std::endl;
            if( fullEihPartials_->getEihEquations( )->getAcceleratingBodyMap( ).count( parameter->getParameterName( ).second.first ) > 0 )
            {
                // If parameter is gravitational parameter, check and create dependency function .
                partialFunctionPair =
                    std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtGravitationalParameter, fullEihPartials_,
                                                 std::placeholders::_1, acceleratedBodyIndex_,
                                                 fullEihPartials_->getEihEquations( )->getAcceleratingBodyMap( ).at( parameter->getParameterName( ).second.first ) ), 1 );
            }
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::ppn_parameter_gamma )
        {
            partialFunctionPair =
                std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtGamma, fullEihPartials_,
                                           std::placeholders::_1, acceleratedBodyIndex_ ), 1 );
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::ppn_parameter_beta )
        {
            partialFunctionPair =
                std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtBeta, fullEihPartials_,
                                           std::placeholders::_1, acceleratedBodyIndex_ ), 1 );
        }
        else
        {
            partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
        }

        return partialFunctionPair;
    }


    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }


    void update( const double currentTime = TUDAT_NAN )
    {
        fullEihPartials_->update( currentTime );
        currentTime_ = currentTime;
    }

    const std::shared_ptr< acceleration_partials::EihEquationsPartials > getFulEihPartials( )
    {
        return fullEihPartials_;
    }


protected:

    const std::shared_ptr< acceleration_partials::EihEquationsPartials > fullEihPartials_;

    std::map< std::string, int > bodyIndices_;

    int acceleratedBodyIndex_;
};


} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_EIHPARTIALS_H
