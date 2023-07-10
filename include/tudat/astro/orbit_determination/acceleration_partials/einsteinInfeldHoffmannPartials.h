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
            numberOfUndergoingBodies_( eihEquations_->getBodiesUndergoingAcceleration( ).size( ) )
    {
        currentTotalAccelerationsWrtPosition_.resize( numberOfUndergoingBodies_ );
        currentTotalAccelerationsWrtVelocity_.resize( numberOfUndergoingBodies_ );
        currentTotalAccelerationsWrtGravitationalParameter_.resize( numberOfUndergoingBodies_ );

        currentSinglePointMassAccelerationsWrtExertingPosition_.resize( numberOfExertingBodies_ );
        currentTotalPointMassAccelerationsWrtPosition_.resize( numberOfExertingBodies_ );
        currentLocalPotentialsWrtPosition_.resize( numberOfExertingBodies_ );
        currentTotalPotentialWrtPosition_.resize( numberOfExertingBodies_ );

        currentTotalScalarTermWrtExertingPosition_.resize( numberOfExertingBodies_ );
        currentTotalScalarTermWrtUndergoingPosition_.resize( numberOfExertingBodies_ );
        currentTotalVectorTermWrtExertingPosition_.resize( numberOfExertingBodies_ );
        currentTotalVectorTermWrtUndergoingPosition_.resize( numberOfExertingBodies_ );

        for( int i = 0; i < numberOfUndergoingBodies_; i++ )
        {
            currentTotalAccelerationsWrtPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalAccelerationsWrtVelocity_[ i ].resize( numberOfExertingBodies_ );
            currentTotalAccelerationsWrtGravitationalParameter_[ i ].resize( numberOfExertingBodies_ );
        }

        for( int i = 0; i < numberOfExertingBodies_; i++ )
        {
            currentSinglePointMassAccelerationsWrtExertingPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalPointMassAccelerationsWrtPosition_[ i ].resize( numberOfExertingBodies_ );
            currentLocalPotentialsWrtPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalPotentialWrtPosition_[ i ].resize( numberOfExertingBodies_ );

            currentTotalScalarTermWrtExertingPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalScalarTermWrtUndergoingPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalVectorTermWrtExertingPosition_[ i ].resize( numberOfExertingBodies_ );
            currentTotalVectorTermWrtUndergoingPosition_[ i ].resize( numberOfExertingBodies_ );
        }


    }

    void update( const double currentTime )
    {
        eihEquations_->update( currentTime );
        for( int i = 0; i < numberOfUndergoingBodies_; i++ )
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

                for( int k = 0; k < 7; k++ )
                {
                    addSingleScalarTermWrtPositionPartial(
                        currentTotalScalarTermWrtExertingPosition_[ i ][ j ], i, j, true, k );
                    addSingleScalarTermWrtPositionPartial(
                        currentTotalScalarTermWrtUndergoingPosition_[ i ][ j ], i, j, false, k );

                    if( k < 3 )
                    {
                        addSingleVectorTermWrtPositionPartial(
                            currentTotalVectorTermWrtExertingPosition_[ i ][ j ], i, j, true, k );
                        addSingleVectorTermWrtPositionPartial(
                            currentTotalVectorTermWrtUndergoingPosition_[ i ][ j ], i, j, false, k );

                    }
                }
            }
        }

        for( int i = 0; i < numberOfUndergoingBodies_; i++ )
        {
            currentTotalAccelerationsWrtPosition_[ i ][ i ].setZero( );
            for( int j = 0; j < numberOfExertingBodies_; j++ )
            {
                if( i != j )
                {
                    currentTotalAccelerationsWrtPosition_[ i ][ j ] =
                        currentSinglePointMassAccelerationsWrtExertingPosition_[ i ][ j ] *
                        ( 1.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getScalarTermMultiplier( i, j ) ) +
                        physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                            eihEquations_->getVectorTermMultiplier( i, j ) * currentLocalPotentialsWrtPosition_.at( i ).at( j );

                    currentTotalAccelerationsWrtPosition_[ i ][ i ] -= currentTotalAccelerationsWrtPosition_[ i ][ j ];

                    currentTotalAccelerationsWrtPosition_[ i ][ j ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                        ( eihEquations_->getSinglePointMassAccelerations( i, j ) * currentTotalScalarTermWrtExertingPosition_.at( i ).at( j ) +
                          eihEquations_->getSingleSourceLocalPotential( i, j ) * currentTotalVectorTermWrtExertingPosition_.at( i ).at( j ) );

                    currentTotalAccelerationsWrtPosition_[ i ][ i ] += physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                       ( eihEquations_->getSinglePointMassAccelerations( i, j ) * currentTotalScalarTermWrtUndergoingPosition_.at( i ).at( j ) +
                         eihEquations_->getSingleSourceLocalPotential( i, j ) * currentTotalVectorTermWrtUndergoingPosition_.at( i ).at( j ) );
                }
            }
        }
    }

    void addSingleScalarTermWrtPositionPartial(
        Eigen::Matrix< double, 1, 3 >& scalarTermWrtPosition, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex )
    {
        switch( termIndex )
        {
        case 0:
            if( wrtExerting )
            {
                scalarTermWrtPosition += currentTotalPotentialWrtPosition_[ bodyUndergoing ][ bodyExerting ];
            }
            else
            {
                scalarTermWrtPosition += currentTotalPotentialWrtPosition_[ bodyUndergoing ][ bodyUndergoing ];
            }
            break;
        case 1:
            if( wrtExerting )
            {
                scalarTermWrtPosition += currentTotalPotentialWrtPosition_[ bodyExerting ][ bodyExerting ];
            }
            else
            {
                scalarTermWrtPosition += currentTotalPotentialWrtPosition_[ bodyExerting ][ bodyUndergoing ];
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
            Eigen::Vector3d normalizedRelativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).normalized( );
            Eigen::Vector3d exertingVelocity = eihEquations_->getVelocity( bodyExerting );
            if( wrtExerting )
            {
                scalarTermWrtPosition += 2.0 * normalizedRelativePosition.dot( exertingVelocity ) * exertingVelocity.transpose( ) *
                    ( Eigen::Matrix3d::Identity( ) - normalizedRelativePosition * normalizedRelativePosition.transpose( ) ) /
                    eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).norm( );
            }
            else
            {
                scalarTermWrtPosition -= 2.0 * normalizedRelativePosition.dot( exertingVelocity ) * exertingVelocity.transpose( ) *
                    ( Eigen::Matrix3d::Identity( ) - normalizedRelativePosition * normalizedRelativePosition.transpose( ) ) /
                    eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).norm( );
            }
            break;
        }
        case 6:
            if( wrtExerting )
            {
                scalarTermWrtPosition += eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) *
                                        currentTotalPointMassAccelerationsWrtPosition_[ bodyExerting ][ bodyExerting ] +
                                        eihEquations_->getTotalPointMassAcceleration( bodyExerting ).transpose( );
            }
            else
            {
                scalarTermWrtPosition += eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) *
                                        currentTotalPointMassAccelerationsWrtPosition_[ bodyExerting ][ bodyUndergoing ] -
                                        eihEquations_->getTotalPointMassAcceleration( bodyExerting ).transpose( );
            }
            break;
        default:
            throw std::runtime_error( "Error when getting EIH scalar term partial w.r.t. positon, index " + std::to_string( termIndex ) + " not allowed." );

        }
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
                scalarTermWrtVelocity += 2.0 * eihEquations_->getVelocity( bodyUndergoing ).transpose( );
            }
            break;
        case 3:
            if( wrtExerting )
            {
                scalarTermWrtVelocity += 2.0 * eihEquations_->getVelocity( bodyExerting ).transpose( );
            }
            break;
        case 4:
            if( wrtExerting )
            {
                scalarTermWrtVelocity += eihEquations_->getVelocity( bodyUndergoing ).transpose( );
            }
            else
            {
                scalarTermWrtVelocity += eihEquations_->getVelocity( bodyExerting ).transpose( );
            }
            break;
        case 5:
        {
            if( wrtExerting )
            {
                Eigen::Vector3d normalizedRelativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).normalized( );
                Eigen::Vector3d exertingVelocity = eihEquations_->getVelocity( bodyExerting );
                scalarTermWrtVelocity += 2.0 * normalizedRelativePosition.dot( exertingVelocity ) * normalizedRelativePosition.transpose( );
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
            Eigen::Vector3d relativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting );


            double sign = ( wrtExerting ? 1.0 : -1.0 );
            vectorTermWrtPosition +=
                sign * eihEquations_->getRelativeVelocity( bodyUndergoing, bodyExerting ) *
                ( eihEquations_->getVelocity( bodyUndergoing ) .transpose( ) * inverseSquareDistance *
                ( Eigen::Matrix3d::Identity( ) - 2.0 * relativePosition * relativePosition.transpose( ) ) );
            break;
        }
        case 1:
        {
            double inverseSquareDistance = eihEquations_->getInverseSquareDistance( bodyUndergoing, bodyExerting );
            Eigen::Vector3d relativePosition = eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting );


            double sign = ( wrtExerting ? 1.0 : -1.0 );
            vectorTermWrtPosition +=
                sign * eihEquations_->getRelativeVelocity( bodyUndergoing, bodyExerting ) *
                ( eihEquations_->getVelocity( bodyExerting ) .transpose( ) * inverseSquareDistance *
                  ( Eigen::Matrix3d::Identity( ) - 2.0 * relativePosition * relativePosition.transpose( ) ) );
            break;
        }
        case 2:
            if( wrtExerting )
            {
                vectorTermWrtPosition += currentTotalPointMassAccelerationsWrtPosition_.at( bodyUndergoing ).at( bodyExerting );
            }
            else
            {
                vectorTermWrtPosition += currentTotalPointMassAccelerationsWrtPosition_.at( bodyUndergoing ).at( bodyUndergoing );
            }
            break;
        default:
            throw std::runtime_error( "Error when getting EIH vector term partial w.r.t. position, index " + std::to_string( termIndex ) + " not allowed." );

        }
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
            double baseMultiplier = -eihEquations_->getLineOfSighSpeed( bodyExerting, bodyUndergoing ) * inverseSquareDistance;
            {
                if( wrtExerting )
                {
                    vectorTermWrtVelocity += -baseMultiplier * Eigen::Matrix3d::Identity( ) + eihEquations_->getRelativeVelocity(
                        bodyUndergoing, bodyExerting ) * eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) * inverseSquareDistance;

                }
                else
                {
                    vectorTermWrtVelocity += baseMultiplier * Eigen::Matrix3d::Identity( );
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
                    vectorTermWrtVelocity -= baseMultiplier * Eigen::Matrix3d::Identity( );

                }
                else
                {
                    vectorTermWrtVelocity += baseMultiplier * Eigen::Matrix3d::Identity( ) + eihEquations_->getRelativeVelocity(
                        bodyUndergoing, bodyExerting ) * eihEquations_->getRelativePositions( bodyUndergoing, bodyExerting ).transpose( ) * inverseSquareDistance;

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




protected:

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations_;

    int numberOfExertingBodies_;

    int numberOfUndergoingBodies_;


    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtPosition_;

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



};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_EIHPARTIALS_H
