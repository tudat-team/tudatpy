/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
#define TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H

#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

inline Eigen::Matrix3d getCrossProductMatrix( const Eigen::Vector3d& vector )
{
    Eigen::Matrix3d crossProductMatrix;
    crossProductMatrix << 0.0, -vector.z( ), vector.y( ), vector.z( ), 0.0, -vector.x( ), -vector.y( ), vector.x( ), 0.0;
    return crossProductMatrix;
}

//! Calculate the Jacobians of position angle and separation distance with respect to both angular positions.
inline std::pair< Eigen::Matrix2d, Eigen::Matrix2d > calculatePositionAngleAndSeparationPartialWrtAngularPositions(
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const Eigen::Vector3d& j2000NorthPoleDirection,
        const bool calculatePositionAngle = true )
{
    if( linkEndStates.size( ) != 3 )
    {
        throw std::runtime_error( "Position-angle partials require states for two transmitters and one receiver." );
    }

    const Eigen::Vector3d firstRelativePosition = linkEndStates.at( 0 ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 );
    const Eigen::Vector3d secondRelativePosition = linkEndStates.at( 1 ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 );

    // This performs all physical singularity checks used by the observation model.
    observation_models::calculatePositionAngleAndSeparation(
            firstRelativePosition, secondRelativePosition, j2000NorthPoleDirection, calculatePositionAngle );

    const Eigen::Vector3d firstLineOfSight = firstRelativePosition.normalized( );
    const Eigen::Vector3d secondLineOfSight = secondRelativePosition.normalized( );
    Eigen::Matrix< double, 2, 3 > firstUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    Eigen::Matrix< double, 2, 3 > secondUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    if( calculatePositionAngle )
    {
        const Eigen::Vector3d northPoleDirection = j2000NorthPoleDirection.normalized( );
        const Eigen::Vector3d unnormalisedEastDirection = northPoleDirection.cross( firstLineOfSight );
        const double eastDirectionNorm = unnormalisedEastDirection.norm( );
        const Eigen::Vector3d eastDirection = unnormalisedEastDirection / eastDirectionNorm;
        const Eigen::Vector3d northDirection = firstLineOfSight.cross( eastDirection );

        const double eastProjection = secondLineOfSight.dot( eastDirection );
        const double northProjection = secondLineOfSight.dot( northDirection );
        const double positionAngleDenominator = eastProjection * eastProjection + northProjection * northProjection;

        const Eigen::Matrix3d eastDirectionPartial = ( Eigen::Matrix3d::Identity( ) - eastDirection * eastDirection.transpose( ) ) *
                getCrossProductMatrix( northPoleDirection ) / eastDirectionNorm;
        const Eigen::Matrix3d northDirectionPartial =
                -getCrossProductMatrix( eastDirection ) + getCrossProductMatrix( firstLineOfSight ) * eastDirectionPartial;

        firstUnitVectorPartial.row( 0 ) = ( northProjection * secondLineOfSight.transpose( ) * eastDirectionPartial -
                                            eastProjection * secondLineOfSight.transpose( ) * northDirectionPartial ) /
                positionAngleDenominator;
        secondUnitVectorPartial.row( 0 ) =
                ( northProjection * eastDirection.transpose( ) - eastProjection * northDirection.transpose( ) ) / positionAngleDenominator;
    }

    const double cosineOfSeparation = firstLineOfSight.dot( secondLineOfSight );
    const double sineOfSeparation = firstLineOfSight.cross( secondLineOfSight ).norm( );
    firstUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * firstLineOfSight - secondLineOfSight ).transpose( ) / sineOfSeparation;
    secondUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * secondLineOfSight - firstLineOfSight ).transpose( ) / sineOfSeparation;

    const double singularityTolerance = 100.0 * std::numeric_limits< double >::epsilon( );
    const auto getUnitVectorPartialWrtAngularPosition = [ singularityTolerance ]( const Eigen::Vector3d& lineOfSight ) {
        const double horizontalNorm = lineOfSight.head< 2 >( ).norm( );
        if( horizontalNorm <= singularityTolerance )
        {
            throw std::runtime_error( "Position-angle partials cannot be formed from right ascension at a pole of the state frame." );
        }

        Eigen::Matrix< double, 3, 2 > partial;
        partial.col( 0 ) << -lineOfSight.y( ), lineOfSight.x( ), 0.0;
        partial.col( 1 ) << -lineOfSight.x( ) * lineOfSight.z( ) / horizontalNorm, -lineOfSight.y( ) * lineOfSight.z( ) / horizontalNorm,
                horizontalNorm;
        return partial;
    };

    return std::make_pair( firstUnitVectorPartial * getUnitVectorPartialWrtAngularPosition( firstLineOfSight ),
                           secondUnitVectorPartial * getUnitVectorPartialWrtAngularPosition( secondLineOfSight ) );
}

//! Scaling object that updates the independent angular-position scalers for both light-time legs.
class PositionAngleAndSeparationPartialScaling : public PositionPartialScaling
{
public:
    PositionAngleAndSeparationPartialScaling( const std::shared_ptr< PositionPartialScaling > firstPartialScaling,
                                              const std::shared_ptr< PositionPartialScaling > secondPartialScaling ):
        firstPartialScaling_( firstPartialScaling ), secondPartialScaling_( secondPartialScaling )
    {}

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd ) override
    {
        if( fixedLinkEnd != observation_models::receiver )
        {
            throw std::runtime_error( "Position-angle and separation observations must be referenced at the receiver." );
        }
        if( linkEndStates.size( ) != 3 || times.size( ) != 3 )
        {
            throw std::runtime_error( "Position-angle partial scaling requires three link-end states and times." );
        }

        firstPartialScaling_->update( { linkEndStates.at( 0 ), linkEndStates.at( 2 ) },
                                      { times.at( 0 ), times.at( 2 ) },
                                      observation_models::receiver,
                                      Eigen::Vector2d::Constant( TUDAT_NAN ) );
        secondPartialScaling_->update( { linkEndStates.at( 1 ), linkEndStates.at( 2 ) },
                                       { times.at( 1 ), times.at( 2 ) },
                                       observation_models::receiver,
                                       Eigen::Vector2d::Constant( TUDAT_NAN ) );
    }

private:
    std::shared_ptr< PositionPartialScaling > firstPartialScaling_;
    std::shared_ptr< PositionPartialScaling > secondPartialScaling_;
};

template< int ObservationSize >
class PositionAngleAndSeparationPartial : public ObservationPartial< ObservationSize >
{
public:
    PositionAngleAndSeparationPartial( const std::shared_ptr< ObservationPartial< 2 > > firstAngularPositionPartial,
                                       const std::shared_ptr< ObservationPartial< 2 > > secondAngularPositionPartial,
                                       const Eigen::Vector3d& j2000NorthPoleDirection,
                                       const int componentIndex ):
        ObservationPartial< ObservationSize >( getParameterIdentifier( firstAngularPositionPartial, secondAngularPositionPartial ) ),
        firstAngularPositionPartial_( firstAngularPositionPartial ), secondAngularPositionPartial_( secondAngularPositionPartial ),
        j2000NorthPoleDirection_( j2000NorthPoleDirection ), componentIndex_( componentIndex )
    {
        if( ( ObservationSize == 2 && componentIndex_ != -1 ) || ( ObservationSize == 1 && componentIndex_ != 0 && componentIndex_ != 1 ) )
        {
            throw std::runtime_error( "Invalid component requested for a position-angle and separation partial." );
        }
    }

    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& =
                    Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) override
    {
        if( linkEndOfFixedTime != observation_models::receiver )
        {
            throw std::runtime_error( "Position-angle and separation partials must be referenced at the receiver." );
        }
        if( states.size( ) != 3 || times.size( ) != 3 )
        {
            throw std::runtime_error( "Position-angle and separation partials require three link-end states and times." );
        }

        const auto angularPositionJacobians =
                calculatePositionAngleAndSeparationPartialWrtAngularPositions( states, j2000NorthPoleDirection_, componentIndex_ != 1 );
        const int firstOutputRow = componentIndex_ < 0 ? 0 : componentIndex_;
        const Eigen::Matrix< double, ObservationSize, 2 > firstTransformation =
                angularPositionJacobians.first.block( firstOutputRow, 0, ObservationSize, 2 );
        const Eigen::Matrix< double, ObservationSize, 2 > secondTransformation =
                angularPositionJacobians.second.block( firstOutputRow, 0, ObservationSize, 2 );

        std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > partials;
        appendTransformedPartials( partials,
                                   firstAngularPositionPartial_,
                                   firstTransformation,
                                   { states.at( 0 ), states.at( 2 ) },
                                   { times.at( 0 ), times.at( 2 ) },
                                   ancillarySettings );
        appendTransformedPartials( partials,
                                   secondAngularPositionPartial_,
                                   secondTransformation,
                                   { states.at( 1 ), states.at( 2 ) },
                                   { times.at( 1 ), times.at( 2 ) },
                                   ancillarySettings );
        return partials;
    }

private:
    static estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier(
            const std::shared_ptr< ObservationPartial< 2 > >& firstPartial,
            const std::shared_ptr< ObservationPartial< 2 > >& secondPartial )
    {
        if( firstPartial != nullptr )
        {
            return firstPartial->getParameterIdentifier( );
        }
        if( secondPartial != nullptr )
        {
            return secondPartial->getParameterIdentifier( );
        }
        throw std::runtime_error( "Cannot create a position-angle partial from two null angular-position partials." );
    }

    static void appendTransformedPartials(
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > >& outputPartials,
            const std::shared_ptr< ObservationPartial< 2 > >& angularPositionPartial,
            const Eigen::Matrix< double, ObservationSize, 2 >& transformation,
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >& ancillarySettings )
    {
        if( angularPositionPartial == nullptr )
        {
            return;
        }
        const auto angularPartials = angularPositionPartial->calculatePartial(
                states, times, observation_models::receiver, ancillarySettings, Eigen::Vector2d::Constant( TUDAT_NAN ) );
        for( const auto& angularPartial : angularPartials )
        {
            outputPartials.push_back( std::make_pair( transformation * angularPartial.first, angularPartial.second ) );
        }
    }

    std::shared_ptr< ObservationPartial< 2 > > firstAngularPositionPartial_;
    std::shared_ptr< ObservationPartial< 2 > > secondAngularPositionPartial_;
    Eigen::Vector3d j2000NorthPoleDirection_;
    int componentIndex_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
