/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/gravitation/polyhedronGravityModel.h"


namespace tudat
{
namespace gravitation
{
//
//void PolyhedronGravitationalAccelerationModel::updateMembers( const double currentTime )
//{
//    if( !( this->currentTime_ == currentTime ) )
//    {
//
//        rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
//
//        subjectPositionFunction_( positionOfBodySubjectToAcceleration_ );
//        sourcePositionFunction_( positionOfBodyExertingAcceleration_ );
//        currentInertialRelativePosition_ =
//                positionOfBodySubjectToAcceleration_ - positionOfBodyExertingAcceleration_ ;
//
//        currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;
//
//        polyhedronCache_->update( currentRelativePosition_ );
//
//        // Compute the current acceleration
//        currentAccelerationInBodyFixedFrame_ = basic_mathematics::calculatePolyhedronGradientOfGravitationalPotential(
//                gravitationalParameterFunction_() / volumeFunction_(),
//                polyhedronCache_->getVerticesCoordinatesRelativeToFieldPoint(),
//                getVerticesDefiningEachFacet_(),
//                getVerticesDefiningEachEdge_(),
//                getFacetDyads_(),
//                getEdgeDyads_(),
//                polyhedronCache_->getPerFacetFactor(),
//                polyhedronCache_->getPerEdgeFactor() );
//
//        currentAcceleration_ = rotationToIntegrationFrame_ * currentAccelerationInBodyFixedFrame_;
//
//        // Compute the current gravitational potential
//        if ( updatePotential_ )
//        {
//            currentPotential_ = basic_mathematics::calculatePolyhedronGravitationalPotential(
//                    gravitationalParameterFunction_( ) / volumeFunction_( ),
//                    polyhedronCache_->getVerticesCoordinatesRelativeToFieldPoint( ),
//                    getVerticesDefiningEachFacet_( ),
//                    getVerticesDefiningEachEdge_( ),
//                    getFacetDyads_( ),
//                    getEdgeDyads_( ),
//                    polyhedronCache_->getPerFacetFactor( ),
//                    polyhedronCache_->getPerEdgeFactor( ) );
//        }
//
//        // Compute the current laplacian
//        if ( updateLaplacianOfPotential_ )
//        {
//            currentLaplacianOfPotential_ = basic_mathematics::calculatePolyhedronLaplacianOfGravitationalPotential(
//                    gravitationalParameterFunction_( ) / volumeFunction_( ),
//                    polyhedronCache_->getPerFacetFactor( ) );
//        }
//    }
//}


} // namespace gravitation

} // namespace tudat
