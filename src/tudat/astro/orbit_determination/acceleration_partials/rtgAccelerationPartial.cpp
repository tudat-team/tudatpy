/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/rtgAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/rtgForceVector.h"

namespace tudat
{

namespace acceleration_partials
{


//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RTGAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case rtg_force_vector: {
                    partialFunction =
                            std::bind( &RTGAccelerationPartial::wrtRTGForceVector,
                                       this,
                                       std::placeholders::_1 );
                    numberOfRows = parameter->getParameterSize( );
                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}


//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RTGAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace tudat::estimatable_parameters;

    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case rtg_force_vector_magnitude: {
                partialFunction =
                        std::bind( &RTGAccelerationPartial::wrtRTGForceVectorMagnitude,
                                   this,
                                   std::placeholders::_1 );
                numberOfRows = parameter->getParameterSize( );

                break;
            }
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}

//! Function for updating common blocks of partial to current state.
void RTGAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        using namespace tudat::basic_mathematics;
        using namespace tudat::linear_algebra;

        rtgAcceleration_->updateMembers( currentTime );

        currentTime_ = currentTime;

    }
}

//! Function to compute the partial w.r.t. reference force vector
void RTGAccelerationPartial::wrtRTGForceVector( Eigen::MatrixXd& partialDerivativeMatrix )
{
    // Compute partial derivative w.r.t. reference force vector
    double partialWrtReferenceForceVector = rtgAcceleration_->getCurrentDecayTerm( ) / rtgAcceleration_->evaluateBodyMassFunction( );
    partialDerivativeMatrix = partialWrtReferenceForceVector * rtgAcceleration_->getCurrentRotationToIntegrationFrameMatrix(  );

}

//! Function to compute the partial w.r.t. magnitude of reference force vector
void RTGAccelerationPartial::wrtRTGForceVectorMagnitude(Eigen::MatrixXd& partialDerivativeMatrix )
{
    // Compute partial derivative w.r.t. magnitude of reference force vector
    Eigen::Vector3d partialWrtReferenceForceMagnitude = rtgAcceleration_->getCurrentDecayTerm( ) /
            rtgAcceleration_->evaluateBodyMassFunction( ) * rtgAcceleration_->getCurrentRotationToIntegrationFrameMatrix(  ) * rtgAcceleration_->getBodyFixedForceUnitVectorAtReferenceEpoch( );
    partialDerivativeMatrix = partialWrtReferenceForceMagnitude;

}

}  // namespace acceleration_partials

}  // namespace tudat