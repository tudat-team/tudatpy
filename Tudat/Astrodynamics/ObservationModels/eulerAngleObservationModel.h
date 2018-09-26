/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EULERANGLEOBSERVATIONMODEL_H
#define TUDAT_EULERANGLEOBSERVATIONMODEL_H

#include <map>

#include <functional>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{

namespace observation_models
{


template< typename ObservationScalarType = double,
          typename TimeType = double >
class EulerAngle313ObservationModel: public ObservationModel< 3, ObservationScalarType, TimeType >
{
public:    
    EulerAngle313ObservationModel(
            const std::function< Eigen::Quaterniond( const TimeType ) > toBodyFixedFrameFunction,
            const std::shared_ptr< ObservationBias< 3 > > observationBiasCalculator = nullptr ):
        ObservationModel< 3, ObservationScalarType, TimeType >( euler_angle_313_observable, observationBiasCalculator ),
      toBodyFixedFrameFunction_( toBodyFixedFrameFunction ){ }

    //! Destructor
    ~EulerAngle313ObservationModel( ){ }

    Eigen::Matrix< ObservationScalarType, 3, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        // Check link end
        if( linkEndAssociatedWithTime != observed_body )
        {
            throw std::runtime_error(
                        "Error when computing euler angle observable, associated link end must be observed_body " );
        }

        // Compute and return state.
        return basic_mathematics::get313EulerAnglesFromQuaternion(
                    toBodyFixedFrameFunction_( time ) ).template cast< ObservationScalarType >( );
    }

    Eigen::Matrix< ObservationScalarType, 3, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        // Check link end
        if( linkEndAssociatedWithTime != observed_body )
        {
            throw std::runtime_error(
                        "Error when computing position observable, associated link end must be observed_body " );
        }

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndTimes.push_back( time );

        linkEndStates.clear( );
        linkEndStates.push_back( Eigen::Matrix< double, 6, 1 >::Constant( TUDAT_NAN ) );

        return basic_mathematics::get313EulerAnglesFromQuaternion(
                    toBodyFixedFrameFunction_( time ) ).template cast< ObservationScalarType >( );
    }

private:

    std::function< Eigen::Quaterniond( const TimeType ) > toBodyFixedFrameFunction_;
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_EULERANGLEOBSERVATIONMODEL_H
