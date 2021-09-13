/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating Z-X-Z Euler angle body orientation angles
/*!
 *  Class for simulating Z-X-Z Euler angle body orientation angles, as generated with e.g. star tracker observations. The
 * observable is a Vector of size 3, with entries [alpha,beta,gamma], where the inertial-to-body-fixed rotation matrix is
 * defined as R=R_z(alpha)R_x(beta)R_z(gamma).
 */
template< typename ObservationScalarType = double,
          typename TimeType = double >
class EulerAngle313ObservationModel: public ObservationModel< 3, ObservationScalarType, TimeType >
{
public:    

    //! Constructor
    /*!
     * Constructor
     * \param toBodyFixedFrameFunction Function that returns the rotation from inertial to body-fixed frame as a function of time
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none)
     */
    EulerAngle313ObservationModel(
            const LinkEnds linkEnds,
            const std::function< Eigen::Quaterniond( const TimeType ) > toBodyFixedFrameFunction,
            const std::shared_ptr< ObservationBias< 3 > > observationBiasCalculator = nullptr ):
        ObservationModel< 3, ObservationScalarType, TimeType >( euler_angle_313_observable, linkEnds, observationBiasCalculator ),
      toBodyFixedFrameFunction_( toBodyFixedFrameFunction ){ }

    //! Destructor
    ~EulerAngle313ObservationModel( ){ }

    //! Function to compute ideal Euler angle observation at given time.
    /*!
     *  This function computes the ideal Euler angle observation at a given time (without biases).
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference. For this observable, the link end states are NOT filled,
     *  as no translational state is used by the observation model.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid (must be observed_body for this derived class)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation: filled with NaN vector by this function
     *  \return Ideal Euler angle position observable.
     */
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

    //! Function that returns the rotation from inertial to body-fixed frame as a function of time
    std::function< Eigen::Quaterniond( const TimeType ) > toBodyFixedFrameFunction_;
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_EULERANGLEOBSERVATIONMODEL_H
