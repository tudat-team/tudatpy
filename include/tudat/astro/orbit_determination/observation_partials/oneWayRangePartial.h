/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYRANGEPARTIAL_H
#define TUDAT_ONEWAYRANGEPARTIAL_H

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to one-way range observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to one-way range observable partial. Implementation is taken
 *  from Moyer(2000) and is separately implemented for fixed receiver and transmitter.
 */
class OneWayRangeScaling: public DirectPositionPartialScaling< 1 >
{
public:

    OneWayRangeScaling( ): DirectPositionPartialScaling< 1 >( observation_models::one_way_range ){ }

    //! Destructor
    ~OneWayRangeScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation = Eigen::VectorXd::Constant( 1, TUDAT_NAN ) );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    virtual Eigen::Matrix< double, 1, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return fixedLinkEndScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    Eigen::Vector1d getLightTimePartialScalingFactor( )
    {
       return ( Eigen::Vector1d( ) << referenceLightTimeCorrectionScaling_ ).finished( );
    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:

    //! Computed scaling factor (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactor_;

    Eigen::Matrix< double, 1, 3 > fixedLinkEndScalingFactor_;

    //! Computed light time correction scaling factor
    double referenceLightTimeCorrectionScaling_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

};

}

}



#endif // TUDAT_ONEWAYRANGEPARTIAL_H
