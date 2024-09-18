/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ANGULARPOSITIONPARTIAL_H
#define TUDAT_ANGULARPOSITIONPARTIAL_H

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
/*!
 * Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver );

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
/*!
 * Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return derivative of (direct geometric) declination w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

//! Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
//! observed object.
/*!
 * Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
 * observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return Derivative of (direct geometric) right ascension and declination w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

//! Derived class for scaling three-dimensional position partial to angular position observable partial
class AngularPositionScaling: public DirectPositionPartialScaling< 2 >
{
public:

    AngularPositionScaling( ): DirectPositionPartialScaling< 2 >( observation_models::angular_position ){ }
    //! Destructor
    ~AngularPositionScaling( ){ }

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
                 const Eigen::VectorXd currentObservation );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 2, 3 > getPositionScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    Eigen::Matrix< double, 2, 3 > getFixedTimePositionScalingFactor(
        const observation_models::LinkEndType linkEndType )
    {
        return scalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    Eigen::Vector2d getLightTimePartialScalingFactor( )
    {
        return referenceLightTimeCorrectionScaling_;
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

    //! Predeclared common scaling factor
    Eigen::Matrix< double, 2, 3 > scalingFactor_;

    //! Computed scaling factor (at receiver)
    Eigen::Matrix< double, 2, 3 > referenceScalingFactor_;

    //! Computed light time correction scaling factor
    Eigen::Vector2d referenceLightTimeCorrectionScaling_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

};


}

}
#endif // ANGULARPOSITIONPARTIAL_H
