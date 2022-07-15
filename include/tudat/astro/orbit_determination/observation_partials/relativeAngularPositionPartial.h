/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVEANGULARPOSITIONPARTIAL_H
#define TUDAT_RELATIVEANGULARPOSITIONPARTIAL_H

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to relative angular position observable partial
class RelativeAngularPositionScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~RelativeAngularPositionScaling( ){ }

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
    Eigen::Matrix< double, 2, 3 > getScalingFactor(
            const observation_models::LinkEndType linkEndType )
    {
        if ( linkEndType == observation_models::transmitter )
        {
            return referenceScalingFactorFirstTransmitter_;
        }
        else if ( linkEndType == observation_models::transmitter2 )
        {
            return - referenceScalingFactorSecondTransmitter_;
        }
        else if ( linkEndType == observation_models::receiver )
        {
            return referenceScalingFactorSecondTransmitter_ - referenceScalingFactorFirstTransmitter_;
        }
        else
        {
            throw std::runtime_error( "Error when geting relative angular position scaling factor, incorrect link end type provided: " +
                                      observation_models::getLinkEndTypeString( linkEndType ) );
        }
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    Eigen::Vector2d getLightTimePartialScalingFactorFirstTransmitter( )
    {
        return referenceLightTimeCorrectionScalingFirstTransmitter_;
    }

    Eigen::Vector2d getLightTimePartialScalingFactorSecondTransmitter( )
    {
        return referenceLightTimeCorrectionScalingSecondTransmitter_;
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

    //! Predeclared common scaling factor for first transmitter
    Eigen::Matrix< double, 2, 3 > scalingFactorFirstTransmitter_;

    //! Computed scaling factor (at receiver) for first transmitter
    Eigen::Matrix< double, 2, 3 > referenceScalingFactorFirstTransmitter_;

    //! Computed light time correction scaling factor for first transmitter
    Eigen::Vector2d referenceLightTimeCorrectionScalingFirstTransmitter_;

    //! Predeclared common scaling factor for second transmitter
    Eigen::Matrix< double, 2, 3 > scalingFactorSecondTransmitter_;

    //! Computed scaling factor (at receiver) for second transmitter
    Eigen::Matrix< double, 2, 3 > referenceScalingFactorSecondTransmitter_;

    //! Computed light time correction scaling factor for second transmitter
    Eigen::Vector2d referenceLightTimeCorrectionScalingSecondTransmitter_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

};

}

}
#endif // RELATIVEANGULARPOSITIONPARTIAL_H
