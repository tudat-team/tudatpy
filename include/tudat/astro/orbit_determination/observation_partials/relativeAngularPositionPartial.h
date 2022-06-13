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

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
/*!
 * Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
 */
Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition2(
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
Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition2(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

//! Function to compute the derivative of (direct geometric) right ascension and declination.
/*!
 * Function to compute the derivative of (direct geometric) right ascension and declination.
 * \param relativeRangeVector Vector from observer to observed object
 * \param isLinkEndReceiver True if the partial is to be computed w.r.t. position of observer, false if it is the observed
 * object
 * \return Derivative of (direct geometric) right ascension and declination.
 */
Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition2(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver );

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

//! Class to compute the partial derivatives of a relative angular position observation partial.
class RelativeAngularPositionPartial: public ObservationPartial< 2 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > RelativeAngularPositionPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param relativeAngularPositionScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    RelativeAngularPositionPartial(
            const std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >& lightTimeCorrectionPartials =
            std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) ):
        ObservationPartial< 2 >( parameterIdentifier ),
        relativeAngularPositionScaler_( relativeAngularPositionScaler ), positionPartialList_( positionPartialList )
    {

        if ( lightTimeCorrectionPartials.size( ) != 0 )
        {
            if ( lightTimeCorrectionPartials.size( ) != 2 )
            {
                throw std::runtime_error( "Error when making apparent distance partials, light time corrections for "
                                          + std::to_string( lightTimeCorrectionPartials.size( ) ) + " links found, instead of 2.");
            }
            else
            {
                lightTimeCorrectionPartialsFirstTransmitter_ = lightTimeCorrectionPartials[ 0 ];
                lightTimeCorrectionPartialsSecondTransmitter_ = lightTimeCorrectionPartials[ 1 ];

                if ( lightTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ) !=
                     lightTimeCorrectionPartialsFunctionsSecondTransmitter_.size( ) )
                {
                    throw std::runtime_error( "Error when making apparent distance partials, number of  light time "
                                              "corrections partials functions not consistent between the first transmitter leg"
                                              "and the second transmitter leg." );
                }

                // Create light time correction partial functions for link between receiver and first transmitter
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                        const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialFirstTransmitter;

                for( unsigned int i = 0; i < lightTimeCorrectionPartialsFirstTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialFirstTransmitter = getLightTimeParameterPartialFunction(
                            parameterIdentifier, lightTimeCorrectionPartialsFirstTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialFirstTransmitter.second != 0 )
                    {
                        lightTimeCorrectionPartialsFunctionsFirstTransmitter_.push_back( lightTimeCorrectionPartialFirstTransmitter.first );
                    }
                }

                // Create light time correction partial functions for link between receiver and second transmitter
                std::pair< std::function< SingleOneWayRangePartialReturnType(
                        const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartialSecondTransmitter;

                for( unsigned int i = 0; i < lightTimeCorrectionPartialsSecondTransmitter_.size( ); i++ )
                {
                    lightTimeCorrectionPartialSecondTransmitter = getLightTimeParameterPartialFunction(
                            parameterIdentifier, lightTimeCorrectionPartialsSecondTransmitter_.at( i ) );
                    if( lightTimeCorrectionPartialSecondTransmitter.second != 0 )
                    {
                        lightTimeCorrectionPartialsFunctionsSecondTransmitter_.push_back( lightTimeCorrectionPartialSecondTransmitter.first );
                    }
                }
            }
        }
    }

    //! Destructor.
    ~RelativeAngularPositionPartial( ){ }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    RelativeAngularPositionPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector2d& currentObservation = Eigen::Vector2d::Constant( TUDAT_NAN ) );

    //! Function to get the number of light-time correction partial functions for receiver - first transmitter leg.
    /*!
     * Number of light-time correction partial functions for receiver - first transmitter leg.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLightTimeCorrectionPartialsFunctionsFirstTransmitter( )
    {
        return lightTimeCorrectionPartialsFunctionsFirstTransmitter_.size( );
    }

    //! Function to get the number of light-time correction partial functions for receiver - second transmitter leg.
    /*!
     * Number of light-time correction partial functions for receiver - second transmitter leg.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLightTimeCorrectionPartialsFunctionsSecondTransmitter( )
    {
        return lightTimeCorrectionPartialsFunctionsSecondTransmitter_.size( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions for first transmitter.
    std::vector< std::function< SingleOneWayRangePartialReturnType( const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lightTimeCorrectionPartialsFunctionsFirstTransmitter_;

    //! List of light-time correction partial functions for second transmitter.
    std::vector< std::function< SingleOneWayRangePartialReturnType( const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lightTimeCorrectionPartialsFunctionsSecondTransmitter_;

    //! List of light-time correction partial objects for first transmitter.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialsFirstTransmitter_;

    //! List of light-time correction partial objects for second transmitter.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialsSecondTransmitter_;

    //! Pre-declare partial for current link end for first transmitter.
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialFirstTransmitter_;

    //! Pre-declare partial for current link end for second transmitter.
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartialSecondTransmitter_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};

}

}
#endif // RELATIVEANGULARPOSITIONPARTIAL_H
