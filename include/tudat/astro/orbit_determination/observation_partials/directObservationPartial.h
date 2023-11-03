/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIRECTOBSERVATIONPARTIAL_H
#define TUDAT_DIRECTOBSERVATIONPARTIAL_H

#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

template< int ObservationSize >
class DirectObservationPartial: public ObservationPartial< ObservationSize >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > ObservationPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > SingleObservationPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleLightTimePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param positionPartialScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    DirectObservationPartial(
            const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
            const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
            lighTimeCorrectionPartials =
            std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) ):
        ObservationPartial< ObservationSize >( parameterIdentifier ), positionPartialScaler_( positionPartialScaler ),
        positionPartialList_( positionPartialList )
    {
        stateEntryIndices_ = observation_models::getSingleLinkStateEntryIndices( positionPartialScaler->getObservableType( ) );

        std::pair< std::function< SingleLightTimePartialReturnType(
                    const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartial;

        // Create light time correction partial functions
        for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
        {
            lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
                        parameterIdentifier, lighTimeCorrectionPartials.at( i ) );
            if( lightTimeCorrectionPartial.second != 0 )
            {
                lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
            }
        }
    }

    //! Destructor.
    ~DirectObservationPartial( ) { }

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
    ObservationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        if( linkEndOfFixedTime != positionPartialScaler_->getCurrentLinkEndType( ) )
        {
            throw std::runtime_error( "Error observation partial and scaling are inconsistent" );
        }

        ObservationPartialReturnType returnPartial;

        // Iterate over all link ends
        for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
             positionPartialIterator_++ )
        {
            int currentIndex = stateEntryIndices_.at( positionPartialIterator_->first );

            currentState_  = states[ currentIndex ];
            currentTime_ = times[ currentIndex ];

            if( positionPartialScaler_->isVelocityScalingNonZero( ) )
            {
                // Scale position partials
                returnPartial.push_back(
                            std::make_pair(
                                positionPartialScaler_->getPositionScalingFactor( positionPartialIterator_->first ) *
                                ( positionPartialIterator_->second->calculatePartialOfPosition(
                                      currentState_ , currentTime_ ) ) +
                                positionPartialScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
                                ( positionPartialIterator_->second->calculatePartialOfVelocity(
                                      currentState_ , currentTime_ ) ), currentTime_ ) );
            }
            else
            {
                returnPartial.push_back(
                            std::make_pair(
                                positionPartialScaler_->getPositionScalingFactor( positionPartialIterator_->first ) *
                                ( positionPartialIterator_->second->calculatePartialOfPosition(
                                      currentState_ , currentTime_ ) ), currentTime_ ) );
            }
        }

        // Add scaled light-time correcion partials.
        for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
        {

            currentLinkTimeCorrectionPartial_ = lighTimeCorrectionPartialsFunctions_.at( i )( states, times );
            returnPartial.push_back(
                        std::make_pair( positionPartialScaler_->getLightTimePartialScalingFactor( ) *
                                        physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartial_.first,
                                        currentLinkTimeCorrectionPartial_.second ) );
        }

        if( useLinkIndependentPartials( ) )
        {
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > additionalPartials =
                    positionPartialScaler_->getLinkIndependentPartials( this->parameterIdentifier_ );
            for( unsigned int i = 0; i < additionalPartials.size( ); i++ )
            {
                returnPartial.push_back( additionalPartials.at( i ) );
            }

        }

        return returnPartial;
    }

    //! Function to get scaling object used for mapping partials of positions to partials of observable
    /*!
     * Function to get scaling object used for mapping partials of positions to partials of observable
     * \return
     */
    std::shared_ptr< PositionPartialScaling > getPositionPartialScaler( )
    {
        return positionPartialScaler_;
    }

    //! Function to get the number of light-time correction partial functions.
    /*!
     * Number of light-time correction partial functions.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
        return static_cast< int >( lighTimeCorrectionPartialsFunctions_.size( ) );
    }

    virtual bool useLinkIndependentPartials( )
    {
        return positionPartialScaler_->useLinkIndependentPartials( );
    }

    observation_models::ObservableType getObservableType( )
    {
        return positionPartialScaler_->getObservableType( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< std::function< SingleLightTimePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctions_;

    //! List of light-time correction partial objects.
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartials_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;

    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > currentLinkTimeCorrectionPartial_;

    std::map< observation_models::LinkEndType, int > stateEntryIndices_;

};


}

}

#endif // directObservationPartial_H
