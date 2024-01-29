/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIFFERENCEDOBSERVATIONPARTIAL_H
#define TUDAT_DIFFERENCEDOBSERVATIONPARTIAL_H

#include <functional>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "tudat/simulation/environment_setup.h"
#include "tudat/astro/observation_models.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to one-way range-rate (differenced) observable partial
class DifferencedObservablePartialScaling: public PositionPartialScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param firstPartialScaling Partial scaling for arc start range observation
     * \param secondPartialScaling Partial scaling for arc end range observation
     */
    DifferencedObservablePartialScaling(
            const std::shared_ptr< PositionPartialScaling > firstPartialScaling,
            const std::shared_ptr< PositionPartialScaling > secondPartialScaling,
            const std::pair< std::vector< int >, std::vector< int > > timeStateIndices,
            const std::function< void( const observation_models::LinkEndType ) > customCheckFunction = nullptr ):
        firstPartialScaling_( firstPartialScaling ), secondPartialScaling_( secondPartialScaling ),
        firstIndices_( timeStateIndices.first ), secondIndices_( timeStateIndices.second ),
        customCheckFunction_( ){ }

    //! Destructor
    ~DifferencedObservablePartialScaling( ){ }


    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    std::shared_ptr< PositionPartialScaling > firstPartialScaling_;

    std::shared_ptr< PositionPartialScaling > secondPartialScaling_;

    std::vector< int > firstIndices_;

    std::vector< int > secondIndices_;

    std::function< void( const observation_models::LinkEndType ) > customCheckFunction_;
};

template< int ObservationSize >
estimatable_parameters::EstimatebleParameterIdentifier getDifferencedPartialParameterIdentifier(
    const std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial,
    const std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial )
{
    if( firstPartial != nullptr )
    {
        return firstPartial->getParameterIdentifier( );
    }
    else if( secondPartial != nullptr )
    {
        return secondPartial->getParameterIdentifier( );
    }
    else
    {
        throw std::runtime_error( "Error when getting differenced parameter partial identified; both partials are null." );
    }
}
//! Class to compute the partial derivatives of a one-way range-rate (differenced) observation
template< int ObservationSize >
class DifferencedObservablePartial: public ObservationPartial< ObservationSize >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param arcStartRangePartial Partial object for arc start range observation
     * \param arcEndRangePartial Partial object for arc end range observation
     */
    DifferencedObservablePartial(
            const std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial,
            const std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial,
            const std::function< double(
                    const observation_models::LinkEndType, const std::vector< Eigen::Vector6d >&,
                    const std::vector< double >&, const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >,
                    const bool ) > scalingFactorFunction,
            const std::pair< std::vector< int >, std::vector< int > >& undifferencedTimeAndStateIndices ):
        ObservationPartial< ObservationSize >( getDifferencedPartialParameterIdentifier< ObservationSize >( firstPartial, secondPartial ) ),
        firstPartial_( firstPartial ),
        secondPartial_( secondPartial ),
        scalingFactorFunction_( scalingFactorFunction ),
        undifferencedTimeAndStateIndices_( undifferencedTimeAndStateIndices )
    {
        if( firstPartial_ != nullptr && secondPartial_ != nullptr )
        {
//            if( firstPartial_->getParameterIdentifier( ) != secondPartial_->getParameterIdentifier( ) )
//            {
//                throw std::runtime_error( "Error when creating differenced observable partial; first and second parameter identifiers are no equal" );
//            }
        }
    }

    //! Destructor
    ~DifferencedObservablePartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation = Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        using namespace observation_partials;

        // Split input times/states for arc start and end ranges
        std::vector< Eigen::Vector6d > firstStates;
        std::vector< double > firstTimes;
        for( unsigned int i = 0; i < undifferencedTimeAndStateIndices_.first.size( ); i++ )
        {
            firstStates.push_back( states.at( undifferencedTimeAndStateIndices_.first.at( i ) ) );
            firstTimes.push_back( times.at( undifferencedTimeAndStateIndices_.first.at( i ) ) );
        }

        std::vector< Eigen::Vector6d > secondStates;
        std::vector< double > secondTimes;
        for( unsigned int i = 0; i < undifferencedTimeAndStateIndices_.second.size( ); i++ )
        {
            secondStates.push_back( states.at( undifferencedTimeAndStateIndices_.second.at( i ) ) );
            secondTimes.push_back( times.at( undifferencedTimeAndStateIndices_.second.at( i ) ) );
        }

        // Obtain arc start and end range partials
        std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > firstPartials;
        if( firstPartial_ != nullptr )
        {
            firstPartials = firstPartial_->calculatePartial( firstStates, firstTimes, linkEndOfFixedTime, ancillarySettings );
        }
        std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > secondPartials;
        if( secondPartial_ != nullptr )
        {
            secondPartials = secondPartial_->calculatePartial( secondStates, secondTimes, linkEndOfFixedTime, ancillarySettings );
        }

        std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > differencedPartials;
        double firstPartialScalingFactor = scalingFactorFunction_(
                linkEndOfFixedTime, states, times, ancillarySettings, true );
        double secondPartialScalingFactor = scalingFactorFunction_(
                linkEndOfFixedTime, states, times, ancillarySettings, false );

        // Scale partials by arc duration
        for( unsigned int i = 0; i < firstPartials.size( ); i++ )
        {
            differencedPartials.push_back(
                        std::make_pair( - firstPartials[ i ].first * firstPartialScalingFactor, firstPartials[ i ].second ) );
        }

        for( unsigned int i = 0; i < secondPartials.size( ); i++ )
        {
            differencedPartials.push_back(
                        std::make_pair( secondPartials[ i ].first * secondPartialScalingFactor, secondPartials[ i ].second ) );
        }

        return differencedPartials;
    }

protected:

    //! Partial object for arc start range observation
    std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial_;

    //! Partial object for arc end range observation
    std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial_;

    const std::function< double(
            const observation_models::LinkEndType, const std::vector< Eigen::Vector6d >&,
            const std::vector< double >&, const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >,
            const bool ) > scalingFactorFunction_;

    const std::pair< std::vector< int >, std::vector< int > > undifferencedTimeAndStateIndices_;
};


}

}


#endif // TUDAT_DIFFERENCEDOBSERVATIONPARTIAL_H
