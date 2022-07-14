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
            const std::pair< std::vector< int >, std::vector< int > > timeStateIndices ):
        firstPartialScaling_( firstPartialScaling ), secondPartialScaling_( secondPartialScaling ),
        firstIndices_( timeStateIndices.first ), secondIndices_( timeStateIndices.second ){ }

    //! Destructor
    ~DifferencedObservablePartialScaling( ){ }


    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation )
    {
        try
        {
            firstPartialScaling_->update(
                        utilities::getVectorEntries( linkEndStates, firstIndices_ ), utilities::getVectorEntries( times, firstIndices_ ),
                        fixedLinkEnd, Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );
            secondPartialScaling_->update(
                        utilities::getVectorEntries( linkEndStates, secondIndices_ ), utilities::getVectorEntries( times, secondIndices_ ),
                        fixedLinkEnd, Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );;
        }
        catch( const std::exception& caughtException )
        {
            std::string exceptionText = std::string( caughtException.what( ) );
            throw std::runtime_error( "Error when computing differenced observation partial scaling, error: " + exceptionText );
        }


    }

    std::shared_ptr< PositionPartialScaling > firstPartialScaling_;

    std::shared_ptr< PositionPartialScaling > secondPartialScaling_;

    std::vector< int > firstIndices_;

    std::vector< int > secondIndices_;
};

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
            const std::function< double( const std::vector< double >&, const observation_models::LinkEndType ) > scalingFactorFunction,
            const std::pair< std::vector< int >, std::vector< int > >& undifferencedTimeAndStateIndices ):
        ObservationPartial< ObservationSize >( firstPartial->getParameterIdentifier( ) ),
        firstPartial_( firstPartial ),
        secondPartial_( secondPartial ),
        scalingFactorFunction_( scalingFactorFunction ),
    undifferencedTimeAndStateIndices_( undifferencedTimeAndStateIndices )
    {
        if( firstPartial_->getParameterIdentifier( ) != secondPartial_->getParameterIdentifier( ) )
        {
            throw std::runtime_error( "Error when creating differenced observable partial; first and second parameter identifiers are no equal" );
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
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) )
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
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > firstPartials =
                firstPartial_->calculatePartial( firstStates, firstTimes, linkEndOfFixedTime );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > secondPartials =
                secondPartial_->calculatePartial( secondStates, secondTimes, linkEndOfFixedTime );
        if( firstPartials.size( ) != secondPartials.size( ) )
        {
            throw std::runtime_error(
                        "Error when calculating differenced observation partials, first and second partial sets are inconsistent(?)" );
        }


        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > differencedPartials;
        double scalingFactor = scalingFactorFunction_( times, linkEndOfFixedTime );

        // Scale partials by arc duration
        for( unsigned int i = 0; i < firstPartials.size( ); i++ )
        {
            differencedPartials.push_back(
                        std::make_pair( -firstPartials[ i ].first * scalingFactor, firstPartials[ i ].second ) );
        }

        for( unsigned int i = 0; i < secondPartials.size( ); i++ )
        {
            differencedPartials.push_back(
                        std::make_pair( secondPartials[ i ].first * scalingFactor, secondPartials[ i ].second ) );
        }

        return differencedPartials;
    }

protected:

    //! Partial object for arc start range observation
    std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial_;

    //! Partial object for arc end range observation
    std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial_;

    std::function< double( const std::vector< double >&, const observation_models::LinkEndType ) > scalingFactorFunction_;

     const std::pair< std::vector< int >, std::vector< int > > undifferencedTimeAndStateIndices_;
};


}

}


#endif // TUDAT_DIFFERENCEDOBSERVATIONPARTIAL_H
