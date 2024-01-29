/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONPARTIALS_H
#define TUDAT_CREATEOBSERVATIONPARTIALS_H

#include <memory>


#include "tudat/simulation/estimation_setup/createDopplerPartials.h"
#include "tudat/simulation/estimation_setup/createNWayRangePartials.h"
#include "tudat/simulation/estimation_setup/createEulerAngleObservationPartials.h"
#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"
#include "tudat/simulation/estimation_setup/createObsevationBiasPartial.h"

namespace tudat
{

namespace observation_partials
{


//! Function to split observation partials and scaling object (produced by observationPartialsAndScaler function) into separate
//! containers
/*!
 *  Function to split observation partials and scaling object (produced by observationPartialsAndScaler function) into separate
 *  containers
 *  \param observationPartialsAndScalers Combined list of observation partials and scaling objects
 *  \param observationPartials List of observation partials, per link ends, and per parameter indices (returned by reference)
 *  \param observationPartialScalers List of position partial scaling objects, per link ends (returned by reference)
 */
template< int ObservationSize >
void splitObservationPartialsAndScalers(
        const std::map< observation_models::LinkEnds,
        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >,
        std::shared_ptr< PositionPartialScaling > > >& observationPartialsAndScalers,
        std::map< observation_models::LinkEnds, std::map< std::pair< int, int >,
        std::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > >& observationPartials,
        std::map< observation_models::LinkEnds, std::shared_ptr< observation_partials::PositionPartialScaling  > >&
        observationPartialScalers )
{
    observationPartials.clear( );
    observationPartialScalers.clear( );
    // Put one-way range partials and scalers in member variables.
    for( typename std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >,
         std::shared_ptr< PositionPartialScaling > > >::const_iterator
         rangePartialPairIterator = observationPartialsAndScalers.begin( ); rangePartialPairIterator != observationPartialsAndScalers.end( );
         rangePartialPairIterator++ )
    {
        observationPartials[ rangePartialPairIterator->first ] = rangePartialPairIterator->second.first;
        observationPartialScalers[ rangePartialPairIterator->first ] = rangePartialPairIterator->second.second;
    }
}


template< typename ParameterType, typename TimeType, int ObservationSize >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > createDifferencedObservablePartials(
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable = false,
        const bool isPartialForConcatenatedObservable = false );

//! Interface class for creating observation partials
/*!
 *  Interface class for creating observation partials. This class is used instead of a single templated free function to
 *  allow ObservationPartial derived classed with different ObservationSize template arguments to be created using the same
 *  interface. This class has template specializations for each value of ObservationSize, and contains a single
 *  createObservationModel function that performs the required operation.
 */
template< int ObservationSize, typename ObservationScalarType, typename TimeType >
class ObservationPartialCreator
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param observationModelList List of observation models, with the link ends of map key, for which partials are to be created
     * \param bodies Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::pair< std::map< std::pair< int, int >,
            std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > > createObservationPartials(
            const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const bool isPartialForDifferencedObservable = false,
            const bool isPartialForConcatenatedObservable = false,
            const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) );
};

//! Interface class for creating observation partials for observables of size 1.
template< typename ObservationScalarType, typename TimeType >
class ObservationPartialCreator< 1, ObservationScalarType, TimeType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param observationModelList List of observation models, with the link ends of map key, for which partials are to be created
     * \param bodies Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::pair< std::map< std::pair< int, int >,
            std::shared_ptr< ObservationPartial< 1 > > >, std::shared_ptr< PositionPartialScaling > > createObservationPartials(
            const std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > observationModel,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const bool isPartialForDifferencedObservable = false,
            const bool isPartialForConcatenatedObservable = false,
                const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
    {
        std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 1 > > >, std::shared_ptr< PositionPartialScaling > > observationPartials;
        switch( observationModel->getObservableType( ) )
        {
        case observation_models::one_way_range:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 1, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::one_way_doppler:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 1, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::two_way_doppler:
            if( isPartialForConcatenatedObservable )
            {
                throw std::runtime_error( "Error when requesting partial creation for 2-way Doppler; concatenated partial not supported" );
            }
//            throw std::runtime_error( "Error, two-way instantaneous Doppler observable currently failing in unit tests, please contact Tudat support" );
            observationPartials = createTwoWayDopplerPartials< ObservationScalarType, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable );
            break;
        case observation_models::n_way_range:
            if( isPartialForConcatenatedObservable )
            {
                throw std::runtime_error( "Error when requesting partial creation for n-way range; concatenated partial not supported" );
            }
            observationPartials = createNWayRangePartials< ObservationScalarType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable );
            break;
        case observation_models::one_way_differenced_range:
            if( isPartialForConcatenatedObservable )
            {
                throw std::runtime_error( "Error when requesting partial creation for 1-way averaged Doppler; concatenated partial not supported" );
            }
            observationPartials = createDifferencedObservablePartials< ObservationScalarType, TimeType, 1 >(
                observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable  );
            break;
        case observation_models::n_way_differenced_range:
        case observation_models::dsn_n_way_averaged_doppler:
            if( isPartialForDifferencedObservable )
            {
                throw std::runtime_error( "Error when requesting partial creation for n-way averaged Doppler; differenced partial not supported" );
            }
            if( isPartialForConcatenatedObservable )
            {
                throw std::runtime_error( "Error when requesting partial creation for n-way averaged Doppler; concatenated partial not supported" );
            }
            observationPartials = createDifferencedObservablePartials< ObservationScalarType, TimeType, 1 >(
                        observationModel, bodies, parametersToEstimate  );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observationModel->getObservableType( ) ) + " of size 1 ";
            throw std::runtime_error( errorMessage );
        }

        return observationPartials;
    }
};

//! Interface class for creating observation partials for observables of size 2.
template< typename ObservationScalarType, typename TimeType >
class ObservationPartialCreator< 2, ObservationScalarType, TimeType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param observationModelList List of observation models, with the link ends of map key, for which partials are to be created
     * \param bodies Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::pair< std::map< std::pair< int, int >,
            std::shared_ptr< ObservationPartial< 2 > > >, std::shared_ptr< PositionPartialScaling > > createObservationPartials(
            const std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > observationModel,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const bool isPartialForDifferencedObservable = false,
            const bool isPartialForConcatenatedObservable = false,
            const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > ( ) )
    {
        std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 2 > > >, std::shared_ptr< PositionPartialScaling > > observationPartials;

        switch( observationModel->getObservableType( ) )
        {
        case observation_models::angular_position:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 2, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::relative_angular_position:
            observationPartials = createDifferencedObservablePartials< ObservationScalarType, TimeType, 2 >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observationModel->getObservableType( ) ) + " of size 2 ";
            throw std::runtime_error( errorMessage );
        }
        return observationPartials;
    }

};

//! Interface class for creating observation partials for observables of size 3.
template< typename ObservationScalarType, typename TimeType >
class ObservationPartialCreator< 3, ObservationScalarType, TimeType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param observationModelList List of observation models, with the link ends of map key, for which partials are to be created
     * \param bodies Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::pair< std::map< std::pair< int, int >,
            std::shared_ptr< ObservationPartial< 3 > > >, std::shared_ptr< PositionPartialScaling > > createObservationPartials(
            const std::shared_ptr< observation_models::ObservationModel< 3, ObservationScalarType, TimeType > > observationModel,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const bool isPartialForDifferencedObservable = false,
            const bool isPartialForConcatenatedObservable = false,
            const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
    {
        std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 3 > > >, std::shared_ptr< PositionPartialScaling > > observationPartials;

        switch( observationModel->getObservableType( ) )
        {
        case observation_models::position_observable:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 3, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::euler_angle_313_observable:
            observationPartials = createEulerAngleObservablePartials< ObservationScalarType >(
                        observationModel->getLinkEnds( ), bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::velocity_observable:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 3, TimeType >(
                        observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        case observation_models::relative_position_observable:
            observationPartials = createSingleLinkObservationPartials< ObservationScalarType, 3, TimeType >(
                    observationModel, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observationModel->getObservableType( ) ) + " of size 3 ";
            throw std::runtime_error( errorMessage );
        }
        return observationPartials;
    }

};


template< typename ObservationScalarType, typename TimeType, int ObservationSize >
std::map< observation_models::LinkEnds,
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > > createObservablePartialsList(
        const std::map< observation_models::LinkEnds,
        std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > > observationModelList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable = false,
        const bool isPartialForConcatenatedObservable = false,
        const std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface =
                std::shared_ptr< propagators::DependentVariablesInterface< TimeType > >( ) )
{
    std::map< observation_models::LinkEnds,
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
    std::shared_ptr< PositionPartialScaling > > > partialsList;

    observation_models::ObservableType observableType = observation_models::undefined_observation_model;
    for( auto it : observationModelList )
    {
        if( observableType == observation_models::undefined_observation_model )
        {
            observableType = it.second->getObservableType( );
        }
        else if( observableType != it.second->getObservableType( ) )
        {
            throw std::runtime_error( "Error when creating differenced observation partials, input models are inconsistent" );
        }
        partialsList[ it.first ] =  ObservationPartialCreator< ObservationSize, ObservationScalarType, TimeType >::createObservationPartials(
                    it.second, bodies, parametersToEstimate, isPartialForDifferencedObservable, isPartialForConcatenatedObservable, dependentVariablesInterface );
    }

    return partialsList;
}


template< int ObservationSize >
class DifferencedObservationPartialCreator
{
public:
    static std::shared_ptr< ObservationPartial< ObservationSize > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial,
            const std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial,
            const observation_models::LinkEnds& linkEnds,
            const simulation_setup::SystemOfBodies& bodies );
};

template< >
class DifferencedObservationPartialCreator< 1 >
{
public:
    static std::shared_ptr< ObservationPartial< 1 > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< 1 > > firstPartial,
            const std::shared_ptr< ObservationPartial< 1 > > secondPartial,
            const observation_models::LinkEnds& linkEnds,
            const simulation_setup::SystemOfBodies& bodies )
    {
        using namespace observation_models;

        std::shared_ptr< ObservationPartial< 1 > > differencedPartial;
        switch( differencedObservableType )
        {
        case one_way_differenced_range:
        {
            if( !isParameterObservationLinkTimeProperty( getDifferencedPartialParameterIdentifier( firstPartial, secondPartial ).first ) )
            {
                if( firstPartial != nullptr )
                {
                    if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( firstPartial ) == nullptr )
                    {
                        throw std::runtime_error( "Error when creating one-way differenced range partial; first input object type is incompatible" );
                    }
                    else if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( firstPartial )->getObservableType( ) != one_way_range )
                    {
                        throw std::runtime_error( "Error when creating one-way differenced range partial; first input observable type is incompatible" );
                    }
                }

                if( secondPartial != nullptr )
                {
                    if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( secondPartial ) == nullptr )
                    {
                        throw std::runtime_error( "Error when creating one-way differenced range partial; second input object type is incompatible" );
                    }
                    else if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( secondPartial )->getObservableType( ) != one_way_range )
                    {
                        throw std::runtime_error( "Error when creating one-way differenced range partial; second input observable type is incompatible" );
                    }
                }
            }

            differencedPartial = std::make_shared< DifferencedObservablePartial< 1 > >(
                        firstPartial, secondPartial, &observation_models::getDifferencedOneWayRangeScalingFactor,
                        getUndifferencedTimeAndStateIndices( one_way_differenced_range, linkEnds.size( ) ) );
            break;
        }
        case n_way_differenced_range:
        {
            if( firstPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( firstPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating n-way differenced range partial; first input object type is incompatible" );
                }
            }

            if( secondPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( secondPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating n-way differenced range partial; second input object type is incompatible" );
                }
            }

            differencedPartial = std::make_shared< DifferencedObservablePartial< 1 > >(
                        firstPartial, secondPartial, &observation_models::getDifferencedNWayRangeScalingFactor,
                        getUndifferencedTimeAndStateIndices( n_way_differenced_range, linkEnds.size( ) ) );
            break;
        }
        case dsn_n_way_averaged_doppler:
        {
            if( firstPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( firstPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating DSN n-way averaged Doppler partial; first input object type is incompatible" );
                }
            }

            if( secondPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( secondPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating DSN n-way averaged Doppler partial; second input object type is incompatible" );
                }
            }

            const std::function< double ( std::vector< FrequencyBands >, double ) > receivedFrequencyFunction =
                    createLinkFrequencyFunction(
                            bodies, linkEnds, observation_models::retransmitter, observation_models::receiver );

            const std::function< double(
                    const observation_models::LinkEndType, const std::vector< Eigen::Vector6d >&,
                    const std::vector< double >&, const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >,
                    const bool ) > scalingFactorFunction =
                            std::bind( observation_models::getDsnNWayAveragedDopplerScalingFactor, receivedFrequencyFunction,
                                       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                                       std::placeholders::_4, std::placeholders::_5 );

            differencedPartial = std::make_shared< DifferencedObservablePartial< 1 > >(
                        firstPartial, secondPartial, scalingFactorFunction,
                        getUndifferencedTimeAndStateIndices( dsn_n_way_averaged_doppler, linkEnds.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial (size 1); observable " + getObservableName( differencedObservableType ) +
                                      " is not differenced. " );

        }
        return differencedPartial;
    }
};

template< >
class DifferencedObservationPartialCreator< 2 >
{
public:
    static std::shared_ptr< ObservationPartial< 2 > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< 2 > > firstPartial,
            const std::shared_ptr< ObservationPartial< 2 > > secondPartial,
            const observation_models::LinkEnds& linkEnds,
            const simulation_setup::SystemOfBodies& bodies )
    {
        using namespace observation_models;

        std::shared_ptr< ObservationPartial< 2 > > differencedPartial;
        switch( differencedObservableType )
        {
        case relative_angular_position:
        {
            if( !isParameterObservationLinkTimeProperty( getDifferencedPartialParameterIdentifier( firstPartial, secondPartial ).first ) )
            {
                if( firstPartial != nullptr )
                {
                    if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( firstPartial ) == nullptr )
                    {
                        throw std::runtime_error( "Error when creating relative angular position partial; first input object type is incompatible" );
                    }
                    else if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( firstPartial )->getObservableType( ) != angular_position )
                    {
                        throw std::runtime_error( "Error when creating relative angular position partial; first input observable type is incompatible" );
                    }
                }

                if( secondPartial != nullptr )
                {
                    if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( secondPartial ) == nullptr )
                    {
                        throw std::runtime_error( "Error when creating relative angular position partial; second input object type is incompatible" );
                    }
                    else if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( secondPartial )->getObservableType( ) != angular_position )
                    {
                        throw std::runtime_error( "Error when creating relative angular position partial; second input observable type is incompatible" );
                    }
                }
            }

            differencedPartial = std::make_shared< DifferencedObservablePartial< 2 > >(
                    firstPartial, secondPartial, &getRelativeAngularPositionScalingFactor,
                    getUndifferencedTimeAndStateIndices( relative_angular_position, linkEnds.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial (size 2); observable " + getObservableName( differencedObservableType ) +
                                      " is not differenced. " );

        }
        return differencedPartial;
    }
};

template< int ObservationSize >
std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
mergeUndifferencedPartialContribution(
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& firstPartialList,
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& secondPartialList )
{
    std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
            mergedPartials;

    std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial;
    std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial;
    for( auto it : firstPartialList )
    {
        firstPartial = it.second;
        if( secondPartialList.count( it.first ) != 0 )
        {
            secondPartial = secondPartialList.at( it.first );
        }
        else
        {
            secondPartial = nullptr;
        }
        mergedPartials[ it.first ] = std::make_pair( firstPartial, secondPartial );
    }

    for( auto it : secondPartialList )
    {
        if( mergedPartials.count( it.first ) == 0 )
        {
            mergedPartials[ it.first ] = std::make_pair( nullptr, it.second );
        }
    }
    return mergedPartials;
}


template< typename ParameterType, typename TimeType, int ObservationSize >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > createDifferencedObservablePartials(
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable,
        const bool isPartialForConcatenatedObservable )
{
    using namespace observation_models;

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
            std::shared_ptr< PositionPartialScaling > > differencedPartialsAndScaling;

    LinkEnds linkEnds = observationModel->getLinkEnds( );
    ObservableType differencedObservableType = observationModel->getObservableType( );
    ObservableType undifferencedObservableType = getUndifferencedObservableType(
                differencedObservableType );

    auto undifferencedObservationModels = UndifferencedObservationModelExtractor< ObservationSize >::extract( observationModel );

    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelFirst =
            undifferencedObservationModels.first;
    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelSecond =
            undifferencedObservationModels.second;

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            firstUndifferencedObservablePartials =
            ObservationPartialCreator<ObservationSize, ParameterType, TimeType >::createObservationPartials(
                undifferencedObservationModelFirst, bodies, parametersToEstimate, true, isPartialForConcatenatedObservable );

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            secondUndifferencedObservablePartials =
            ObservationPartialCreator<ObservationSize, ParameterType, TimeType >::createObservationPartials(
                undifferencedObservationModelSecond, bodies, parametersToEstimate, true, isPartialForConcatenatedObservable );

    std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
            mergedPartials = mergeUndifferencedPartialContribution(
                firstUndifferencedObservablePartials.first, secondUndifferencedObservablePartials.first );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > differencedObservationPartialList;

    // Iterate over all one-way range partials and create one-way range rate partial from them.
    for( auto it : mergedPartials )
    {
        // Create range rate partial.
        differencedObservationPartialList[ it.first ] =
                DifferencedObservationPartialCreator< ObservationSize >::createDifferencedObservationPartial(
                    differencedObservableType,
                    it.second.first,
                    it.second.second,
                    linkEnds,
                    bodies );
    }


    // Create bias partials
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =  parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        std::shared_ptr< ObservationPartial< ObservationSize > > currentDifferencedObservationPartial;
        if( !isParameterObservationLinkTimeProperty( parameterIterator->second->getParameterName( ).first ) &&
            isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first ) )
        {
            currentDifferencedObservationPartial = createObservationPartialWrtLinkProperty< ObservationSize >(
                        linkEnds, undifferencedObservableType, parameterIterator->second, bodies, isPartialForDifferencedObservable, isPartialForConcatenatedObservable );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDifferencedObservationPartial != nullptr )
        {
            // Add partial to the list.
            std::pair< double, double > currentPair = std::pair< int, int >( parameterIterator->first,
                                                                             parameterIterator->second->getParameterSize( ) );
            differencedObservationPartialList[ currentPair ] = currentDifferencedObservationPartial;
        }
    }

    differencedPartialsAndScaling = std::make_pair(
                differencedObservationPartialList, ObservationPartialScalingCreator< ObservationSize >::
                template createDifferencedPositionPartialScalingObject< ParameterType, TimeType >(
                    differencedObservableType, firstUndifferencedObservablePartials.second,
                    secondUndifferencedObservablePartials.second, bodies ) );
    return differencedPartialsAndScaling;


}



}

}


#endif // TUDAT_CREATEOBSERVATIONPARTIALS_H
