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


#include "tudat/simulation/estimation_setup/createRelativeAngularPositionPartials.h"
#include "tudat/simulation/estimation_setup/createDopplerPartials.h"
#include "tudat/simulation/estimation_setup/createNWayRangePartials.h"
#include "tudat/simulation/estimation_setup/createEulerAngleObservationPartials.h"
#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/simulation/estimation_setup/createDifferencedObservablePartials.h"

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
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    std::shared_ptr< ObservationPartial< ObservationSize > > >,
    std::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::map< observation_models::LinkEnds,
            std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
            observationModelList,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate );
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
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    std::shared_ptr< ObservationPartial< 1 > > >,
    std::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::map< observation_models::LinkEnds,
            std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > >
            observationModelList,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 1 > > >,
                std::shared_ptr< PositionPartialScaling > > > observationPartialList;
        switch( observableType )
        {
        case observation_models::one_way_range:
            observationPartialList = createSingleLinkObservationPartialsList< ObservationScalarType, TimeType, 1 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::one_way_doppler:
            observationPartialList = createSingleLinkObservationPartialsList< ObservationScalarType, TimeType, 1 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::two_way_doppler:
            observationPartialList = createTwoWayDopplerPartials< ObservationScalarType, TimeType >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::one_way_differenced_range:
            observationPartialList = createDifferencedObservablePartialsList< ObservationScalarType, TimeType, 1 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::n_way_range:
            observationPartialList = createNWayRangePartials< ObservationScalarType >(
                        observationModelList, bodies, parametersToEstimate );

            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observableType ) + " of size 1 ";
            throw std::runtime_error( errorMessage );
        }

        return observationPartialList;
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
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    std::shared_ptr< ObservationPartial< 2 > > >,
    std::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::map< observation_models::LinkEnds,
            std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > > observationModelList,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 2 > > >,
                std::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::angular_position:
            observationPartialList = createSingleLinkObservationPartialsList< ObservationScalarType, TimeType, 2 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::relative_angular_position:
            observationPartialList = createRelativeAngularPositionPartials< ObservationScalarType >(
                    utilities::createVectorFromMapKeys( observationModelList ), bodies, parametersToEstimate,
                    getLightTimeCorrectionsList( observationModelList ) );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observableType ) + " of size 2 ";
            throw std::runtime_error( errorMessage );
        }
        return observationPartialList;
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
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    std::shared_ptr< ObservationPartial< 3 > > >,
    std::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::map< observation_models::LinkEnds,
            std::shared_ptr< observation_models::ObservationModel< 3, ObservationScalarType, TimeType > > >
            observationModelList,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 3 > > >,
                std::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::position_observable:
            observationPartialList = createSingleLinkObservationPartialsList< ObservationScalarType, TimeType, 3 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        case observation_models::euler_angle_313_observable:
            observationPartialList = createEulerAngleObservablePartials< ObservationScalarType >(
                        utilities::createVectorFromMapKeys( observationModelList ), bodies, parametersToEstimate );
            break;

        case observation_models::velocity_observable:
            observationPartialList = createSingleLinkObservationPartialsList< ObservationScalarType, TimeType, 3 >(
                        observationModelList, bodies, parametersToEstimate );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    std::to_string( observableType ) + " of size 3 ";
            throw std::runtime_error( errorMessage );        }
        return observationPartialList;
    }

};


}

}


#endif // TUDAT_CREATEOBSERVATIONPARTIALS_H
