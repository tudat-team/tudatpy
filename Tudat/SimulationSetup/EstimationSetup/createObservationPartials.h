/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"

#include "Tudat/SimulationSetup/EstimationSetup/createAngularPositionPartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"

namespace tudat
{

namespace observation_partials
{

//! Typedef for list of light time corrections for a list of link ends
typedef std::map< observation_models::LinkEnds,
std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >
PerLinkEndPerLightTimeSolutionCorrections;


//! Function to retrieve a list of light-time corrections per link end from a list of observation models.
/*!
 *  Function to retrieve a list of light-time corrections per link end from a list of observation models.
 *  \param observationModels Map of observation models (may not be of mixed type) with LinkEnds of observable as map key
 *  \return Map of light-time corrections, with associated link ends as key.
 */
template< typename ObservationScalarType, typename TimeType, int ObservationSize  >
PerLinkEndPerLightTimeSolutionCorrections getLightTimeCorrectionsList(
        const std::map< observation_models::LinkEnds, boost::shared_ptr< observation_models::ObservationModel<
        ObservationSize, ObservationScalarType, TimeType> > > observationModels )
{
    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionsList;
    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;

    // Retrieve type of observable
    observation_models::ObservableType observableType = observationModels.begin( )->second->getObservableType( );

    // Iterate over link ends
    for( typename  std::map< observation_models::LinkEnds, boost::shared_ptr< observation_models::ObservationModel<
         ObservationSize, ObservationScalarType, TimeType> > >::const_iterator
         observationModelIterator = observationModels.begin( );
         observationModelIterator != observationModels.end( ); observationModelIterator++ )
    {
        // Clear list, for current link ends.
        currentLightTimeCorrections.clear( );

        if( observationModelIterator->second->getObservableType( ) != observableType )
        {
            throw std::runtime_error( "Error when making grouped light time correction list, observable type is not constant" );
        }
        else
        {
            // Check type of observable
            switch( observableType )
            {
            case observation_models::oneWayRange:
            {
                boost::shared_ptr< observation_models::OneWayRangeObservationModel
                        < ObservationScalarType, TimeType> > oneWayRangeModel =
                        boost::dynamic_pointer_cast< observation_models::OneWayRangeObservationModel
                        < ObservationScalarType, TimeType> >
                        ( observationModelIterator->second );
                currentLightTimeCorrections.push_back(
                            oneWayRangeModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
                break;
            }
            case observation_models::angular_position:
            {
                boost::shared_ptr< observation_models::AngularPositionObservationModel
                        < ObservationScalarType, TimeType> > angularPositionModel =
                        boost::dynamic_pointer_cast< observation_models::AngularPositionObservationModel
                        < ObservationScalarType, TimeType> >
                        ( observationModelIterator->second );
                currentLightTimeCorrections.push_back(
                            angularPositionModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
                break;
            }
            case observation_models::position_observable:
            {
                break;
            }
            default:
                std::string errorMessage =
                        "Error in light time correction list creation, observable type " +
                        boost::lexical_cast< std::string >( observableType ) + " not recognized.";
                throw std::runtime_error( errorMessage );
            }

            // Add light-time correctionsfor current link ends
            lightTimeCorrectionsList[ observationModelIterator->first ] = currentLightTimeCorrections;
        }

    }
    return lightTimeCorrectionsList;
}

template< int ObservationSize >
void splitObservationPartialsAndScalers(
        const std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservationSize > > >,
        boost::shared_ptr< PositionPartialScaling > > >& observationPartialsAndScalers,
        std::map< observation_models::LinkEnds, std::map< std::pair< int, int >, boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > >& observationPartials,
        std::map< observation_models::LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > >& observationPartialScalers )
{
    observationPartials.clear( );
    observationPartialScalers.clear( );
    // Put one-way range partials and scalers in member variables.
    for( typename std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservationSize > > >,
         boost::shared_ptr< PositionPartialScaling > > >::const_iterator
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
template< int ObservationSize, typename ParameterType >
class ObservationPartialCreator
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< ObservationSize > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) );
};

//! Interface class for creating observation partials for observables of size 1.
template< typename ParameterType >
class ObservationPartialCreator< 1, ParameterType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< 1 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                boost::shared_ptr< ObservationPartial< 1 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;
        switch( observableType )
        {
        case observation_models::oneWayRange:
            observationPartialList = createOneWayRangePartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate, lightTimeCorrections );
            break;

        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    boost::lexical_cast< std::string >( observableType ) + " of size 1 ";
            throw std::runtime_error( errorMessage );
        }

        return observationPartialList;
    }
};

//! Interface class for creating observation partials for observables of size 2.
template< typename ParameterType >
class ObservationPartialCreator< 2, ParameterType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< 2 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                boost::shared_ptr< ObservationPartial< 2 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::angular_position:
            observationPartialList = createAngularPositionPartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate, lightTimeCorrections );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    boost::lexical_cast< std::string >( observableType ) + " of size 2 ";
            throw std::runtime_error( errorMessage );
        }
        return observationPartialList;
    }

};

//! Interface class for creating observation partials for observables of size 3.
template< typename ParameterType >
class ObservationPartialCreator< 3, ParameterType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< 3 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                boost::shared_ptr< ObservationPartial< 3 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::position_observable:
            observationPartialList = createPositionObservablePartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate );
            break;
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    boost::lexical_cast< std::string >( observableType ) + " of size 3 ";
            throw std::runtime_error( errorMessage );        }
        return observationPartialList;
    }

};

//! Interface class for creating observation partials for observables of size 3.
template< typename ParameterType >
class ObservationPartialCreator< 6, ParameterType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< 6 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                boost::shared_ptr< ObservationPartial< 6 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    boost::lexical_cast< std::string >( observableType ) + " of size 6 ";
            throw std::runtime_error( errorMessage );
        }
        return observationPartialList;
    }

};

//! Interface class for creating observation partials for observables of dynamic sizwe.
template< typename ParameterType >
class ObservationPartialCreator< Eigen::Dynamic, ParameterType >
{
public:

    //! Function to create a list of observation partial objects, and associated scaling objects
    /*!
     * Function to create a list of observation partial objects, and associated scaling objects
     * \param observableType Type of observable for which partials are to be created
     * \param linkEnds List of link ends of observable which partials are to be created
     * \param bodyMap Map of body objects that comprises the environment
     * \param parametersToEstimate Parameters for which partial derivatives are to be computed
     * \param lightTimeCorrections List of light-time corrections (per link ends set) for that are used in computation of
     * the observables
     * \return Map with list of observation partials. Key is associated link ends. Value is a list of observation partial
     * objects, one for each parameter w.r.t. which the observation partial is non-zero (in general). The format is a pair
     * with:
     * first: map with pair (parameter start entry in parameter vector and parameter size) and associated observation partial
     * object as value of map
     * second: PositionPartialScaling object associated with all partials of single LinkEnds.
     */
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
    boost::shared_ptr< ObservationPartial< Eigen::Dynamic > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >,
                boost::shared_ptr< ObservationPartial< Eigen::Dynamic > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        default:
            std::string errorMessage =
                    "Error when making observation partial set, could not recognize observable " +
                    boost::lexical_cast< std::string >( observableType ) + " of size dynamic ";
            throw std::runtime_error( errorMessage );        }
        return observationPartialList;
    }

};

}

}


#endif // TUDAT_CREATEOBSERVATIONPARTIALS_H
