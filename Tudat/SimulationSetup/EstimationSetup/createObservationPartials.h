#ifndef CREATEOBSERVATIONPARTIALS_H
#define CREATEOBSERVATIONPARTIALS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/angularPositionObservationModel.h"

#include "Tudat/SimulationSetup/EstimationSetup/createAngularPositionPartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"

namespace tudat
{

namespace observation_partials
{

typedef std::map< observation_models::LinkEnds, std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > > PerLinkEndPerLightTimeSolutionCorrections;


template< typename ObservationScalarType, typename TimeType, typename StateScalarType, int ObservationSize  >
PerLinkEndPerLightTimeSolutionCorrections getLightTimeCorrectionsList(
        const std::map< observation_models::LinkEnds, boost::shared_ptr< observation_models::ObservationModel<
        ObservationSize, ObservationScalarType, TimeType, StateScalarType > > > observationModels )
{
    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionsList;
    observation_models::ObservableType observableType = observationModels.begin( )->second->getObservableType( );

    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;

    for( typename  std::map< observation_models::LinkEnds, boost::shared_ptr< observation_models::ObservationModel<
         ObservationSize, ObservationScalarType, TimeType, StateScalarType > > >::const_iterator observationModelIterator = observationModels.begin( );
         observationModelIterator != observationModels.end( ); observationModelIterator++ )
    {
        currentLightTimeCorrections.clear( );

        if( observationModelIterator->second->getObservableType( ) != observableType )
        {
            std::cerr<<"Error when making grouped light time correction list, observable type is not constant"<<std::endl;
        }
        else
        {
            switch( observableType )
            {
            case observation_models::oneWayRange:
            {
                boost::shared_ptr< observation_models::OneWayRangeObservationModel< ObservationScalarType, TimeType, StateScalarType > > oneWayRangeModel =
                        boost::dynamic_pointer_cast< observation_models::OneWayRangeObservationModel< ObservationScalarType, TimeType, StateScalarType > >
                        ( observationModelIterator->second );
                currentLightTimeCorrections.push_back( oneWayRangeModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
                break;
            }
            case observation_models::angular_position:
            {
                boost::shared_ptr< observation_models::AngularPositionObservationModel< ObservationScalarType, TimeType, StateScalarType > > angularPositionModel =
                        boost::dynamic_pointer_cast< observation_models::AngularPositionObservationModel< ObservationScalarType, TimeType, StateScalarType > >
                        ( observationModelIterator->second );
                currentLightTimeCorrections.push_back( angularPositionModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
                break;
            }
            case observation_models::position_observable:
            {
                break;
            }
            default:
                std::cerr<<"Error in light time correction list creation, observable type "<<observableType<<" not recognized"<<std::endl;
            }
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

template< int ObservationSize, typename ParameterType >
class ObservationPartialCreator
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservationSize > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections = PerLinkEndPerLightTimeSolutionCorrections( ) );
};

template< typename ParameterType >
class ObservationPartialCreator< 1, ParameterType >
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections = PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;
        switch( observableType )
        {
        case observation_models::oneWayRange:
            observationPartialList = createOneWayRangePartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate, lightTimeCorrections );
            break;

        default:
            std::cerr<<"Error when making observation partial set, could not recognize observable "<<observableType<<std::endl;
        }

        return observationPartialList;
    }
};

template< typename ParameterType >
class ObservationPartialCreator< 2, ParameterType >
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 2 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections = PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 2 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::angular_position:
            observationPartialList = createAngularPositionPartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate, lightTimeCorrections );
            break;
        default:
            std::cerr<<"Error when making observation partial set, could not recognize observable "<<observableType<<std::endl;
        }
        return observationPartialList;
    }

};

template< typename ParameterType >
class ObservationPartialCreator< 3, ParameterType >
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 3 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections = PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 3 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        case observation_models::position_observable:
            observationPartialList = createPositionObservablePartials< ParameterType >(
                        linkEnds, bodyMap, parametersToEstimate );
            break;
        default:
            std::cerr<<"Error when making observation partial set of 3 dimensions, could not recognize observable "<<observableType<<std::endl;
        }
        return observationPartialList;
    }

};

template< typename ParameterType >
class ObservationPartialCreator< 6, ParameterType >
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 6 > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections =
            PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 6 > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        default:
            std::cerr<<"Error when making observation partial set, could not recognize observable "<<observableType<<std::endl;
        }
        return observationPartialList;
    }

};

template< typename ParameterType >
class ObservationPartialCreator< Eigen::Dynamic, ParameterType >
{
public:
    static std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< Eigen::Dynamic > > >,
    boost::shared_ptr< PositionPartialScaling > > > createObservationPartials(
            const observation_models::ObservableType observableType,
            const std::vector< observation_models::LinkEnds >& linkEnds,
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrections = PerLinkEndPerLightTimeSolutionCorrections( ) )
    {
        std::map< observation_models::LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< Eigen::Dynamic > > >,
                boost::shared_ptr< PositionPartialScaling > > > observationPartialList;

        switch( observableType )
        {
        default:
            std::cerr<<"Error when making observation partial set, could not recognize observable "<<observableType<<std::endl;
        }
        return observationPartialList;
    }

};

}

}


#endif // CREATEOBSERVATIONPARTIALS_H
