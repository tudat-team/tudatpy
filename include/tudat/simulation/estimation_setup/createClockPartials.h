#ifndef CREATECLOCKPARTIALS_H
#define CREATECLOCKPARTIALS_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/orbit_determination/observation_partials/clockParameterPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters//polynomialClockCorrections.h"


namespace tudat
{

namespace observation_partials
{

std::map< int, double > getTimingPartialMultipliers( const observation_models::ObservableType observable );

template< typename ParameterType, int ObservationSize >
std::map< int, std::shared_ptr< TimingPartial > > createTimingPartialWrtClockProperty(
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::vector< std::shared_ptr< observation_models::ObservationBias< ObservationSize > > > clockInducedBiases )
{
    std::map< int, std::shared_ptr< TimingPartial > > timingPartials;

    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        if( ( linkEndIterator->second.bodyName_ == parameterToEstimate->getParameterName( ).second.first ) &&
                ( linkEndIterator->second.stationName_ == parameterToEstimate->getParameterName( ).second.second ) )
        {
            for( unsigned int j = 0; j < clockInducedBiases.size( ); j++ )
            {
                std::shared_ptr< observation_models::ClockInducedRangeBias< 1 > > clockInducedBias =
                        std::dynamic_pointer_cast< observation_models::ClockInducedRangeBias< 1 > >( clockInducedBiases.at( j ) );
                if( clockInducedBias == nullptr )
                {
                    throw std::runtime_error( "Error when making clock partials, clock-induced bias not found." );
                }
                std::shared_ptr< system_models::TimingSystem > timingSystem = clockInducedBias->getTimingSystem( );
                if( timingSystem != nullptr )
                {
                    switch ( parameterToEstimate->getParameterName( ).first )
                    {
                        case estimatable_parameters::global_polynomial_clock_corrections:
                        {
                            std::shared_ptr<estimatable_parameters::GlobalPolynomialClockCorrections> globalPolynomialClockCorrections =
                                    std::dynamic_pointer_cast<estimatable_parameters::GlobalPolynomialClockCorrections>(
                                            parameterToEstimate );

                            if( globalPolynomialClockCorrections->getLinkEndId( ) == clockInducedBias->getLinkEndId( ) )
                            {
                                if ( timingSystem.get( ) != globalPolynomialClockCorrections->getTimingSystem( ).get( ))
                                {
                                    std::cout<<timingSystem.get( )<<" "<<globalPolynomialClockCorrections->getTimingSystem( ).get( )<<std::endl;
                                    throw std::runtime_error(
                                             "Error when making range partial w.r.t. polynimial global clock variation parameter of " +
                                            parameterToEstimate->getParameterName( ).second.first  + ", " +
                                            parameterToEstimate->getParameterName( ).second.second +
                                            ", pointers are inconsistent." );
                                }

                                std::vector<int> linkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                                        observableType, linkEndIterator->first, linkEnds.size( ));
                                for ( unsigned int i = 0; i < linkEndIndices.size( ); i++ ) {
                                    timingPartials[ linkEndIndices.at( i ) ] =
                                            std::make_shared<TimingPartialWrtGlobalPolynomialCorrectionCoefficients>(
                                                    timingSystem, globalPolynomialClockCorrections->getCorrectionPowers( ));
                                }
                            }
                            break;
                        }
                        case estimatable_parameters::arc_wise_polynomial_clock_corrections:
                        {
                            std::shared_ptr<estimatable_parameters::MultiArcClockCorrections> arcWisePolynomialClockCorrections =
                                    std::dynamic_pointer_cast<estimatable_parameters::MultiArcClockCorrections>(
                                            parameterToEstimate );

                            if( arcWisePolynomialClockCorrections->getLinkEndId( ) == clockInducedBias->getLinkEndId( ) )
                            {
                                if ( timingSystem.get( ) != arcWisePolynomialClockCorrections->getTimingSystem( ).get( ))
                                {

                                    throw std::runtime_error(
                                            "Error when making range partial w.r.t. polynimial arcwise clock variation parameter of " +
                                            parameterToEstimate->getParameterName( ).second.first  + ", " +
                                            parameterToEstimate->getParameterName( ).second.second +
                                            ", pointers are inconsistent" );
                                }

                                std::vector<int> linkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                                        observableType, linkEndIterator->first, linkEnds.size( ));
                                for ( unsigned int i = 0; i < linkEndIndices.size( ); i++ ) {
                                    timingPartials[ linkEndIndices.at( i ) ] =
                                            std::make_shared<TimingPartialWrtArcwisePolynomialCorrectionCoefficients>(
                                                    timingSystem, arcWisePolynomialClockCorrections->getCorrectionPowers( ),
                                                    arcWisePolynomialClockCorrections->getArcIndices( ));
                                }
                            }

                            break;
                        }
                        default:
                            std::cerr << "Error when making range partial w.r.t. clock parameter, parameter " <<
                                      parameterToEstimate->getParameterName( ).first << " not found" << std::endl;
                    }
                }
                else
                {
                    throw std::runtime_error( "Error, timing system not found when creating clock partials" );
                }
            }
        }
    }
    return timingPartials;
}

template< typename ParameterType, int ObservationSize >
class ObservationPartialWrtClockCreator
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
    static std::shared_ptr< ObservationPartial< ObservationSize > > createPartialWrtClockProperty(
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType,
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationBias< ObservationSize > > > clockInducedBiases )
    {
        return nullptr;
    }
};

template< typename ParameterType >
class ObservationPartialWrtClockCreator< ParameterType, 1 >
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
    static std::shared_ptr< ObservationPartial< 1 > > createPartialWrtClockProperty(
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType,
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationBias< 1 > > > clockInducedBiases )
    {

        std::shared_ptr< ObservationPartial< 1 > > observationPartial = nullptr;
        switch( observableType )
        {
            case observation_models::one_way_range:
            {
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialList = createTimingPartialWrtClockProperty(
                        linkEnds, observableType, parameterToEstimate, clockInducedBiases );
                if( timingPartialList.size( ) > 0 )
                {
                    observationPartial = std::make_shared< RangePartialWrtClockParameter >(
                            parameterToEstimate->getParameterName( ), timingPartialList, getTimingPartialMultipliers( observableType ) );

                }
                break;
            }
            case observation_models::n_way_range:
            {
                std::map< int, std::shared_ptr< TimingPartial > > timingPartialList = createTimingPartialWrtClockProperty(
                    linkEnds, observableType, parameterToEstimate, clockInducedBiases );
                if( timingPartialList.size( ) > 0 )
                {
                    observationPartial = std::make_shared< RangePartialWrtClockParameter >(
                        parameterToEstimate->getParameterName( ), timingPartialList, getTimingPartialMultipliers( observableType ) );

                }
                break;
            }
            default:
                break;
        }

        return observationPartial;
    }
};






}

}

#endif // CREATECLOCKPARTIALS_H
