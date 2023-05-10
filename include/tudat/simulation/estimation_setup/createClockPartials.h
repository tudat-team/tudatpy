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

template< typename ParameterType >
std::map< int, std::shared_ptr< TimingPartial > > createTimingPartialWrtClockProperty(
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< observation_models::ObservationBias< 1 > > clockInducedBiases )
{
    std::map< int, std::shared_ptr< TimingPartial > > timingPartials;

    for( observation_models::LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        if( ( linkEndIterator->second.bodyName_ == parameterToEstimate->getParameterName( ).second.first ) &&
                ( linkEndIterator->second.stationName_ == parameterToEstimate->getParameterName( ).second.second ) )
        {
            std::shared_ptr< system_models::TimingSystem > timingSystem =
                    simulation_setup::getTimingSystem( parameterToEstimate->getParameterName( ).second, clockInducedBiases );
            if( timingSystem != nullptr )
            {
                switch( parameterToEstimate->getParameterName( ).first )
                {
                case estimatable_parameters::global_polynomial_clock_corrections:
                {
                    std::shared_ptr< estimatable_parameters::GlobalPolynomialClockCorrections > globalPolynomialClockCorrections =
                            std::dynamic_pointer_cast< estimatable_parameters::GlobalPolynomialClockCorrections >( parameterToEstimate );

                    std::cout<<"Clock powers ";
                    for( unsigned int i = 0; i < globalPolynomialClockCorrections->getCorrectionPowers( ).size( ); i++ )
                    {
                        std::cout<<globalPolynomialClockCorrections->getCorrectionPowers( ).at( i )<<" ";
                    }
                    std::cout<<std::endl;

                    if( timingSystem.get( ) != globalPolynomialClockCorrections->getTimingSystem( ).get( ) )
                    {
                        std::cerr<<"Error when making range partial w.r.t. polynimial global clock variation parameter of "<<
                                   parameterToEstimate->getParameterName( ).second.first<<" "<<parameterToEstimate->getParameterName( ).second.second<<
                                   ", pointers are inconsistent "<<timingSystem.get( )<<" "<<globalPolynomialClockCorrections->getTimingSystem( ).get( )<<std::endl;                        }

                    std::vector< int > linkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                                observableType, linkEndIterator->first, linkEnds.size( ) );
                    for( unsigned int i = 0; i < linkEndIndices.size( ); i++ )
                    {
                        timingPartials[ linkEndIndices.at( i ) ] =
                                std::make_shared< TimingPartialWrtGlobalPolynomialCorrectionCoefficients >(
                                    timingSystem, globalPolynomialClockCorrections->getCorrectionPowers( ) );
                    }

                    break;
                }
                case estimatable_parameters::arc_wise_polynomial_clock_corrections:
                {
                    std::shared_ptr< estimatable_parameters::MultiArcClockCorrections > arcWisePolynomialClockCorrections =
                            std::dynamic_pointer_cast< estimatable_parameters::MultiArcClockCorrections >( parameterToEstimate );

                    if( timingSystem.get( ) != arcWisePolynomialClockCorrections->getTimingSystem( ).get( ) )
                    {
                        std::cerr<<"Error when making range partial w.r.t. polynimial global clock variation parameter of "<<
                                   parameterToEstimate->getParameterName( ).second.first<<" "<<parameterToEstimate->getParameterName( ).second.second
                                <<", pointers are inconsistent"<<timingSystem.get( )<<" "<<arcWisePolynomialClockCorrections->getTimingSystem( ).get( )<<std::endl;                        }

                    std::vector< int > linkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                                observableType, linkEndIterator->first, linkEnds.size( ) );
                    for( unsigned int i = 0; i < linkEndIndices.size( ); i++ )
                    {
                        timingPartials[ linkEndIndices.at( i ) ] =
                                std::make_shared< TimingPartialWrtArcwisePolynomialCorrectionCoefficients >(
                                    timingSystem, arcWisePolynomialClockCorrections->getCorrectionPowers( ),
                                    arcWisePolynomialClockCorrections->getArcIndices( ) );
                    }

                    break;
                }
                default:
                    std::cerr<<"Error when making range partial w.r.t. clock parameter, parameter "<<
                               parameterToEstimate->getParameterName( ).first<<" not found"<<std::endl;
                }
            }
        }
    }
    return timingPartials;
}

template< typename ParameterType >
std::shared_ptr< ObservationPartial< 1 > > createPartialWrtClockProperty(
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< observation_models::ObservationBias< 1 > > clockInducedBiases )
{
    std::map< int, std::shared_ptr< TimingPartial > > timingPartialList = createTimingPartialWrtClockProperty(
                linkEnds, observableType, parameterToEstimate, clockInducedBiases );
    std::shared_ptr< ObservationPartial< 1 > > observationPartial;
    if( timingPartialList.size( ) > 0 )
    {
        observationPartial = std::make_shared< RangePartialWrtClockParameter >(
                    parameterToEstimate->getParameterName( ), timingPartialList, getTimingPartialMultipliers( observableType ) );
    }

    return observationPartial;
}

}

}

#endif // CREATECLOCKPARTIALS_H
