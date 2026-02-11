#ifndef TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H
#define TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <iostream>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/propagators/relativisticTimeStateDerivative.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"

namespace tudat
{
namespace simulation_setup
{

template< typename StateScalarType = double, typename TimeType = double >
class DirectRelativisticTimeConverterSettings
{
public:

    DirectRelativisticTimeConverterSettings(
        const std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > >& baryCentricToBodyCentricConversionSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& numericalIntegrationSettings,
        const std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >& bodyCentricToTopocentricConversionSettings =
            std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > >( ) )
        : baryCentricToBodyCentricConversionSettings_( baryCentricToBodyCentricConversionSettings ),
          bodyCentricToTopocentricConversionSettings_( bodyCentricToTopocentricConversionSettings ),
          numericalIntegrationSettings_( numericalIntegrationSettings )
    {
        using namespace propagators;

        if( !( baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) == first_order_barycentric_to_bodycentric ||
               baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) == second_order_barycentric_to_bodycentric ) )
        {
            std::cerr<<"[ERROR] Invalid type for barycentric-to-bodycentric conversion: "
                     << baryCentricToBodyCentricConversionSettings_->getRelativisticStateDerivativeType( ) <<std::endl;
        }

        for( unsigned int i = 0; i < bodyCentricToTopocentricConversionSettings_.size( ); ++i )
        {
            auto topocentricSetting = bodyCentricToTopocentricConversionSettings_.at( i );
            if( topocentricSetting->getRelativisticStateDerivativeType( ) != first_order_bodycentric_to_topocentric )
            {
                std::cerr<<"[ERROR] Invalid type for topocentric setting for station "
                         << topocentricSetting->getReferencePointId( ).second << ": "
                         << topocentricSetting->getRelativisticStateDerivativeType( ) <<std::endl;
            }

            if( topocentricSetting->getReferencePointId( ).first !=
                baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first )
            {
                std::cerr<<"[ERROR] Station "
                         << topocentricSetting->getReferencePointId( ).second
                         << " is on body "
                         << topocentricSetting->getReferencePointId( ).first
                         << ", but barycentric setting uses body "
                         << baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first
                         << std::endl;
            }
        }

        associatedBody_ = baryCentricToBodyCentricConversionSettings_->getReferencePointId( ).first;
    }

    std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > getBaryCentricToBodyCentricConversionSettings( ) const
    {
        return baryCentricToBodyCentricConversionSettings_;
    }

    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > > getBodyCentricToTopocentricConversionSettings( ) const
    {
        return bodyCentricToTopocentricConversionSettings_;
    }

    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getNumericalIntegrationSettings( ) const
    {
        return numericalIntegrationSettings_;
    }

    std::string getAssociatedBody( ) const
    {
        return associatedBody_;
    }

private:

    std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > baryCentricToBodyCentricConversionSettings_;

    std::vector< std::shared_ptr< propagators::RelativisticTimeStatePropagatorSettings< StateScalarType, TimeType > > > bodyCentricToTopocentricConversionSettings_;

    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > numericalIntegrationSettings_;

    std::string associatedBody_;
};

//! Setup a single relativistic time converter (barycentric → bodycentric → topocentric)
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverter(
    const std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > >& conversionSettings,
    const SystemOfBodies& bodies );

//! Setup multiple relativistic time converters (one per associated body)
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverters(
    const SystemOfBodies& bodies,
    const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > > >& converterSettings );

} // namespace simulation_setup
} // namespace tudat

#endif // TUDAT_CREATE_RELATIVISTIC_TIME_CONVERTER_H
