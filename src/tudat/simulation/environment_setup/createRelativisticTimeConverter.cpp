#include <cmath>

#include "tudat/basics/utilities.h"
#include "tudat/astro/propagators/integrateEquations.h"
#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"


namespace tudat
{

namespace simulation_setup
{

template< typename StateScalarType, typename TimeType >
void setRelativisticTimeConverter(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > >& conversionSettings,
        const SystemOfBodies& bodyMap )
{
    const auto barycentricSettings = conversionSettings->getBaryCentricToBodyCentricConversionSettings( );
    const auto topocentricConversions = conversionSettings->getBodyCentricToTopocentricConversionSettings( );

    const TimeType initialTime = barycentricSettings->getInitialTime( );
    if( !std::isfinite( static_cast< double >( initialTime ) ) )
    {
        throw std::runtime_error( "Error in setRelativisticTimeConverter: initial propagation time is not finite." );
    }

    const auto baseIntegratorSettings = conversionSettings->getNumericalIntegrationSettings( );
    if( baseIntegratorSettings == nullptr )
    {
        throw std::runtime_error( "Error in setRelativisticTimeConverter: no integrator settings provided." );
    }

    barycentricSettings->getOutputSettings( )->setIntegratedResult( true );

    auto barycentricIntegratorSettings = baseIntegratorSettings->clone( );
    barycentricIntegratorSettings->initialTimeDeprecated_ = initialTime;
    barycentricSettings->setIntegratorSettings( barycentricIntegratorSettings );

    std::cerr << "[relativistic_time] creating barycentric simulator for "
              << barycentricSettings->getReferencePointId( ).first << std::endl;
    propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > barycentricSimulator(
        bodyMap,
        barycentricSettings,
        true );
    std::cerr << "[relativistic_time] barycentric simulator created" << std::endl;

    for( const auto& topocentricSettings : topocentricConversions )
    {
        if ( topocentricSettings->getRelativisticStateDerivativeType( ) != propagators::first_order_bodycentric_to_topocentric )
        {
            throw std::runtime_error(
                "Error in setRelativisticTimeConverter: inconsistent derivative type for topocentric conversion." );
        }

        const TimeType topocentricInitialTime = topocentricSettings->getInitialTime( );
        if( !std::isfinite( static_cast< double >( topocentricInitialTime ) ) )
        {
            throw std::runtime_error( "Error in setRelativisticTimeConverter: topocentric propagation time is not finite." );
        }

        topocentricSettings->getOutputSettings( )->setIntegratedResult( true );

        auto topocentricIntegratorSettings = baseIntegratorSettings->clone( );
        topocentricIntegratorSettings->initialTimeDeprecated_ = topocentricInitialTime;
        topocentricSettings->setIntegratorSettings( topocentricIntegratorSettings );

        std::cerr << "[relativistic_time] creating topocentric simulator for "
                  << topocentricSettings->getReferencePointId( ).second << std::endl;
        propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > topocentricSimulator(
            bodyMap,
            topocentricSettings,
            true );
        std::cerr << "[relativistic_time] topocentric simulator created" << std::endl;
    }
}

template< typename StateScalarType, typename TimeType >
void setRelativisticTimeConverters(
        const SystemOfBodies& bodyMap,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > > >& converterSettings )
{
    for ( const auto& converterEntry : converterSettings )
    {
        setRelativisticTimeConverter( converterEntry.second, bodyMap );
    }
}

} // namespace simulation_setup

} // namespace tudat

namespace tudat
{
namespace simulation_setup
{

template void setRelativisticTimeConverter<double, double>(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > >& conversionSettings,
        const SystemOfBodies& bodyMap );

template void setRelativisticTimeConverters<double, double>(
        const SystemOfBodies& bodyMap,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > > >& converterSettings );

template void setRelativisticTimeConverter<double, Time>(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< double, Time > >& conversionSettings,
        const SystemOfBodies& bodyMap );

template void setRelativisticTimeConverters<double, Time>(
        const SystemOfBodies& bodyMap,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< double, Time > > >& converterSettings );

} // namespace simulation_setup
} // namespace tudat
