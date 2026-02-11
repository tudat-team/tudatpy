#include <algorithm>
#include <sstream>
#include <typeinfo>
#include <vector>

#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"

namespace tudat
{

namespace simulation_setup
{

namespace detail
{

template< typename StateScalarType, typename TimeType, typename InitialStateParameterType >
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartialsImpl(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< StateScalarType, TimeType > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parametersToEstimate )
{
    orbit_determination::StateDerivativePartialsMap relativisticStatePartials;

    const auto bodiesToEstimate =
            estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate( parametersToEstimate );

    if( bodiesToEstimate.count( propagators::proper_time ) == 0 )
    {
        return relativisticStatePartials;
    }

    const auto& referencePointsToEstimate = bodiesToEstimate.at( propagators::proper_time );
    relativisticStatePartials.resize( referencePointsToEstimate.size( ) );

    for( const auto& stateDerivativeModel : stateDerivativeModels )
    {
        auto relativisticTimeModel =
                std::dynamic_pointer_cast< RelativisticTimeStateDerivative< StateScalarType, TimeType > >(
                        stateDerivativeModel );
        if( relativisticTimeModel == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating relativistic time state derivative partials, model type is not relativistic time derivative." );
        }

        auto stateDerivativePartial =
                createRelativisticTimeStateDerivativePartial< StateScalarType, TimeType, InitialStateParameterType >(
                        relativisticTimeModel, bodies, parametersToEstimate );

        const auto referencePoint = relativisticTimeModel->getReferencePoint( );
        auto referenceIterator = std::find( referencePointsToEstimate.begin( ), referencePointsToEstimate.end( ), referencePoint );
        if( referenceIterator == referencePointsToEstimate.end( ) )
        {
            std::ostringstream errorMessage;
            errorMessage << "Error when creating relativistic time state derivative partials. Requested reference point ("
                         << referencePoint.first << ", " << referencePoint.second
                         << ") is not part of the estimated proper-time states.";
            throw std::runtime_error( errorMessage.str( ) );
        }

        const int parameterIndex = static_cast< int >( std::distance( referencePointsToEstimate.begin( ), referenceIterator ) );
        relativisticStatePartials[ parameterIndex ].push_back( stateDerivativePartial );
    }

    return relativisticStatePartials;
}

template< typename StateScalarType, typename TimeType, typename InitialStateParameterType >
void throwUnsupportedCombinationError( )
{
    std::ostringstream errorMessage;
    errorMessage << "Error when creating relativistic time state derivative partials. The combination of state scalar type ("
                 << typeid( StateScalarType ).name( ) << "), time type (" << typeid( TimeType ).name( )
                 << "), and parameter scalar type (" << typeid( InitialStateParameterType ).name( )
                 << ") is currently not supported.";
    throw std::runtime_error( errorMessage.str( ) );
}

}  // namespace detail

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< double, double, double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, double > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate )
{
    return detail::createRelativisticTimeStateDerivativePartialsImpl( stateDerivativeModels, bodies, parametersToEstimate );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< double, double, long double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, double > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< double, double, long double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< double, Time, double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, Time > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< double, Time, double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< double, Time, long double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, Time > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< double, Time, long double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< long double, double, double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, double > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< long double, double, double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< long double, double, long double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, double > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< long double, double, long double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< long double, Time, double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, Time > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< long double, Time, double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

template<>
orbit_determination::StateDerivativePartialsMap createRelativisticTimeStateDerivativePartials< long double, Time, long double >(
        const std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, Time > > >& stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate )
{
    detail::throwUnsupportedCombinationError< long double, Time, long double >( );
    return orbit_determination::StateDerivativePartialsMap( );
}

}  // namespace simulation_setup

}  // namespace tudat
