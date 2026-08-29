/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace tudat
{

namespace propagators
{

namespace
{
template< typename TimeType, typename StateScalarType, int StateSize >
std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, StateSize, 1 > > >
createConfiguredStateInterpolator( const std::map< TimeType, Eigen::Matrix< StateScalarType, StateSize, 1 > >& stateMap,
                                   const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    if( interpolatorSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating integrated state interpolator, settings are null." );
    }

    return interpolators::createOneDimensionalInterpolator( stateMap, interpolatorSettings );
}
}  // namespace

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > > createStateInterpolator(
        const std::map< double, Eigen::Matrix< double, 6, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< double, double, 6 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > > createStateInterpolator(
        const std::map< double, Eigen::Matrix< long double, 6, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< double, long double, 6 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > > createStateInterpolator(
        const std::map< Time, Eigen::Matrix< long double, 6, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< Time, long double, 6 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > > createStateInterpolator(
        const std::map< Time, Eigen::Matrix< double, 6, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< Time, double, 6 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< double, Eigen::Matrix< double, 7, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< double, double, 7 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< long double, 7, 1 > >& stateMap,
                                   const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< double, long double, 7 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< Time, Eigen::Matrix< double, 7, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< Time, double, 7 >( stateMap, interpolatorSettings );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< Time, Eigen::Matrix< long double, 7, 1 > >& stateMap,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    return createConfiguredStateInterpolator< Time, long double, 7 >( stateMap, interpolatorSettings );
}

}  // namespace propagators

}  // namespace tudat
