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


//! Function to create an interpolator for the new translational state of a body.
template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< double, 6, 1 > >& stateMap )
{

    return std::make_shared<
            interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );

}

//! Function to create an interpolator for the new translational state of a body.
template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return std::make_shared<
            interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}

//! Function to create an interpolator for the new translational state of a body.
template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > >
createStateInterpolator( const std::map< Time, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return std::make_shared<
            interpolators::LagrangeInterpolator< Time, Eigen::Matrix< long double, 6, 1 >, long double > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}


//! Function to create an interpolator for the new translational state of a body.
template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > >
createStateInterpolator( const std::map< Time, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return std::make_shared<
            interpolators::LagrangeInterpolator< Time, Eigen::Matrix< double, 6, 1 >, long double > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}


template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 7, 1 > > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}

template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}

template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< Time, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< double, 7, 1 >, long double > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}

template< >
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< Time, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< long double, 7, 1 >, long double > >(
                stateMap, 6, interpolators::huntingAlgorithm, interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary );
}

} // namespace propagators

} // namespace tudat
