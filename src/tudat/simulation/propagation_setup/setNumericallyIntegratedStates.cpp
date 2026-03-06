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
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/ephemerides/timeEphemerisDirectFromMetric.h"
#include "tudat/astro/ephemerides/timeEphemerisWithFirstOrderDirectConversion.h"
#include <stdexcept>

namespace tudat
{

namespace propagators
{



std::pair<
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >,
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >
> createRelativisticTimeInterpolators(std::map< double, double >& originalToTargetTimeMap)
{
    std::map< double, double > targetToOriginalTimeMap;
    for( auto mapIterator = originalToTargetTimeMap.begin(); mapIterator != originalToTargetTimeMap.end(); mapIterator++ )
    {
        targetToOriginalTimeMap[ mapIterator->first + mapIterator->second ] = -mapIterator->second;
    }
    return std::make_pair(
        std::make_shared< interpolators::LinearInterpolator< double, double > >( originalToTargetTimeMap ),
        std::make_shared< interpolators::LinearInterpolator< double, double > >( targetToOriginalTimeMap )
    );
}

std::pair<
    std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >,
    std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >
> createRelativisticTimeInterpolators(std::map< Time, double >& originalToTargetTimeMap)
{
    std::map< Time, double > targetToOriginalTimeMap;
    for( auto mapIterator = originalToTargetTimeMap.begin(); mapIterator != originalToTargetTimeMap.end(); mapIterator++ )
    {
        targetToOriginalTimeMap[ mapIterator->first + mapIterator->second ] = -mapIterator->second;
    }
    return std::make_pair(
        std::make_shared< interpolators::LinearInterpolator< Time, double > >( originalToTargetTimeMap ),
        std::make_shared< interpolators::LinearInterpolator< Time, double > >( targetToOriginalTimeMap )
    );
}

//! Function to create an interpolator for the new translational state of a body.
template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > > createStateInterpolator(
        const std::map< double, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

//! Function to create an interpolator for the new translational state of a body.
template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > > createStateInterpolator(
        const std::map< double, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

//! Function to create an interpolator for the new translational state of a body.
template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > > createStateInterpolator(
        const std::map< Time, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< long double, 6, 1 >, long double > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

//! Function to create an interpolator for the new translational state of a body.
template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > > createStateInterpolator(
        const std::map< Time, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< double, 6, 1 >, long double > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< double, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 7, 1 > > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< Time, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< double, 7, 1 >, long double > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}

template<>
std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 7, 1 > > > createRotationalStateInterpolator(
        const std::map< Time, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return std::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< long double, 7, 1 >, long double > >(
            stateMap,
            6,
            interpolators::huntingAlgorithm,
            interpolators::lagrange_cubic_spline_boundary_interpolation,
            interpolators::throw_exception_at_boundary );
}



template< >
void resetIntegratedDirectFromMetricTimeEphemeris< double, double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    if( startIndexAndSize.second != 1 )
    {
        throw std::runtime_error(
                    "Error when resetting integrated time ephemeris, found requested size " +
                    std::to_string( startIndexAndSize.second ) );
    }
    std::map< double, double > floatingPointValueNumericalSolution;

    for( typename std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator solutionIterator =
         numericalSolution.begin( ); solutionIterator != numericalSolution.end( ); solutionIterator++ )
    {
        floatingPointValueNumericalSolution[ solutionIterator->first ] =  solutionIterator->second( startIndexAndSize.first );
    }

    std::pair< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > > timeInterpolators =
            createRelativisticTimeInterpolators( floatingPointValueNumericalSolution );

    if( bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) == nullptr )
    {
        std::shared_ptr< TimeEphemeris > newTimeEphemeris = std::make_shared< TimeEphemerisDirectFromMetric >(
                    referencePointIdentifier.first );
        bodies.getBody( referencePointIdentifier.first )->setTimeScaleConverter( newTimeEphemeris );
    }

    std::shared_ptr< TimeEphemerisDirectFromMetric > timeEphemeris =
            std::dynamic_pointer_cast< TimeEphemerisDirectFromMetric >(
                bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) );
    if( timeEphemeris == nullptr )
    {
        throw std::runtime_error(
                    "Error when resetting integrated direct-from-metric time ephemeris, no TimeEphemerisDirectFromMetric object found." );
    }

    timeEphemeris->resetGlobalToProperTimeInterpolators(
                timeInterpolators.first, timeInterpolators.second, referencePointIdentifier.second );
}

template< >
void resetIntegratedDirectFromMetricTimeEphemeris< Time, double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< Time, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    if( startIndexAndSize.second != 1 )
    {
        throw std::runtime_error(
                    "Error when resetting integrated time ephemeris, found requested size " +
                    std::to_string( startIndexAndSize.second ) );
    }

    std::map< Time, double > floatingPointValueNumericalSolution;
    for( const auto& solutionEntry : numericalSolution )
    {
        floatingPointValueNumericalSolution[ solutionEntry.first ] = solutionEntry.second( startIndexAndSize.first );
    }

    std::pair< std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > > > timeInterpolators =
            createRelativisticTimeInterpolators( floatingPointValueNumericalSolution );

    if( bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) == nullptr )
    {
        std::shared_ptr< TimeEphemeris > newTimeEphemeris = std::make_shared< TimeEphemerisDirectFromMetric >(
                    referencePointIdentifier.first );
        bodies.getBody( referencePointIdentifier.first )->setTimeScaleConverter( newTimeEphemeris );
    }

    std::shared_ptr< TimeEphemerisDirectFromMetric > timeEphemeris =
            std::dynamic_pointer_cast< TimeEphemerisDirectFromMetric >(
                bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) );
    if( timeEphemeris == nullptr )
    {
        throw std::runtime_error(
                    "Error when resetting integrated direct-from-metric time ephemeris, no TimeEphemerisDirectFromMetric object found." );
    }

    timeEphemeris->resetGlobalToProperTimeInterpolators(
                timeInterpolators.first, timeInterpolators.second, referencePointIdentifier.second );
}

template< >
void resetIntegratedDirectFromMetricTimeEphemeris< double, long double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    throw std::runtime_error( "Cannot reset integrated direct from metric time ephemeris for (double, long double)" );
}

template< >
void resetIntegratedDirectFromMetricTimeEphemeris< Time, long double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    throw std::runtime_error( "Cannot reset integrated direct from metric time ephemeris for (Time, long double)" );
}


template< >
void resetIntegratedPostNewtonianTimeEphemeris< double, double >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    // resetIntegratedPostNewtonianTimeEphemeris
    if( startIndexAndSize.second != 1 )
    {
        throw std::runtime_error(
                    "Error when resetting integrated time ephemeris, found requested size " +
                    std::to_string( startIndexAndSize.second ) );
    }
    std::map< double, double > floatingPointValueNumericalSolution;

    for( typename std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator solutionIterator =
         numericalSolution.begin( ); solutionIterator != numericalSolution.end( ); solutionIterator++ )
    {
        floatingPointValueNumericalSolution[ solutionIterator->first ] =  solutionIterator->second( startIndexAndSize.first );
    }

    std::pair< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > > timeInterpolators =
            createRelativisticTimeInterpolators( floatingPointValueNumericalSolution );

    if( bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) == nullptr )
    {

        std::shared_ptr< TimeEphemeris > newTimeEphemeris = std::make_shared< TimeEphemerisWithFirstOrderDirectConversion >(
                    TimeEphemerisFromPostNewtonianExpansion::TimeDifferenceInterpolator( ),
                    TimeEphemerisFromPostNewtonianExpansion::TimeDifferenceInterpolator( ), referencePointIdentifier.first,
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >, bodies.getBody( referencePointIdentifier.first ), std::placeholders::_1 ) );
        bodies.getBody( referencePointIdentifier.first )->setTimeScaleConverter( newTimeEphemeris );
    }

    std::shared_ptr< TimeEphemerisFromPostNewtonianExpansion > timeEphemeris =
            std::dynamic_pointer_cast< TimeEphemerisFromPostNewtonianExpansion >(
                bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) );
    if( timeEphemeris == nullptr )
    {
        throw std::runtime_error(
                    "Error when resetting integrated post-Newtonian time ephemeris, no TimeEphemerisFromPostNewtonianExpansion object found." );
    }

    if( referencePointIdentifier.second == "" )
    {
        timeEphemeris->resetBarycentricToBodycentricInterpolators( timeInterpolators.first, timeInterpolators.second );
    }
    else
    {
        if( timeEphemeris->doesReferencePointTopocentricConverterExist( referencePointIdentifier.second ) )
        {
            timeEphemeris->resetBodycentricToTopocentricInterpolators(
                        timeInterpolators.first, timeInterpolators.second, referencePointIdentifier.second  );
        }
        else
        {
            std::shared_ptr< ground_stations::GroundStationState > groundStationState =
                    std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( referencePointIdentifier.first ) )->getGroundStation(
                        referencePointIdentifier.second )->getNominalStationState( );
            timeEphemeris->resetBodycentricToTopocentricInterpolators(
                        timeInterpolators.first, timeInterpolators.second, referencePointIdentifier.second,
                        std::bind( &ground_stations::GroundStationState::getNominalCartesianPosition, groundStationState ) );
        }
    }
}


template< >
void resetIntegratedPostNewtonianTimeEphemeris< Time, double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< Time, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    if( startIndexAndSize.second != 1 )
    {
        throw std::runtime_error(
                    "Error when resetting integrated time ephemeris, found requested size " +
                    std::to_string( startIndexAndSize.second ) );
    }

    std::map< Time, double > floatingPointValueNumericalSolution;
    for( const auto& solutionEntry : numericalSolution )
    {
        floatingPointValueNumericalSolution[ solutionEntry.first ] = solutionEntry.second( startIndexAndSize.first );
    }

    std::pair< std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >,
            std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > > > timeInterpolators =
            createRelativisticTimeInterpolators( floatingPointValueNumericalSolution );

    if( bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) == nullptr )
    {
        std::shared_ptr< TimeEphemeris > newTimeEphemeris = std::make_shared< TimeEphemerisWithFirstOrderDirectConversion >(
                    TimeEphemerisFromPostNewtonianExpansion::TimeDifferenceInterpolator( ),
                    TimeEphemerisFromPostNewtonianExpansion::TimeDifferenceInterpolator( ),
                    referencePointIdentifier.first,
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, Time >,
                               bodies.getBody( referencePointIdentifier.first ),
                               std::placeholders::_1 ) );
        bodies.getBody( referencePointIdentifier.first )->setTimeScaleConverter( newTimeEphemeris );
    }

    std::shared_ptr< TimeEphemerisFromPostNewtonianExpansion > timeEphemeris =
            std::dynamic_pointer_cast< TimeEphemerisFromPostNewtonianExpansion >(
                bodies.getBody( referencePointIdentifier.first )->getTimeScaleConverter( ) );
    if( timeEphemeris == nullptr )
    {
        throw std::runtime_error(
                    "Error when resetting integrated post-Newtonian time ephemeris, no TimeEphemerisFromPostNewtonianExpansion object found." );
    }

    if( referencePointIdentifier.second == "" )
    {
        timeEphemeris->resetBarycentricToBodycentricInterpolators( timeInterpolators.first, timeInterpolators.second );
    }
    else
    {
        if( timeEphemeris->doesReferencePointTopocentricConverterExist( referencePointIdentifier.second ) )
        {
            timeEphemeris->resetBodycentricToTopocentricInterpolators(
                        timeInterpolators.first, timeInterpolators.second, referencePointIdentifier.second );
        }
        else
        {
            std::shared_ptr< ground_stations::GroundStationState > groundStationState =
                    std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( referencePointIdentifier.first ) )->getGroundStation(
                        referencePointIdentifier.second )->getNominalStationState( );
            timeEphemeris->resetBodycentricToTopocentricInterpolators(
                        timeInterpolators.first,
                        timeInterpolators.second,
                        referencePointIdentifier.second,
                        std::bind( &ground_stations::GroundStationState::getNominalCartesianPosition, groundStationState ) );
        }
    }
}

template< >
void resetIntegratedPostNewtonianTimeEphemeris< double, long double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    throw std::runtime_error( "Cannot reset integrated time ephemeris for (double, long double)" );
}

template< >
void resetIntegratedPostNewtonianTimeEphemeris< Time, long double >(
        const simulation_setup::SystemOfBodies& bodies, const std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >& numericalSolution,
        const std::pair< std::string, std::string > referencePointIdentifier,
        const std::pair< int, int >& startIndexAndSize )
{
    throw std::runtime_error( "Cannot reset integrated time ephemeris for (Time, long double)" );
}

}  // namespace propagators

}  // namespace tudat
