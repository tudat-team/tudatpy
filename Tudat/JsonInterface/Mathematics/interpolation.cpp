/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace interpolators
{

// InterpolatorSettings

//! Create a `json` object from a shared pointer to a `InterpolatorSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< InterpolatorSettings >& interpolatorSettings )
{
    if ( ! interpolatorSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Interpolation::Interpolator;

    const InterpolatorTypes interpolatorType = interpolatorSettings->getInterpolatorType( );
    jsonObject[ K::type ] = interpolatorType;
    jsonObject[ K::lookupScheme ] = interpolatorSettings->getSelectedLookupScheme( );
    jsonObject[ K::useLongDoubleTimeStep ] = interpolatorSettings->getUseLongDoubleTimeStep( );
    jsonObject[ K::boundaryHandling ] = interpolatorSettings->getBoundaryHandling( );

    switch ( interpolatorType )
    {
    case linear_interpolator:
    case cubic_spline_interpolator:
    case hermite_spline_interpolator:
    case piecewise_constant_interpolator:
        return;
    case lagrange_interpolator:
    {
        std::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                std::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        assertNonnullptrPointer( lagrangeInterpolatorSettings );
        jsonObject[ K::order ] = lagrangeInterpolatorSettings->getInterpolatorOrder( );
        jsonObject[ K::lagrangeBoundaryHandling ] = lagrangeInterpolatorSettings->getLagrangeBoundaryHandling( );
        return;
    }
    default:
        handleUnimplementedEnumValue( interpolatorType, oneDimensionalInterpolatorTypes,
                                      unsupportedOneDimensionalInterpolatorTypes );
    }
}

//! Create a shared pointer to a `InterpolatorSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< InterpolatorSettings >& interpolatorSettings )
{
    using namespace json_interface;
    using K = Keys::Interpolation::Interpolator;

    // Get interpolator type
    const InterpolatorTypes interpolatorType =
            getValue< InterpolatorTypes >( jsonObject, K::type );

    switch ( interpolatorType )
    {
    case linear_interpolator:
    case cubic_spline_interpolator:
    case hermite_spline_interpolator:
    case piecewise_constant_interpolator:
    {
        InterpolatorSettings defaults( linear_interpolator );
        interpolatorSettings = std::make_shared< InterpolatorSettings >(
                    interpolatorType,
                    getValue( jsonObject, K::lookupScheme, defaults.getSelectedLookupScheme( ) ),
                    getValue( jsonObject, K::useLongDoubleTimeStep, defaults.getUseLongDoubleTimeStep( ) ),
                    getValue( jsonObject, K::boundaryHandling, defaults.getBoundaryHandling( ) ) );
        return;
    }
    case lagrange_interpolator:
    {
        LagrangeInterpolatorSettings defaults( 0 );
        interpolatorSettings = std::make_shared< LagrangeInterpolatorSettings >(
                    getValue< double >( jsonObject, K::order ),
                    getValue( jsonObject, K::useLongDoubleTimeStep, defaults.getUseLongDoubleTimeStep( ) ),
                    getValue( jsonObject, K::lookupScheme, defaults.getSelectedLookupScheme( ) ),
                    getValue( jsonObject, K::lagrangeBoundaryHandling, defaults.getLagrangeBoundaryHandling( ) ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( interpolatorType, oneDimensionalInterpolatorTypes,
                                      unsupportedOneDimensionalInterpolatorTypes );
    }
}

} // namespace interpolators


namespace simulation_setup
{

// ModelInterpolationSettings

//! Create a `json` object from a shared pointer to a `ModelInterpolationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    if ( ! modelInterpolationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Interpolation::ModelInterpolation;

    jsonObject[ K::initialTime ] = modelInterpolationSettings->initialTime_;
    jsonObject[ K::finalTime ] = modelInterpolationSettings->finalTime_;
    jsonObject[ K::timeStep ] = modelInterpolationSettings->timeStep_;
    jsonObject[ K::interpolator ] = modelInterpolationSettings->interpolatorSettings_;
}

//! Create a shared pointer to a `ModelInterpolationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ModelInterpolationSettings >& modelInterpolationSettings )
{
    using namespace json_interface;
    using K = Keys::Interpolation::ModelInterpolation;

    ModelInterpolationSettings defaults;
    modelInterpolationSettings = std::make_shared< ModelInterpolationSettings >(
                getValue( jsonObject, K::initialTime, defaults.initialTime_ ),
                getValue( jsonObject, K::finalTime, defaults.finalTime_ ),
                getValue( jsonObject, K::timeStep, defaults.timeStep_ ),
                getValue( jsonObject, K::interpolator, defaults.interpolatorSettings_ ) );
}

} // namespace simulation_setup

} // namespace tudat
