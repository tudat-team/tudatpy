/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/interface/json/environment/gravityFieldVariation.h"

#include "tudat/interface/json/math/interpolation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldVariationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GravityFieldVariationSettings >& variationSettings )
{
    if ( ! variationSettings )
    {
        return;
    }
    using namespace gravitation;
    using namespace json_interface;
    using K = Keys::Body::GravityFieldVariation;

    const BodyDeformationTypes bodyDeformationType = variationSettings->getBodyDeformationType( );
    jsonObject[ K::bodyDeformationType ] = bodyDeformationType;
    jsonObject[ K::modelInterpolation ] = variationSettings->getInterpolatorSettings( );

    switch ( bodyDeformationType )
    {
    case basic_solid_body:
    {
        std::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > basicSolidBodySettings =
                std::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( variationSettings );
        assertNonnullptrPointer( basicSolidBodySettings );
        jsonObject[ K::deformingBodies ] = basicSolidBodySettings->getDeformingBodies( );
        jsonObject[ K::loveNumbers ] = basicSolidBodySettings->getLoveNumbers( );
        return;
    }
    case tabulated_variation:
    {
        std::shared_ptr< TabulatedGravityFieldVariationSettings > tabulatedSettings =
                std::dynamic_pointer_cast< TabulatedGravityFieldVariationSettings >( variationSettings );
        assertNonnullptrPointer( tabulatedSettings );
        jsonObject[ K::cosineCoefficientCorrections ] = tabulatedSettings->getCosineCoefficientCorrections( );
        jsonObject[ K::sineCoefficientCorrections ] = tabulatedSettings->getSineCoefficientCorrections( );
        jsonObject[ K::minimumDegree ] = tabulatedSettings->getMinimumDegree( );
        jsonObject[ K::minimumOrder ] = tabulatedSettings->getMinimumOrder( );
        return;
    }
    default:
        handleUnimplementedEnumValue( bodyDeformationType, bodyDeformationTypes, unsupportedBodyDeformationTypes );
    }
}

//! Create a shared pointer to a `GravityFieldVariationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GravityFieldVariationSettings >& variationSettings )
{
    using namespace gravitation;
    using namespace interpolators;
    using namespace json_interface;
    using K = Keys::Body::GravityFieldVariation;

    // Body deformation type
    const BodyDeformationTypes bodyDeformationType =
            getValue< BodyDeformationTypes >( jsonObject, K::bodyDeformationType );

    switch ( bodyDeformationType ) {
    case basic_solid_body:
    {
        throw std::runtime_error( "Error, BasicSolidBodyGravityFieldVariationSettings not yet enabled in JSON interface." );
//        BasicSolidBodyGravityFieldVariationSettings defaults( { }, { } );
//        variationSettings = std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
//                    getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ),
//                    getValue< std::vector< std::vector< std::complex< double > > > >( jsonObject, K::loveNumbers ),
//                    getValue( jsonObject, K::modelInterpolation, defaults.getInterpolatorSettings( ) ) );
        return;
    }
    case tabulated_variation:
        variationSettings = std::make_shared< TabulatedGravityFieldVariationSettings >(
                    getValue< std::map< double, Eigen::MatrixXd > >( jsonObject, K::cosineCoefficientCorrections ),
                    getValue< std::map< double, Eigen::MatrixXd > >( jsonObject, K::sineCoefficientCorrections ),
                    getValue< int >( jsonObject, K::minimumDegree ),
                    getValue< int >( jsonObject, K::minimumOrder ),
                    getValue< std::shared_ptr< InterpolatorSettings > >( jsonObject, K::interpolator ) );
        return;
    default:
        handleUnimplementedEnumValue( bodyDeformationType, bodyDeformationTypes, unsupportedBodyDeformationTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
