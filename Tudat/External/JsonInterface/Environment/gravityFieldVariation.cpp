/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "gravityFieldVariation.h"

#include "Tudat/External/JsonInterface/Mathematics/interpolator.h"

namespace tudat
{

namespace gravitation
{

//! Convert `BodyDeformationTypes` to `json`.
void to_json( json& jsonObject, const BodyDeformationTypes& bodyDeformationType )
{
    jsonObject = json_interface::stringFromEnum( bodyDeformationType, bodyDeformationTypes );
}

//! Convert `json` to `BodyDeformationTypes`.
void from_json( const json& jsonObject, BodyDeformationTypes& bodyDeformationType )
{
    bodyDeformationType = json_interface::enumFromString( jsonObject.get< std::string >( ), bodyDeformationTypes );
}

} // namespace gravitation


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldVariationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< GravityFieldVariationSettings >& variationSettings )
{
    if ( variationSettings )
    {
        using namespace json_interface;
        using K = Keys::Body::GravityFieldVariation;

        // Common parameters
        jsonObject[ K::bodyDeformationType ] = variationSettings->getBodyDeformationType( );
        jsonObject[ K::modelInterpolation ] = variationSettings->getInterpolatorSettings( );

        /// BasicSolidBodyGravityFieldVariationSettings
        boost::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > basicSolidBodySettings =
                boost::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( variationSettings );
        if ( basicSolidBodySettings )
        {
            jsonObject[ K::deformingBodies ] = basicSolidBodySettings->getDeformingBodies( );
            jsonObject[ K::loveNumbers ] = basicSolidBodySettings->getLoveNumbers( );
            jsonObject[ K::referenceRadius ] = basicSolidBodySettings->getBodyReferenceRadius( );
            return;
        }

        /// TabulatedGravityFieldVariationSettings
        boost::shared_ptr< TabulatedGravityFieldVariationSettings > tabulatedSettings =
                boost::dynamic_pointer_cast< TabulatedGravityFieldVariationSettings >( variationSettings );
        if ( tabulatedSettings )
        {
            jsonObject[ K::cosineCoefficientCorrections ] = tabulatedSettings->getCosineCoefficientCorrections( );
            jsonObject[ K::sineCoefficientCorrections ] = tabulatedSettings->getSineCoefficientCorrections( );
            jsonObject[ K::minimumDegree ] = tabulatedSettings->getMinimumDegree( );
            jsonObject[ K::minimumOrder ] = tabulatedSettings->getMinimumOrder( );
            return;
        }
    }
}

//! Create a shared pointer to a `GravityFieldVariationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< GravityFieldVariationSettings >& variationSettings )
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
        BasicSolidBodyGravityFieldVariationSettings defaults( { }, { }, TUDAT_NAN );
        variationSettings = boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                    getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ),
                    getValue< std::vector< std::vector< std::complex< double > > > >( jsonObject, K::loveNumbers ),
                    getNumeric< double >( jsonObject, K::referenceRadius ),
                    getValue< boost::shared_ptr< ModelInterpolationSettings > >( jsonObject, K::modelInterpolation ) );
        return;
    }
    case tabulated_variation:
        variationSettings = boost::make_shared< TabulatedGravityFieldVariationSettings >(
                    getValue< std::map< double, Eigen::MatrixXd > >( jsonObject, K::cosineCoefficientCorrections ),
                    getValue< std::map< double, Eigen::MatrixXd > >( jsonObject, K::sineCoefficientCorrections ),
                    getNumeric< int >( jsonObject, K::minimumDegree ),
                    getNumeric< int >( jsonObject, K::minimumOrder ),
                    getValue< boost::shared_ptr< ModelInterpolationSettings > >(
                        jsonObject, K::modelInterpolation )->interpolatorSettings_ );
        return;
    default:
        throw std::runtime_error( stringFromEnum( bodyDeformationType, bodyDeformationTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
