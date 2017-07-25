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
        using Keys = Keys::Body::GravityFieldVariation;

        // Common parameters
        jsonObject[ Keys::bodyDeformationType ] = variationSettings->getBodyDeformationType( );
        jsonObject[ Keys::modelInterpolation ] = variationSettings->getInterpolatorSettings( );

        /// BasicSolidBodyGravityFieldVariationSettings
        boost::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > basicSolidBodySettings =
                boost::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( variationSettings );
        if ( basicSolidBodySettings )
        {
            jsonObject[ Keys::deformingBodies ] = basicSolidBodySettings->getDeformingBodies( );
            jsonObject[ Keys::loveNumbers ] = basicSolidBodySettings->getLoveNumbers( );
            jsonObject[ Keys::referenceRadius ] = basicSolidBodySettings->getBodyReferenceRadius( );
            return;
        }

        /// TabulatedGravityFieldVariationSettings
        boost::shared_ptr< TabulatedGravityFieldVariationSettings > tabulatedSettings =
                boost::dynamic_pointer_cast< TabulatedGravityFieldVariationSettings >( variationSettings );
        if ( tabulatedSettings )
        {
            jsonObject[ Keys::cosineCoefficientCorrections ] = tabulatedSettings->getCosineCoefficientCorrections( );
            jsonObject[ Keys::sineCoefficientCorrections ] = tabulatedSettings->getSineCoefficientCorrections( );
            jsonObject[ Keys::minimumDegree ] = tabulatedSettings->getMinimumDegree( );
            jsonObject[ Keys::minimumOrder ] = tabulatedSettings->getMinimumOrder( );
            return;
        }
    }
}

/*
//! Convert `json` to `InterpolatorSettings` shared pointer.
void from_json( const json& jsonObject, boost::shared_ptr< GravityFieldVariationSettings >& variationSettings )
{
    variationSettings = json_interface::createGravityFieldVariationSettings( jsonObject );
}
*/

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `GravityFieldVariationSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::GravityFieldVariationSettings > createGravityFieldVariationSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace gravitation;
    using namespace interpolators;
    using namespace simulation_setup;
    using Keys = Keys::Body::GravityFieldVariation;

    // Body deformation type
    const BodyDeformationTypes bodyDeformationType =
            getValue< BodyDeformationTypes >( settings, keyTree + Keys::bodyDeformationType );

    switch ( bodyDeformationType ) {
    case basic_solid_body:
    {
        BasicSolidBodyGravityFieldVariationSettings defaults( { }, { }, TUDAT_NAN );
        return boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                    getValue< std::vector< std::string > >( settings, keyTree + Keys::deformingBodies ),
                    getValue< std::vector< std::vector< std::complex< double > > > >(
                        settings, keyTree + Keys::loveNumbers ),
                    getNumeric< double >( settings, keyTree + Keys::referenceRadius ),
                    createModelInterpolationSettings(
                        settings, keyTree + Keys::modelInterpolation, defaults.getInterpolatorSettings( ) ) );
    }
    case tabulated_variation:
        return boost::make_shared< TabulatedGravityFieldVariationSettings >(
                    getValue< std::map< double, Eigen::MatrixXd > >(
                        settings, keyTree + Keys::cosineCoefficientCorrections ),
                    getValue< std::map< double, Eigen::MatrixXd > >(
                        settings, keyTree + Keys::sineCoefficientCorrections ),
                    getNumeric< int >( settings, keyTree + Keys::minimumDegree ),
                    getNumeric< int >( settings, keyTree + Keys::minimumOrder ),
                    createModelInterpolationSettings(
                        settings, keyTree + Keys::modelInterpolation )->interpolatorSettings_ );
    default:
        throw std::runtime_error( stringFromEnum( bodyDeformationType, bodyDeformationTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
