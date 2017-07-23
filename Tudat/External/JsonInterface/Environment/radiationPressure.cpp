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

#include "radiationPressure.h"

namespace tudat
{

namespace simulation_setup
{

//! Convert `RadiationPressureType` to `json`.
void to_json( json& jsonObject, const RadiationPressureType& radiationPressureType )
{
    jsonObject = json_interface::stringFromEnum( radiationPressureType, radiationPressureTypes );
}

//! Convert `json` to `RadiationPressureType`.
void from_json( const json& jsonObject, RadiationPressureType& radiationPressureType )
{
    radiationPressureType = json_interface::enumFromString( jsonObject.get< std::string >( ), radiationPressureTypes );
}

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings )
{
    if ( radiationPressureInterfaceSettings )
    {
        using namespace json_interface;
        using Keys = Keys::Body::RadiationPressure;

        // Get type
        jsonObject[ Keys::type ] = radiationPressureInterfaceSettings->getRadiationPressureType( );

        /// CannonBallRadiationPressureInterfaceSettings
        boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallRadiationPressureInterfaceSettings =
                boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        if ( cannonBallRadiationPressureInterfaceSettings )
        {
            jsonObject[ Keys::referenceArea ] = cannonBallRadiationPressureInterfaceSettings->getArea( );
            jsonObject[ Keys::radiationPressureCoefficient ] =
                    cannonBallRadiationPressureInterfaceSettings->getRadiationPressureCoefficient( );
            jsonObject[ Keys::ocultingBodies ] =
                    cannonBallRadiationPressureInterfaceSettings->getOccultingBodies( );
            return;
        }
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `RadiationPressureInterfaceSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::RadiationPressureInterfaceSettings > createRadiationPressureInterfaceSettings(
        const json& settings, const std::string& sourceBodyName, const KeyTree& keyTree, const double& fallbackArea )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::RadiationPressure;

    // Get radiation pressure coefficient type (cannonBall by default)
    const RadiationPressureType radiationPressureType =
            getValue( settings, keyTree + sourceBodyName + Keys::type, cannon_ball );

    switch ( radiationPressureType )
    {
    case cannon_ball:
    {
        CannonBallRadiationPressureInterfaceSettings defaults( "", TUDAT_NAN, TUDAT_NAN );
        return boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    sourceBodyName,
                    getNumeric( settings, keyTree + sourceBodyName + Keys::referenceArea, fallbackArea ),
                    getValue< double >( settings, keyTree + sourceBodyName + Keys::radiationPressureCoefficient ),
                    getValue( settings, keyTree + sourceBodyName + Keys::ocultingBodies,
                              defaults.getOccultingBodies( ) ) );
    }
    default:
        throw std::runtime_error( stringFromEnum( radiationPressureType, radiationPressureTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
