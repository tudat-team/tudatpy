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

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to get a string representing a 'named identification' of an acceleration type
std::string getAccelerationModelName( const AvailableAcceleration accelerationType )
{
    std::string accelerationName;
    switch( accelerationType )
    {
    case central_gravity:
        accelerationName = "central gravity ";
        break;
    case aerodynamic:
        accelerationName = "aerodynamic ";
        break;
    case cannon_ball_radiation_pressure:
        accelerationName = "cannonball radiation pressure ";
        break;
    case spherical_harmonic_gravity:
        accelerationName = "spherical harmonic gravity ";
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationName = "mutual spherical harmonic gravity ";
        break;
    case third_body_central_gravity:
        accelerationName = "third-body central gravity ";
        break;
    case third_body_spherical_harmonic_gravity:
        accelerationName = "third-body spherical harmonic gravity ";
        break;
    case third_body_mutual_spherical_harmonic_gravity:
        accelerationName = "third-body mutual spherical harmonic gravity ";
        break;
    case thrust_acceleration:
        accelerationName = "thrust ";
        break;
    case relativistic_correction_acceleration:
        accelerationName  = "relativistic correction ";
        break;
    case empirical_acceleration:
        accelerationName  = "empirical correction ";
        break;
    case direct_tidal_dissipation_acceleration:
        accelerationName  = "direct tidal dissipation ";
        break;
    default:
        std::string errorMessage = "Error, acceleration type " +
                std::to_string( accelerationType ) +
                "not found when retrieving acceleration name ";
        throw std::runtime_error( errorMessage );
    }
    return accelerationName;
}

//! Function to identify the derived class type of an acceleration model.
AvailableAcceleration getAccelerationModelType(
        const boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
        accelerationModel )
{
    using namespace tudat::aerodynamics;
    using namespace tudat::electro_magnetism;
    using namespace tudat::gravitation;

    // Nominal type is undefined
    AvailableAcceleration accelerationType = undefined_acceleration;

    // Check for each accelerarion mdoel type implemented as AvailableAcceleration.
    if( boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                accelerationModel ) != NULL )
    {
        accelerationType = central_gravity;
    }
    else if( boost::dynamic_pointer_cast< CannonBallRadiationPressureAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = cannon_ball_radiation_pressure;
    }
    else if( boost::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = third_body_central_gravity;
    }
    else if( boost::dynamic_pointer_cast< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != NULL )
    {
        accelerationType = third_body_spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != NULL )
    {
        accelerationType = third_body_mutual_spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != NULL  )
    {
        accelerationType = spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel ) != NULL )
    {
        accelerationType = mutual_spherical_harmonic_gravity;
    }
    else if( boost::dynamic_pointer_cast< AerodynamicAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = aerodynamic;
    }
    else if( boost::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                 accelerationModel ) != NULL )
    {
        accelerationType = thrust_acceleration;
    }
    else if( boost::dynamic_pointer_cast< relativity::RelativisticAccelerationCorrection >(
                 accelerationModel ) != NULL )
    {
        accelerationType = relativistic_correction_acceleration;
    }
    else if( boost::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( accelerationModel ) != NULL )
    {
        accelerationType = empirical_acceleration;
    }
    else if( boost::dynamic_pointer_cast<  gravitation::DirectTidalDissipationAcceleration >( accelerationModel ) != NULL )
    {
        accelerationType = direct_tidal_dissipation_acceleration;
    }
    else
    {
        throw std::runtime_error(
                    "Error, acceleration model not identified when getting acceleration type." );
    }

    // Return identified type.
    return accelerationType;

}

//! Function to identify the type of a mass rate model.
AvailableMassRateModels getMassRateModelType(
        const boost::shared_ptr< MassRateModel > massRateModel )
{
    // Nominal type is undefined
    AvailableMassRateModels massRateType = undefined_mass_rate_model;

    // Check for each mass rate mdoel type implemented as AvailableMassRateModels.
    if( boost::dynamic_pointer_cast< basic_astrodynamics::CustomMassRateModel >(
                massRateModel ) != NULL )
    {
        massRateType = custom_mass_rate_model;
    }
    else if( boost::dynamic_pointer_cast< propulsion::FromThrustMassRateModel >( massRateModel ) != NULL )
    {
        massRateType = from_thrust_mass_rate_model;
    }
    else
    {
        throw std::runtime_error(
                    "Error, mass rate model not identified when getting mass rate model type." );
    }
    return massRateType;
}

//! Function to get all acceleration models of a given type from a list of models
std::vector< boost::shared_ptr< AccelerationModel3d > > getAccelerationModelsOfType(
        const std::vector< boost::shared_ptr< AccelerationModel3d > >& fullList,
        const AvailableAcceleration modelType )
{
    std::vector< boost::shared_ptr< AccelerationModel3d > > accelerationList;
    for( unsigned int i = 0; i < fullList.size( ); i++ )
    {
        if( getAccelerationModelType( fullList.at( i ) ) == modelType )
        {
            accelerationList.push_back( fullList.at( i  ) );
        }
    }
    return accelerationList;
}

//! Function to check whether an acceleration type is a direct gravitational acceleration
bool isAccelerationDirectGravitational( const AvailableAcceleration accelerationType )
{
    bool accelerationIsDirectGravity = 0;
    if( ( accelerationType == central_gravity ) ||
            ( accelerationType == spherical_harmonic_gravity ) ||
            ( accelerationType == mutual_spherical_harmonic_gravity ) )
    {
        accelerationIsDirectGravity = 1;
    }

    return accelerationIsDirectGravity;
}

//! Function to get the third-body counterpart of a direct gravitational acceleration type
AvailableAcceleration getAssociatedThirdBodyAcceleration( const AvailableAcceleration accelerationType )
{

    AvailableAcceleration thirdBodyAccelerationType;
    if( !isAccelerationDirectGravitational( accelerationType ) )
    {
        std::string errorMessage = "Error when getting third-body gravity type, requested type: " +
                std::to_string( accelerationType ) + " is not a direct gravity acceleration";
        throw std::runtime_error( errorMessage );
    }
    else if( accelerationType == central_gravity )
    {
        thirdBodyAccelerationType = third_body_central_gravity;
    }
    else if( accelerationType == spherical_harmonic_gravity )
    {
        thirdBodyAccelerationType = third_body_spherical_harmonic_gravity;
    }
    else if( accelerationType == mutual_spherical_harmonic_gravity )
    {
        thirdBodyAccelerationType = third_body_mutual_spherical_harmonic_gravity;
    }
    else
    {        std::string errorMessage = "Error when getting thirdbody gravity type, requested type: " +
                std::to_string( accelerationType ) + " is not recognized.";
        throw std::runtime_error( errorMessage );

    }
    return thirdBodyAccelerationType;
}

} // namespace basic_astrodynamics

} // namespace tudat
