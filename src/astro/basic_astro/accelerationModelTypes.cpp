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

#include "tudat/astro/basic_astro/accelerationModelTypes.h"

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
    case point_mass_gravity:
        accelerationName = "central gravity ";
        break;
    case aerodynamic:
        accelerationName = "aerodynamic ";
        break;
    case spherical_harmonic_gravity:
        accelerationName = "spherical harmonic gravity ";
        break;
    case mutual_spherical_harmonic_gravity:
        accelerationName = "mutual spherical harmonic gravity ";
        break;
    case polyhedron_gravity:
        accelerationName = "polyhedron gravity ";
        break;
    case ring_gravity:
        accelerationName = "ring gravity ";
        break;
    case third_body_point_mass_gravity:
        accelerationName = "third-body central gravity ";
        break;
    case third_body_spherical_harmonic_gravity:
        accelerationName = "third-body spherical harmonic gravity ";
        break;
    case third_body_mutual_spherical_harmonic_gravity:
        accelerationName = "third-body mutual spherical harmonic gravity ";
        break;
    case third_body_polyhedron_gravity:
        accelerationName = "third-body polyhedron gravity ";
        break;
    case third_body_ring_gravity:
        accelerationName = "third-body ring gravity ";
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
    case direct_tidal_dissipation_in_central_body_acceleration:
        accelerationName  = "direct tidal dissipation in central body ";
        break;
    case direct_tidal_dissipation_in_orbiting_body_acceleration:
        accelerationName  = "direct tidal dissipation in orbiting body ";
        break;
    case radiation_pressure:
        accelerationName  = "radiation pressure acceleration";
        break;
    case cannon_ball_radiation_pressure:
        accelerationName  = "cannonball radiation pressure acceleration (deprecated)";
        break;
    case momentum_wheel_desaturation_acceleration:
        accelerationName  = "momentum wheen desaturation acceleration ";
        break;
    case yarkovsky_acceleration:
        accelerationName = "yarkovsky acceleration";
        break;
    case custom_acceleration:
        accelerationName = "custom acceleration";
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
        const std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
        accelerationModel )
{
    using namespace tudat::aerodynamics;
    using namespace tudat::electromagnetism;
    using namespace tudat::gravitation;

    // Nominal type is undefined
    AvailableAcceleration accelerationType = undefined_acceleration;

    // Check for each accelerarion mdoel type implemented as AvailableAcceleration.
    if( std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                accelerationModel ) != nullptr )
    {
        accelerationType = point_mass_gravity;
    }
    else if( std::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = third_body_point_mass_gravity;
    }
    else if( std::dynamic_pointer_cast< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = third_body_spherical_harmonic_gravity;
    }
    else if( std::dynamic_pointer_cast< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = third_body_mutual_spherical_harmonic_gravity;
    }
    else if( std::dynamic_pointer_cast< ThirdBodyPolyhedronGravitationalAccelerationModel >( accelerationModel ) != nullptr )
    {
        accelerationType = third_body_polyhedron_gravity;
    }
    else if( std::dynamic_pointer_cast< ThirdBodyRingGravitationalAccelerationModel >( accelerationModel ) != nullptr )
    {
        accelerationType = third_body_ring_gravity;
    }
    else if( std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                 accelerationModel ) != nullptr  )
    {
        accelerationType = spherical_harmonic_gravity;
    }
    else if( std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel ) != nullptr )
    {
        accelerationType = mutual_spherical_harmonic_gravity;
    }
    else if( std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >( accelerationModel ) != nullptr  )
    {
        accelerationType = polyhedron_gravity;
    }
    else if( std::dynamic_pointer_cast< RingGravitationalAccelerationModel >( accelerationModel ) != nullptr  )
    {
        accelerationType = ring_gravity;
    }
    else if( std::dynamic_pointer_cast< AerodynamicAcceleration >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = aerodynamic;
    }
    else if( std::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = thrust_acceleration;
    }
    else if( std::dynamic_pointer_cast< propulsion::MomentumWheelDesaturationThrustAcceleration >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = momentum_wheel_desaturation_acceleration;
    }
    else if( std::dynamic_pointer_cast< relativity::RelativisticAccelerationCorrection >(
                 accelerationModel ) != nullptr )
    {
        accelerationType = relativistic_correction_acceleration;
    }
    else if( std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( accelerationModel ) != nullptr )
    {
        accelerationType = empirical_acceleration;
    }
    else if( std::dynamic_pointer_cast<  gravitation::DirectTidalDissipationAcceleration >( accelerationModel ) != nullptr )
    {
        std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > dissipationAcceleration =
             std::dynamic_pointer_cast<  gravitation::DirectTidalDissipationAcceleration >( accelerationModel );
        if( dissipationAcceleration->getModelTideOnPlanet( ) )
        {
            accelerationType = direct_tidal_dissipation_in_central_body_acceleration;
        }
        else
        {
            accelerationType = direct_tidal_dissipation_in_orbiting_body_acceleration;
        }
    }
    else if( std::dynamic_pointer_cast< RadiationPressureAcceleration >( accelerationModel ) != NULL )
    {
        accelerationType = radiation_pressure;
    }
    else if (std::dynamic_pointer_cast<YarkovskyAcceleration>(accelerationModel) != nullptr)
    {
        accelerationType = yarkovsky_acceleration;
    }
    else if( std::dynamic_pointer_cast< CustomAccelerationModel >( accelerationModel ) != nullptr )
    {
        accelerationType = custom_acceleration;
    }
    else if( std::dynamic_pointer_cast< relativity::EinsteinInfeldHoffmannAcceleration >( accelerationModel ) != nullptr )
    {
        accelerationType = einstein_infeld_hoffmann_acceleration;
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
        const std::shared_ptr< MassRateModel > massRateModel )
{
    // Nominal type is undefined
    AvailableMassRateModels massRateType = undefined_mass_rate_model;

    // Check for each mass rate mdoel type implemented as AvailableMassRateModels.
    if( std::dynamic_pointer_cast< basic_astrodynamics::CustomMassRateModel >(
                massRateModel ) != nullptr )
    {
        massRateType = custom_mass_rate_model;
    }
    else if( std::dynamic_pointer_cast< propulsion::FromThrustMassRateModel >( massRateModel ) != nullptr )
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
std::vector< std::shared_ptr< AccelerationModel3d > > getAccelerationModelsOfType(
        const std::vector< std::shared_ptr< AccelerationModel3d > >& fullList,
        const AvailableAcceleration modelType )
{
    std::vector< std::shared_ptr< AccelerationModel3d > > accelerationList;
    for( unsigned int i = 0; i < fullList.size( ); i++ )
    {
        if( getAccelerationModelType( fullList.at( i ) ) == modelType )
        {
            accelerationList.push_back( fullList.at( i  ) );
        }
        else if( std::dynamic_pointer_cast< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration >( fullList.at( i ) ) != nullptr &&
            modelType == cannon_ball_radiation_pressure )
        {
            std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration  > radiationPressureAcceleration =
                std::dynamic_pointer_cast< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration >( fullList.at( i ) );
            if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( radiationPressureAcceleration->getTargetModel( ) ) != nullptr )
            {
                accelerationList.push_back( fullList.at( i ) );
            }
        }
    }
    return accelerationList;
}


//! Function to check whether an acceleration type is a direct gravitational acceleration
bool isAccelerationDirectGravitational( const AvailableAcceleration accelerationType )
{
    bool accelerationIsDirectGravity = 0;
    if( accelerationType == point_mass_gravity || accelerationType == spherical_harmonic_gravity ||
        accelerationType == mutual_spherical_harmonic_gravity || accelerationType == polyhedron_gravity ||
        accelerationType == ring_gravity )
    {
        accelerationIsDirectGravity = 1;
    }

    return accelerationIsDirectGravity;
}

//! Function to check whether an acceleration type is a third-body gravitational acceleration
bool isAccelerationFromThirdBody( const AvailableAcceleration accelerationType )
{
    bool accelerationIsFromThirdBody = false;
    if( accelerationType == third_body_point_mass_gravity || accelerationType == third_body_spherical_harmonic_gravity ||
        accelerationType == third_body_mutual_spherical_harmonic_gravity ||
        accelerationType == third_body_polyhedron_gravity || accelerationType == third_body_ring_gravity )
    {
        accelerationIsFromThirdBody = true;
    }

    return accelerationIsFromThirdBody;
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
    else if( accelerationType == point_mass_gravity )
    {
        thirdBodyAccelerationType = third_body_point_mass_gravity;
    }
    else if( accelerationType == spherical_harmonic_gravity )
    {
        thirdBodyAccelerationType = third_body_spherical_harmonic_gravity;
    }
    else if( accelerationType == mutual_spherical_harmonic_gravity )
    {
        thirdBodyAccelerationType = third_body_mutual_spherical_harmonic_gravity;
    }
    else if( accelerationType == polyhedron_gravity )
    {
        thirdBodyAccelerationType = third_body_polyhedron_gravity;
    }
    else if( accelerationType == ring_gravity )
    {
        thirdBodyAccelerationType = third_body_ring_gravity;
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
