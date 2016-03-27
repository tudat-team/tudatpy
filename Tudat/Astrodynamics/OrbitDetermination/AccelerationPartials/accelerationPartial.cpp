
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{


//AccelerationDependenceSettings::AccelerationDependenceSettings( int flag )
//{
//    if( flag == 0 )
//    {
//        accelerationDependence[ central_gravity ][ gravitational_parameter ] = true;
//        accelerationDependence[ central_gravity ][ constant_rotation_rate ] = false;
//        accelerationDependence[ central_gravity ][ rotation_pole_position ] = false;
//        accelerationDependence[ central_gravity ][ constant_drag_coefficient ] = false;
//        accelerationDependence[ central_gravity ][ exponential_atmosphere_scale_height ] = false;
//        accelerationDependence[ central_gravity ][ radiation_pressure_coefficient ] = false;
//        accelerationDependence[ central_gravity ][ spherical_harmonics_cosine_coefficient_block ] = true;
//        accelerationDependence[ central_gravity ][ spherical_harmonics_sine_coefficient_block ] = true;
//        accelerationDependence[ central_gravity ][ ground_station_position ] = false;

//        accelerationDependence[ constant_drag_aerodynamic ][ gravitational_parameter ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ constant_rotation_rate ] = true;
//        accelerationDependence[ constant_drag_aerodynamic ][ rotation_pole_position ] = true;
//        accelerationDependence[ constant_drag_aerodynamic ][ constant_drag_coefficient ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ exponential_atmosphere_scale_height ] = true;
//        accelerationDependence[ constant_drag_aerodynamic ][ radiation_pressure_coefficient ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ spherical_harmonics_cosine_coefficient_block ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ spherical_harmonics_sine_coefficient_block ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ ground_station_position ] = false;

//        accelerationDependence[ cannon_ball_radiation_pressure ][ gravitational_parameter ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ constant_rotation_rate ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ rotation_pole_position ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ constant_drag_coefficient ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ exponential_atmosphere_scale_height ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ radiation_pressure_coefficient ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ spherical_harmonics_cosine_coefficient_block ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ spherical_harmonics_sine_coefficient_block ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ ground_station_position ] = false;

//    }
//    else if( flag == 1 )
//    {
//        accelerationDependence[ central_gravity ][ gravitational_parameter ] = false;
//        accelerationDependence[ central_gravity ][ constant_rotation_rate ] = false;
//        accelerationDependence[ central_gravity ][ rotation_pole_position ] = false;
//        accelerationDependence[ central_gravity ][ constant_drag_coefficient ] = false;
//        accelerationDependence[ central_gravity ][ exponential_atmosphere_scale_height ] = false;
//        accelerationDependence[ central_gravity ][ radiation_pressure_coefficient ] = false;
//        accelerationDependence[ central_gravity ][ spherical_harmonics_cosine_coefficient_block ] = false;
//        accelerationDependence[ central_gravity ][ spherical_harmonics_sine_coefficient_block ] = false;
//        accelerationDependence[ central_gravity ][ ground_station_position ] = false;

//        accelerationDependence[ constant_drag_aerodynamic ][ gravitational_parameter ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ constant_rotation_rate ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ rotation_pole_position ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ constant_drag_coefficient ] = true;
//        accelerationDependence[ constant_drag_aerodynamic ][ exponential_atmosphere_scale_height ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ radiation_pressure_coefficient ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ spherical_harmonics_cosine_coefficient_block ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ spherical_harmonics_sine_coefficient_block ] = false;
//        accelerationDependence[ constant_drag_aerodynamic ][ ground_station_position ] = false;

//        accelerationDependence[ cannon_ball_radiation_pressure ][ gravitational_parameter ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ constant_rotation_rate ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ rotation_pole_position ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ constant_drag_coefficient ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ exponential_atmosphere_scale_height ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ radiation_pressure_coefficient ] = true;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ spherical_harmonics_cosine_coefficient_block ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ spherical_harmonics_sine_coefficient_block ] = false;
//        accelerationDependence[ cannon_ball_radiation_pressure ][ ground_station_position ] = false;
//    }
//}




//bool isAccelerationDependentOnParameter( const EstimatebleParameterIdentifier parameterIdentifier,
//                                         const basic_astrodynamics::AvailableAcceleration accelerationType,
//                                         const std::string& acceleratingBody, const std::string& acceleratedBody )
//{
//    bool isDependencyTrue = false;

//    if( parameterIdentifier.second.first == acceleratingBody )
//    {
//        isDependencyTrue = defaultAccelerationDependenceExertingBody.accelerationDependence
//                [ accelerationType ][ parameterIdentifier.first ];
//    }
//    else if( parameterIdentifier.second.first == acceleratedBody )
//    {
//        isDependencyTrue = defaultAccelerationDependenceUndergoingBody.accelerationDependence
//                [ accelerationType ][ parameterIdentifier.first ];
//    }

//    return isDependencyTrue;
//}


}

}

}
