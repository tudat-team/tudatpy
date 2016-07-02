/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Aerodynamics/standardAtmosphere.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace aerodynamics
{

//! Constructor, sets objects and functions from which relevant environment and state variables
//! are retrieved.
FlightConditions::FlightConditions(
        const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const boost::function< double( const Eigen::Vector3d ) > altitudeFunction,
        const boost::function< basic_mathematics::Vector6d( ) > stateOfVehicle,
        const boost::function< basic_mathematics::Vector6d( ) > stateOfCentralBody,
        const boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
        transformationToCentralBodyFrame,
        const boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
        aerodynamicAngleCalculator ):
    atmosphereModel_( atmosphereModel ),
    altitudeFunction_( altitudeFunction ),
    stateOfVehicle_( stateOfVehicle ),
    stateOfCentralBody_( stateOfCentralBody ),
    transformationToCentralBodyFrame_( transformationToCentralBodyFrame ),
    aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
    aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),currentAltitude_( TUDAT_NAN ),
    currentLatitude_( TUDAT_NAN ), currentLongitude_( TUDAT_NAN ), currentTime_( TUDAT_NAN )
{
    updateLatitudeAndLongitude_ = 0;

    if( boost::dynamic_pointer_cast< aerodynamics::StandardAtmosphere >( atmosphereModel_ ) ==
            NULL )
    {
        throw std::runtime_error( "Error when making flight conditions, no atmosphere is found" );
    }

    if( updateLatitudeAndLongitude_ && aerodynamicAngleCalculator_== NULL )
    {
        throw std::runtime_error( "Error when making flight conditions, angles are to be updated, but no calculator is set" );
    }
}

//! Function to set custom dependency of aerodynamic coefficients
void FlightConditions::setAerodynamicCoefficientsIndependentVariableFunction(
        const AerodynamicCoefficientsIndependentVariables independentVariable,
        const boost::function< double( ) > coefficientDependency )
{
    if( ( independentVariable == mach_number_dependent ) ||
            ( independentVariable == angle_of_attack_dependent ) ||
            ( independentVariable == angle_of_sideslip_dependent ) )
    {
        throw std::runtime_error(
                    std::string( "Error when setting aerodynamic coefficient function dependency, value of parameter " ) +
                    boost::lexical_cast< std::string >( independentVariable ) +
                    std::string(", will not  be used." ) );
    }
    else
    {
        customCoefficientDependencies_[ independentVariable ] = coefficientDependency;
    }
}

//! Function to update all flight conditions.
void FlightConditions::updateConditions( const double currentTime )
{
    currentTime_ = currentTime;

    // Calculate state of vehicle in global frame and corotating frame.
    currentBodyCenteredState_ = stateOfVehicle_( ) - stateOfCentralBody_( );
    currentBodyCenteredPseudoBodyFixedState_ = transformationToCentralBodyFrame_(
                currentBodyCenteredState_ );

    // Calculate altitute and airspeed of vehicle.
    currentAltitude_ =
            altitudeFunction_( currentBodyCenteredPseudoBodyFixedState_.segment( 0, 3 ) );
    currentAirspeed_ = currentBodyCenteredPseudoBodyFixedState_.segment( 3, 3 ).norm( );

    // Update aerodynamic/geometric angles.
    if( aerodynamicAngleCalculator_!= NULL )
    {
        aerodynamicAngleCalculator_->update( );
    }

    // Update latitude and longitude (if required)
    if( updateLatitudeAndLongitude_ )
    {
        currentLatitude_ = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::latitude_angle );
        currentLongitude_ = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::longitude_angle );
    }

    // Update density
    currentDensity_ = atmosphereModel_->getDensity( currentAltitude_, currentLongitude_,
                                                    currentLatitude_, currentTime_ );

    // Calculate independent variables for aerodynamic coefficients.
    std::vector< double > aerodynamicCoefficientIndependentVariables;
    for( unsigned int i = 0; i < aerodynamicCoefficientInterface_->
         getNumberOfIndependentVariables( ); i++ )
    {
        switch( aerodynamicCoefficientInterface_->getIndependentVariableName( i ) )
        {
        //Calculate Mach number if needed.
        case mach_number_dependent:
            aerodynamicCoefficientIndependentVariables.push_back(
                        aerodynamics::computeMachNumber(
                        currentAirspeed_, atmosphereModel_->getSpeedOfSound(
                            currentAltitude_, currentLongitude_, currentLatitude_, currentTime_ ) ) );
            break;
        //Get angle of attack if needed.
        case angle_of_attack_dependent:

            if( aerodynamicAngleCalculator_== NULL )
            {
                throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of attack" );
            }
            aerodynamicCoefficientIndependentVariables.push_back(
                        aerodynamicAngleCalculator_->getAerodynamicAngle(
                            reference_frames::angle_of_attack ) );
            break;
        //Get angle of sideslip if needed.
        case angle_of_sideslip_dependent:
            if( aerodynamicAngleCalculator_== NULL )
            {
                throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of sideslip" );
            }
            aerodynamicCoefficientIndependentVariables.push_back(
                        aerodynamicAngleCalculator_->getAerodynamicAngle(
                            reference_frames::angle_of_sideslip ) );
            break;
        default:
            if( customCoefficientDependencies_.count(
                        aerodynamicCoefficientInterface_->getIndependentVariableName( i ) ) == 0 )
            {
                throw std::runtime_error( "Error, did not recognize aerodynamic coefficient dependency "
                                          + boost::lexical_cast< std::string >(
                                              aerodynamicCoefficientInterface_->getIndependentVariableName( i ) ) );
            }
            else
            {
                aerodynamicCoefficientIndependentVariables.push_back(
                            customCoefficientDependencies_.at(
                                aerodynamicCoefficientInterface_->getIndependentVariableName( i ) )( ) );
            }
        }
    }

    // Update aerodynamic coefficients.
    aerodynamicCoefficientInterface_->updateCurrentCoefficients(
                aerodynamicCoefficientIndependentVariables );


}


}

}
