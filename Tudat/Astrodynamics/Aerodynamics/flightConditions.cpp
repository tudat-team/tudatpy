/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

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
        const boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
        aerodynamicAngleCalculator ):
    atmosphereModel_( atmosphereModel ),
    altitudeFunction_( altitudeFunction ),
    aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
    aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),currentAltitude_( TUDAT_NAN ),
    currentLatitude_( TUDAT_NAN ), currentLongitude_( TUDAT_NAN ), currentTime_( TUDAT_NAN )
{
    updateLatitudeAndLongitude_ = 0;

    bodyCenteredPseudoBodyFixedStateFunction_ = boost::bind(
                &reference_frames::AerodynamicAngleCalculator::getCurrentBodyFixedState, aerodynamicAngleCalculator_ );

    if( boost::dynamic_pointer_cast< aerodynamics::StandardAtmosphere >( atmosphereModel_ ) == NULL )
    {
        updateLatitudeAndLongitude_ = 1;
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
    if( !( currentTime == currentTime_ ) )
    {
        currentTime_ = currentTime;

        // Update aerodynamic angles (but not angles w.r.t. body-fixed frame).
        if( aerodynamicAngleCalculator_!= NULL )
        {
            aerodynamicAngleCalculator_->update( currentTime, false );
        }

        // Calculate state of vehicle in global frame and corotating frame.
        currentBodyCenteredPseudoBodyFixedState_ = bodyCenteredPseudoBodyFixedStateFunction_( );

        // Calculate altitute and airspeed of vehicle.
        currentAltitude_ =
                altitudeFunction_( currentBodyCenteredPseudoBodyFixedState_.segment( 0, 3 ) );
        currentAirspeed_ = currentBodyCenteredPseudoBodyFixedState_.segment( 3, 3 ).norm( );



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

        updateAerodynamicCoefficientInput( );

        // Update angles from aerodynamic to body-fixed frame (if relevant).
        if( aerodynamicAngleCalculator_!= NULL )
        {
            aerodynamicAngleCalculator_->update( currentTime, true );
            updateAerodynamicCoefficientInput( );
        }

        // Update aerodynamic coefficients.
        aerodynamicCoefficientInterface_->updateCurrentCoefficients(
                    aerodynamicCoefficientIndependentVariables_ );
    }
}

void FlightConditions::updateAerodynamicCoefficientInput( )
{
    aerodynamicCoefficientIndependentVariables_.clear( );

    // Calculate independent variables for aerodynamic coefficients.
    for( unsigned int i = 0; i < aerodynamicCoefficientInterface_->
         getNumberOfIndependentVariables( ); i++ )
    {
        switch( aerodynamicCoefficientInterface_->getIndependentVariableName( i ) )
        {
        //Calculate Mach number if needed.
        case mach_number_dependent:
            aerodynamicCoefficientIndependentVariables_.push_back(
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
            aerodynamicCoefficientIndependentVariables_.push_back(
                        aerodynamicAngleCalculator_->getAerodynamicAngle(
                            reference_frames::angle_of_attack ) );
            break;
        //Get angle of sideslip if needed.
        case angle_of_sideslip_dependent:
            if( aerodynamicAngleCalculator_== NULL )
            {
                throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of sideslip" );
            }
            aerodynamicCoefficientIndependentVariables_.push_back(
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
                aerodynamicCoefficientIndependentVariables_.push_back(
                            customCoefficientDependencies_.at(
                                aerodynamicCoefficientInterface_->getIndependentVariableName( i ) )( ) );
            }
        }
    }
}

//! Function to set the angle of attack to trimmed conditions.
boost::shared_ptr< TrimOrientationCalculator > setTrimmedConditions(
        const boost::shared_ptr< FlightConditions > flightConditions )
{
    // Create trim object.
    boost::shared_ptr< TrimOrientationCalculator > trimOrientation =
            boost::make_shared< TrimOrientationCalculator >(
                flightConditions->getAerodynamicCoefficientInterface( ) );

    // Create angle-of-attack function from trim object.
    boost::function< std::vector< double >( ) > untrimmedIndependentVariablesFunction =
            boost::bind( &FlightConditions::getAerodynamicCoefficientIndependentVariables,
                         flightConditions );
    flightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::bind( &TrimOrientationCalculator::findTrimAngleOfAttackFromFunction, trimOrientation,
                             untrimmedIndependentVariablesFunction ) );

    return trimOrientation;
}

} // namespace aerodynamics

} // namespace tudat
