/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/aerodynamics/standardAtmosphere.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace aerodynamics
{


//! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
FlightConditions::FlightConditions( const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
                  const std::shared_ptr< reference_frames::AerodynamicAngleCalculator >
                  aerodynamicAngleCalculator ):
    shapeModel_( shapeModel ),
    aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),
    currentTime_( TUDAT_NAN )
{
    // Link body-state function.
    bodyCenteredPseudoBodyFixedStateFunction_ = std::bind(
                &reference_frames::AerodynamicAngleCalculator::getCurrentAirspeedBasedBodyFixedState, aerodynamicAngleCalculator_ );

    // Check if given body shape is an oblate spheroid and set geodetic latitude function if so
    if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ) != nullptr )
    {
        geodeticLatitudeFunction_ = std::bind(
                    &basic_astrodynamics::OblateSpheroidBodyShapeModel::getGeodeticLatitude,
                    std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ),
                    std::placeholders::_1, 1.0E-4 );
    }
    scalarFlightConditions_.resize( 12 );
    isScalarFlightConditionComputed_ = allScalarFlightConditionsUncomputed;
}

//! Function to update all flight conditions.
void FlightConditions::updateConditions( const double currentTime )
{

    if( !( currentTime == currentTime_ ) )
    {
        currentTime_ = currentTime;

        // Update aerodynamic angles (but not angles w.r.t. body-fixed frame).
        if( aerodynamicAngleCalculator_!= nullptr )
        {
            aerodynamicAngleCalculator_->update( currentTime, false );
        }

        // Calculate state of vehicle in global frame and corotating frame.
        currentBodyCenteredAirspeedBasedBodyFixedState_ = bodyCenteredPseudoBodyFixedStateFunction_( );

        // Update angles from aerodynamic to body-fixed frame (if relevant).
        if( aerodynamicAngleCalculator_!= nullptr )
        {
            aerodynamicAngleCalculator_->update( currentTime, true );
        }

    }
}


//! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
AtmosphericFlightConditions::AtmosphericFlightConditions(
        const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
        const std::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
        const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const std::function< double( const std::string& ) > controlSurfaceDeflectionFunction ):
    FlightConditions( shapeModel, aerodynamicAngleCalculator ),
    atmosphereModel_( atmosphereModel ),
    aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
    controlSurfaceDeflectionFunction_( controlSurfaceDeflectionFunction )
{

    if(  aerodynamicAngleCalculator_== nullptr )
    {
        throw std::runtime_error( "Error when making flight conditions, angles are to be updated, but no calculator is set" );
    }
}

//! Function to set custom dependency of aerodynamic coefficients
void AtmosphericFlightConditions::setAerodynamicCoefficientsIndependentVariableFunction(
        const AerodynamicCoefficientsIndependentVariables independentVariable,
        const std::function< double( ) > coefficientDependency )
{
    if( ( independentVariable == mach_number_dependent ) ||
            ( independentVariable == angle_of_attack_dependent ) ||
            ( independentVariable == angle_of_sideslip_dependent )||
            ( independentVariable == altitude_dependent ) ||
            ( independentVariable == time_dependent ) )
    {
        throw std::runtime_error(
                    std::string( "Error when setting aerodynamic coefficient function dependency, value of parameter " ) +
                    std::to_string( independentVariable ) +
                    std::string(", will not  be used." ) );
    }
    else
    {
        customCoefficientDependencies_[ independentVariable ] = coefficientDependency;
    }
}

//! Function to update all flight conditions.
void AtmosphericFlightConditions::updateConditions( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        currentTime_ = currentTime;

        // Update aerodynamic angles (but not angles w.r.t. body-fixed frame).
        if( aerodynamicAngleCalculator_!= nullptr )
        {
            aerodynamicAngleCalculator_->update( currentTime, false );
        }

        // Calculate state of vehicle in global frame and corotating frame.
        currentBodyCenteredAirspeedBasedBodyFixedState_ = bodyCenteredPseudoBodyFixedStateFunction_( );

        updateAerodynamicCoefficientInput( );

        // Update angles from aerodynamic to body-fixed frame (if relevant).
        if( aerodynamicAngleCalculator_!= nullptr )
        {
            aerodynamicAngleCalculator_->update( currentTime, true );
            updateAerodynamicCoefficientInput( );
        }

        // Update aerodynamic coefficients.
        if( aerodynamicCoefficientInterface_ != nullptr )
        {
            aerodynamicCoefficientInterface_->updateFullCurrentCoefficients(
                        aerodynamicCoefficientIndependentVariables_, controlSurfaceAerodynamicCoefficientIndependentVariables_,
                        currentTime_ );
        }
    }
}


void AtmosphericFlightConditions::updateAtmosphereInput( )
{
    if( ( isScalarFlightConditionComputed_.at( latitude_flight_condition ) == 0 ||
          isScalarFlightConditionComputed_.at( longitude_flight_condition ) == 0 ) )
    {
        computeLatitudeAndLongitude( );
    }

    if( isScalarFlightConditionComputed_.at( altitude_flight_condition ) == 0 )
    {
        computeAltitude( );
    }
}


//! Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
double AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariable(
        const AerodynamicCoefficientsIndependentVariables independentVariableType,
        const std::string& secondaryIdentifier )
{
    double currentIndependentVariable;
    switch( independentVariableType )
    {
    //Calculate Mach number if needed.
    case mach_number_dependent:
        currentIndependentVariable = getCurrentMachNumber( );
        break;
        //Get angle of attack if needed.
    case angle_of_attack_dependent:
        if( aerodynamicAngleCalculator_== nullptr )
        {
            throw std::runtime_error( "Error, aerodynamic angle calculator is nullptr, but require angle of attack" );
        }
        currentIndependentVariable = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::angle_of_attack );
        break;
        //Get angle of sideslip if needed.
    case angle_of_sideslip_dependent:
        if( aerodynamicAngleCalculator_== nullptr )
        {
            throw std::runtime_error( "Error, aerodynamic angle calculator is nullptr, but require angle of sideslip" );
        }
        currentIndependentVariable = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::angle_of_sideslip );
        break;
    case altitude_dependent:
        currentIndependentVariable = getCurrentAltitude( );
        break;
    case time_dependent:
        currentIndependentVariable = currentTime_;
        break;
    case control_surface_deflection_dependent:
    {
        try
        {
            currentIndependentVariable = controlSurfaceDeflectionFunction_( secondaryIdentifier );
        }
        catch( std::runtime_error const& )
        {
            throw std::runtime_error( "Error, control surface " + secondaryIdentifier + "not recognized when updating coefficients" );
        }
        break;
    }
    default:
        if( customCoefficientDependencies_.count( independentVariableType ) == 0 )
        {
            throw std::runtime_error( "Error, did not recognize aerodynamic coefficient dependency "
                                      + std::to_string( independentVariableType ) );
        }
        else
        {
            currentIndependentVariable = customCoefficientDependencies_.at( independentVariableType )( );
        }
    }

    return currentIndependentVariable;
}

//! Function to update the independent variables of the aerodynamic coefficient interface
void AtmosphericFlightConditions::updateAerodynamicCoefficientInput( )
{

    aerodynamicCoefficientIndependentVariables_.clear( );
    if( aerodynamicCoefficientInterface_ == nullptr )
    {
        throw std::runtime_error( "Error when calcualting aerodynamic coefficient input, no coefficient interface is defined" );
    }
    // Calculate independent variables for aerodynamic coefficients.
    for( unsigned int i = 0; i < aerodynamicCoefficientInterface_->getNumberOfIndependentVariables( ); i++ )
    {
        aerodynamicCoefficientIndependentVariables_.push_back(
                    getAerodynamicCoefficientIndependentVariable(
                        aerodynamicCoefficientInterface_->getIndependentVariableName( i ) ) );
    }

    controlSurfaceAerodynamicCoefficientIndependentVariables_.clear( );
    for( unsigned int i = 0; i < aerodynamicCoefficientInterface_->getNumberOfControlSurfaces( ); i++ )
    {
        std::string currentControlSurface = aerodynamicCoefficientInterface_->getControlSurfaceName( i );
        for( unsigned int j = 0; j < aerodynamicCoefficientInterface_->getNumberOfControlSurfaceIndependentVariables( currentControlSurface ); j++ )
        {
            controlSurfaceAerodynamicCoefficientIndependentVariables_[ currentControlSurface ].push_back(
                        getAerodynamicCoefficientIndependentVariable(
                            aerodynamicCoefficientInterface_->getControlSurfaceIndependentVariableName(
                                currentControlSurface, j ), currentControlSurface ) );
        }
    }
}

} // namespace aerodynamics

} // namespace tudat
