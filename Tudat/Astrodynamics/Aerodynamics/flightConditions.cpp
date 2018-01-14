/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace aerodynamics
{

//! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
FlightConditions::FlightConditions(
        const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
        const boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const boost::function< double( const std::string& ) > controlSurfaceDeflectionFunction ):
    atmosphereModel_( atmosphereModel ),
    shapeModel_( shapeModel ),
    aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
    aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),
    controlSurfaceDeflectionFunction_( controlSurfaceDeflectionFunction ),
    currentTime_( TUDAT_NAN )
{
    // Check if given body shape is an oblate spheroid and set geodetic latitude function if so
    if( boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ) != NULL )
    {
        geodeticLatitudeFunction_ = boost::bind(
                    &basic_astrodynamics::OblateSpheroidBodyShapeModel::getGeodeticLatitude,
                    boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ),
                    _1, 1.0E-4 );
    }

    // Link body-state function.
    bodyCenteredPseudoBodyFixedStateFunction_ = boost::bind(
                &reference_frames::AerodynamicAngleCalculator::getCurrentAirspeedBasedBodyFixedState, aerodynamicAngleCalculator_ );

    // Check if atmosphere requires latitude and longitude update.
    if( boost::dynamic_pointer_cast< aerodynamics::StandardAtmosphere >( atmosphereModel_ ) == NULL )
    {
        updateLatitudeAndLongitudeForAtmosphere_ = 1;
    }
    else
    {
        updateLatitudeAndLongitudeForAtmosphere_ = 0;
    }
    isLatitudeAndLongitudeSet_ = 0;


    if( updateLatitudeAndLongitudeForAtmosphere_ && aerodynamicAngleCalculator_== NULL )
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
            ( independentVariable == angle_of_sideslip_dependent )||
            ( independentVariable == altitude_dependent ) )
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
        currentBodyCenteredAirspeedBasedBodyFixedState_ = bodyCenteredPseudoBodyFixedStateFunction_( );

        updateAerodynamicCoefficientInput( );

        // Update angles from aerodynamic to body-fixed frame (if relevant).
        if( aerodynamicAngleCalculator_!= NULL )
        {
            aerodynamicAngleCalculator_->update( currentTime, true );
            updateAerodynamicCoefficientInput( );
        }

        // Update aerodynamic coefficients.
        aerodynamicCoefficientInterface_->updateFullCurrentCoefficients(
                    aerodynamicCoefficientIndependentVariables_, controlSurfaceAerodynamicCoefficientIndependentVariables_ );
    }
}

//! Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
double FlightConditions::getAerodynamicCoefficientIndependentVariable(
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

        if( aerodynamicAngleCalculator_== NULL )
        {
            throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of attack" );
        }
        currentIndependentVariable = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::angle_of_attack );
        break;
        //Get angle of sideslip if needed.
    case angle_of_sideslip_dependent:
        if( aerodynamicAngleCalculator_== NULL )
        {
            throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of sideslip" );
        }
        currentIndependentVariable = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::angle_of_sideslip );
        break;
    case altitude_dependent:
        if( aerodynamicAngleCalculator_== NULL )
        {
            throw std::runtime_error( "Error, aerodynamic angle calculator is null, but require angle of sideslip" );
        }
        currentIndependentVariable = getCurrentAltitude( );
        break;
    case control_surface_deflection_dependent:
    {
        try
        {
            currentIndependentVariable = controlSurfaceDeflectionFunction_( secondaryIdentifier );
        }
        catch( std::runtime_error )
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
void FlightConditions::updateAerodynamicCoefficientInput( )
{
    aerodynamicCoefficientIndependentVariables_.clear( );
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
    boost::function< std::map< std::string, std::vector< double > >( ) > untrimmedControlSurfaceIndependentVariablesFunction =
            boost::bind( &FlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables,
                         flightConditions );
    flightConditions->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::bind( &TrimOrientationCalculator::findTrimAngleOfAttackFromFunction, trimOrientation,
                             untrimmedIndependentVariablesFunction, untrimmedControlSurfaceIndependentVariablesFunction ) );

    return trimOrientation;
}

} // namespace aerodynamics

} // namespace tudat
