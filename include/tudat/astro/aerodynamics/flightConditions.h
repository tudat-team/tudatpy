/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FLIGHTCONDITIONS_H
#define TUDAT_FLIGHTCONDITIONS_H

#include <vector>

#include <functional>
#include <boost/bind/bind.hpp>

#include "tudat/astro/aerodynamics/trimOrientation.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/basics/basicTypedefs.h"

using namespace boost::placeholders;

namespace tudat
{

namespace aerodynamics
{

//! Class for calculating aerodynamic flight characteristics of a vehicle during numerical
//! integration, in the absence of an atmosphere.
/*!
 *  Class for calculating aerodynamic flight characteristics of a vehicle during numerical
 *  integration, in the absence of an atmosphere. Class is used to ensure that dependent variables such as altitude, etc.
 *  are only calculated once during each numerical integration step. The get functions of this class are linked to the various
 *  models in the code that subsequently require these values. In the case of atmospheric flight, the AtmosphericFlightConditions
 *  derived class should be used.
 */
class FlightConditions
{
protected:

    //! List of variables that can be computed by flight condition
    enum FlightConditionVariables
    {
        altitude_flight_condition = 0,
        density_flight_condition = 1,
        pressure_flight_condition = 2,
        temperature_flight_condition = 3,
        latitude_flight_condition = 4,
        longitude_flight_condition = 5,
        mach_number_flight_condition = 6,
        speed_of_sound_flight_condition = 7,
        airspeed_flight_condition = 8,
        geodetic_latitude_condition = 9,
        dynamic_pressure_condition = 10,
        aerodynamic_heat_rate = 11
    };

public:

    //! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param shapeModel Model describing the shape of the body w.r.t. which the flight is taking place.
     *  \param aerodynamicAngleCalculator Object from which the aerodynamic/trajectory angles
     *  of the vehicle are calculated.
     */
    FlightConditions( const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
                      const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            std::shared_ptr< reference_frames::AerodynamicAngleCalculator >( ) );

    //! Destructor
    virtual ~FlightConditions( ){ }

    //! Function to update all flight conditions.
    /*!
     *  Function to update all flight conditions (altitude, etc.) to
     *  current state of vehicle and central body.
     *  \param currentTime Time to which conditions are to be updated.
     */
    virtual void updateConditions( const double currentTime );

    //! Function to retrieve (and compute if necessary) the current altitude
    /*!
     * Function to retrieve (and compute if necessary) the current altitude
     * \return Current altitude
     */
    double getCurrentAltitude( )
    {
        if( isScalarFlightConditionComputed_.at( altitude_flight_condition ) == 0 )
        {
            computeAltitude( );
        }
        return scalarFlightConditions_.at( altitude_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current longitude
    /*!
     * Function to retrieve (and compute if necessary) the current longitude
     * \return Current longitude
     */
    double getCurrentLongitude( )
    {
        if( isScalarFlightConditionComputed_.at( longitude_flight_condition ) == 0 )
        {
            computeLatitudeAndLongitude( );
        }
        return scalarFlightConditions_.at( longitude_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current geodetic latitude
    /*!
     * Function to retrieve (and compute if necessary) the current geodetic latitude
     * \return Current geodetic latitude
     */
    double getCurrentGeodeticLatitude( )
    {
        if( isScalarFlightConditionComputed_.at( geodetic_latitude_condition ) == 0 )
        {
            computeGeodeticLatitude( );
        }
        return scalarFlightConditions_.at( geodetic_latitude_condition );
    }

    //! Function to return the current time of the AtmosphericFlightConditions
    /*!
     *  Function to return the current time of the AtmosphericFlightConditions.
     *  \return Current time of the AtmosphericFlightConditions
     */
    double getCurrentTime( )
    {
        return currentTime_;
    }

    //! Function to (re)set aerodynamic angle calculator object
    /*!
     *  Function to (re)set aerodynamic angle calculator object
     *  \param aerodynamicAngleCalculator Aerodynamic angle calculator object to set.
     */
    void setAerodynamicAngleCalculator(
            const std::shared_ptr< reference_frames::AerodynamicAngleCalculator >
            aerodynamicAngleCalculator )
    {
        aerodynamicAngleCalculator_ = aerodynamicAngleCalculator;
        bodyCenteredPseudoBodyFixedStateFunction_ = std::bind(
                    &reference_frames::AerodynamicAngleCalculator::getCurrentAirspeedBasedBodyFixedState, aerodynamicAngleCalculator_ );
    }

    //! Function to return aerodynamic angle calculator object
    /*!
     *  Function to return aerodynamic angle calculator object
     *  \return Aerodynamic angle calculator object
     */
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator >
    getAerodynamicAngleCalculator( )
    {
        return aerodynamicAngleCalculator_;
    }

    //! Function to reset the current time of the flight conditions.
    /*!
     *  Function to reset the current time of the flight conditions. This function is typically sused to set the current time
     *  to NaN, indicating the need to recompute all quantities for the next time computation.
     * \param currentTime
     */
    virtual void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        isScalarFlightConditionComputed_ = allScalarFlightConditionsUncomputed;
        aerodynamicAngleCalculator_->resetCurrentTime( currentTime_ );
    }

    //! Function to return current central body-fixed state of vehicle.
    /*!
     *  Function to return central body-fixed state of vehicle.
     *  \return Current central body-fixed state of vehicle.
     */
    Eigen::Vector6d getCurrentBodyCenteredBodyFixedState( )
    {
        return currentBodyCenteredAirspeedBasedBodyFixedState_;
    }

protected:

    //! Function to compute and set the current latitude and longitude
    void computeLatitudeAndLongitude( )
    {
        scalarFlightConditions_[ latitude_flight_condition ] = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::latitude_angle );
        isScalarFlightConditionComputed_[ latitude_flight_condition ] = true;

        scalarFlightConditions_[ longitude_flight_condition ] = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::longitude_angle );
        isScalarFlightConditionComputed_[ longitude_flight_condition ] = true;
    }

    //! Function to compute and set the current altitude
    void computeAltitude( )
    {
        scalarFlightConditions_[ altitude_flight_condition ] =
                shapeModel_->getAltitude( currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 0, 3 ) );
        isScalarFlightConditionComputed_[ altitude_flight_condition ] = true;
    }

    //! Function to compute and set the current geodetic latitude.
    void computeGeodeticLatitude( )
    {
        if( !( geodeticLatitudeFunction_ == nullptr ) )
        {
            scalarFlightConditions_[ geodetic_latitude_condition ] = geodeticLatitudeFunction_(
                        currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 0, 3 ) );
            isScalarFlightConditionComputed_[ geodetic_latitude_condition ] = true;
        }
        else
        {
            if( isScalarFlightConditionComputed_.at( latitude_flight_condition ) == 0 )
            {
                computeLatitudeAndLongitude( );
            }
            scalarFlightConditions_[ geodetic_latitude_condition ] = scalarFlightConditions_[ latitude_flight_condition ] ;
            isScalarFlightConditionComputed_[ geodetic_latitude_condition ] = true;
        }
    }

    //! Name of central body (i.e. body with the atmosphere)
    std::string centralBody_;

    //! Model describing the shape of the body w.r.t. which the flight is taking place.
    const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel_;

    //! Object from which the aerodynamic/trajectory angles of the vehicle are calculated.
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    //! Function to return the current state of the vehicle in a body-fixed frame.
    std::function< Eigen::Vector6d( ) > bodyCenteredPseudoBodyFixedStateFunction_;

    //! Current state of vehicle in base frame for Body objects.
    Eigen::Vector6d currentBodyCenteredState_;

    //! Current state of vehicle in body-fixed frame.
    Eigen::Vector6d currentBodyCenteredAirspeedBasedBodyFixedState_;

    //! Current time of propagation.
    double currentTime_;

    //! List of atmospheric/flight properties computed at current time step.
    std::vector< double > scalarFlightConditions_;

    std::vector< bool > isScalarFlightConditionComputed_;

    const std::vector< bool > allScalarFlightConditionsUncomputed = std::vector< bool >( 12, false );

    //! Function from which to compute the geodetic latitude as function of body-fixed position (empty if equal to
    //! geographic latitude).
    std::function< double( const Eigen::Vector3d& ) > geodeticLatitudeFunction_;

};

//! Class for calculating aerodynamic flight characteristics of a vehicle during numerical
//! integration.
/*!
 *  Class for calculating aerodynamic flight characteristics of a vehicle during numerical
 *  integration. Class is used to ensure that dependent variables such as density, altitude, etc.
 *  are only calculated once during each numerical integration step. The get functions of this class
 *  are linked to the various models in the code that subsequently require these values.
 */
class AtmosphericFlightConditions: public FlightConditions
{

public:

    //! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param atmosphereModel Atmosphere model of atmosphere through which vehicle is flying
     *  \param shapeModel Model describing the shape of the body w.r.t. which the flight is taking place.
     *  \param aerodynamicCoefficientInterface Class from which the aerodynamic (force and moment)
     *  coefficients are retrieved
     *  \param aerodynamicAngleCalculator Object from which the aerodynamic/trajectory angles
     *  of the vehicle are calculated.
     *  \param controlSurfaceDeflectionFunction Function returning control surface deflection, with input the control
     *  surface identifier.
     */
    AtmosphericFlightConditions( const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
                                 const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel,
                                 const std::shared_ptr< AerodynamicCoefficientInterface >
                                 aerodynamicCoefficientInterface,
                                 const std::shared_ptr< reference_frames::AerodynamicAngleCalculator >
                                 aerodynamicAngleCalculator,
                                 const std::function< double( const std::string& )> controlSurfaceDeflectionFunction =
            std::function< double( const std::string& )>( ) );


    //! Function to update all flight conditions.
    /*!
     *  Function to update all flight conditions (altitude, density, force coefficients) to
     *  current state of vehicle and central body.
     *  \param currentTime Time to which conditions are to be updated.
     */
    void updateConditions( const double currentTime );

    //! Function to retrieve (and compute if necessary) the current freestream density
    /*!
     * Function to retrieve (and compute if necessary) the current freestream density
     * \return Current freestream density
     */
    double getCurrentDensity( )
    {
        if( isScalarFlightConditionComputed_.at( density_flight_condition ) == 0 )
        {
            computeDensity( );
        }
        return scalarFlightConditions_.at( density_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current freestream temperature
    /*!
     * Function to retrieve (and compute if necessary) the current freestream temperature
     * \return Current freestream temperature
     */
    double getCurrentFreestreamTemperature( )
    {
        if( isScalarFlightConditionComputed_.at( temperature_flight_condition ) == 0 )
        {
            computeTemperature( );
        }
        return scalarFlightConditions_.at( temperature_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current freestream dynamic pressure
    /*!
     * Function to retrieve (and compute if necessary) the current freestream dynamic pressure
     * \return Current freestream dynamic pressure
     */
    double getCurrentDynamicPressure( )
    {
        if( isScalarFlightConditionComputed_.at( dynamic_pressure_condition ) == 0 )
        {
            computeDynamicPressure( );
        }
        return scalarFlightConditions_.at( dynamic_pressure_condition );
    }

    //! Function to retrieve (and compute if necessary) the current aerodynamic heat rate
    /*!
     * Function to retrieve (and compute if necessary) the current aerodynamic heat rate
     * \return Current aerodynamic heat rate
     */
    double getCurrentAerodynamicHeatRate( )
    {
        if( isScalarFlightConditionComputed_.at( aerodynamic_heat_rate ) == 0 )
        {
            computeAerodynamicHeatRate( );
        }
        return scalarFlightConditions_.at( aerodynamic_heat_rate );
    }

    //! Function to retrieve (and compute if necessary) the current freestream pressure
    /*!
     * Function to retrieve (and compute if necessary) the current freestream pressure
     * \return Current freestream dynamic pressure
     */
    double getCurrentPressure( )
    {
        if( isScalarFlightConditionComputed_.at( pressure_flight_condition ) == 0 )
        {
            computeFreestreamPressure( );
        }
        return scalarFlightConditions_.at( pressure_flight_condition );
    }

    /*!
     * Function to retrieve (and compute if necessary) the current airspeed
     * \return Current airspeed
     */
    double getCurrentAirspeed( )
    {
        if( isScalarFlightConditionComputed_.at( airspeed_flight_condition ) == 0 )
        {
            computeAirspeed( );
        }
        return scalarFlightConditions_.at( airspeed_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current speed of sound
    /*!
     * Function to retrieve (and compute if necessary) the current speed of sound
     * \return Current speed of sound
     */
    double getCurrentSpeedOfSound( )
    {
        if( isScalarFlightConditionComputed_.at( speed_of_sound_flight_condition ) == 0 )
        {
            computeSpeedOfSound( );
        }
        return scalarFlightConditions_.at( speed_of_sound_flight_condition );
    }

    //! Function to retrieve (and compute if necessary) the current Mach number
    /*!
     * Function to retrieve (and compute if necessary) the current Mach number
     * \return Current Mach number
     */
    double getCurrentMachNumber( )
    {
        if( isScalarFlightConditionComputed_.at( mach_number_flight_condition ) == 0 )
        {
            computeMachNumber( );
        }
        return scalarFlightConditions_.at( mach_number_flight_condition );
    }

    //! Function to return atmosphere model object
    /*!
     *  Function to return atmosphere model object
     *  \return Atmosphere model object
     */
    std::shared_ptr< aerodynamics::AtmosphereModel > getAtmosphereModel( ) const
    {
        return atmosphereModel_;
    }

    //! Function to set custom dependency of aerodynamic coefficients
    /*!
     * Function to set custom dependency of aerodynamic coefficients. If needed, the
     * AerodynamicCoefficientsIndependentVariables enum may be expanded to include e.g. control surface deflections, the
     * values of which will then be retrieved from the function set here
     * \param independentVariable Identifier of independent variable
     * \param coefficientDependency Function returning the current value of the independent variable.
     */
    void setAerodynamicCoefficientsIndependentVariableFunction(
            const AerodynamicCoefficientsIndependentVariables independentVariable,
            const std::function< double( ) > coefficientDependency );

    //! Function to return current central body-fixed velocity of vehicle.
    /*!
     *  Function to return central body-fixed velocity of vehicle.
     *  \return Current central body-fixed velocity of vehicle.
     */
    Eigen::Vector3d getCurrentAirspeedBasedVelocity( )
    {
        return currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 3, 3 );
    }

    //! Function to return object from which the aerodynamic coefficients are obtained.
    /*!
     *  Function to return object from which the aerodynamic coefficients are obtained.
     *  \return Object from which the aerodynamic coefficients are obtained.
     */
    std::shared_ptr< AerodynamicCoefficientInterface > getAerodynamicCoefficientInterface( )
    {
        return aerodynamicCoefficientInterface_;
    }

    //! Function to return list of independent variables of the aerodynamic coefficient interface
    /*!
     *  Function to return list of independent variables of the aerodynamic coefficient interface
     *  \return List of independent variables of the aerodynamic coefficient interface
     */
    std::vector< double > getAerodynamicCoefficientIndependentVariables( )
    {
        if( aerodynamicCoefficientIndependentVariables_.size( ) !=
                aerodynamicCoefficientInterface_->getNumberOfIndependentVariables( ) )
        {
            updateAerodynamicCoefficientInput( );
        }

        return aerodynamicCoefficientIndependentVariables_;
    }

    //! Function to return list of independent variables of the control surface aerodynamic coefficient interface
    /*!
     *  Function to return list of independent variables of the control surface aerodynamic coefficient interface
     *  \return List of independent variables of the control surface aerodynamic coefficient interface, with map key
     *  the control surface identifiers.
     */
    std::map< std::string, std::vector< double > > getControlSurfaceAerodynamicCoefficientIndependentVariables( )
    {
        if( controlSurfaceAerodynamicCoefficientIndependentVariables_.size( ) !=
                aerodynamicCoefficientInterface_->getNumberOfControlSurfaces( ) )
        {
            updateAerodynamicCoefficientInput( );
        }

        return controlSurfaceAerodynamicCoefficientIndependentVariables_;
    }

    //! Function to reset the current time of the flight conditions.
    /*!
     *  Function to reset the current time of the flight conditions. This function is typically sused to set the current time
     *  to NaN, indicating the need to recompute all quantities for the next time computation.
     * \param currentTime
     */
    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;

        isScalarFlightConditionComputed_ = allScalarFlightConditionsUncomputed;
        aerodynamicAngleCalculator_->resetCurrentTime( currentTime_ );
        aerodynamicCoefficientIndependentVariables_.clear( );
        controlSurfaceAerodynamicCoefficientIndependentVariables_.clear( );
    }

private:

    //! Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
    /*!
     * Function to (compute and) retrieve the value of an independent variable of aerodynamic coefficients
     * \param independentVariableType Identifier for type of independent variable
     * \param secondaryIdentifier String used as secondary identifier of independent variable (e.g. control surface name).
     * \return Current value of requested independent variable.
     */
    double getAerodynamicCoefficientIndependentVariable(
            const AerodynamicCoefficientsIndependentVariables independentVariableType,
            const std::string& secondaryIdentifier = "" );

    //! Function to update input to atmosphere model (altitude, as well as latitude and longitude if needed).
    void updateAtmosphereInput( );

    //! Function to compute and set the current freestream density
    void computeDensity( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ density_flight_condition ] =
                atmosphereModel_->getDensity(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
        isScalarFlightConditionComputed_[ density_flight_condition ] = true;
    }

    //! Function to compute and set the current freestream temperature
    void computeTemperature( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ temperature_flight_condition ] =
                atmosphereModel_->getTemperature(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
        isScalarFlightConditionComputed_[ temperature_flight_condition ] = true;
    }

    //! Function to compute and set the current freestream pressure.
    void computeFreestreamPressure( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ pressure_flight_condition ] =
                atmosphereModel_->getPressure(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
        isScalarFlightConditionComputed_[ pressure_flight_condition ] = true;
    }


    //! Function to compute and set the current speed of sound
    void computeSpeedOfSound( )
    {
        updateAtmosphereInput( );
        scalarFlightConditions_[ speed_of_sound_flight_condition ]  =
                atmosphereModel_->getSpeedOfSound(
                    scalarFlightConditions_.at( altitude_flight_condition ),
                    scalarFlightConditions_.at( longitude_flight_condition ),
                    scalarFlightConditions_.at( latitude_flight_condition ), currentTime_ );
        isScalarFlightConditionComputed_[ speed_of_sound_flight_condition ] = true;
    }

    //! Function to compute and set the current airspeed
    void computeAirspeed( )
    {
        scalarFlightConditions_[ airspeed_flight_condition ] = currentBodyCenteredAirspeedBasedBodyFixedState_.segment( 3, 3 ).norm( );
        isScalarFlightConditionComputed_[ airspeed_flight_condition ] = true;
    }

    //! Function to compute and set the current freestream dynamic pressure.
    void computeDynamicPressure( )
    {
        double currentAirspeed = getCurrentAirspeed( );
        scalarFlightConditions_[ dynamic_pressure_condition ] = 0.5 *
                getCurrentDensity( ) * currentAirspeed * currentAirspeed;
        isScalarFlightConditionComputed_[ dynamic_pressure_condition ] = true;
    }

    //! Function to compute and set the current aerodynamic heat rate.
    void computeAerodynamicHeatRate( )
    {
        double currentAirspeed = getCurrentAirspeed( );
        scalarFlightConditions_[ aerodynamic_heat_rate ] = 0.5 *
                getCurrentDensity( ) * currentAirspeed * currentAirspeed * currentAirspeed;
        isScalarFlightConditionComputed_[ aerodynamic_heat_rate ] = true;
    }

    //! Function to compute and set the current Mach number.
    void computeMachNumber( )
    {
        scalarFlightConditions_[ mach_number_flight_condition ] =
                getCurrentAirspeed( ) / getCurrentSpeedOfSound( );
        isScalarFlightConditionComputed_[ mach_number_flight_condition ] = true;
    }

    //! Function to update the independent variables of the aerodynamic coefficient interface
    void updateAerodynamicCoefficientInput( );


    //! Atmosphere model of atmosphere through which vehicle is flying
    std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel_;

    //! Object from which the aerodynamic coefficients are obtained.
    std::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! Function returning control surface deflection, with input the control surface identifier.
    std::function< double( const std::string& ) > controlSurfaceDeflectionFunction_;


    //! List of custom functions for aerodynamic coefficient dependencies.
    std::map< AerodynamicCoefficientsIndependentVariables, std::function< double( ) > > customCoefficientDependencies_;


    //! Current list of independent variables of the aerodynamic coefficient interface
    std::vector< double > aerodynamicCoefficientIndependentVariables_;

    //! List of independent variables of the control surface aerodynamic coefficient interface, with map key
    //! control surface identifiers.
    std::map< std::string, std::vector< double > > controlSurfaceAerodynamicCoefficientIndependentVariables_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_FLIGHTCONDITIONS_H
