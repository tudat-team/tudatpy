/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace aerodynamics
{


//! Class for calculating aerodynamic flight characteristics of a vehicle during numerical
//! integration.
/*!
 *  Class for calculating aerodynamic flight characteristics of a vehicle during numerical
 *  integration. Class is used to ensure that dependent variables such as density, altitude, etc.
 *  are only calculated once during each numerical integration step. The get functions of this class
 *  are linked to the various models in the code that subsequently require these values.
 */
class FlightConditions
{
public:

    //! Constructor, sets objects and functions from which relevant environment and state variables
    //! are retrieved.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param atmosphereModel Atmosphere model of atmosphere through which vehicle is flying
     *  \param altitudeFunction Function returning the altitude of the vehicle as a function of
     *  its body-fixed position.
     *  \param aerodynamicCoefficientInterface Class from which the aerodynamic (force and moment)
     *  coefficients are retrieved
     *  \param stateOfVehicle Function returning the current state of the vehicle
     *  (in the global frame)
     *  \param stateOfCentralBody Function returning the current state of the central body
     *  (in the global frame)
     *  \param transformationToCentralBodyFrame Function transforming the inertial body-centered to
     *  the body-centered, body-fixed (co-rotating) frame.
     *  \param aerodynamicCoefficientInterface Object from which the aerodynamic coefficients
     *  are obtained.
     *  \param aerodynamicAngleCalculator Object from which the aerodynamic/trajectory angles
     *  of the vehicle are calculated.
     */
    FlightConditions( const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
                      const boost::function< double( const Eigen::Vector3d ) > altitudeFunction,
                      const boost::function< basic_mathematics::Vector6d( ) > stateOfVehicle,
                      const boost::function< basic_mathematics::Vector6d( ) > stateOfCentralBody,
                      const boost::function< basic_mathematics::Vector6d(
                          const basic_mathematics::Vector6d& ) >
                      transformationToCentralBodyFrame,
                      const boost::shared_ptr< AerodynamicCoefficientInterface >
                      aerodynamicCoefficientInterface,
                      const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
                      aerodynamicAngleCalculator =
            boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >( ) );

    //! Function to update all flight conditions.
    /*!
     *  Function to update all flight conditions (altitude, density, force coefficients) to
     *  current state of vehicle and central body.
     *  \param currentTime Time to which conditions are to be updated.
     */
    void updateConditions( const double currentTime );

    //! Function to return altitude
    /*!
     *  Function to return altitude that was set by previous call of updateConditions function.
     *  \return Current altitude
     */
    double getCurrentAltitude( ) const
    {
        return currentAltitude_;
    }

    //! Function to return density
    /*!
     *  Function to return density that was set by previous call of updateConditions or
     *  updateDensity function.
     *  \return Current density
     */
    double getCurrentDensity( ) const
    {
        return currentDensity_;
    }

    //! Function to return airspeed
    /*!
     *  Function to return airspeed that was set by previous call of updateConditions.
     *  \return Current airspeed
     */
    double getCurrentAirspeed( ) const
    {
        return currentAirspeed_;
    }

    //! Function to return speed of sound
    /*!
     *  Function to return speed of from quantities computed by previous call of updateConditions.
     *  \return Current speed of sound.
     */
    double getCurrentSpeedOfSound( )
    {
        return atmosphereModel_->getSpeedOfSound(
                    currentAltitude_, currentLongitude_, currentLatitude_, currentTime_ );
    }

    //! Function to return the current time of the FlightConditions
    /*!
     *  Function to return the current time of the FlightConditions.
     *  \return Current time of the FlightConditions
     */
    double getCurrentTime( ) const
    {
        return currentTime_;
    }
    //! Function to return atmosphere model object
    /*!
     *  Function to return atmosphere model object
     *  \return Atmosphere model object
     */
    boost::shared_ptr< aerodynamics::AtmosphereModel > getAtmosphereModel( ) const
    {
        return atmosphereModel_;
    }

    //! Function to (re)set aerodynamic angle calculator object
    /*!
     *  Function to (re)set aerodynamic angle calculator object
     *  \param aerodynamicAngleCalculator Aerodynamic angle calculator object to set.
     */
    void setAerodynamicAngleCalculator(
            const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
            aerodynamicAngleCalculator )
    {
        aerodynamicAngleCalculator_ = aerodynamicAngleCalculator;
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
            const boost::function< double( ) > coefficientDependency );

    //! Function to return current central body-fixed state of vehicle.
    /*!
     *  Function to return central body-fixed state of vehicle.
     *  \return Current central body-fixed state of vehicle.
     */
    basic_mathematics::Vector6d getCurrentBodyCenteredBodyFixedState( )
    {
        return currentBodyCenteredPseudoBodyFixedState_;

    }

    //! Function to return aerodynamic angle calculator object
    /*!
     *  Function to return aerodynamic angle calculator object
     *  \return Aerodynamic angle calculator object
     */
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
    getAerodynamicAngleCalculator( )
    {
        return aerodynamicAngleCalculator_;
    }


private:

    //! Name of central body (i.e. body with the atmosphere)
    std::string centralBody_;

    //! Atmosphere model of atmosphere through which vehicle is flying
    boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel_;

    //! Function returning the altitude of the vehicle as a function of its body-fixed position.
    const boost::function< double( const Eigen::Vector3d ) > altitudeFunction_;

    //! Function returning the current state of the vehicle (in the global frame)
    boost::function< basic_mathematics::Vector6d( ) > stateOfVehicle_;

    //! Function returning the current state of the central body (in the global frame)
    boost::function< basic_mathematics::Vector6d( ) > stateOfCentralBody_;

    //! Function transforming the inertial body-centered to the body-centered, body-fixed
    //! co-rotating) frame.
    boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
    transformationToCentralBodyFrame_;

    //! Object from which the aerodynamic coefficients are obtained.
    boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! Object from which the aerodynamic/trajectory angles of the vehicle are calculated.
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    //! Current state of vehicle in base frame for Body objects.
    basic_mathematics::Vector6d currentBodyCenteredState_;

    //! Current state of vehicle in body-fixed frame.
    basic_mathematics::Vector6d currentBodyCenteredPseudoBodyFixedState_;

    //! Current density at vehicle's position.
    double currentDensity_;

    //! Current airspeed at vehicle's position.
    double currentAirspeed_;

    //! Current altitude of vehicle above central body's shapeModel_
    double currentAltitude_;

    //! Current latitude of vehicle above central body.
    double currentLatitude_;

    //! Current longitude of vehicle above central body.
    double currentLongitude_;

    //! Current time of propagation.
    double currentTime_;

    //! List of custom functions for aerodynamic coefficient dependencies.
    std::map< AerodynamicCoefficientsIndependentVariables, boost::function< double( ) > > customCoefficientDependencies_;

    //! Boolean setting whether latitude and longitude are to be updated by updateConditions().
    bool updateLatitudeAndLongitude_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_FLIGHTCONDITIONS_H
