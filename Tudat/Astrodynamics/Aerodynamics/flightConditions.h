#ifndef FLIGHTCONDITIONS_H
#define FLIGHTCONDITIONS_H

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


//! Class for calculating aerodynamic flight characteristics of a vehicle during numerical integration.
/*!
 *  Class for calculating aerodynamic flight characteristics of a vehicle during numerical integration. Class is used to ensure that dependent variables
 *  such as density, altitude, etc. are only calculated once during each numerical integration step. The get functions of this class are
 *  linked to the various models in the code that subsequently require these values.
 */
class FlightConditions
{
public:

    //! Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables are retrieved.
     *  \param centralBody Name of central body (i.e. body with the atmosphere)
     *  \param atmosphereModel Atmosphere model of atmosphere through which vehicle is flying
     *  \param altitudeFunction Function returning the altitude of the vehicle as a function of
     *  its body-fixed position.
     *  \param aerodynamicCoefficientInterface Class from which the aerodynamic (force and moment) coefficients are retrieved
     *  \param stateOfVehicle Function returning the current state of the vehicle (in the global frame)
     *  \param stateOfCentralBody Function returning the current state of the central body (in the global frame)
     *  \param transformationToCentralBodyFrame Function transforming the inertial body-centered to
     *  the body-centered, body-fixed (co-rotating) frame.
     *  \param currentTimeFunction Function returning the current time.
     *  \param aerodynamicCoefficientInterface Object from which the aerodynamic coefficients
     *  are obtained.
     *  \param aerodynamicAngleCalculator Object from which the aerodynamic/trajectory angles
     *  of the vehicle are calculated.
     */
    FlightConditions( const std::string& centralBody,
                      const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
                      const boost::function< double( const Eigen::Vector3d ) > altitudeFunction,
                      const boost::function< basic_mathematics::Vector6d( ) > stateOfVehicle,
                      const boost::function< basic_mathematics::Vector6d( ) > stateOfCentralBody,
                      const boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
                      transformationToCentralBodyFrame,
                      const boost::function< double( ) > currentTimeFunction,
                      const boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
                      const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
                      aerodynamicAngleCalculator = boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >( ) );

    //! Function to update all flight conditions.
    /*!
     *  Function to update all flight conditions (altitude, density, force coefficients) to current state of vehicle and central body.
     */
    void updateConditions( );

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
     *  Function to return density that was set by previous call of updateConditions or updateDensity function.
     *  \return Current altitude
     */
    double getCurrentDensity( ) const
    {
        return currentDensity_;
    }

    double getCurrentAirspeed( ) const
    {
        return currentAirspeed_;
    }


    //! Function to return central body name
    /*!
     *  Function to return central body name
     *  \return Name of central body
     */
    std::string getCentralBodyName( ) const
    {
        return centralBody_;
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

    void setAerodynamicAngleCalculator(
            const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator )
    {
        aerodynamicAngleCalculator_ = aerodynamicAngleCalculator;
    }

    basic_mathematics::Vector6d getCurrentBodyCenteredBodyFixedState( )
    {
        return currentBodyCenteredPseudoBodyFixedState_;

    }

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

    //! Function transforming the inertial body-centered to the body-centered, body-fixed (co-rotating) frame.
    boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
    transformationToCentralBodyFrame_;

    //! Function returning the current time.
    boost::function< double( ) > currentTimeFunction_;

    //! Object from which the aerodynamic coefficients are obtained.
    boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! Object from which the aerodynamic/trajectory angles of the vehicle are calculated.
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    basic_mathematics::Vector6d currentBodyCenteredState_;

    basic_mathematics::Vector6d currentBodyCenteredPseudoBodyFixedState_;
    //! Current density at vehicle's position.
    double currentDensity_;

    double currentAirspeed_;

    //! Current altitude of vehicle above central body's shapeModel_
    double currentAltitude_;

    double currentLatitude_;

    double currentLongitude_;

    double currentTime_;

    bool updateLatitudeAndLongitude_;
};

}

}
#endif // FLIGHTCONDITIONS_H
