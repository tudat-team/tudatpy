
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Aerodynamics/standardAtmosphere.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace aerodynamics
{


FlightConditions::FlightConditions(
        const std::string& centralBody,
        const boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel,
        const boost::function< double( const Eigen::Vector3d ) > altitudeFunction,
        const boost::function< basic_mathematics::Vector6d( ) > stateOfVehicle,
        const boost::function< basic_mathematics::Vector6d( ) > stateOfCentralBody,
        const boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
        transformationToCentralBodyFrame,
        const boost::function< double( ) > currentTimeFunction,
        const boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator >
        aerodynamicAngleCalculator ):
    centralBody_( centralBody ), atmosphereModel_( atmosphereModel ),
    altitudeFunction_( altitudeFunction ),
    stateOfVehicle_( stateOfVehicle ),
    stateOfCentralBody_( stateOfCentralBody ),
    transformationToCentralBodyFrame_( transformationToCentralBodyFrame ),
    currentTimeFunction_( currentTimeFunction ),
    aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
    aerodynamicAngleCalculator_( aerodynamicAngleCalculator )
{
    updateLatitudeAndLongitude_ = 0;

    if( boost::dynamic_pointer_cast< aerodynamics::StandardAtmosphere >( atmosphereModel_ ) == NULL )
    {
        throw( "" );
    }

    if( updateLatitudeAndLongitude_ && aerodynamicAngleCalculator_== NULL )
    {
        throw( "" );
    }
}

void FlightConditions::updateConditions(  )
{
    currentTime_ = currentTimeFunction_( );
    currentBodyCenteredState_ = stateOfVehicle_( ) - stateOfCentralBody_( );
    currentBodyCenteredPseudoBodyFixedState_ = transformationToCentralBodyFrame_(
                currentBodyCenteredState_ );
    currentAltitude_ = altitudeFunction_( currentBodyCenteredPseudoBodyFixedState_.segment( 0, 3 ) );
    currentAirspeed_ = currentBodyCenteredPseudoBodyFixedState_.segment( 3, 3 ).norm( );

    if( aerodynamicAngleCalculator_!= NULL )
    {
        aerodynamicAngleCalculator_->update( );
    }

    if( updateLatitudeAndLongitude_ )
    {
        currentLatitude_ = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::latitude_angle );
        currentLongitude_ = aerodynamicAngleCalculator_->getAerodynamicAngle(
                    reference_frames::longitude_angle );
    }

    std::vector< double > aerodynamicCoefficientIndependentVariables;
    for( unsigned int i = 0; i < aerodynamicCoefficientInterface_->
         getNumberOfIndependentVariables( ); i++ )
    {
        switch( aerodynamicCoefficientInterface_->getIndependentVariableName( i ) )
        {
        case mach_number_dependent:
            aerodynamicCoefficientIndependentVariables.push_back(
                        currentAirspeed_ / atmosphereModel_->getSpeedOfSound(
                            currentAltitude_, currentLongitude_, currentLatitude_, currentTime_ ) );
            break;
        case angle_of_attack_dependent:
            if( aerodynamicAngleCalculator_== NULL )
            {
                throw( "" );
            }
            aerodynamicCoefficientIndependentVariables.push_back(
                        aerodynamicAngleCalculator_->getAerodynamicAngle( reference_frames::angle_of_attack ) );
            break;
        case angle_of_sideslip_dependent:
            if( aerodynamicAngleCalculator_== NULL )
            {
                throw( "" );
            }
            aerodynamicCoefficientIndependentVariables.push_back(
                        aerodynamicAngleCalculator_->getAerodynamicAngle( reference_frames::angle_of_sideslip ) );
        default:
            throw( "" );
        }
    }

    currentDensity_ = atmosphereModel_->getDensity( currentAltitude_, currentLongitude_,
                                                    currentLatitude_, currentTime_ );

}


}

}
