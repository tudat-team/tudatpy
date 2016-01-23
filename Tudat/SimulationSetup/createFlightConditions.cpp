#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/SimulationSetup/accelerationModelTypes.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const boost::function< double( ) > angleOfAttackFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > angleOfSideslipFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > bankAngleFunction =
        boost::lambda::constant ( 0.0 ) )
{
    if( centralBody->getAtmosphereModel( ) == NULL )
    {
        throw( "" );
    }

    if( centralBody->getShapeModel( ) == NULL )
    {
        throw( "" );
    }

    if( centralBody->getRotationalEphemeris( ) == NULL )
    {
        throw( "" );
    }

    if( bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) == NULL )
    {
        throw( "" );
    }

    boost::function< double( const Eigen::Vector3d ) > altitudeFunction =
            boost::bind( &basic_astrodynamics::BodyShapeModel::getAltitude,
                         centralBody->getShapeModel( ), _1 );

    boost::function< Eigen::Quaterniond( ) > rotationToFrameFunction =
            boost::bind( &Body::getCurrentRotationToLocalFrame, centralBody );
    boost::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction =
            boost::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, centralBody );
    boost::function< basic_mathematics::Vector6d( const basic_mathematics::Vector6d& ) >
            transformationToCentralBodyFrame =
            boost::bind(
                static_cast< basic_mathematics::Vector6d(&)(
                    const basic_mathematics::Vector6d&,
                    const boost::function< Eigen::Quaterniond( ) >,
                    const boost::function< Eigen::Matrix3d( ) > ) >( &ephemerides::transformStateToFrame ),
                _1, rotationToFrameFunction,
                rotationMatrixToFrameDerivativeFunction );


    boost::shared_ptr< aerodynamics::FlightConditions > flightConditions =
            boost::make_shared< aerodynamics::FlightConditions >(
                centralBody->getAtmosphereModel( ), altitudeFunction,
                boost::bind( &Body::getState, bodyWithFlightConditions ),
                boost::bind( &Body::getState, centralBody ),
                transformationToCentralBodyFrame,
                boost::bind( &Body::getCurrentTime, centralBody ),
                bodyWithFlightConditions->getAerodynamicCoefficientInterface( ) );

    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
            boost::make_shared< reference_frames::AerodynamicAngleCalculator >(
                boost::bind( &aerodynamics::FlightConditions::getCurrentBodyCenteredBodyFixedState,
                             flightConditions ),
                angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction );

    flightConditions->setAerodynamicAngleCalculator( aerodynamicAngleCalculator );

    return flightConditions;


}

}

}
