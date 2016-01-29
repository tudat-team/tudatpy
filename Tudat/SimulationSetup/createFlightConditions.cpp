#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > createAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    boost::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface;

    switch( coefficientSettings->getAerodynamicCoefficientType( ) )
    {
    case constant_aerodynamic_coefficients:
    {
        boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantCoefficientSettings =
                boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >(
                    coefficientSettings );

        if( constantCoefficientSettings == NULL )
        {
            std::cerr<<"Error, expected constant aerodynamic coefficients for body "<<body<<std::endl;
        }
        else
        {
            coefficientInterface = createConstantCoefficientAerodynamicCoefficientInterface(
                        constantCoefficientSettings->getConstantForceCoefficient( ),
                        constantCoefficientSettings->getConstantMomentCoefficient( ),
                        constantCoefficientSettings->getReferenceLength( ),
                        constantCoefficientSettings->getReferenceArea( ),
                        constantCoefficientSettings->getReferenceLength( ),
                        constantCoefficientSettings->getMomentReferencePoint( ),
                        constantCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                        constantCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
        }
        break;
    }
    case tabulated_coefficients:
    {
        int numberOfDimensions = coefficientSettings->getIndependentVariableNames( ).size( );

        switch( numberOfDimensions )
        {
        case 1:
        {
            boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 1 > > tabulatedCoefficientSettings =
                    boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 1 > >(
                        coefficientSettings );

            if( tabulatedCoefficientSettings == NULL )
            {
                std::cerr<<"Error, expected tabulated aerodynamic coefficients of size 1 for body "<<body<<std::endl;
            }
            else
            {
                coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 1 >(
                            tabulatedCoefficientSettings->getIndependentVariables( ),
                            tabulatedCoefficientSettings->getForceCoefficients( ),
                            tabulatedCoefficientSettings->getMomentCoefficients( ),
                            tabulatedCoefficientSettings->getIndependentVariableNames( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getReferenceArea( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getMomentReferencePoint( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
            }
        }
        case 2:
        {

            boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 2 > > tabulatedCoefficientSettings =
                    boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 2 > >(
                        coefficientSettings );

            if( tabulatedCoefficientSettings == NULL )
            {
                std::cerr<<"Error, expected tabulated aerodynamic coefficients of size 2 for body "<<body<<std::endl;
            }
            else
            {
                coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 2 >(
                            tabulatedCoefficientSettings->getIndependentVariables( ),
                            tabulatedCoefficientSettings->getForceCoefficients( ),
                            tabulatedCoefficientSettings->getMomentCoefficients( ),
                            tabulatedCoefficientSettings->getIndependentVariableNames( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getReferenceArea( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getMomentReferencePoint( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
            }
        }
        case 3:
        {

            boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< 3 > > tabulatedCoefficientSettings =
                    boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< 3 > >(
                        coefficientSettings );

            if( tabulatedCoefficientSettings == NULL )
            {
                std::cerr<<"Error, expected tabulated aerodynamic coefficients of size 3 for body "<<body<<std::endl;
            }
            else
            {
                coefficientInterface = createTabulatedCoefficientAerodynamicCoefficientInterface< 3 >(
                            tabulatedCoefficientSettings->getIndependentVariables( ),
                            tabulatedCoefficientSettings->getForceCoefficients( ),
                            tabulatedCoefficientSettings->getMomentCoefficients( ),
                            tabulatedCoefficientSettings->getIndependentVariableNames( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getReferenceArea( ),
                            tabulatedCoefficientSettings->getReferenceLength( ),
                            tabulatedCoefficientSettings->getMomentReferencePoint( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                            tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
            }
        }
        default:
            std::cerr<<"Error when making tabulated aerodynamic coefficient interface, "<<
                       numberOfDimensions<<" dimensions not yet implemented"<<std::endl;
        }

    }
    default:
        std::cerr<<"Error, did not recognize aerodynamic coefficient settings for "<<body<<std::endl;
    }
    return coefficientInterface;
}

boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const boost::function< double( ) > angleOfAttackFunction,
        const boost::function< double( ) > angleOfSideslipFunction,
        const boost::function< double( ) > bankAngleFunction )
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
