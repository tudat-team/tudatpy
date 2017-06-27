/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< aerodynamics::AerodynamicTorque > createAerodynamicTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    // Check existence of required environment models
    if( bodyUndergoingTorque->getAerodynamicCoefficientInterface( ) == NULL )
    {
        throw std::runtime_error( "Error when making aerodynamic torque, body " +
                                  nameOfBodyUndergoingTorque +
                                  "has no aerodynamic coefficients." );
    }

    if( bodyExertingTorque->getAtmosphereModel( ) == NULL )
    {
        throw std::runtime_error(  "Error when making aerodynamic torque, central body " +
                                   nameOfBodyExertingTorque + " has no atmosphere model.");
    }

    if( bodyExertingTorque->getShapeModel( ) == NULL )
    {
        throw std::runtime_error( "Error when making aerodynamic torque, central body " +
                                  nameOfBodyExertingTorque + " has no shape model." );
    }

    // Retrieve flight conditions; create object if not yet extant.
    boost::shared_ptr< aerodynamics::FlightConditions > bodyFlightConditions =
            bodyUndergoingTorque->getFlightConditions( );

    if( bodyFlightConditions == NULL )
    {
        bodyUndergoingTorque->setFlightConditions(
                    createFlightConditions( bodyUndergoingTorque,
                                            bodyExertingTorque,
                                            nameOfBodyUndergoingTorque,
                                            nameOfBodyExertingTorque ) );
        bodyFlightConditions = bodyUndergoingTorque->getFlightConditions( );
    }

    // Retrieve frame in which aerodynamic coefficients are defined.
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > aerodynamicCoefficients =
            bodyUndergoingTorque->getAerodynamicCoefficientInterface( );
    reference_frames::AerodynamicsReferenceFrames torqueFrame;
    if( aerodynamicCoefficients->getAreCoefficientsInAerodynamicFrame( ) )
    {
        torqueFrame = reference_frames::aerodynamic_frame;
    }
    else
    {
        torqueFrame = reference_frames::body_frame;
    }

    // Create function to transform from frame of aerodynamic coefficienrs to that of propagation.
    boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationFunction(
                bodyFlightConditions->getAerodynamicAngleCalculator( ),
                torqueFrame,
                boost::bind( &Body::getCurrentRotationToGlobalFrame, bodyExertingTorque ),
                reference_frames::body_frame );


    boost::function< Eigen::Vector3d( ) > coefficientFunction =
            boost::bind( &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                         aerodynamicCoefficients );
    boost::function< Eigen::Vector3d( ) > coefficientInPropagationFrameFunction =
            boost::bind( &reference_frames::transformVectorFunctionFromVectorFunctions,
                         coefficientFunction, toPropagationFrameTransformation );

    // Create torque model.
    return boost::make_shared< aerodynamics::AerodynamicTorque >(
                coefficientInPropagationFrameFunction,
                boost::bind( &aerodynamics::FlightConditions::getCurrentDensity, bodyFlightConditions ),
                boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed, bodyFlightConditions ),
                boost::bind( &aerodynamics::AerodynamicCoefficientInterface::getReferenceArea, aerodynamicCoefficients ),
                boost::bind( &aerodynamics::AerodynamicCoefficientInterface::getReferenceLengths, aerodynamicCoefficients ),
                aerodynamicCoefficients->getAreCoefficientsInNegativeAxisDirection( ) );
}

boost::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > createSecondDegreeGravitationalTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    boost::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction =
            boost::bind( &simulation_setup::Body::getPosition, bodyUndergoingTorque );
    boost::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction =
            boost::bind( &simulation_setup::Body::getPosition, bodyExertingTorque );


    boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel = bodyExertingTorque->getGravityFieldModel( );

    boost::function< double( ) > gravitationalParameterOfAttractingBodyFunction;
    if( gravityFieldModel ==  NULL )
    {
        std::cerr<<"Error when making second degree gravitational torque, "<<nameOfBodyExertingTorque<<" does not possess a gravity field"<<std::endl;
    }
    else
    {
        gravitationalParameterOfAttractingBodyFunction = boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                                      gravityFieldModel );
    }

    boost::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction  =
            boost::bind( &simulation_setup::Body::getBodyInertiaTensor, bodyUndergoingTorque );

    boost::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction =
            boost::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyUndergoingTorque );

    return boost::make_shared<gravitation::SecondDegreeGravitationalTorqueModel >(
                positionOfBodySubjectToTorqueFunction, gravitationalParameterOfAttractingBodyFunction, inertiaTensorOfRotatingBodyFunction,
                positionOfBodyExertingTorqueFunction, rotationToBodyFixedFrameFunction );

}

//! Function to create torque model object.
boost::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const boost::shared_ptr< TorqueSettings > torqueSettings,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel;

    switch( torqueSettings->torqueType_ )
    {
    case basic_astrodynamics::second_order_gravitational_torque:
    {
        torqueModel = createSecondDegreeGravitationalTorqueModel(
                    bodyUndergoingTorque, bodyExertingTorque, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque );
        break;
    }
    case basic_astrodynamics::aerodynamic_torque:
    {
        torqueModel = createAerodynamicTorqueModel(
                    bodyUndergoingTorque, bodyExertingTorque, nameOfBodyUndergoingTorque, nameOfBodyExertingTorque );
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize type " + boost::lexical_cast< std::string >( torqueSettings->torqueType_ ) +
                    " when making torque model" );
    }

    return torqueModel;
}


}

}

