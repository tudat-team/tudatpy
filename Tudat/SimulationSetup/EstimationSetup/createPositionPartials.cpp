#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/SimulationSetup/EstimationSetup/createPositionPartials.h"


namespace tudat
{

namespace observation_partials
{


//! Function to create partial object(s) of rotation matrix wrt a (double) parameter.
boost::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterToEstimate )
{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    boost::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    boost::shared_ptr< Body > currentBody = bodyMap.at( parameterToEstimate->getParameterName( ).second.first );

    // Check for which rotation model parameter the partial object is to be created.
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::constant_rotation_rate:

        if( boost::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >(
                    currentBody->getRotationalEphemeris( ) ) == NULL )
        {
            throw std::runtime_error( "Warning, body's rotation model is not simple when making position w.r.t. constant rtoation rate partial" ) ;
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = boost::make_shared< RotationMatrixPartialWrtConstantRotationRate >(
                    boost::dynamic_pointer_cast< SimpleRotationalEphemeris>( currentBody->getRotationalEphemeris( ) ) );
        break;
    default:
        std::string errorMessage = "Warning, rotation matrix partial not implemented for parameter " +
                boost::lexical_cast< std::string >( parameterToEstimate->getParameterName( ).first );
        throw std::runtime_error( errorMessage );

        break;
    }

    return rotationMatrixPartial;


}

//! Function to create partial object(s) of rotation matrix wrt a (vector) parameter.
boost::shared_ptr< RotationMatrixPartial > createRotationMatrixPartialsWrtParameter(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate )

{
    using namespace simulation_setup;
    using namespace ephemerides;

    // Declare return object.
    boost::shared_ptr< RotationMatrixPartial >  rotationMatrixPartial;

    // Get body for rotation of which partial is to be created.
    boost::shared_ptr< Body > currentBody = bodyMap.at( parameterToEstimate->getParameterName( ).second.first );

    // Check for which rotation model parameter the partial object is to be created.
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::rotation_pole_position:


        if( boost::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >(
                    currentBody->getRotationalEphemeris( ) ) == NULL )
        {
            std::string errorMessage = "Warning, body's rotation model is not simple when making position w.r.t. pole position partial";
            throw std::runtime_error( errorMessage );
        }

        // Create rotation matrix partial object
        rotationMatrixPartial = boost::make_shared< RotationMatrixPartialWrtPoleOrientation >(
                    boost::dynamic_pointer_cast< SimpleRotationalEphemeris>( currentBody->getRotationalEphemeris( ) ) );
        break;

    default:
        std::string errorMessage = "Warning, rotation matrix partial not implemented for parameter " +
                boost::lexical_cast< std::string >( parameterToEstimate->getParameterName( ).first );
        throw std::runtime_error( errorMessage );
        break;
    }

    return rotationMatrixPartial;

}


}

}
