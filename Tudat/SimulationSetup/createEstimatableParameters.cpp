#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/SimulationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace simulation_setup
{


using namespace simulation_setup;
using namespace ephemerides;
using namespace gravitation;
using namespace estimatable_parameters;

boost::shared_ptr< EstimatableParameter< double > > createDoubleParameterToEstimate(
        const boost::shared_ptr< EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
{
    boost::shared_ptr< EstimatableParameter< double > > doubleParameterToEstimate;

    if( isDoubleParameter( doubleParameterName->parameterType_.first ) != true )
    {
        std::string errorMessage = "Error when requesting to make double parameter " +
                boost::lexical_cast< std::string >( doubleParameterName->parameterType_.first ) + " of " +
                boost::lexical_cast< std::string >( doubleParameterName->parameterType_.second.first ) +
                ", parameter is not a double parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        std::string currentBodyName = doubleParameterName->parameterType_.second.first;
        boost::shared_ptr< Body > currentBody;

        if( ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Error when creating parameters to estimate, body " +
                    boost::lexical_cast< std::string >( currentBodyName ) +
                    "  not in body map " +
                    boost::lexical_cast< std::string >( doubleParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( currentBodyName != "" )
        {
            currentBody = bodyMap.at( currentBodyName );
        }

        switch( doubleParameterName->parameterType_.first )
        {
        case gravitational_parameter:
        {
            if( currentBody->getGravityFieldModel( )== NULL )
            {
                std::string errorMessage = "Error, body " +
                        boost::lexical_cast< std::string >( currentBodyName ) +
                        " has no gravity field, cannot estimate gravitational parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                boost::shared_ptr< GravityFieldModel > gravityFieldModel = currentBody->getGravityFieldModel( );
                doubleParameterToEstimate = boost::make_shared< GravitationalParameter >
                        ( gravityFieldModel, currentBodyName );
            }
            break;
        }
        case radiation_pressure_coefficient:
        {
            if( currentBody->getRadiationPressureInterfaces( ).size( ) == 0 )
            {
                std::string errorMessage = "Error, no radiation pressure interfaces found in body " +
                        boost::lexical_cast< std::string >( currentBodyName) +
                        " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else if( currentBody->getRadiationPressureInterfaces( ).size( ) > 1 )
            {
                std::string errorMessage = "Error, multiple radiation pressure interfaces found in body " +
                        boost::lexical_cast< std::string >( currentBodyName) +
                        " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = boost::make_shared< RadiationPressureCoefficient >(
                            currentBody->getRadiationPressureInterfaces( ).begin( )->second,
                            currentBodyName );
            }
            break;
        }
        default:
            throw std::runtime_error( "Warning, this double parameter has not yet been implemented when making parameters" );
            break;
        }
    }

    return doubleParameterToEstimate;
}

boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const boost::shared_ptr< EstimatableParameterSettings >& vectorParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
{
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > vectorParameterToEstimate;

    if( isDoubleParameter( vectorParameterName->parameterType_.first ) != false )
    {
        std::string errorMessage = "Error when requesting to make vector parameter " +
                boost::lexical_cast< std::string >( vectorParameterName->parameterType_.first ) +
                " of  " + boost::lexical_cast< std::string >( vectorParameterName->parameterType_.second.first ) +
                ", parameter is not a vector parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        std::string currentBodyName = vectorParameterName->parameterType_.second.first;
        boost::shared_ptr< Body > currentBody;

        if( ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Warning when creating parameters to estimate, body " +
                    boost::lexical_cast< std::string >( currentBodyName ) +
                    "not in body map " +
                    boost::lexical_cast< std::string >( vectorParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( ( currentBodyName != "" ) )
        {
            currentBody = bodyMap.at( currentBodyName );
        }

        switch( vectorParameterName->parameterType_.first )
        {
        default:
            std::string errorMessage = "Warning, this vector parameter (" +
                    boost::lexical_cast< std::string >( vectorParameterName->parameterType_.first ) +
                    ") has not yet been implemented when making parameters";
            throw std::runtime_error( errorMessage );

            break;
        }
    }

    return vectorParameterToEstimate;
}

}

}
