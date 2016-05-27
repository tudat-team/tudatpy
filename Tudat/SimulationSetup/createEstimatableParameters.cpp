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
        std::cerr<<"Error when requesting to make double parameter "<<doubleParameterName->parameterType_.first<<" of "<<
                   doubleParameterName->parameterType_.second.first<<", parameter is not a double parameter "<<std::endl;
    }
    else
    {
        std::string currentBodyName = doubleParameterName->parameterType_.second.first;
        boost::shared_ptr< Body > currentBody;

        if( ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::cerr<<"Warning when creating parameters to estimate, body "<<currentBodyName;
            std::cerr<<" not in body map "<<doubleParameterName->parameterType_.first<<std::endl;
            currentBody = boost::shared_ptr< Body >( );
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
                std::cerr<<"Warning, body "<<currentBodyName<<" has no gravity field";
                std::cerr<<", cannot estimate gravitational parameter."<<std::endl;
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
                std::cerr<<"Error, no radiation pressure interfaces found in body "<<currentBodyName<<" when making Cr parameter"<<std::endl;
            }
            else if( currentBody->getRadiationPressureInterfaces( ).size( ) > 1 )
            {
                std::cerr<<"Error, multiple radiation pressure interfaces found in body "<<currentBodyName<<" when making Cr parameter"<<std::endl;
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
            std::cerr<<"Warning, this double parameter has not yet been implemented when making parameters"<<std::endl;
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
        std::cerr<<"Error when requesting to make vector parameter "<<vectorParameterName->parameterType_.first<<" of "<<
                   vectorParameterName->parameterType_.second.first<<", parameter is not a vector parameter "<<std::endl;
    }
    else
    {
        std::string currentBodyName = vectorParameterName->parameterType_.second.first;
        boost::shared_ptr< Body > currentBody;

        if( ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::cerr<<"Warning when creating parameters to estimate, body "<<currentBodyName;
            std::cerr<<" not in body map "<<vectorParameterName->parameterType_.first<<std::endl;
        }
        else if( ( currentBodyName != "" ) )
        {
            currentBody = bodyMap.at( currentBodyName );
        }

        switch( vectorParameterName->parameterType_.first )
        {
        default:
            std::cerr<<"Warning, this vector parameter ("<<vectorParameterName->parameterType_.first<<") has not yet been implemented when making parameters"<<std::endl;
            break;
        }
    }

    return vectorParameterToEstimate;
}

}

}
