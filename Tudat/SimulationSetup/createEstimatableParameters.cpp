#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicSineCoefficients.h"
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

//! Function to create an interface object for estimating a parameter defined by a single double value
boost::shared_ptr< EstimatableParameter< double > > createDoubleParameterToEstimate(
        const boost::shared_ptr< EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
{
    boost::shared_ptr< EstimatableParameter< double > > doubleParameterToEstimate;

    // Check input consistency.
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
        // Check if body associated with parameter exists.
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

        // Identify parameter type.
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
        case constant_rotation_rate:
        {
            if( boost::dynamic_pointer_cast< SimpleRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == NULL )
            {
                std::string errorMessage = "Warning, no simple rotational ephemeris present in body " + currentBodyName +
                        " when making constant rotation rate parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = boost::make_shared< RotationRate >(
                            boost::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >
                            ( currentBody->getRotationalEphemeris( ) ), currentBodyName );
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

    // Check input consistency.
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
        // Check if body associated with parameter exists.
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

        // Identify parameter type.
        switch( vectorParameterName->parameterType_.first )
        {
        case rotation_pole_position:
            if( boost::dynamic_pointer_cast< SimpleRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == NULL )
            {
                std::string errorMessage = "Warning, no simple rotational ephemeris present in body " + currentBodyName +
                        " when making constant rotation orientation parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                vectorParameterToEstimate = boost::make_shared< ConstantRotationalOrientation >
                        ( boost::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >
                          ( currentBody->getRotationalEphemeris( ) ), currentBodyName );

            }
            break;
        case spherical_harmonics_cosine_coefficient_block:
        {
            boost::shared_ptr< GravityFieldModel > gravityField = currentBody->getGravityFieldModel( );
            boost::shared_ptr< SphericalHarmonicsGravityField > shGravityField =
                    boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >( gravityField );
            if( shGravityField == NULL )
            {
                std::string errorMessage = "Error, requested spherical harmonic cosine coefficient block parameter of " +
                        boost::lexical_cast< std::string >( vectorParameterName->parameterType_.second.first ) +
                        ", but body does not have a spherical harmonic gravity field.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Check if spherical harmonic gravity field is static or time-dependent; set associated
                // functions accordingly
                boost::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentShField =
                        boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( shGravityField );

                boost::function< Eigen::MatrixXd( ) > getCosineCoefficientsFunction;
                boost::function< void( Eigen::MatrixXd ) > setCosineCoefficientsFunction;

                if( timeDependentShField == NULL )
                {
                    getCosineCoefficientsFunction = boost::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                                                 shGravityField );
                    setCosineCoefficientsFunction = boost::bind( &SphericalHarmonicsGravityField::setCosineCoefficients,
                                                                 shGravityField, _1 );
                }
                else
                {
                    getCosineCoefficientsFunction = boost::bind(
                                &TimeDependentSphericalHarmonicsGravityField::getNominalCosineCoefficients,
                                                                 timeDependentShField );
                    setCosineCoefficientsFunction = boost::bind(
                                &TimeDependentSphericalHarmonicsGravityField::setNominalCosineCoefficients,
                                                                 timeDependentShField, _1 );
                }

                // Create cosine coefficients estimation object.
                boost::shared_ptr< SphericalHarmonicEstimatableParameterSettings > blockParameterSettings =
                        boost::dynamic_pointer_cast< SphericalHarmonicEstimatableParameterSettings >( vectorParameterName );
                if( blockParameterSettings != NULL )
                {
                    vectorParameterToEstimate = boost::make_shared< SphericalHarmonicsCosineCoefficients >(
                                getCosineCoefficientsFunction,
                                setCosineCoefficientsFunction,
                                blockParameterSettings->blockIndices_,
                                vectorParameterName->parameterType_.second.first );
                }
                else
                {
                    throw std::runtime_error( "Error, expected SphericalHarmonicEstimatableParameterSettings for cosine coefficients" );
                }
            }
            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            boost::shared_ptr< GravityFieldModel > gravityField = currentBody->getGravityFieldModel( );
            boost::shared_ptr< SphericalHarmonicsGravityField > shGravityField =
                    boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >( gravityField );
            if( shGravityField == NULL )
            {
                std::string errorMessage = "Error, requested spherical harmonic sine coefficient block parameter of " +
                        boost::lexical_cast< std::string >( vectorParameterName->parameterType_.second.first ) +
                        ", but body does not have a spherical harmonic gravity field.";
                throw std::runtime_error( errorMessage );

            }
            else
            {
                boost::shared_ptr< SphericalHarmonicEstimatableParameterSettings > blockParameterSettings =
                        boost::dynamic_pointer_cast< SphericalHarmonicEstimatableParameterSettings >( vectorParameterName );

                // Check if spherical harmonic gravity field is static or time-dependent; set associated
                // functions accordingly
                boost::function< Eigen::MatrixXd( ) > getSineCoefficientsFunction;
                boost::function< void( Eigen::MatrixXd ) > setSineCoefficientsFunction;
                boost::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentShField =
                        boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( shGravityField );

                if( timeDependentShField == NULL )
                {
                    getSineCoefficientsFunction = boost::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                                               shGravityField );
                    setSineCoefficientsFunction = boost::bind( &SphericalHarmonicsGravityField::setSineCoefficients,
                                                               shGravityField, _1 );
                }
                else
                {
                    getSineCoefficientsFunction = boost::bind(
                                &TimeDependentSphericalHarmonicsGravityField::getNominalSineCoefficients,
                                                               timeDependentShField );
                    setSineCoefficientsFunction = boost::bind(
                                &TimeDependentSphericalHarmonicsGravityField::setNominalSineCoefficients,
                                                               timeDependentShField, _1 );
                }

                // Create sine coefficients estimation object.
                if( blockParameterSettings != NULL )
                {
                    vectorParameterToEstimate = boost::make_shared< SphericalHarmonicsSineCoefficients >(
                                getSineCoefficientsFunction,
                                setSineCoefficientsFunction,
                                blockParameterSettings->blockIndices_,
                                vectorParameterName->parameterType_.second.first );
                }
                else
                {
                    throw std::runtime_error( "Error, expected SphericalHarmonicEstimatableParameterSettings for sine coefficients" );
                }
            }

            break;
        }
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
