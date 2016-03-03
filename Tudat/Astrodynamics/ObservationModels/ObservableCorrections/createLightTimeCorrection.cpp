
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/Relativity/Metrics/metric.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/secondOrderRelativisticLightTimeCorrection.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/massMultipoleLightTimeCorrection.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/troposphericDelayFunctions.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/ionosphericDelayFunctions.h"

namespace tudat
{

namespace observation_models
{

boost::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings >& correctionSettings,
        const NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver )
{

    using namespace tudat::ephemerides;
    using namespace tudat::bodies;
    using namespace tudat::gravitation;
    using namespace tudat::relativity;

    boost::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    switch( correctionSettings->getCorrectionType( ) )
    {
    case marini_murray_troposperic:
    {
        if( boost::dynamic_pointer_cast< TroposphericCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::string bodyWithAtmosphere =
                    boost::dynamic_pointer_cast< TroposphericCorrectionSettings >( correctionSettings )->getBodyWithAtmosphere( );

            boost::shared_ptr< OpticalTroposphericDelayFunction > delayFunction;

            if( bodyWithAtmosphere == transmitter.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<transmitter.first<<" is not a celestial body when making marini murray correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > transmitterCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) );
                    delayFunction = boost::make_shared< MariniMurrayDelayFunction >(
                                transmitterCelestialBody->getGroundStation( transmitter.second ) );
                    lightTimeCorrection = boost::make_shared< TroposphericDelayCorrectionInterface >(
                                delayFunction, boost::bind(
                                    &PointingAnglesCalculator::calculateElevationAngle,
                                    transmitterCelestialBody->getGroundStation(
                                        transmitter.second )->getPointingAnglesCalculator( ), _1, _2 ), 0 );
                }

            }
            else if( bodyWithAtmosphere == receiver.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<receiver.first<<" is not a celestial body when making marini murray correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > receiverCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) );
                    delayFunction = boost::make_shared< MariniMurrayDelayFunction >(
                                receiverCelestialBody->getGroundStation( receiver.second ) );
                    lightTimeCorrection = boost::make_shared< TroposphericDelayCorrectionInterface >(
                                delayFunction, boost::bind( &PointingAnglesCalculator::calculateElevationAngle,
                                                            receiverCelestialBody->getGroundStation(
                                                                receiver.second )->getPointingAnglesCalculator( ), _1, _2 ), 1 );
                }

            }
            else
            {
                std::cerr<<"Error, atmosphere that is to delay signal does not belong to transmitter or receiver when making marini murray correction"<<std::endl;
            }
        }
        else
        {
            std::cerr<<"Error, correction settings type (marini murray) does not coincide with data type."<<std::endl;
        }
        break;
    }

    case fcul_troposphericDelay:
    {
        if( boost::dynamic_pointer_cast< TroposphericCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::string bodyWithAtmosphere =
                    boost::dynamic_pointer_cast< TroposphericCorrectionSettings >( correctionSettings )->getBodyWithAtmosphere( );

            boost::shared_ptr< OpticalTroposphericDelayFunction > delayFunction;

            if( bodyWithAtmosphere == transmitter.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<transmitter.first<<" is not a celestial body when making fcul correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > transmitterCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) );
                    delayFunction = boost::make_shared< FculaDelayFunction >(
                                transmitterCelestialBody->getGroundStation( transmitter.second ) );
                    lightTimeCorrection = boost::make_shared< TroposphericDelayCorrectionInterface >(
                                delayFunction, boost::bind( &PointingAnglesCalculator::calculateElevationAngle,
                                                            transmitterCelestialBody->getGroundStation(
                                                                transmitter.second )->getPointingAnglesCalculator( ), _1, _2 ), 0 );
                }

            }
            else if( bodyWithAtmosphere == receiver.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<receiver.first<<" is not a celestial body when making fcul correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > receiverCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) );
                    delayFunction = boost::make_shared< FculaDelayFunction >(
                                receiverCelestialBody->getGroundStation( receiver.second ) );
                    lightTimeCorrection = boost::make_shared< TroposphericDelayCorrectionInterface >(
                                delayFunction, boost::bind( &PointingAnglesCalculator::calculateElevationAngle,
                                                            receiverCelestialBody->getGroundStation(
                                                                receiver.second )->getPointingAnglesCalculator( ), _1, _2 ), 1 );
                }

            }
            else
            {
                std::cerr<<"Error, atmosphere that is to delay signal does not belong to transmitter or receiver when making fcul correction"<<std::endl;
            }
        }
        else
        {
            std::cerr<<"Error, correction settings type (fcul) does not coincide with data type."<<std::endl;
        }
        break;
    }
    case ionex_ionosphericDelay:
    {
        if( boost::dynamic_pointer_cast< IonosphericCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::string bodyWithAtmosphere =
                    boost::dynamic_pointer_cast< IonosphericCorrectionSettings >( correctionSettings )->getBodyWithAtmosphere( );

            boost::shared_ptr< IonosphericDelayFunction > delayFunction;

            if( bodyWithAtmosphere == transmitter.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<transmitter.first<<" is not a celestial body when making ionex correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > transmitterCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( transmitter.first ) );

                    boost::shared_ptr< OblateSpheroidBodyShapeModel > bodyShape =
                            boost::dynamic_pointer_cast< OblateSpheroidBodyShapeModel >( transmitterCelestialBody->getShapeModel( ) );

                    delayFunction = boost::make_shared< IonexDelayFunction >(
                                transmitterCelestialBody->getGroundStation( transmitter.second ),
                                bodyShape );
                    lightTimeCorrection = boost::make_shared< IonosphericDelayCorrectionInterface >(
                                delayFunction,
                                boost::bind( &PointingAnglesCalculator::calculateElevationAngle,
                                             transmitterCelestialBody->getGroundStation( transmitter.second )->getPointingAnglesCalculator( ),
                                             _1, _2 ),
                                boost::bind( &PointingAnglesCalculator::calculationAzimuthAngle,
                                             transmitterCelestialBody->getGroundStation( transmitter.second )->getPointingAnglesCalculator( ),
                                             _1, _2 ),
                                0 );
                }

            }
            else if( bodyWithAtmosphere == receiver.first )
            {
                if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) ) == NULL )
                {
                    std::cerr<<"Error, body "<<receiver.first<<" is not a celestial body when making ionex correction"<<std::endl;
                }
                else
                {
                    boost::shared_ptr< CelestialBody > receiverCelestialBody =
                            boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( receiver.first ) );
                    boost::shared_ptr< OblateSpheroidBodyShapeModel > bodyShape =
                            boost::dynamic_pointer_cast< OblateSpheroidBodyShapeModel >( receiverCelestialBody->getShapeModel( ) );

                    delayFunction = boost::make_shared< IonexDelayFunction >(
                                receiverCelestialBody->getGroundStation( receiver.second ),
                                bodyShape );
                    lightTimeCorrection = boost::make_shared< IonosphericDelayCorrectionInterface >(
                                delayFunction,
                                boost::bind( &PointingAnglesCalculator::calculateElevationAngle,
                                             receiverCelestialBody->getGroundStation( receiver.second )->getPointingAnglesCalculator( ),
                                             _1, _2 ),
                                boost::bind( &PointingAnglesCalculator::calculationAzimuthAngle,
                                             receiverCelestialBody->getGroundStation( receiver.second )->getPointingAnglesCalculator( ),
                                             _1, _2 ),
                                1 );
                }

            }
            else
            {
                std::cerr<<"Error, atmosphere that is to delay signal does not belong to transmitter or receiver when making ionex correction"<<std::endl;
            }
        }
        else
        {
            std::cerr<<"Error, correction settings type (ionex) does not coincide with data type."<<std::endl;
        }
        break;
    }

    case first_order_relativistic:
    {
        if( boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making 1st order relativistic light time correction, could not find body "<<perturbingBodies[ i ]<<std::endl;
                }
                else
                {
                    perturbingBodyStateFunctions.push_back( boost::bind( &Body::getStateInBaseFrameFromEphemeris,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );
                    if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                            getGravityFieldModel( ) == NULL )
                    {
                        std::cerr<<"Error when making 1st order time correction, body "<<perturbingBodies[ i ]<<" has no gravity field "<<std::endl;

                    }
                    else
                    {
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                 boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                                                 getGravityFieldModel( ) ) );
                    }
                }
            }

            lightTimeCorrection = boost::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ) );

        }
        else
        {
            std::cerr<<"Error, correction settings type (1st order relativistic) does not coincide with data type."<<std::endl;
        }

        break;
    }

    case second_order_relativistic:
    {
        if( boost::dynamic_pointer_cast< SecondOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< SecondOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making 2nd order relativistic light time correction, could not find body "<<perturbingBodies[ i ]<<std::endl;
                }
                else
                {
                    perturbingBodyStateFunctions.push_back( boost::bind( &Body::getStateInBaseFrameFromEphemeris,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );
                    if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                            getGravityFieldModel( ) == NULL )
                    {
                        std::cerr<<"Error when making 2nd order time correction, body "<<perturbingBodies[ i ]<<" has no gravity field "<<std::endl;

                    }
                    else
                    {
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                 boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                                                 getGravityFieldModel( ) ) );
                    }
                }
            }

            lightTimeCorrection = boost::make_shared< SecondOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ),
                        boost::bind( &PPNParameterSet::getParameterBeta, ppnParameterSet ),
                        boost::bind( &PPNParameterSet::getParameterDelta, ppnParameterSet ) );

        }
        else
        {
            std::cerr<<"Error, correction settings type (2nd order relativistic) does not coincide with data type."<<std::endl;
        }

        break;
    }
    case j2_relativistic:
    {
        if( boost::dynamic_pointer_cast< J2RelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< J2RelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;
            std::vector< boost::function< Eigen::Vector3d( const double ) > > perturbingBodyRotationAxisFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyJ2Coefficients;
            std::vector< boost::function< double( ) > > perturbingBodyReferenceRadii;

            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making j2 light time correction, could not find body "<<perturbingBodies[ i ]<<std::endl;
                }
                else
                {
                    perturbingBodyStateFunctions.push_back( boost::bind( &Body::getStateInBaseFrameFromEphemeris,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );



                    if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                            getGravityFieldModel( ) == NULL )
                    {
                        std::cerr<<"Error when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has no gravity field "<<std::endl;

                    }
                    else
                    {
                        boost::shared_ptr< gravitation::GravityFieldModel > bodyGravityField =
                                boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                                getGravityFieldModel( );
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter, bodyGravityField ) );
                        if( boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodyGravityField ) == NULL )
                        {
                            std::cerr<<"Error when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has no spherical harmonic gravity field "<<std::endl;
                        }
                        else
                        {
                            if( boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodyGravityField ) != NULL )
                            {
                                std::cerr<<"Warning when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has "<<
                                           "a time-dependent gravity field; will only use static part for light time correction"<<std::endl;
                            }
                            boost::shared_ptr< SphericalHarmonicsGravityField > bodySphericalHarmonicField =
                                    boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodyGravityField );
                            perturbingBodyReferenceRadii.push_back(
                                        boost::bind( &SphericalHarmonicsGravityField::getReferenceRadius, bodySphericalHarmonicField ) );
                            perturbingBodyJ2Coefficients.push_back(
                                        boost::bind( &SphericalHarmonicsGravityField::getUnnormalizedCosineCoefficient, bodySphericalHarmonicField,
                                                     2, 0 ) );

                            if( bodyMap.at( perturbingBodies[ i ] )->getRotationalEphemeris( ) == NULL )
                            {
                                std::cerr<<"Error when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has no rotation model"<<std::endl;

                            }
                            else
                            {
                                boost::shared_ptr< RotationalEphemeris > bodyRotationModel = bodyMap.at(
                                            perturbingBodies[ i ] )->getRotationalEphemeris( );
                                perturbingBodyRotationAxisFunctions.push_back(
                                            boost::bind( &RotationalEphemeris::getRotationalVelocityVectorInBaseFrame, bodyRotationModel, _1 ) );

                                if( bodyRotationModel->getBaseFrameOrientation( ) != "ECLIPJ2000" )
                                {
                                    std::cerr<<"Warning when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has "<<
                                               "rotation mdodel with "<<bodyRotationModel->getTargetFrameOrientation( )<<" base frame"<<std::endl;
                                }

                                if( bodyRotationModel->getTargetFrameOrientation( ) != bodySphericalHarmonicField->getFixedReferenceFrame( ) )
                                {
                                    std::cerr<<"Warning when making j2 light time correction, body "<<perturbingBodies[ i ]<<" has "<<
                                               "rotation mdodel with "<<bodySphericalHarmonicField->getFixedReferenceFrame( )<<" as fixed frame, "<<
                                               "but rotation model has "<<bodyRotationModel->getTargetFrameOrientation( )<<" as target frame "<<std::endl;
                                }

                            }
                        }
                    }
                }

            }

            lightTimeCorrection = boost::make_shared< J2RelativisticLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodyRotationAxisFunctions,
                        perturbingBodyJ2Coefficients, perturbingBodyJ2Coefficients, perturbingBodies,
                        boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ) );

        }
        else
        {
            std::cerr<<"Error, correction settings type (2nd order relativistic) does not coincide with data type."<<std::endl;
        }
        break;
    }
    default:
    {
        std::cerr<<"Error, light time correction type not recognized."<<std::endl;
        break;
    }

    }
    return lightTimeCorrection;
}

boost::shared_ptr< LightTimeDerivativeCorrection > createLightTimeDerivativeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver )
{
    using namespace tudat::ephemerides;
    using namespace tudat::bodies;
    using namespace tudat::gravitation;
    using namespace tudat::relativity;

    boost::shared_ptr< LightTimeDerivativeCorrection > lightTimeCorrectionDerivative;

    switch( correctionSettings->getCorrectionType( ) )
    {
    case first_order_relativistic:
    {
        if( boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making 1st order relativistic light time correction, could not find body "<<perturbingBodies[ i ]<<std::endl;
                }
                else
                {
                    perturbingBodyStateFunctions.push_back( boost::bind( &Body::getStateInBaseFrameFromEphemeris,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );
                    if( boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                            getGravityFieldModel( ) == NULL )
                    {
                        std::cerr<<"Error when making 1st order time correction, body "<<perturbingBodies[ i ]<<" has no gravity field "<<std::endl;

                    }
                    else
                    {
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                 boost::dynamic_pointer_cast< CelestialBody >( bodyMap.at( perturbingBodies[ i ] ) )->
                                                 getGravityFieldModel( ) ) );
                    }
                }
            }

            lightTimeCorrectionDerivative = boost::make_shared< FirstOrderLightTimeCorrectionDerivativeCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        boost::bind( &PPNParameterSet::getParameterGamma, ppnParameterSet ) );

        }
        else
        {
            std::cerr<<"Error, correction settings type (1st order relativistic) does not coincide with data type."<<std::endl;
        }

        break;
    }
    default:
        std::cerr<<"Error, cannot create light time correction derivative for type "<<correctionSettings->getCorrectionType( )<<std::endl;
    }

    return lightTimeCorrectionDerivative;
}

}

}
