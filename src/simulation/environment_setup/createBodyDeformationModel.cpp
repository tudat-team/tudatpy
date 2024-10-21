#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{

namespace simulation_setup
{


std::shared_ptr< basic_astrodynamics::Iers2010EarthDeformation > createDefaultEarthIers2010DeformationModel(
        const std::shared_ptr< ephemerides::Ephemeris > earthEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > lunarEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > solarEphemeris,
        const std::shared_ptr< ephemerides::RotationalEphemeris > earthRotation,
        const std::function< double( ) > gravitionalParametersOfEarth,
        const std::function< double( ) > gravitionalParametersOfMoon,
        const std::function< double( ) > gravitionalParametersOfSun,
        const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction )
{
    using namespace tudat::ephemerides;
    using namespace tudat::basic_astrodynamics;

    std::vector< std::function< double( ) > > gravitationalParameters;
    gravitationalParameters.resize( 2 );
    gravitationalParameters[ 0 ] = gravitionalParametersOfMoon;
    gravitationalParameters[ 1 ] = gravitionalParametersOfSun;

    std::vector< std::function< Eigen::Vector6d( const double ) > > ephemerides;
    ephemerides.resize( 2 );
    ephemerides[ 0 ] = std::bind( &Ephemeris::getCartesianState, lunarEphemeris, std::placeholders::_1 );
    ephemerides[ 1 ] = std::bind( &Ephemeris::getCartesianState, solarEphemeris, std::placeholders::_1 );

    double equatorialRadius = 6378136.6;

    using namespace iers_2010_parameters;

    std::map< int, std::pair< double, double > > nominalDisplacementLoveNumbers;
    nominalDisplacementLoveNumbers[ 2 ] = std::make_pair( PRINCIPAL_DEGREE_TWO_LOVE_NUMBER, PRINCIPAL_DEGREE_TWO_SHIDA_NUMBER );
    nominalDisplacementLoveNumbers[ 3 ] = std::make_pair( PRINCIPAL_DEGREE_THREE_LOVE_NUMBER, PRINCIPAL_DEGREE_THREE_SHIDA_NUMBER );

    std::vector< double > latitudeTerms;
    latitudeTerms.resize( 2 );
    latitudeTerms[ 0 ] = DEGREE_TWO_LATITUDE_LOVE_NUMBER;
    latitudeTerms[ 1 ] = DEGREE_TWO_LATITUDE_SHIDA_NUMBER;

    std::vector< bool > areTermsCalculated;
    areTermsCalculated.resize( 6 );
    for( int i = 0; i < 6 ; i++ )
    {
        areTermsCalculated[ i ] = 1;
    }

    std::vector< double > correctionNumbers;
    correctionNumbers.resize( 6 );
    correctionNumbers[ 2 ] = IMAGINARY_DEGREE_TWO_DIURNAL_LOVE_NUMBER;
    correctionNumbers[ 3 ] = IMAGINARY_DEGREE_TWO_DIURNAL_SHIDA_NUMBER;
    correctionNumbers[ 4 ] = IMAGINARY_DEGREE_TWO_SEMIDIURNAL_LOVE_NUMBER;
    correctionNumbers[ 5 ] = IMAGINARY_DEGREE_TWO_SEMIDIURNAL_SHIDA_NUMBER;
    correctionNumbers[ 0 ] = DEGREE_TWO_DIURNAL_TOROIDAL_LOVE_NUMBER;
    correctionNumbers[ 1 ] = DEGREE_TWO_SEMIDIURNAL_TOROIDAL_LOVE_NUMBER;

    std::string longPeriodFile = paths::getEarthDeformationDataFilesPath( ) +
            "/longPeriodDisplacementFrequencyDependence.txt";
    std::string diurnalFile = paths::getEarthDeformationDataFilesPath( ) +
            "/diurnalDisplacementFrequencyDependence2.txt";

    std::shared_ptr< Iers2010EarthDeformation > deformationModel = std::make_shared< Iers2010EarthDeformation >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              ephemerides,
              std::bind( &RotationalEphemeris::getRotationToTargetFrame, earthRotation, std::placeholders::_1 ),
              gravitionalParametersOfEarth, gravitationalParameters,
              equatorialRadius, nominalDisplacementLoveNumbers, latitudeTerms, areTermsCalculated,
              correctionNumbers, diurnalFile, longPeriodFile, doodsonArgumentFunction );

    return deformationModel;

}

std::shared_ptr< basic_astrodynamics::BodyDeformationModel > createBodyDeformationModel(
        const std::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const SystemOfBodies& bodyMap )
{
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    std::shared_ptr< BodyDeformationModel > bodyDeformationModel;

    switch( bodyDeformationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    {
        std::shared_ptr< BasicSolidBodyDeformationSettings > basicSolidBodyDeformationSettings =
                std::dynamic_pointer_cast< BasicSolidBodyDeformationSettings >( bodyDeformationSettings );

        if( basicSolidBodyDeformationSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body deformation model, expected basic solid body settings for " + body );
        }
        else
        {
            std::vector< std::string > deformingBodies = basicSolidBodyDeformationSettings->getDeformingBodies( );
            std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyEphemerides;
            std::vector< std::function< double( ) > > gravitionalParametersOfDeformingBodies;

            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                if( bodyMap.count( deformingBodies.at( i ) ) == 0 )
                {
                    throw std::runtime_error( "Error when making basic solid body deformation model, deforming body not found: " +
                                              deformingBodies.at( i ) );
                }

                deformingBodyEphemerides.push_back(
                            std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodyMap.at( deformingBodies.at( i ) ), std::placeholders::_1 ) );

                std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel =
                        bodyMap.at( deformingBodies.at( i ) )->getGravityFieldModel( );
                if( gravityFieldModel == nullptr )
                {
                    throw std::runtime_error( "Error, no gravity field model of " + deformingBodies.at( i ) +
                                              " found when making basic body deformation of " + body );
                }
                gravitionalParametersOfDeformingBodies.push_back(
                            std::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel ) );
            }

            std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel = bodyMap.at( body )->getGravityFieldModel( );
            if( gravityFieldModel == NULL )
            {
                throw std::runtime_error("Error, no gravity field model of " + body +
                                         " found when making basic body deformation of " + body );
            }

            std::function< double( ) > gravitionalParameterOfDeformedBody =
                    std::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel );


            double deformationReferenceRadius = basicSolidBodyDeformationSettings->getBodyReferenceRadius( );
            if( deformationReferenceRadius != deformationReferenceRadius )
            {
                std::shared_ptr< BodyShapeModel > bodyShapeModel =
                        bodyMap.at( body )->getShapeModel( );
                if( bodyShapeModel == nullptr )
                {
                    throw std::runtime_error("Error, when making basic body deformation of " + body + ", no reference radius, and no shape model specified" );
                }
                else
                {
                    deformationReferenceRadius = bodyShapeModel->getAverageRadius( );
                }
            }

            bodyDeformationModel = std::make_shared< basic_astrodynamics::BasicTidalBodyDeformation >(
                        std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodyMap.at( body ),
                                   std::placeholders::_1 ),
                        deformingBodyEphemerides,
                        std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame, bodyMap.at( body )->getRotationalEphemeris( ),
                                   std::placeholders::_1 ),
                        gravitionalParameterOfDeformedBody,
                        gravitionalParametersOfDeformingBodies,
                        deformationReferenceRadius,
                        basicSolidBodyDeformationSettings->getDisplacementLoveNumbers( ) );

        }
        break;
    }
    case iers_2010:
    {
        if( body != "Earth" )
        {
            std::cerr<<"Warning, you should use IERS 2010 deformation model only for Earth, now adding to body"<<std::endl;
        }
        if( bodyMap.count( "Earth" ) == 0 )
        {
            throw std::runtime_error( "Error when making IERS 2010 deformation model, Earth not found" );
        }
        if( bodyMap.count( "Moon" ) == 0 )
        {
            throw std::runtime_error( "Error when making IERS 2010 deformation model, Moon not found" );
        }
        if( bodyMap.count( "Sun" ) == 0 )
        {
            throw std::runtime_error( "Error when making IERS 2010 deformation model, Sun not found" );
        }

        bodyDeformationModel = createDefaultEarthIers2010DeformationModel(
                    bodyMap.at( "Earth" )->getEphemeris( ),
                    bodyMap.at( "Moon" )->getEphemeris( ),
                    bodyMap.at( "Sun" )->getEphemeris( ),
                    bodyMap.at( "Earth" )->getRotationalEphemeris( ),
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, bodyMap.at( "Earth" )->getGravityFieldModel( ) ),
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, bodyMap.at( "Moon" )->getGravityFieldModel( ) ),
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, bodyMap.at( "Sun" )->getGravityFieldModel( ) ) );
        break;
    }
    case pole_tide:
    {
        if( body != "Earth" )
        {
            std::cerr<<"Warning, you should use pole tide deformation model only for Earth, now adding to body"<<std::endl;
        }
        if( bodyMap.count( body ) == 0 )
        {
            throw std::runtime_error( "Error when making pole tide deformation model, " + body  + " not found" );
        }
        else if( std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( bodyMap.at( body )->getRotationalEphemeris( ))  == nullptr )
        {
             throw std::runtime_error( "Error when making pole tide deformation model, " + body  + ", model requires GCRS<->ITRS rotation model" );
        }
        else
        {
            bodyDeformationModel = std::make_shared< basic_astrodynamics::PoleTideDeformation >(
                std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( bodyMap.at( body )->getRotationalEphemeris( ) )
                    ->getAnglesCalculator( )->getPolarMotionCalculator( ) );
        }
        break;
    }
    case ocean_tide:
    {
        std::shared_ptr< OceanTideBodyDeformationSettings > oceanTideBodyDeformationSettings =
            std::dynamic_pointer_cast< OceanTideBodyDeformationSettings >( bodyDeformationSettings );

        if( oceanTideBodyDeformationSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body deformation model, expected ocean tide body settings for " + body );
        }

        if( body != "Earth" )
        {
            std::cerr<<"Warning, you should use ocean tide deformation model only for Earth, now adding to body"<<std::endl;
        }

        if( bodyMap.count( body ) == 0 )
        {
            throw std::runtime_error( "Error when making pole tide deformation model, " + body  + " not found" );
        }
        else
        {
            bodyDeformationModel = std::make_shared< basic_astrodynamics::OceanTideEarthDeformation >(
                oceanTideBodyDeformationSettings->getBlqFiles( ) );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize body deformation settings type " +
                                  std::to_string( bodyDeformationSettings->getBodyDeformationType( ) ) +
                                  " of body " + body );
    }

    return bodyDeformationModel;
}

}

}
