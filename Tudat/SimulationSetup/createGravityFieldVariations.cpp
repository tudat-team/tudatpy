#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/tabulatedGravityFieldVariations.h"
#include "Tudat/SimulationSetup/createGravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< TabulatedGravityFieldVariationSettings > createTabulatedGravityFieldVariationSettings(
        const boost::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVariations,
        const double startTime,
        const double endTime,
        const double interpolationStep,
        const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections;
    std::map< double, Eigen::MatrixXd > sineCoefficientCorrections;
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > currentCorrections;
    
    double currentTime = startTime - interpolationStep;
    
    while( currentTime < endTime )
    {
        currentTime += interpolationStep;
        currentCorrections = gravityFieldVariations->calculateSphericalHarmonicsCorrections( currentTime );
        
        cosineCoefficientCorrections[ currentTime ] = currentCorrections.first;
        sineCoefficientCorrections[ currentTime ] = currentCorrections.second;
        
    }
    
    return boost::make_shared< TabulatedGravityFieldVariationSettings >(
                cosineCoefficientCorrections, sineCoefficientCorrections,
                gravityFieldVariations->getMinimumDegree( ), gravityFieldVariations->getMinimumOrder( ), interpolatorSettings );
}

boost::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const boost::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const NamedBodyMap bodyMap )
{
    using namespace tudat::gravitation;
    
    boost::shared_ptr< GravityFieldVariations > gravityFieldVariationModel;
    
    switch( gravityFieldVariationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    {
        boost::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > basicSolidBodyGravityVariationSettings =
                boost::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( gravityFieldVariationSettings );
        if( basicSolidBodyGravityVariationSettings == NULL )
        {
            std::cerr<<"Error, expected basic solid body gravity field settings for "<<body<<std::endl;
        }
        else
        {
            std::vector< std::string > deformingBodies = basicSolidBodyGravityVariationSettings->getDeformingBodies( );
            boost::function< basic_mathematics::Vector6d( const double ) > deformedBodyStateFunction;
            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > deformingBodyStateFunctions;
            std::vector< boost::function< double( ) > > gravitionalParametersOfDeformingBodies;
            
            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                if( bodyMap.count( deformingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making basic solid body gravity field variation, deforming body not found: "<<deformingBodies[ i ]<<std::endl;
                }
                
                if( gravityFieldVariationSettings->getInterpolateVariation( ) )
                {
                    deformingBodyStateFunctions.push_back(
                                boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( deformingBodies[ i ] ), _1 ) );
                }
                else
                {
                    deformingBodyStateFunctions.push_back(
                                boost::bind( &Body::getState, bodyMap.at( deformingBodies[ i ] ) ) );
                }

                if( bodyMap.at( deformingBodies[ i ] )->getGravityFieldModel( ) == NULL )
                {
                    std::cerr<<"Error, could not find gravity field model in body "<<deformingBodies[ i ]<<" when making basic sh variation for body "<<body<<std::endl;
                }
                else
                {
                    gravitionalParametersOfDeformingBodies.push_back(
                                boost::bind( &GravityFieldModel::getGravitationalParameter,
                                             bodyMap.at( deformingBodies[ i ] )->getGravityFieldModel( ) ) );
                }
            }
            if( gravityFieldVariationSettings->getInterpolateVariation( ) )
            {
                deformedBodyStateFunction = boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( body ), _1 );
            }
            else
            {
                deformedBodyStateFunction = boost::bind( &Body::getState, bodyMap.at( body ) );
                
            }
            boost::function< double( ) > gravitionalParameterOfDeformedBody =
                    boost::bind( &GravityFieldModel::getGravitationalParameter, bodyMap.at( body )->getGravityFieldModel( ) );
            
            gravityFieldVariationModel = boost::make_shared< BasicSolidBodyTideGravityFieldVariations >(
                        deformedBodyStateFunction,
                        boost::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                                     bodyMap.at( body )->getRotationalEphemeris( ), _1, basic_astrodynamics::JULIAN_DAY_ON_J2000  ),
                        deformingBodyStateFunctions,
                        basicSolidBodyGravityVariationSettings->getBodyReferenceRadius( ),
                        gravitionalParameterOfDeformedBody,
                        gravitionalParametersOfDeformingBodies,
                        basicSolidBodyGravityVariationSettings->getLoveNumbers( ),
                        deformingBodies );
        }
        break;
    }    
    case tabulated_variation:
    {
        boost::shared_ptr< TabulatedGravityFieldVariationSettings > tabulatedGravityFieldVariationSettings =
                boost::dynamic_pointer_cast< TabulatedGravityFieldVariationSettings >( gravityFieldVariationSettings );
        if( tabulatedGravityFieldVariationSettings == NULL )
        {
            std::cerr<<"Error, expected tabulated gravity field variation settings for "<<body<<std::endl;
        }
        else
        {
            gravityFieldVariationModel = boost::make_shared< TabulatedGravityFieldVariations >
                    (  tabulatedGravityFieldVariationSettings->getCosineCoefficientCorrections( ),
                       tabulatedGravityFieldVariationSettings->getSineCoefficientCorrections( ),
                       tabulatedGravityFieldVariationSettings->getMinimumDegree( ),
                       tabulatedGravityFieldVariationSettings->getMinimumOrder( ),
                       tabulatedGravityFieldVariationSettings->getInterpolatorSettings( ) );
        }
        break;
    }
    default:
    {
        std::cerr<<"Error, this case, "<<gravityFieldVariationSettings->getBodyDeformationType( )<<" not implemented for gravity field variations"<<std::endl;
    }
        
    }
    
    if( gravityFieldVariationModel == NULL )
    {
        std::cerr<<"Model IS NULL"<<std::endl;
    }
    
    return gravityFieldVariationModel;
    
}

}

}
