#include "tudat/astro/orbit_determination/metric_partials/solarSystemMetricPartials.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <Eigen/Geometry>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

std::vector< Eigen::Matrix< double, 4, 4 > > calculateSolarSystemMetricPartialFromPotentialPartials(
        const Eigen::VectorXd& scalarPotentialPartials,
        const Eigen::Matrix< double, 3, Eigen::Dynamic >& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet )
{
    std::vector< Eigen::Matrix< double, 4, 4 > > metricPartials;

    if( scalarPotentialPartials.rows( ) != vectorPotentialPartials.cols( ) )
    {
        std::cerr << "Error when calculating solar system metric partial, sizes are inconsistent" << std::endl;
        return metricPartials;
    }

    metricPartials.resize( scalarPotentialPartials.rows( ) );

    const double currentBeta = currentPpnParameterSet->getParameterBeta( );
    const double currentGamma = currentPpnParameterSet->getParameterGamma( );

    const double timeTimeMultiplier = 2.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT -
            4.0 * currentBeta * physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT * totalScalarPotential;
    const double spaceTimeMultiplier = -2.0 * ( currentGamma + 1.0 ) * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT;
    const Eigen::Matrix3d spaceSpaceMultiplier = 2.0 * currentGamma * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            Eigen::Matrix3d::Identity( );

    for( int i = 0; i < scalarPotentialPartials.rows( ); ++i )
    {
        metricPartials[ i ].setZero( );
        metricPartials[ i ]( 0, 0 ) = timeTimeMultiplier * scalarPotentialPartials( i );

        Eigen::Vector3d spaceTimeContribution =
                spaceTimeMultiplier * vectorPotentialPartials.block( 0, i, 3, 1 );
        metricPartials[ i ].block( 1, 0, 3, 1 ) = spaceTimeContribution;
        metricPartials[ i ].block( 0, 1, 1, 3 ) = spaceTimeContribution.transpose( );

        metricPartials[ i ].block( 1, 1, 3, 3 ) = spaceSpaceMultiplier * scalarPotentialPartials( i );
    }

    return metricPartials;
}

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialFromPotentialPartials(
        const double scalarPotentialPartials,
        const Eigen::Vector3d& vectorPotentialPartials,
        const double totalScalarPotential,
        const std::shared_ptr< relativity::PPNParameterSet > currentPpnParameterSet )
{
    Eigen::VectorXd singleScalarPartial( 1 );
    singleScalarPartial << scalarPotentialPartials;

    Eigen::Matrix< double, 3, 1 > singleVectorPartial = vectorPotentialPartials;

    return calculateSolarSystemMetricPartialFromPotentialPartials(
                singleScalarPartial,
                singleVectorPartial,
                totalScalarPotential,
                currentPpnParameterSet ).at( 0 );
}

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< SphericalHarmonicWrapper > shWrapper,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices )
{
    ( void ) shWrapper;
    ( void ) bodyIndex;
    ( void ) coefficientBlockIndices;
    throw std::runtime_error( "Solar system metric partials wrt spherical harmonic coefficients require SolarSystemMetric context." );
}

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const int maximumDegree,
        const int maximumOrder,
        const int minimumDegree,
        const int minimumOrder )
{
    std::map< int, std::pair< int, int > > coefficientBlocks;
    const auto blockVector = estimatable_parameters::getSphericalHarmonicBlockIndices(
            minimumDegree, minimumOrder, maximumDegree, maximumOrder );
    for( const auto& block : blockVector )
    {
        coefficientBlocks[ block.first ] = block;
    }

    return calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
                solarSystemMetric,
                bodyIndex,
                coefficientBlocks );
}

std::pair< std::vector< Eigen::Matrix< double, 4, 4 > >, std::vector< Eigen::Matrix< double, 4, 4 > > >
calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const int bodyIndex,
        const std::map< int, std::pair< int, int > >& coefficientBlockIndices )
{
    std::vector< Eigen::Matrix< double, 4, 4 > > cosineCoefficientMetricPartials;
    std::vector< Eigen::Matrix< double, 4, 4 > > sineCoefficientMetricPartials;

    if( solarSystemMetric->getBodySphericalHarmonicGravityWrappers( ).count( bodyIndex ) == 0 )
    {
        std::cerr << "Error cannot calculate metric partials wrt sh coefficients for body "
                  << solarSystemMetric->getBodyList( ).at( bodyIndex ) << std::endl;
        return std::make_pair( cosineCoefficientMetricPartials, sineCoefficientMetricPartials );
    }

    if( coefficientBlockIndices.empty( ) )
    {
        return std::make_pair( cosineCoefficientMetricPartials, sineCoefficientMetricPartials );
    }

    int totalNumberOfCoefficients = 0;
    for( const auto& block : coefficientBlockIndices )
    {
        totalNumberOfCoefficients += block.second.second;
    }

    std::shared_ptr< SphericalHarmonicWrapper > shWrapper =
            solarSystemMetric->getBodySphericalHarmonicGravityWrappers( ).at( bodyIndex );

    Eigen::VectorXd cosinePotentialPartials = Eigen::VectorXd::Zero( totalNumberOfCoefficients );
    std::vector< double > sinePotentialPartialsVector;
    sinePotentialPartialsVector.reserve( totalNumberOfCoefficients );

    basic_mathematics::LegendreCache* legendreCache = shWrapper->getLegendreCache( );
    if( legendreCache == nullptr )
    {
        throw std::runtime_error( "Legendre cache pointer not set in SolarSystemMetricPartial." );
    }

    Eigen::Vector3d bodyRelativePosition = solarSystemMetric->getCurrentBodyRelativePosition( bodyIndex );
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                shWrapper->getRotationToBodyFixedFrame( ) * bodyRelativePosition );

    const double sineOfLatitude = std::sin( mathematical_constants::PI / 2.0 - sphericalPosition( 1 ) );
    legendreCache->update( sineOfLatitude );

    const double bodyDistance = sphericalPosition( 0 );
    const double bodyReferenceRadius = shWrapper->getSphericalHarmonicGravityField( )->getReferenceRadius( );
    const double baseFactor = shWrapper->getSphericalHarmonicGravityField( )->getGravitationalParameter( ) / bodyDistance;
    const double radiusRatio = bodyReferenceRadius / bodyDistance;

    int cosineCounter = 0;

    for( const auto& block : coefficientBlockIndices )
    {
        const int degree = block.first;
        const int startOrder = block.second.first;
        const int numberOfOrders = block.second.second;

        const double premultiplier = baseFactor * std::pow( radiusRatio, degree );

        for( int order = startOrder; order < startOrder + numberOfOrders; ++order )
        {
            const double legendreValue = basic_mathematics::computeGeodesyLegendrePolynomialFromCache(
                        degree, order, *legendreCache );
            const double cosineLongitude = std::cos( order * sphericalPosition( 2 ) );
            cosinePotentialPartials( cosineCounter ) = premultiplier * legendreValue * cosineLongitude;
            ++cosineCounter;

            if( order > 0 )
            {
                const double sineLongitude = std::sin( order * sphericalPosition( 2 ) );
                sinePotentialPartialsVector.push_back( premultiplier * legendreValue * sineLongitude );
            }
        }
    }

    Eigen::VectorXd sinePotentialPartials = Eigen::VectorXd::Zero( static_cast< int >( sinePotentialPartialsVector.size( ) ) );
    for( int i = 0; i < static_cast< int >( sinePotentialPartialsVector.size( ) ); ++i )
    {
        sinePotentialPartials( i ) = sinePotentialPartialsVector[ i ];
    }

    const Eigen::Vector3d bodyVelocity = solarSystemMetric->getCurrentBodyStates( bodyIndex ).segment( 3, 3 );

    Eigen::Matrix< double, 3, Eigen::Dynamic > cosineVectorPotentialPartials =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, cosinePotentialPartials.rows( ) );
    for( int i = 0; i < cosinePotentialPartials.rows( ); ++i )
    {
        cosineVectorPotentialPartials.block( 0, i, 3, 1 ) = cosinePotentialPartials( i ) * bodyVelocity;
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > sineVectorPotentialPartials =
            Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, sinePotentialPartials.rows( ) );
    for( int i = 0; i < sinePotentialPartials.rows( ); ++i )
    {
        sineVectorPotentialPartials.block( 0, i, 3, 1 ) = sinePotentialPartials( i ) * bodyVelocity;
    }

    cosineCoefficientMetricPartials = calculateSolarSystemMetricPartialFromPotentialPartials(
                cosinePotentialPartials,
                cosineVectorPotentialPartials,
                solarSystemMetric->getCurrentScalarPotential( ),
                solarSystemMetric->getPpnParameterSet( ) );

    sineCoefficientMetricPartials = calculateSolarSystemMetricPartialFromPotentialPartials(
                sinePotentialPartials,
                sineVectorPotentialPartials,
                solarSystemMetric->getCurrentScalarPotential( ),
                solarSystemMetric->getPpnParameterSet( ) );

    return std::make_pair( cosineCoefficientMetricPartials, sineCoefficientMetricPartials );
}

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterGamma(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric )
{
    Eigen::Matrix< double, 4, 4 > currentPartial = Eigen::Matrix< double, 4, 4 >::Zero( );
    currentPartial.block( 0, 1, 1, 3 ) =
            -physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * solarSystemMetric->getCurrentVectorPotential( ).transpose( );
    currentPartial.block( 1, 0, 3, 1 ) =
            -physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * solarSystemMetric->getCurrentVectorPotential( );
    currentPartial.block( 1, 1, 3, 3 ) = 2.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            solarSystemMetric->getCurrentScalarPotential( ) * Eigen::Matrix3d::Identity( );
    return currentPartial;
}

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtPpnParameterBeta(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric )
{
    Eigen::Matrix< double, 4, 4 > currentPartial = Eigen::Matrix< double, 4, 4 >::Zero( );
    currentPartial( 0, 0 ) = -2.0 * physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT *
            solarSystemMetric->getCurrentScalarPotentialSquaredTerms( );
    return currentPartial;
}

Eigen::Matrix< double, 4, 4 > calculateSolarSystemMetricPartialWrtBodyGravitationalParameter(
        const int bodyIndex, const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric )
{
    const double normalizedCurrentBodyScalarPotential =
            solarSystemMetric->getCurrentSingleBodyScalarPotential( bodyIndex ) /
            solarSystemMetric->getCurrentGravitationalParameter( bodyIndex );

    return calculateSolarSystemMetricPartialFromPotentialPartials(
                normalizedCurrentBodyScalarPotential,
                normalizedCurrentBodyScalarPotential * solarSystemMetric->getCurrentBodyStates( bodyIndex ).segment( 3, 3 ),
                solarSystemMetric->getCurrentScalarPotential( ),
                solarSystemMetric->getPpnParameterSet( ) );
}

SolarSystemMetricPartial::SolarSystemMetricPartial(
        const std::shared_ptr< relativity::SolarSystemMetric > solarSystemMetric,
        const std::pair< std::string, std::string >& evaluationPointName,
        const std::function< Eigen::Quaterniond( const double ) > rotationFunctionFromBodyFixedFrame,
        const std::map< int, std::shared_ptr< SphericalHarmonicPartialWrapper > >& sphericalHarmonicPartialWrappers ):
    MetricPartial( ),
    solarSystemMetric_( solarSystemMetric ),
    evaluationPointName_( evaluationPointName ),
    sphericalHarmonicPartialWrappers_( sphericalHarmonicPartialWrappers ),
    bodyList_( solarSystemMetric_->getBodyList( ) ),
    rotationFunctionFromBodyFixedFrame_( rotationFunctionFromBodyFixedFrame )
{
    currentPartialWrtEvaluationPointPosition_.assign( 6, Eigen::Matrix< double, 4, 4 >::Zero( ) );
    currentPartialsWrtMetricBodyStates_.resize( bodyList_.size( ) );
    for( auto& bodyPartials : currentPartialsWrtMetricBodyStates_ )
    {
        bodyPartials.assign( 6, Eigen::Matrix< double, 4, 4 >::Zero( ) );
    }
}

bool SolarSystemMetricPartial::isMetricPartialWrtTranslationalStateNonNull( const std::string& bodyName )
{
    if( bodyName == evaluationPointName_.first )
    {
        return true;
    }

    return std::find( bodyList_.begin( ), bodyList_.end( ), bodyName ) != bodyList_.end( );
}

std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
SolarSystemMetricPartial::getDerivativeFunctionWrtStateOfIntegratedBody(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > partialFunction =
            std::make_pair( std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >( ), 0 );

    if( integratedStateType == propagators::translational_state )
    {
        if( !stateReferencePoint.second.empty( ) )
        {
            std::cerr << "Error when getting solar system metric partial function, cannot parse translational state of ("
                      << stateReferencePoint.first << "," << stateReferencePoint.second << ")." << std::endl;
        }
        else if( stateReferencePoint.first == evaluationPointName_.first )
        {
            partialFunction = std::make_pair(
                        [ this ]( ){ return getPartialWrtReferencePointPosition( ); }, 6 );
        }
        else
        {
            auto iterator = std::find( bodyList_.begin( ), bodyList_.end( ), stateReferencePoint.first );
            if( iterator != bodyList_.end( ) )
            {
                const int currentBodyIndex = static_cast< int >( std::distance( bodyList_.begin( ), iterator ) );
                partialFunction = std::make_pair(
                            [ this, currentBodyIndex ]( ){ return getCurrentPartialWrtBodyState( currentBodyIndex ); }, 6 );
            }
        }
    }

    return partialFunction;
}

std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > SolarSystemMetricPartial::getParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< Eigen::Matrix< double, 4, 4 >( ) >, int > partialFunctionPair =
            std::make_pair( std::function< Eigen::Matrix< double, 4, 4 >( ) >( ), 0 );

    switch( parameter->getParameterName( ).first )
    {
    case estimatable_parameters::gravitational_parameter:
    {
        auto bodyFindIterator = std::find( bodyList_.begin( ), bodyList_.end( ), parameter->getParameterName( ).second.first );
        if( bodyFindIterator != bodyList_.end( ) )
        {
            const int bodyIndex = static_cast< int >( std::distance( bodyList_.begin( ), bodyFindIterator ) );
            doubleParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                    [ this, bodyIndex ]( )
                    {
                        return calculateSolarSystemMetricPartialWrtBodyGravitationalParameter( bodyIndex, solarSystemMetric_ );
                    };

            partialFunctionPair = std::make_pair(
                        [ this, identifier = parameter->getParameterName( ) ]( )
                        {
                            return getDoubleParameterPartial( identifier );
                        },
                        1 );
        }
        break;
    }
    case estimatable_parameters::ppn_parameter_gamma:
        doubleParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                [ this ]( ) { return calculateSolarSystemMetricPartialWrtPpnParameterGamma( solarSystemMetric_ ); };
        partialFunctionPair = std::make_pair(
                    [ this, identifier = parameter->getParameterName( ) ]( )
                    {
                        return getDoubleParameterPartial( identifier );
                    },
                    1 );
        break;
    case estimatable_parameters::ppn_parameter_beta:
        doubleParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                [ this ]( ) { return calculateSolarSystemMetricPartialWrtPpnParameterBeta( solarSystemMetric_ ); };
        partialFunctionPair = std::make_pair(
                    [ this, identifier = parameter->getParameterName( ) ]( )
                    {
                        return getDoubleParameterPartial( identifier );
                    },
                    1 );
        break;
    default:
        break;
    }

    return partialFunctionPair;
}

std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int >
SolarSystemMetricPartial::getParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::pair< std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >, int > partialFunctionPair =
            std::make_pair( std::function< std::vector< Eigen::Matrix< double, 4, 4 > >( ) >( ), 0 );

    switch( parameter->getParameterName( ).first )
    {
    case estimatable_parameters::spherical_harmonics_cosine_coefficient_block:
    {
        auto bodyFindIterator = std::find( bodyList_.begin( ), bodyList_.end( ), parameter->getParameterName( ).second.first );
        if( bodyFindIterator != bodyList_.end( ) )
        {
            const int bodyIndex = static_cast< int >( std::distance( bodyList_.begin( ), bodyFindIterator ) );
            estimatable_parameters::EstimatebleParameterIdentifier correspondingSineParameterId = parameter->getParameterName( );
            correspondingSineParameterId.first = estimatable_parameters::spherical_harmonics_sine_coefficient_block;

            auto cosineCoefficientsParameter =
                    std::dynamic_pointer_cast< estimatable_parameters::SphericalHarmonicsCosineCoefficients >( parameter );
            vectorParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                    [ this,
                      bodyIndex,
                      parameterId = parameter->getParameterName( ),
                      blockIndices = cosineCoefficientsParameter->getBlockIndices( ) ]( )
                    {
                        std::map< int, std::pair< int, int > > blockMap;
                        for( const auto& block : blockIndices )
                        {
                            blockMap[ block.first ] = block;
                        }
                        return getUpdatedSphericalHarmonicPartials( bodyIndex, parameterId, blockMap );
                    };

            vectorParameterPartialFunctions_.erase( correspondingSineParameterId );

            partialFunctionPair = std::make_pair(
                        [ this, identifier = parameter->getParameterName( ) ]( )
                        {
                            return getVectorParameterPartial( identifier );
                        },
                        parameter->getParameterSize( ) );
        }
        break;
    }
    case estimatable_parameters::spherical_harmonics_sine_coefficient_block:
    {
        auto bodyFindIterator = std::find( bodyList_.begin( ), bodyList_.end( ), parameter->getParameterName( ).second.first );
        if( bodyFindIterator != bodyList_.end( ) )
        {
            const int bodyIndex = static_cast< int >( std::distance( bodyList_.begin( ), bodyFindIterator ) );

            estimatable_parameters::EstimatebleParameterIdentifier correspondingCosineParameterId = parameter->getParameterName( );
            correspondingCosineParameterId.first = estimatable_parameters::spherical_harmonics_cosine_coefficient_block;

            if( vectorParameterPartialFunctions_.count( correspondingCosineParameterId ) == 0 )
            {
                auto sineCoefficientsParameter =
                        std::dynamic_pointer_cast< estimatable_parameters::SphericalHarmonicsSineCoefficients >( parameter );
                vectorParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                        [ this,
                          bodyIndex,
                          parameterId = parameter->getParameterName( ),
                          blockIndices = sineCoefficientsParameter->getBlockIndices( ) ]( )
                        {
                            std::map< int, std::pair< int, int > > blockMap;
                            for( const auto& block : blockIndices )
                            {
                                blockMap[ block.first ] = block;
                            }
                            return getUpdatedSphericalHarmonicPartials( bodyIndex, parameterId, blockMap );
                        };
            }

            partialFunctionPair = std::make_pair(
                        [ this, identifier = parameter->getParameterName( ) ]( )
                        {
                            return getVectorParameterPartial( identifier );
                        },
                        parameter->getParameterSize( ) );
        }
        break;
    }
    case estimatable_parameters::ground_station_position:
        if( parameter->getParameterName( ).second == evaluationPointName_ )
        {
            vectorParameterPartialFunctions_[ parameter->getParameterName( ) ] =
                    [ this ]( )
                    {
                        return calculateSolarSystemMetricPartialWrtGroundStationPosition( );
                    };
            partialFunctionPair = std::make_pair(
                        [ this, identifier = parameter->getParameterName( ) ]( )
                        {
                            return getVectorParameterPartial( identifier );
                        },
                        3 );
        }
        break;
    default:
        break;
    }

    return partialFunctionPair;
}

void SolarSystemMetricPartial::update( )
{
    const double currentBeta = solarSystemMetric_->getPpnParameterSet( )->getParameterBeta( );
    const double currentGamma = solarSystemMetric_->getPpnParameterSet( )->getParameterGamma( );

    std::fill( currentPartialWrtEvaluationPointPosition_.begin( ), currentPartialWrtEvaluationPointPosition_.end( ),
               Eigen::Matrix< double, 4, 4 >::Zero( ) );

    for( auto& bodyPartials : currentPartialsWrtMetricBodyStates_ )
    {
        std::fill( bodyPartials.begin( ), bodyPartials.end( ), Eigen::Matrix< double, 4, 4 >::Zero( ) );
    }

    Eigen::Vector3d totalScalarPotentialGradient = Eigen::Vector3d::Zero( );

    for( unsigned int i = 0; i < bodyList_.size( ); ++i )
    {
        const double currentScalarPotential = solarSystemMetric_->getCurrentSingleBodyScalarPotential( static_cast< int >( i ) );
        const Eigen::Vector3d currentScalarPotentialGradient = solarSystemMetric_->getCurrentScalarPotentialGradient( static_cast< int >( i ) );
        const Eigen::Vector3d currentBodyVelocity = solarSystemMetric_->getCurrentBodyStates( static_cast< int >( i ) ).segment( 3, 3 );

        const std::vector< Eigen::Matrix< double, 4, 4 > > currentBodyPositionPartials =
                calculateSolarSystemMetricPartialFromPotentialPartials(
                    -currentScalarPotentialGradient,
                    -( currentBodyVelocity * currentScalarPotentialGradient.transpose( ) ),
                    solarSystemMetric_->getCurrentScalarPotential( ),
                    solarSystemMetric_->getPpnParameterSet( ) );

        totalScalarPotentialGradient += currentScalarPotentialGradient;

        for( int j = 0; j < 3; ++j )
        {
            currentPartialsWrtMetricBodyStates_[ i ][ j ] = currentBodyPositionPartials.at( j );

            currentPartialsWrtMetricBodyStates_[ i ][ j + 3 ].setZero( );
            currentPartialsWrtMetricBodyStates_[ i ][ j + 3 ]( 0, j + 1 ) =
                    -2.0 * ( currentGamma + 1.0 ) * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * currentScalarPotential;
            currentPartialsWrtMetricBodyStates_[ i ][ j + 3 ]( j + 1, 0 ) =
                    currentPartialsWrtMetricBodyStates_[ i ][ j + 3 ]( 0, j + 1 );

            currentPartialWrtEvaluationPointPosition_[ j ] -= currentPartialsWrtMetricBodyStates_[ i ][ j ];
        }
    }

    for( int j = 0; j < 3; ++j )
    {
        currentPartialWrtEvaluationPointPosition_[ j ]( 0, 0 ) =
                2.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * totalScalarPotentialGradient( j ) *
                ( 1.0 - 2.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * currentBeta *
                  solarSystemMetric_->getCurrentScalarPotential( ) );
    }

    updateParameterPartials( );
}

Eigen::Matrix< double, 4, 4 > SolarSystemMetricPartial::wrtScaledTime( )
{
    Eigen::Matrix< double, 4, 4 > timePartial = Eigen::Matrix< double, 4, 4 >::Zero( );

    for( unsigned int i = 0; i < bodyList_.size( ); ++i )
    {
        const Eigen::Matrix< double, 6, 1 > currentBodyState = solarSystemMetric_->getCurrentBodyStates( static_cast< int >( i ) );

        for( int j = 0; j < 3; ++j )
        {
            timePartial += currentPartialsWrtMetricBodyStates_[ i ][ j ] * currentBodyState( 3 + j );
            timePartial += currentPartialsWrtMetricBodyStates_[ i ][ j + 3 ] *
                    solarSystemMetric_->getCurrentBodyAccelerations( static_cast< int >( i ) )( j );
        }
    }

    const double inherentScalarPotentialTimePartial = solarSystemMetric_->getCurrentInherentScalarPotentialTimePartial( );
    const Eigen::Vector3d inherentVectorPotentialTimePartial =
            solarSystemMetric_->getCurrentInherentVectorPotentialTimePartial( );

    timePartial( 0, 0 ) += 2.0 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * inherentScalarPotentialTimePartial *
            ( 1.0 - 4.0 * solarSystemMetric_->getPpnParameterSet( )->getParameterBeta( ) *
              solarSystemMetric_->getCurrentScalarPotential( ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT );

    timePartial.block( 1, 1, 3, 3 ) += 2.0 * solarSystemMetric_->getPpnParameterSet( )->getParameterGamma( ) *
            Eigen::Matrix3d::Identity( ) * inherentScalarPotentialTimePartial * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;

    const double vectorScaling = -( solarSystemMetric_->getPpnParameterSet( )->getParameterGamma( ) + 1.0 ) *
            physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT;
    timePartial.block( 1, 0, 3, 1 ) += vectorScaling * inherentVectorPotentialTimePartial;
    timePartial.block( 0, 1, 1, 3 ) += vectorScaling * inherentVectorPotentialTimePartial.transpose( );

    timePartial /= physical_constants::SPEED_OF_LIGHT;

    return timePartial;
}

std::vector< Eigen::Matrix< double, 4, 4 > > SolarSystemMetricPartial::getCurrentPartialWrtBodyState( const int bodyIndex )
{
    return currentPartialsWrtMetricBodyStates_.at( bodyIndex );
}

std::vector< Eigen::Matrix< double, 4, 4 > > SolarSystemMetricPartial::getPartialWrtReferencePointPosition( )
{
    return currentPartialWrtEvaluationPointPosition_;
}

Eigen::Matrix< double, 4, 4 > SolarSystemMetricPartial::getDoubleParameterPartial(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    auto iterator = currentDoubleParameterPartials_.find( parameterIdentifier );
    if( iterator != currentDoubleParameterPartials_.end( ) )
    {
        return iterator->second;
    }

    std::cerr << "Error, requested solar system metric partial " << parameterIdentifier.first << " "
              << parameterIdentifier.second.first << " " << parameterIdentifier.second.second
              << ", but partial not found." << std::endl;
    return Eigen::Matrix< double, 4, 4 >::Zero( );
}

std::vector< Eigen::Matrix< double, 4, 4 > > SolarSystemMetricPartial::getVectorParameterPartial(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    auto iterator = currentVectorParameterPartials_.find( parameterIdentifier );
    if( iterator != currentVectorParameterPartials_.end( ) )
    {
        return iterator->second;
    }

    std::cerr << "Error, requested solar system metric partial " << parameterIdentifier.first << " "
              << parameterIdentifier.second.first << " " << parameterIdentifier.second.second
              << ", but partial not found." << std::endl;
    return std::vector< Eigen::Matrix< double, 4, 4 > >( );
}

std::vector< Eigen::Matrix< double, 4, 4 > > SolarSystemMetricPartial::calculateSolarSystemMetricPartialWrtGroundStationPosition( )
{
    const Eigen::Matrix3d currentRotationFromBodyFixedFrame =
            rotationFunctionFromBodyFixedFrame_( solarSystemMetric_->getCurrentTime( ) ).toRotationMatrix( );
    return getMetricPartialWrtGroundStationPosition( currentPartialWrtEvaluationPointPosition_, currentRotationFromBodyFixedFrame );
}

void SolarSystemMetricPartial::updateParameterPartials( )
{
    for( const auto& parameterFunction : doubleParameterPartialFunctions_ )
    {
        currentDoubleParameterPartials_[ parameterFunction.first ] = parameterFunction.second( );
    }

    for( const auto& parameterFunction : vectorParameterPartialFunctions_ )
    {
        currentVectorParameterPartials_[ parameterFunction.first ] = parameterFunction.second( );
    }
}

std::vector< Eigen::Matrix< double, 4, 4 > > SolarSystemMetricPartial::getUpdatedSphericalHarmonicPartials(
        const int bodyIndex,
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const std::map< int, std::pair< int, int > > coefficientBlockIndices )
{
    const auto coefficientPartials = calculateSolarSystemMetricPartialWrtSphericalHarmonicCoefficients(
                solarSystemMetric_, bodyIndex, coefficientBlockIndices );

    std::vector< Eigen::Matrix< double, 4, 4 > > partialsToReturn;
    estimatable_parameters::EstimatebleParameterIdentifier dualParameter = parameterId;

    if( parameterId.first == estimatable_parameters::spherical_harmonics_cosine_coefficient_block )
    {
        dualParameter.first = estimatable_parameters::spherical_harmonics_sine_coefficient_block;
        partialsToReturn = coefficientPartials.first;
        currentVectorParameterPartials_[ dualParameter ] = coefficientPartials.second;
    }
    else if( parameterId.first == estimatable_parameters::spherical_harmonics_sine_coefficient_block )
    {
        dualParameter.first = estimatable_parameters::spherical_harmonics_cosine_coefficient_block;
        partialsToReturn = coefficientPartials.second;
        currentVectorParameterPartials_[ dualParameter ] = coefficientPartials.first;
    }
    else
    {
        throw std::runtime_error( "Error when updating metric SH partials, input parameter is inconsistent" );
    }

    return partialsToReturn;
}

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat
