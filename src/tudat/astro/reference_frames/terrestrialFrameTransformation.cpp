#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include <Eigen/LU>
#include <cmath>
#include <stdexcept>
#include <optional>
namespace tudat
{
namespace reference_frames
{
namespace
{
void validateEpoch( const TerrestrialFrameEvaluationEpoch& epoch )
{
    if( !std::isfinite( epoch.decimalYear ) || !std::isfinite( epoch.decimalYearRateWithRespectToInputSeconds ) ||
        epoch.decimalYearRateWithRespectToInputSeconds <= 0 )
        throw std::runtime_error( "Terrestrial transformation requires explicit finite decimal year and positive dYear/dInputSecond" );
}
Eigen::Matrix3d crossMatrix( const Eigen::Vector3d& vector )
{
    Eigen::Matrix3d matrix;
    matrix << 0, -vector.z( ), vector.y( ), vector.z( ), 0, -vector.x( ), -vector.y( ), vector.x( ), 0;
    return matrix;
}
void validateAffine( const TerrestrialFrameAffineTransformation& affine )
{
    if( !affine.translation.allFinite( ) || !affine.matrix.allFinite( ) || !affine.translationRate.allFinite( ) ||
        !affine.matrixRate.allFinite( ) )
        throw std::runtime_error( "Nonfinite terrestrial affine transformation" );
}
}  // namespace
Eigen::Vector3d TerrestrialFrameAffineTransformation::transformPosition( const Eigen::Vector3d& position ) const
{
    validateAffine( *this );
    if( !position.allFinite( ) ) throw std::runtime_error( "Nonfinite terrestrial position" );
    return translation + matrix * position;
}
Eigen::Vector6d TerrestrialFrameAffineTransformation::transformState( const Eigen::Vector6d& state ) const
{
    validateAffine( *this );
    if( !state.allFinite( ) ) throw std::runtime_error( "Nonfinite terrestrial state" );
    Eigen::Vector6d result;
    result.head< 3 >( ) = translation + matrix * state.head< 3 >( );
    result.tail< 3 >( ) = translationRate + matrixRate * state.head< 3 >( ) + matrix * state.tail< 3 >( );
    return result;
}
TerrestrialFrameAffineTransformation TerrestrialFrameAffineTransformation::inverse( ) const
{
    validateAffine( *this );
    Eigen::FullPivLU< Eigen::Matrix3d > decomposition( matrix );
    if( !decomposition.isInvertible( ) ) throw std::runtime_error( "Singular terrestrial-frame affine matrix" );
    TerrestrialFrameAffineTransformation result;
    result.matrix = decomposition.inverse( );
    result.matrixRate = -result.matrix * matrixRate * result.matrix;
    result.translation = -result.matrix * translation;
    result.translationRate = -result.matrixRate * translation - result.matrix * translationRate;
    validateAffine( result );
    return result;
}
TerrestrialFrameAffineTransformation TerrestrialFrameAffineTransformation::composeAfter(
        const TerrestrialFrameAffineTransformation& before ) const
{
    validateAffine( *this );
    validateAffine( before );
    TerrestrialFrameAffineTransformation result;
    result.matrix = matrix * before.matrix;
    result.translation = translation + matrix * before.translation;
    result.matrixRate = matrixRate * before.matrix + matrix * before.matrixRate;
    result.translationRate = translationRate + matrixRate * before.translation + matrix * before.translationRate;
    validateAffine( result );
    return result;
}
namespace
{
TerrestrialFrameAffineTransformation evaluateHelmertCoordinates( const input_output::TerrestrialFrameTransformationParameters& p,
                                                                 double decimalYear,
                                                                 std::optional< double > rate )
{
    input_output::validateTerrestrialFrameTransformationParameters( p );
    if( !std::isfinite( decimalYear ) ) throw std::runtime_error( "Nonfinite terrestrial-frame decimal year" );
    if( p.rotationConvention == input_output::HelmertRotationConvention::unspecified )
        throw std::runtime_error( "Explicit Helmert rotation convention required" );
    const double sign = p.rotationConvention == input_output::HelmertRotationConvention::positionVector ? 1.0 : -1.0;
    const double elapsed = p.parameterRateYearSeconds * ( decimalYear - p.referenceDecimalYear );
    // Convert normalized SI rates back to the published parameter-year difference, then chain its physical differential.
    const double rateFactor = rate ? p.parameterRateYearSeconds * *rate : 0.0;
    TerrestrialFrameAffineTransformation result;
    result.translation = p.translationMetres + elapsed * p.translationRateMetresPerSecond;
    result.matrix = ( 1 + p.scale + elapsed * p.scaleRatePerSecond ) * Eigen::Matrix3d::Identity( ) +
            sign * crossMatrix( p.rotationRadians + elapsed * p.rotationRateRadiansPerSecond );
    result.translationRate = rateFactor * p.translationRateMetresPerSecond;
    result.matrixRate =
            rateFactor * ( p.scaleRatePerSecond * Eigen::Matrix3d::Identity( ) + sign * crossMatrix( p.rotationRateRadiansPerSecond ) );
    validateAffine( result );
    return result;
}
}  // namespace

TerrestrialFrameAffineTransformation evaluateLinearizedHelmert( const input_output::TerrestrialFrameTransformationParameters& p,
                                                                const TerrestrialFrameEvaluationEpoch& epoch )
{
    validateEpoch( epoch );
    return evaluateHelmertCoordinates( p, epoch.decimalYear, epoch.decimalYearRateWithRespectToInputSeconds );
}

input_output::TerrestrialFrameTransformationData getConventionalItrfTransformationData( ItrfTransformationCatalogue catalogue )
{
    switch( catalogue )
    {
        case ItrfTransformationCatalogue::ignItrf2014:
            return input_output::readOfficialItrf2014TransformationData( );
        case ItrfTransformationCatalogue::ignItrf2020:
            return input_output::readOfficialItrf2020TransformationData( );
        default:
            throw std::runtime_error( "Unknown conventional ITRF parameter catalogue" );
    }
}
TerrestrialFrameTransformationModel::TerrestrialFrameTransformationModel( const input_output::TerrestrialFrameTransformationData& data,
                                                                          Clock clock ): data_( data ), clock_( std::move( clock ) )
{
    input_output::validateTerrestrialFrameTransformationData( data_ );
    for( const auto& entry : data_.records )
        if( entry.second.rotationConvention == input_output::HelmertRotationConvention::unspecified )
            throw std::runtime_error( "Terrestrial-frame model requires a declared rotation convention" );
}
TerrestrialFrameAffineTransformation TerrestrialFrameTransformationModel::getTransformation(
        const std::string& source,
        const std::string& target,
        const TerrestrialFrameEvaluationEpoch& epoch ) const
{
    validateEpoch( epoch );
    if( source.empty( ) || target.empty( ) ) throw std::runtime_error( "Empty terrestrial-frame realization" );
    // Release provenance is not a new datum; aligned updates retain the same realization.
    if( source == target ) return {};
    TerrestrialFrameAffineTransformation fromHub, toHub;
    if( source != data_.referenceFrame.realization ) toHub = evaluateLinearizedHelmert( data_.getParameters( source ), epoch ).inverse( );
    if( target != data_.referenceFrame.realization ) fromHub = evaluateLinearizedHelmert( data_.getParameters( target ), epoch );
    return fromHub.composeAfter( toHub );
}
TerrestrialFrameAffineTransformation TerrestrialFrameTransformationModel::getTransformationAtTime( const std::string& source,
                                                                                                   const std::string& target,
                                                                                                   const Time& inputTime ) const
{
    if( !clock_ )
        throw std::runtime_error( "Native terrestrial-frame time evaluation requires the model-owned explicit decimal-year clock" );
    return getTransformation( source, target, clock_( inputTime ) );
}
Eigen::Vector3d TerrestrialFrameTransformationModel::transformPosition( const Eigen::Vector3d& position,
                                                                        const std::string& source,
                                                                        const std::string& target,
                                                                        const TerrestrialFrameEvaluationEpoch& epoch ) const
{
    return getTransformation( source, target, epoch ).transformPosition( position );
}
Eigen::Vector3d TerrestrialFrameTransformationModel::transformPositionAtDecimalYear( const Eigen::Vector3d& position,
                                                                                     const std::string& source,
                                                                                     const std::string& target,
                                                                                     double decimalYear ) const
{
    if( !std::isfinite( decimalYear ) || source.empty( ) || target.empty( ) )
        throw std::runtime_error( "Invalid value-only terrestrial-frame coordinates" );
    if( source == target )
    {
        if( !position.allFinite( ) ) throw std::runtime_error( "Nonfinite terrestrial position" );
        return position;
    }
    // These internal affine objects carry only value coordinates. Their unused
    // derivative entries never escape this position-only API.
    TerrestrialFrameAffineTransformation fromHub, toHub;
    if( source != data_.referenceFrame.realization )
        toHub = evaluateHelmertCoordinates( data_.getParameters( source ), decimalYear, std::nullopt ).inverse( );
    if( target != data_.referenceFrame.realization )
        fromHub = evaluateHelmertCoordinates( data_.getParameters( target ), decimalYear, std::nullopt );
    return fromHub.composeAfter( toHub ).transformPosition( position );
}
Eigen::Vector6d TerrestrialFrameTransformationModel::transformState( const Eigen::Vector6d& state,
                                                                     const std::string& source,
                                                                     const std::string& target,
                                                                     const TerrestrialFrameEvaluationEpoch& epoch ) const
{
    return getTransformation( source, target, epoch ).transformState( state );
}
Eigen::Vector3d TerrestrialFrameTransformationModel::transformPositionAtTime( const Eigen::Vector3d& position,
                                                                              const std::string& source,
                                                                              const std::string& target,
                                                                              const Time& inputTime ) const
{
    return getTransformationAtTime( source, target, inputTime ).transformPosition( position );
}
Eigen::Vector6d TerrestrialFrameTransformationModel::transformStateAtTime( const Eigen::Vector6d& state,
                                                                           const std::string& source,
                                                                           const std::string& target,
                                                                           const Time& inputTime ) const
{
    return getTransformationAtTime( source, target, inputTime ).transformState( state );
}
TerrestrialFrameEvaluationEpoch getLegacyItrfEvaluationEpoch( double secondsSinceJ2000 )
{
    if( !std::isfinite( secondsSinceJ2000 ) ) throw std::runtime_error( "Nonfinite legacy ITRF epoch" );
    return { 2000 + secondsSinceJ2000 / 31557600.0, 1 / 31557600.0 };
}
}  // namespace reference_frames
}  // namespace tudat
