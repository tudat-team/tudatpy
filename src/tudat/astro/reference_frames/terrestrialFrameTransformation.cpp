#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include "tudat/math/basic/linearAlgebra.h"
#include <Eigen/LU>
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace reference_frames
{
namespace
{

//! Require a finite published year; state evaluation also needs a positive years-per-second rate.
void validateCoordinates( double decimalYear, const std::optional< double >& yearRate )
{
    if( !std::isfinite( decimalYear ) || ( yearRate && ( !std::isfinite( *yearRate ) || *yearRate <= 0 ) ) )
    {
        throw std::runtime_error( "Terrestrial transformation requires finite decimal year and positive dYear/dInputSecond for state" );
    }
}

//! Reject nonfinite affine values and derivatives before applying, inverting or composing them.
void validateAffine( const TerrestrialFrameAffineTransformation& affine )
{
    if( !affine.translation.allFinite( ) || !affine.matrix.allFinite( ) || !affine.translationRate.allFinite( ) ||
        !affine.matrixRate.allFinite( ) )
    {
        throw std::runtime_error( "Nonfinite terrestrial affine transformation" );
    }
}

//! Form an affine map from normalized SI parameters; absent rates remain unused by position evaluation.
TerrestrialFrameAffineTransformation evaluateHelmertCoordinates( const input_output::TerrestrialFrameTransformationParameters& p,
                                                                 double decimalYear,
                                                                 std::optional< double > rate )
{
    // The declared rotation convention sets the sign; elapsed is expressed in the published rate-year unit.
    const double sign = p.rotationConvention == input_output::HelmertRotationConvention::positionVector ? 1.0 : -1.0;

    const double elapsed = p.parameterRateYearSeconds * ( decimalYear - p.referenceDecimalYear );

    // Chain SI-normalized parameter rates to the caller's input seconds only when a state was requested.
    const double rateFactor = rate ? p.parameterRateYearSeconds * *rate : 0.0;

    TerrestrialFrameAffineTransformation result;

    result.translation = p.translationMetres + elapsed * p.translationRateMetresPerSecond;
    result.matrix = ( 1 + p.scale + elapsed * p.scaleRatePerSecond ) * Eigen::Matrix3d::Identity( ) +
            sign * linear_algebra::getCrossProductMatrix( p.rotationRadians + elapsed * p.rotationRateRadiansPerSecond );
    result.translationRate = rateFactor * p.translationRateMetresPerSecond;
    result.matrixRate = rateFactor *
            ( p.scaleRatePerSecond * Eigen::Matrix3d::Identity( ) +
              sign * linear_algebra::getCrossProductMatrix( p.rotationRateRadiansPerSecond ) );
    validateAffine( result );
    return result;
}

}  // namespace

//! Apply metre-valued translation and the complete dimensionless affine matrix.
Eigen::Vector3d TerrestrialFrameAffineTransformation::transformPosition( const Eigen::Vector3d& position ) const
{
    validateAffine( *this );
    if( !position.allFinite( ) )
    {
        throw std::runtime_error( "Nonfinite terrestrial position" );
    }
    return translation + matrix * position;
}

//! Apply position and parameter-rate contributions to velocity in the same input seconds.
Eigen::Vector6d TerrestrialFrameAffineTransformation::transformState( const Eigen::Vector6d& state ) const
{
    validateAffine( *this );
    if( !state.allFinite( ) )
    {
        throw std::runtime_error( "Nonfinite terrestrial state" );
    }

    // Cartesian position is followed by its derivative in the caller's input time variable.
    Eigen::Vector6d result;

    result.head< 3 >( ) = translation + matrix * state.head< 3 >( );
    result.tail< 3 >( ) = translationRate + matrixRate * state.head< 3 >( ) + matrix * state.tail< 3 >( );
    return result;
}

//! Invert and differentiate the full matrix, retaining scale rather than assuming an orthogonal rotation.
TerrestrialFrameAffineTransformation TerrestrialFrameAffineTransformation::inverse( ) const
{
    validateAffine( *this );

    // A general solve is required because the affine matrix includes scale.
    Eigen::FullPivLU< Eigen::Matrix3d > decomposition( matrix );

    if( !decomposition.isInvertible( ) )
    {
        throw std::runtime_error( "Singular terrestrial-frame affine matrix" );
    }

    // Differentiate the inverse matrix before forming inverse translation and its rate.
    TerrestrialFrameAffineTransformation result;

    result.matrix = decomposition.inverse( );
    result.matrixRate = -result.matrix * matrixRate * result.matrix;
    result.translation = -result.matrix * translation;
    result.translationRate = -result.matrixRate * translation - result.matrix * translationRate;
    validateAffine( result );
    return result;
}

//! Compose the supplied map first and include all matrix and translation product-rule terms.
TerrestrialFrameAffineTransformation TerrestrialFrameAffineTransformation::composeAfter(
        const TerrestrialFrameAffineTransformation& before ) const
{
    validateAffine( *this );
    validateAffine( before );

    // Values and derivatives use the same input seconds for both maps.
    TerrestrialFrameAffineTransformation result;

    result.matrix = matrix * before.matrix;
    result.translation = translation + matrix * before.translation;
    result.matrixRate = matrixRate * before.matrix + matrix * before.matrixRate;
    result.translationRate = translationRate + matrixRate * before.translation + matrix * before.translationRate;
    validateAffine( result );
    return result;
}

//! Evaluate one published parameter record with a caller-supplied year and physical year rate.
TerrestrialFrameAffineTransformation evaluateLinearizedHelmert( const input_output::TerrestrialFrameTransformationParameters& p,
                                                                double decimalYear,
                                                                double decimalYearRateWithRespectToInputSeconds )
{
    validateCoordinates( decimalYear, decimalYearRateWithRespectToInputSeconds );
    input_output::validateTerrestrialFrameTransformationParameters( p );
    if( p.rotationConvention == input_output::HelmertRotationConvention::unspecified )
    {
        throw std::runtime_error( "Explicit Helmert rotation convention required" );
    }
    return evaluateHelmertCoordinates( p, decimalYear, decimalYearRateWithRespectToInputSeconds );
}

//! Store the validated catalogue and its declared rotation conventions.
TerrestrialFrameTransformationModel::TerrestrialFrameTransformationModel( const input_output::TerrestrialFrameTransformationData& data ):
    data_( data )
{
    input_output::validateTerrestrialFrameTransformationData( data_ );
    for( const auto& entry : data_.records )
    {
        if( entry.second.rotationConvention == input_output::HelmertRotationConvention::unspecified )
        {
            throw std::runtime_error( "Terrestrial-frame model requires a declared rotation convention" );
        }
    }
}

//! Share the source-to-hub-to-target calculation between position and state evaluation.
TerrestrialFrameAffineTransformation TerrestrialFrameTransformationModel::evaluateTransformation( const std::string& source,
                                                                                                  const std::string& target,
                                                                                                  double decimalYear,
                                                                                                  std::optional< double > yearRate ) const
{
    validateCoordinates( decimalYear, yearRate );
    if( source.empty( ) || target.empty( ) )
    {
        throw std::runtime_error( "Empty terrestrial-frame realization" );
    }
    if( source == target )
    {
        return {};
    }

    // The catalogue's reference realization is the common hub for both directed maps.
    TerrestrialFrameAffineTransformation fromHub, toHub;

    if( source != data_.referenceFrame.realization )
    {
        toHub = evaluateHelmertCoordinates( data_.getParameters( source ), decimalYear, yearRate ).inverse( );
    }
    if( target != data_.referenceFrame.realization )
    {
        fromHub = evaluateHelmertCoordinates( data_.getParameters( target ), decimalYear, yearRate );
    }
    return fromHub.composeAfter( toHub );
}

//! Return affine values and derivatives using the explicit years-per-input-second mapping.
TerrestrialFrameAffineTransformation TerrestrialFrameTransformationModel::getTransformation(
        const std::string& source,
        const std::string& target,
        double decimalYear,
        double decimalYearRateWithRespectToInputSeconds ) const
{
    return evaluateTransformation( source, target, decimalYear, decimalYearRateWithRespectToInputSeconds );
}

//! Transform metre-valued position at the supplied decimal year.
Eigen::Vector3d TerrestrialFrameTransformationModel::transformPosition( const Eigen::Vector3d& position,
                                                                        const std::string& source,
                                                                        const std::string& target,
                                                                        double decimalYear ) const
{
    return evaluateTransformation( source, target, decimalYear, std::nullopt ).transformPosition( position );
}

//! Transform Cartesian position and velocity with the caller's explicit year rate.
Eigen::Vector6d TerrestrialFrameTransformationModel::transformState( const Eigen::Vector6d& state,
                                                                     const std::string& source,
                                                                     const std::string& target,
                                                                     double decimalYear,
                                                                     double decimalYearRateWithRespectToInputSeconds ) const
{
    return getTransformation( source, target, decimalYear, decimalYearRateWithRespectToInputSeconds ).transformState( state );
}

}  // namespace reference_frames
}  // namespace tudat
