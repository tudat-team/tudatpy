#ifndef TUDAT_READ_TERRESTRIAL_FRAME_TRANSFORMATION_DATA_H
#define TUDAT_READ_TERRESTRIAL_FRAME_TRANSFORMATION_DATA_H

#include <istream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace tudat
{
namespace input_output
{
struct TerrestrialFrameIdentity {
    std::string realization;
    //! Catalogue release provenance; an aligned release update is not a new datum.
    std::string release;
};
enum class HelmertRotationConvention { unspecified, positionVector, coordinateFrame };
struct TerrestrialFrameTransformationSource {
    std::string source = "<stream>";
    std::string sourceRelease;
    std::string sourceSha256;
};
//! Neutral SI records. No time mapping, propagation, table selection or frame composition.
struct TerrestrialFrameTransformationParameters {
    TerrestrialFrameIdentity sourceFrame, targetFrame;
    double referenceDecimalYear = 0;
    Eigen::Vector3d translationMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRadians = Eigen::Vector3d::Zero( );
    double scale = 0;
    Eigen::Vector3d translationRateMetresPerSecond = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRateRadiansPerSecond = Eigen::Vector3d::Zero( );
    double scaleRatePerSecond = 0;
    //! Unit conversion only. Decimal-year propagation multiplies SI rates by this factor.
    double parameterRateYearSeconds = 31557600.0;
    HelmertRotationConvention rotationConvention = HelmertRotationConvention::unspecified;
    std::string originalTranslationUnits = "mm", originalRotationUnits = "mas", originalScaleUnits = "ppb";
    std::string originalRateTimeUnits = "Julian year (365.25 days)";
    TerrestrialFrameTransformationSource provenance;
    std::string originalParameterLine, originalRateLine;
    std::size_t parameterSourceLine = 0, rateSourceLine = 0;
};
struct TerrestrialFrameTransformationData {
    TerrestrialFrameIdentity referenceFrame;
    TerrestrialFrameTransformationSource provenance;
    std::vector< std::string > sourceHeaders;
    std::map< std::string, TerrestrialFrameTransformationParameters > records;
    const TerrestrialFrameTransformationParameters& getParameters( const std::string& targetRealization ) const;
};
void validateTerrestrialFrameTransformationParameters( const TerrestrialFrameTransformationParameters& parameters );
void validateTerrestrialFrameTransformationData( const TerrestrialFrameTransformationData& data );
//! Strict published IGN table reader, including declared unit/order/position-vector equation checks.
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( std::istream& stream,
                                                                           const TerrestrialFrameTransformationSource& source = {} );
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( const std::string& filePath,
                                                                           const TerrestrialFrameTransformationSource& source = {} );
//! Pinned raw official input, parsed by the same neutral reader. Model selection belongs to the consuming layer.
TerrestrialFrameTransformationData readOfficialItrf2014TransformationData( );
TerrestrialFrameTransformationData readOfficialItrf2020TransformationData( );
}  // namespace input_output
}  // namespace tudat
#endif
