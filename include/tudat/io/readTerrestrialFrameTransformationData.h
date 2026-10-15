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

//! Identifier of a terrestrial reference frame (for instance ITRF2014).
struct TerrestrialFrameIdentity {
    //! Name of the frame realization.
    std::string realization;

    //! Catalogue release provenance; an aligned release update is not a new datum.
    std::string release;
};

//! Sign convention of the rotation parameters of a Helmert transformation.
enum class HelmertRotationConvention { unspecified, positionVector, coordinateFrame };

//! Origin of a table of terrestrial frame transformation parameters.
struct TerrestrialFrameTransformationSource {
    //! Name of the file (or other source), its release, and its SHA-256 checksum.
    std::string source = "<stream>";
    std::string sourceRelease;
    std::string sourceSha256;
};

//! Parameters of a Helmert transformation between two terrestrial reference frames, in SI units, as given in the table.
struct TerrestrialFrameTransformationParameters {
    //! Frames from which and to which the transformation is applied.
    TerrestrialFrameIdentity sourceFrame, targetFrame;

    //! Epoch (decimal year) at which the parameters below apply.
    double referenceDecimalYear = 0;

    //! Translation (metres), rotation (radians) and scale (dimensionless) at the reference epoch.
    Eigen::Vector3d translationMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRadians = Eigen::Vector3d::Zero( );
    double scale = 0;

    //! Rates of the translation, rotation and scale, per second.
    Eigen::Vector3d translationRateMetresPerSecond = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRateRadiansPerSecond = Eigen::Vector3d::Zero( );
    double scaleRatePerSecond = 0;

    //! Unit conversion only. Decimal-year propagation multiplies SI rates by this factor.
    double parameterRateYearSeconds = 31557600.0;

    //! Sign convention of the rotation parameters.
    HelmertRotationConvention rotationConvention = HelmertRotationConvention::unspecified;

    //! Units of the parameters and rates in the table.
    std::string originalTranslationUnits = "mm", originalRotationUnits = "mas", originalScaleUnits = "ppb";
    std::string originalRateTimeUnits = "Julian year (365.25 days)";

    //! Origin of the table, and the text and line numbers of the parameter and rate lines.
    TerrestrialFrameTransformationSource provenance;
    std::string originalParameterLine, originalRateLine;
    std::size_t parameterSourceLine = 0, rateSourceLine = 0;
};

//! Table of Helmert transformation parameters from one reference frame to several other frames.
struct TerrestrialFrameTransformationData {
    //! Frame from which all transformations in the table start, and the origin of the table.
    TerrestrialFrameIdentity referenceFrame;
    TerrestrialFrameTransformationSource provenance;

    //! Header lines of the table.
    std::vector< std::string > sourceHeaders;

    //! Transformation parameters, with the name of the target frame realization as key.
    std::map< std::string, TerrestrialFrameTransformationParameters > records;

    //! Parameters of the transformation to the given target frame realization; throws if it is not in the table.
    const TerrestrialFrameTransformationParameters& getParameters( const std::string& targetRealization ) const;
};

void validateTerrestrialFrameTransformationParameters( const TerrestrialFrameTransformationParameters& parameters );
void validateTerrestrialFrameTransformationData( const TerrestrialFrameTransformationData& data );

//! Strict published IGN table reader, including declared unit/order/position-vector equation checks.
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( std::istream& stream,
                                                                           const TerrestrialFrameTransformationSource& source = {} );
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( const std::string& filePath,
                                                                           const TerrestrialFrameTransformationSource& source = {} );

//! Official IGN tables from the station_locations resource directory, parsed by the same reader.
TerrestrialFrameTransformationData readOfficialItrf2014TransformationData( );
TerrestrialFrameTransformationData readOfficialItrf2020TransformationData( );
}  // namespace input_output
}  // namespace tudat
#endif
