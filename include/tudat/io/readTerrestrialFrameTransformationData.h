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
enum class HelmertRotationConvention { positionVector, coordinateFrame };

//! Parameters of a Helmert transformation between two terrestrial reference frames, in SI units.
struct TerrestrialFrameTransformationParameters {
    //! Frames from which and to which the transformation is applied.
    TerrestrialFrameIdentity sourceFrame, targetFrame;

    //! Epoch (decimal year) at which the parameters below apply.
    double referenceDecimalYear = 0;

    //! Translation (metres), rotation (radians) and scale (dimensionless) at the reference epoch.
    Eigen::Vector3d translationMetres = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRadians = Eigen::Vector3d::Zero( );
    double scale = 0;

    //! Rates of the translation, rotation and scale, per second (of a year of 365.25 days).
    Eigen::Vector3d translationRateMetresPerSecond = Eigen::Vector3d::Zero( );
    Eigen::Vector3d rotationRateRadiansPerSecond = Eigen::Vector3d::Zero( );
    double scaleRatePerSecond = 0;

    //! Sign convention of the rotation parameters.
    HelmertRotationConvention rotationConvention = HelmertRotationConvention::positionVector;
};

//! Table of Helmert transformation parameters from one reference frame to several other frames.
struct TerrestrialFrameTransformationData {
    //! Frame from which all transformations in the table start, and the name of the file of the table.
    TerrestrialFrameIdentity referenceFrame;
    std::string source;

    //! Transformation parameters, with the name of the target frame realization as key.
    std::map< std::string, TerrestrialFrameTransformationParameters > records;

    //! Parameters of the transformation to the given target frame realization; throws if it is not in the table.
    const TerrestrialFrameTransformationParameters& getParameters( const std::string& targetRealization ) const;
};

//! Read an IGN table of transformation parameters (for instance Transfo-ITRF2020_TRFs.txt), with position-vector rotations.
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( std::istream& stream, const std::string& source = "<stream>" );
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( const std::string& filePath );

//! Official IGN tables from the station_locations resource directory.
TerrestrialFrameTransformationData readOfficialItrf2014TransformationData( );
TerrestrialFrameTransformationData readOfficialItrf2020TransformationData( );
}  // namespace input_output
}  // namespace tudat
#endif
