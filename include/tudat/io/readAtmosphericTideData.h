#ifndef TUDAT_READ_ATMOSPHERIC_TIDE_DATA_H
#define TUDAT_READ_ATMOSPHERIC_TIDE_DATA_H

#include <istream>
#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace tudat
{
namespace input_output
{
enum class AtmosphericLoadingReferenceOrigin { unspecified, solidEarthCentre, wholeEarthCentreOfMass };
enum class AtmosphericTideLatitudeConvention { unspecified, geocentric, geodetic };

//! Geometric interpretation of provider component labels, separate from the grid lookup latitude convention.
enum class AtmosphericTideDisplacementBasis { unspecified, geocentricRadialNorthEast, geodeticUpNorthEast };
//! Explicitly identified conventional source profile. A model name or filename alone does not establish this profile.
enum class AtmosphericTideSourceProfile { unspecified, iers2010Rp03NoIbFarrell };

//! Neutral source metadata. No origin or latitude convention is inferred from a filename.
struct AtmosphericTideMetadata {
    std::string source;
    std::vector< std::string > sourceHeaders;
    std::string modelName;
    std::string greensFunction;
    AtmosphericLoadingReferenceOrigin referenceOrigin = AtmosphericLoadingReferenceOrigin::unspecified;
    AtmosphericTideLatitudeConvention latitudeConvention = AtmosphericTideLatitudeConvention::unspecified;
    //! Optional expected dimensions for source formats with no count header; zero means infer the provided extent.
    std::size_t expectedLongitudeCount = 0;
    std::size_t expectedLatitudeCount = 0;
    AtmosphericTideDisplacementBasis displacementBasis = AtmosphericTideDisplacementBasis::unspecified;
    AtmosphericTideSourceProfile sourceProfile = AtmosphericTideSourceProfile::unspecified;
    std::string originalCoefficientUnits = "mm";
    std::string originalCoefficientOrder = "radial,north,east; cosS1,sinS1,cosS2,sinS2";
};

//! Regular provider grid; axes and both seam endpoints are retained unchanged, without interpolation.
struct AtmosphericTideGrid {
    AtmosphericTideMetadata metadata;
    std::vector< double > longitudesDegrees;
    std::vector< double > latitudesDegrees;
    //! Longitude-major, latitude-fast. Rows radial/north/east, columns cosS1/sinS1/cosS2/sinS2, METRES.
    std::vector< Eigen::Matrix< double, 3, 4 > > coefficientsMetres;
};
struct AtmosphericTideStation {
    std::string stationIdentifier;
    std::string providerStationIdentifierLine;
    Eigen::Matrix< double, 3, 4 > coefficientsMetres = Eigen::Matrix< double, 3, 4 >::Zero( );
    AtmosphericTideMetadata metadata;
    bool hasProviderCoordinates = false;
    Eigen::Vector2d providerLongitudeLatitudeDegrees = Eigen::Vector2d::Zero( );
};
struct AtmosphericTideStationData {
    AtmosphericTideMetadata metadata;
    std::map< std::string, AtmosphericTideStation > stations;
    const AtmosphericTideStation& getStation( const std::string& identifier ) const;
};

void validateAtmosphericTideGrid( const AtmosphericTideGrid& grid );
void validateAtmosphericTideStationData( const AtmosphericTideStationData& data );

//! Read Bernese 14-column RP03-format grids. Millimetres become metres; headers/model/origin/units remain metadata.
AtmosphericTideGrid readAtmosphericTideGrid( std::istream& stream, const AtmosphericTideMetadata& metadata = {} );
AtmosphericTideGrid readAtmosphericTideGrid( const std::string& path, const AtmosphericTideMetadata& metadata = {} );
//! Read Bernese ATL station blocks: identifier then three rows of four signed millimetre coefficients.
AtmosphericTideStationData readAtmosphericTideStationData( std::istream& stream, const AtmosphericTideMetadata& metadata = {} );
AtmosphericTideStationData readAtmosphericTideStationData( const std::string& path, const AtmosphericTideMetadata& metadata = {} );
}  // namespace input_output
}  // namespace tudat
#endif
