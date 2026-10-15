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
//! Origin of the reference frame of atmospheric tide loading coefficients.
enum class AtmosphericLoadingReferenceOrigin { unspecified, solidEarthCentre, wholeEarthCentreOfMass };

//! Type of latitude used to look up a station in an atmospheric tide grid.
enum class AtmosphericTideLatitudeConvention { unspecified, geocentric, geodetic };

//! Directions of the radial (or up), north and east components of the coefficients in the file.
enum class AtmosphericTideDisplacementBasis { unspecified, geocentricRadialNorthEast, geodeticUpNorthEast };

//! Known product whose conventions are set explicitly by the user (a model name or file name alone does not select it).
enum class AtmosphericTideSourceProfile { unspecified, iers2010Rp03NoIbFarrell };

//! Description of an atmospheric tide loading file. The origin and latitude convention are never derived from the file name.
struct AtmosphericTideMetadata {
    //! Name of the file (or other source) and its header lines.
    std::string source;
    std::vector< std::string > sourceHeaders;

    //! Name of the atmospheric tide model and of the Green's functions used to compute the loading.
    std::string modelName;
    std::string greensFunction;

    //! Origin of the reference frame, and the latitude convention of the grid.
    AtmosphericLoadingReferenceOrigin referenceOrigin = AtmosphericLoadingReferenceOrigin::unspecified;
    AtmosphericTideLatitudeConvention latitudeConvention = AtmosphericTideLatitudeConvention::unspecified;

    //! Optional expected dimensions for source formats with no count header; zero means infer the provided extent.
    std::size_t expectedLongitudeCount = 0;
    std::size_t expectedLatitudeCount = 0;

    //! Directions of the displacement components, and the known product (if any) that these settings describe.
    AtmosphericTideDisplacementBasis displacementBasis = AtmosphericTideDisplacementBasis::unspecified;
    AtmosphericTideSourceProfile sourceProfile = AtmosphericTideSourceProfile::unspecified;

    //! Units and order of the coefficients in the file.
    std::string originalCoefficientUnits = "mm";
    std::string originalCoefficientOrder = "radial,north,east; cosS1,sinS1,cosS2,sinS2";
};

//! Regular grid of atmospheric tide loading coefficients, as in the file (including both copies of the 0/360 degree seam).
struct AtmosphericTideGrid {
    //! Description of the file.
    AtmosphericTideMetadata metadata;

    //! Longitudes and latitudes (degrees) of the grid points.
    std::vector< double > longitudesDegrees;
    std::vector< double > latitudesDegrees;

    //! Longitude-major, latitude-fast. Rows radial/north/east, columns cosS1/sinS1/cosS2/sinS2, METRES.
    std::vector< Eigen::Matrix< double, 3, 4 > > coefficientsMetres;
};

//! Atmospheric tide loading coefficients of a single station.
struct AtmosphericTideStation {
    //! Station name, and the line of the file that identifies the station.
    std::string stationIdentifier;
    std::string providerStationIdentifierLine;

    //! Coefficients (metres): rows radial, north, east; columns cos S1, sin S1, cos S2, sin S2.
    Eigen::Matrix< double, 3, 4 > coefficientsMetres = Eigen::Matrix< double, 3, 4 >::Zero( );

    //! Description of the file.
    AtmosphericTideMetadata metadata;

    //! Longitude and latitude (degrees) of the station, if the file gives them.
    bool hasProviderCoordinates = false;
    Eigen::Vector2d providerLongitudeLatitudeDegrees = Eigen::Vector2d::Zero( );
};

//! Atmospheric tide loading coefficients of a set of stations, as read from a file.
struct AtmosphericTideStationData {
    //! Description of the file.
    AtmosphericTideMetadata metadata;

    //! Stations, with their name as key.
    std::map< std::string, AtmosphericTideStation > stations;

    //! Return the exact provider identifier; absent stations raise std::out_of_range.
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
