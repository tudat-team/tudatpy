#ifndef TUDAT_READ_OCEAN_LOADING_DATA_H
#define TUDAT_READ_OCEAN_LOADING_DATA_H

#include <array>
#include <istream>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! Origin convention explicitly declared by an ocean loading provider.
enum class OceanLoadingReferenceOrigin { unspecified, solidEarthCentre, wholeEarthCentreOfMass };

//! BLQ coefficients in provider units and ordering; independent of dynamics and estimation.
struct OceanLoadingStation {
    std::string stationIdentifier;

    //! Rows: up/radial, west, south. Columns: M2,S2,N2,K2,K1,O1,P1,Q1,Mf,Mm,Ssa.
    Eigen::Matrix< double, 3, 11 > amplitudesMetres = Eigen::Matrix< double, 3, 11 >::Zero( );

    //! Positive Greenwich phase lags in degrees, in the same row/column order as amplitudes.
    Eigen::Matrix< double, 3, 11 > phaseLagsDegrees = Eigen::Matrix< double, 3, 11 >::Zero( );
    std::string source;
    std::vector< std::string > sourceHeaders;
    std::string oceanTideModel;
    std::string greensFunction;
    OceanLoadingReferenceOrigin referenceOrigin = OceanLoadingReferenceOrigin::unspecified;

    //! Provider lon/lat in degrees and height in metres, without reinterpretation of the source values.
    Eigen::Vector3d providerLongitudeLatitudeHeight = Eigen::Vector3d::Zero( );
    bool hasProviderCoordinates = false;
};

//! A BLQ data set with source provenance and validated station records.
struct OceanLoadingData {
    std::string source;
    std::vector< std::string > sourceHeaders;
    std::map< std::string, OceanLoadingStation > stations;

    //! Exact identifier lookup; throws std::out_of_range if the requested station is absent.
    const OceanLoadingStation& getStation( const std::string& stationIdentifier ) const;
};

//! Canonical BLQ column ordering; aliases MF, MM and SSA in headers are case insensitive.
const std::array< std::string, 11 >& getOceanLoadingConstituentNames( );

//! Validate header conventions, finite coefficients, nonnegative amplitudes and station identities.
//! Constructors may explicitly allow an empty station set for later append; readers and append sources reject it.
void validateOceanLoadingData( const OceanLoadingData& data, bool allowEmpty = false );

//! Append validated station records and provider headers, rejecting duplicates before changing the destination.
//! An empty destination is allowed; each appended source must contain valid station records.
void appendOceanLoadingData( OceanLoadingData& destination, const OceanLoadingData& source );

//! Read a BLQ stream, preserving header provenance and U/W/S metre/degree conventions.
/*!
 * Variable numbers of comment lines and blank lines, CRLF and scientific notation are supported.
 * Each station must have exactly six rows of eleven values. Duplicate identifiers, incomplete
 * blocks, nonfinite values and explicitly incompatible column ordering or units are rejected.
 * CMC:NO means solid-Earth centre (CE); CMC:YES means the joint Earth/ocean mass centre (CM).
 * An absent origin declaration remains unspecified.
 */
OceanLoadingData readOceanLoadingData( std::istream& stream, const std::string& source = "<stream>" );

//! Filename adapter for the stream parser. No input is required during subsequent model evaluation.
OceanLoadingData readOceanLoadingData( const std::string& fileName );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READ_OCEAN_LOADING_DATA_H
