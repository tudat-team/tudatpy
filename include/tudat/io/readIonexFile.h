#ifndef TUDAT_READIONEXFILE_H
#define TUDAT_READIONEXFILE_H

#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>
#include <iostream>

namespace tudat
{

namespace input_output
{
// https://cddis.nasa.gov/archive/gnss/products/ionex
//! Struct to store a parsed IONEX TEC map and metadata
struct IonexTecMap
{
public:
    IonexTecMap( ) = default;

    // Data grid
    std::vector< double > epochs;
    std::vector< double > latitudes;
    std::vector< double > longitudes;
    std::map< double, Eigen::MatrixXd > tecMaps;

    // Metadata from file
    double epochStart = -1.0;
    double epochEnd = -1.0;

    double latMin = -90.0;
    double latMax = 90.0;
    double dLat = 0.0;

    double lonMin = -180.0;
    double lonMax = 180.0;
    double dLon = 0.0;

    double hgtMin = 450.0;
    double hgtMax = 450.0;
    double dHgt = 0.0;
    double referenceIonosphereHeight_ = 450.0e3;
    
    std::string solutionType;  // e.g., "RAP", "FIN"
    std::string samplingInterval; // e.g., "01H", "02H"
    std::string contentType; // GIM, ROT

    //! Validate consistency
    void validate( ) const
    {
        for ( const auto& entry : tecMaps )
        {
            const Eigen::MatrixXd& tec = entry.second;
            if ( static_cast< std::size_t >( tec.rows( ) ) != latitudes.size( ) ||
                 static_cast< std::size_t >( tec.cols( ) ) != longitudes.size( ) )
            {
                throw std::runtime_error( "IONEX: Data matrix size inconsistent with lat/lon grid." );
            }
        }
        if ( epochs.size( ) != tecMaps.size( ) )
        {
            throw std::runtime_error( "IONEX: Mismatch between epochs and stored TEC maps." );
        }
    }

    //! Print metadata to console
    void printMetadata( ) const
    {
        std::cout << std::fixed;
        std::cout << "----------------------------------------\n";
        std::cout << "Parsed IONEX Metadata Summary:\n";
        std::cout << "  EPOCH OF FIRST MAP: " << epochStart << " [s since J2000]\n";
        std::cout << "  EPOCH OF LAST  MAP: " << epochEnd << " [s since J2000]\n";

        std::cout << "  LAT grid: [" << latMax << ", " << latMin << "], dLat = " << dLat << " ("
                  << latitudes.size( ) << " points)\n";
        std::cout << "  LON grid: [" << lonMin << ", " << lonMax << "], dLon = " << dLon << " ("
                  << longitudes.size( ) << " points)\n";
        std::cout << "  HGT: " << hgtMin << " to " << hgtMax << " km, step = " << dHgt << "\n";

        std::cout << "  Sampling: " << samplingInterval << ", Content: " << contentType
                  << ", Solution: " << solutionType << "\n";
        std::cout << "  # Epochs parsed: " << epochs.size( ) << "\n";
        std::cout << "----------------------------------------\n";
    }
};

//! Read single IONEX file and fill TEC data
void readIonexFile( const std::string& filePath, IonexTecMap& data );

//! Read and merge multiple IONEX files
void readIonexFiles( const std::vector< std::string >& filePaths, IonexTecMap& data );

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READIONEXFILE_H
