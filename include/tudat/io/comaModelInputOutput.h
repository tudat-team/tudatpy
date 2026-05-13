/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_COMA_MODEL_INPUT_OUTPUT_H
#define TUDAT_COMA_MODEL_INPUT_OUTPUT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/time_period.hpp>
#include <boost/filesystem.hpp>
#include <boost/variant.hpp>

#include "tudat/basics/identityElements.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/interpolators/interpolator.h"

namespace tudat
{

namespace simulation_setup
{

// ---- coefficient storage: rows = coeffs, cols = {C,S} ----
using StokesBlock = Eigen::Matrix< double, Eigen::Dynamic, 2, Eigen::RowMajor >;

struct FileMeta
{
    // Epochs in seconds since J2000 (or your preferred convention)
    double start_epoch{};
    double end_epoch{};
    std::string source_tag; // filename / description
};

// Degree-major mapping (n = 0..nmax, m = 0..n)
inline std::size_t nm_to_index_deg_major( int n, int m )
{
    // index(n,m) = n(n+1)/2 + m
    if(n < 0 || m < 0 || m > n) throw std::out_of_range( "nm_to_index: invalid (n,m)" );
    return static_cast< std::size_t >(n * ( n + 1 ) / 2 + m);
}

inline int index_to_n_deg_major( std::size_t k )
{
    // invert k = n(n+1)/2 + m  => n = floor((sqrt(8k+1)-1)/2)
    const double nd = ( std::sqrt( 8.0 * static_cast< double >(k) + 1.0 ) - 1.0 ) * 0.5;
    return static_cast< int >(std::floor( nd + 1e-12 ));
}

inline std::pair< int, int > index_to_nm_deg_major( std::size_t k )
{
    const int n = index_to_n_deg_major( k );
    const std::size_t base = static_cast< std::size_t >(n * ( n + 1 ) / 2);
    const int m = static_cast< int >(k - base);
    return { n, m };
}

// ============= Core Data Models (Pure data storage) =============

// ---- Stokes Dataset ----
class ComaStokesDataset
{
public:
    struct FileMeta
    {
        double start_epoch{};
        double end_epoch{};
        std::string source_tag;
        double referenceRadius{};
    };

    // Factory method
    static ComaStokesDataset create(
            std::vector< FileMeta > files,
            std::vector< double > radii,
            std::vector< double > solLongs,
            int nmax,
            bool computeReducedCoeffs = true )
    {
        if(files.empty( ) || radii.empty( ) || solLongs.empty( ) || nmax < 0)
            throw std::runtime_error( "StokesDataset: invalid metadata." );

        // Filter radii to only include those <= reference radius from all files
        // Find the minimum reference radius across all files
        double minReferenceRadius = std::numeric_limits< double >::max( );
        for(const auto& file: files)
        {
            if(file.referenceRadius > 1e-10) // Only consider valid reference radii
            {
                minReferenceRadius = std::min( minReferenceRadius, file.referenceRadius );
            }
        }

        // If no valid reference radius found, use all radii (backward compatibility)
        std::vector< double > validRadii;
        if(minReferenceRadius == std::numeric_limits< double >::max( ) || minReferenceRadius < 1e-10)
        {
            // No valid reference radius found - use all radii (backward compatibility for tests)
            validRadii = radii;
        }
        else
        {
            // Filter valid radii (those <= minimum reference radius)
            std::vector< double > discardedRadii;
            for(const auto& radius: radii)
            {
                if(radius <= minReferenceRadius)
                {
                    validRadii.push_back( radius );
                }
                else
                {
                    discardedRadii.push_back( radius );
                }
            }

            // Warn about discarded radii
            if(!discardedRadii.empty( ))
            {
                std::cerr << "Warning: The following radii exceed the reference radius ("
                        << minReferenceRadius << " m) and will be discarded from Stokes coefficient computation: ";
                for(size_t i = 0; i < discardedRadii.size( ); ++i)
                {
                    std::cerr << discardedRadii[ i ] << " m";
                    if(i < discardedRadii.size( ) - 1) std::cerr << ", ";
                }
                std::cerr << ". For radii beyond the reference radius, reduced coefficients and decay terms will be used during runtime." <<
                        std::endl;
            }

            // Ensure reference radius is included in validRadii
            bool referenceRadiusPresent = false;
            for(const auto& radius: validRadii)
            {
                if(std::abs( radius - minReferenceRadius ) < 1e-10)
                {
                    referenceRadiusPresent = true;
                    break;
                }
            }
            if(validRadii.empty( ))
            {
                throw std::runtime_error( "StokesDataset: No valid radii found (all radii exceed reference radius)." );
            }
            if(!referenceRadiusPresent)
            {
                std::cerr << "Warning: Reference radius (" << minReferenceRadius
                        << " m) was not present in the provided radii vector. Adding it automatically for proper interpolation." <<
                        std::endl;
                validRadii.push_back( minReferenceRadius );
                std::sort( validRadii.begin( ), validRadii.end( ) );
            }
        }

        ComaStokesDataset g;
        g.files_ = std::move( files );
        g.radii_ = std::move( validRadii );
        g.lons_ = std::move( solLongs );
        g.nmax_ = nmax;
        g.hasReducedCoeffs_ = computeReducedCoeffs;

        g.n_files_ = g.files_.size( );
        g.n_radii_ = g.radii_.size( );
        g.n_lons_ = g.lons_.size( );
        g.n_coeffs_ = static_cast< std::size_t >(( nmax + 1 ) * ( nmax + 2 ) / 2);

        const std::size_t totalRows = g.n_files_ * g.n_radii_ * g.n_lons_ * g.n_coeffs_;
        g.data_.setZero( static_cast< Eigen::Index >(totalRows), 2 );

        // Only allocate reduced coefficient storage if needed
        if( computeReducedCoeffs )
        {
            const std::size_t reducedTotalRows = g.n_files_ * g.n_lons_ * g.n_coeffs_;
            g.reduced_data_.setZero( static_cast< Eigen::Index >(reducedTotalRows), 2 );
        }

        return g;
    }

    // Pure data accessors
    std::size_t nFiles( ) const
    {
        return n_files_;
    }

    std::size_t nRadii( ) const
    {
        return n_radii_;
    }

    std::size_t nLongitudes( ) const
    {
        return n_lons_;
    }

    std::size_t nCoeffs( ) const
    {
        return n_coeffs_;
    }

    int nmax( ) const
    {
        return nmax_;
    }

    const std::vector< double >& radii( ) const
    {
        return radii_;
    }

    const std::vector< double >& lons( ) const
    {
        return lons_;
    }

    const std::vector< FileMeta >& files( ) const
    {
        return files_;
    }

    bool hasReducedCoeffs( ) const
    {
        return hasReducedCoeffs_;
    }

    // Convenience accessor for reference radius
    double getReferenceRadius( std::size_t f = 0 ) const
    {
        if(f >= n_files_) throw std::out_of_range( "File index out of range" );
        return files_[ f ].referenceRadius;
    }

    // Block access
    auto block( std::size_t f, std::size_t r, std::size_t l )
    {
        const std::size_t start = startRow_( f, r, l );
        return data_.block( static_cast< Eigen::Index >(start),
                            0,
                            static_cast< Eigen::Index >(n_coeffs_),
                            2 );
    }

    auto block( std::size_t f, std::size_t r, std::size_t l ) const
    {
        const std::size_t start = startRow_( f, r, l );
        return data_.block( static_cast< Eigen::Index >(start),
                            0,
                            static_cast< Eigen::Index >(n_coeffs_),
                            2 );
    }

    // Coefficient matrices
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
    getCoefficientMatrices( std::size_t f, std::size_t r, std::size_t l ) const
    {
        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( nmax_ + 1, nmax_ + 1 );
        Eigen::MatrixXd sine = Eigen::MatrixXd::Zero( nmax_ + 1, nmax_ + 1 );

        const auto blk = block( f, r, l );
        for(std::size_t k = 0; k < n_coeffs_; ++k)
        {
            auto nm = index_to_nm_deg_major( k );
            int n = nm.first;
            int m = nm.second;
            cosine( n, m ) = blk( static_cast< Eigen::Index >(k), 0 );
            sine( n, m ) = blk( static_cast< Eigen::Index >(k), 1 );
        }
        return { cosine, sine };
    }

    // Single coefficient access
    void setCoeff( std::size_t f,
                   std::size_t r,
                   std::size_t l,
                   int n,
                   int m,
                   double C,
                   double S )
    {
        const std::size_t k = nm_to_index_deg_major( n, m );
        if(k >= n_coeffs_) throw std::out_of_range( "setCoeff: (n,m) exceeds nmax" );
        const std::size_t row = startRow_( f, r, l ) + k;
        data_( static_cast< Eigen::Index >(row), 0 ) = C;
        data_( static_cast< Eigen::Index >(row), 1 ) = S;
    }

    std::pair< double, double > getCoeff( std::size_t f,
                                          std::size_t r,
                                          std::size_t l,
                                          int n,
                                          int m ) const
    {
        const std::size_t k = nm_to_index_deg_major( n, m );
        if(k >= n_coeffs_) throw std::out_of_range( "getCoeff: (n,m) exceeds nmax" );
        const std::size_t row = startRow_( f, r, l ) + k;
        return { data_( static_cast< Eigen::Index >(row), 0 ),
                 data_( static_cast< Eigen::Index >(row), 1 ) };
    }

    const StokesBlock& data( ) const
    {
        return data_;
    }

    StokesBlock& data( )
    {
        return data_;
    }

    // Reduced coefficient access methods (no radius dimension)
    auto reducedBlock( std::size_t f, std::size_t l )
    {
        const std::size_t start = reducedStartRow_( f, l );
        return reduced_data_.block( static_cast< Eigen::Index >(start),
                                    0,
                                    static_cast< Eigen::Index >(n_coeffs_),
                                    2 );
    }

    auto reducedBlock( std::size_t f, std::size_t l ) const
    {
        const std::size_t start = reducedStartRow_( f, l );
        return reduced_data_.block( static_cast< Eigen::Index >(start),
                                    0,
                                    static_cast< Eigen::Index >(n_coeffs_),
                                    2 );
    }

    void setReducedCoeff( std::size_t f, std::size_t l, int n, int m, double C, double S )
    {
        const std::size_t k = nm_to_index_deg_major( n, m );
        if(k >= n_coeffs_) throw std::out_of_range( "setReducedCoeff: (n,m) exceeds nmax" );
        const std::size_t row = reducedStartRow_( f, l ) + k;
        reduced_data_( static_cast< Eigen::Index >(row), 0 ) = C;
        reduced_data_( static_cast< Eigen::Index >(row), 1 ) = S;
    }

    std::pair< double, double > getReducedCoeff( std::size_t f, std::size_t l, int n, int m ) const
    {
        const std::size_t k = nm_to_index_deg_major( n, m );
        if(k >= n_coeffs_) throw std::out_of_range( "getReducedCoeff: (n,m) exceeds nmax" );
        const std::size_t row = reducedStartRow_( f, l ) + k;
        return { reduced_data_( static_cast< Eigen::Index >(row), 0 ),
                 reduced_data_( static_cast< Eigen::Index >(row), 1 ) };
    }

    const StokesBlock& reducedData( ) const
    {
        return reduced_data_;
    }

    StokesBlock& reducedData( )
    {
        return reduced_data_;
    }

private:
    std::size_t startRow_( std::size_t f, std::size_t r, std::size_t l ) const
    {
        if(f >= n_files_ || r >= n_radii_ || l >= n_lons_)
            throw std::out_of_range( "StokesDataset: index out of range." );
        const std::size_t cell = ( ( f * n_radii_ ) + r ) * n_lons_ + l;
        return cell * n_coeffs_;
    }

    std::size_t reducedStartRow_( std::size_t f, std::size_t l ) const
    {
        if(f >= n_files_ || l >= n_lons_)
            throw std::out_of_range( "StokesDataset: reduced index out of range." );
        const std::size_t cell = f * n_lons_ + l;
        return cell * n_coeffs_;
    }

    std::vector< FileMeta > files_;
    std::vector< double > radii_;
    std::vector< double > lons_;
    int nmax_{};
    std::size_t n_files_{}, n_radii_{}, n_lons_{}, n_coeffs_{};
    bool hasReducedCoeffs_{ true }; // Flag indicating if reduced coefficients are computed/stored
    StokesBlock data_;
    StokesBlock reduced_data_; // For reduced coefficients (n_files × n_lons × n_coeffs)
};

// ---- Poly Dataset ----
class ComaPolyDataset
{
public:
    struct FileMeta
    {
        double referenceRadius{};
        Eigen::ArrayXd powersInvRadius;
        std::vector< std::pair< double, double > > timePeriods;
        int maxDegreeSH{};
        Eigen::Index numRadialTerms{};
        Eigen::Index numIntervals{};
        std::string sourcePath;
        std::string startDateString; // Store start date string for reference
        std::string endDateString; // Store end date string for reference
    };

    // Simple accessors
    std::size_t getNumFiles( ) const
    {
        return numPolyCoefFiles_;
    }

    const FileMeta& getFileMeta( std::size_t f ) const
    {
        boundsCheck_( f );
        return fileMeta_[ f ];
    }

    const Eigen::MatrixXd& getPolyCoefficients( std::size_t f ) const
    {
        boundsCheck_( f );
        return polyCoefficients_[ f ];
    }

    const Eigen::ArrayXXi& getSHDegreeAndOrderIndices( std::size_t f ) const
    {
        boundsCheck_( f );
        return SHDegreeAndOrderIndices_[ f ];
    }

    // Convenience accessors
    double getReferenceRadius( std::size_t f ) const
    {
        boundsCheck_( f );
        return fileMeta_[ f ].referenceRadius;
    }

    const Eigen::ArrayXd& getPowersInvRadius( std::size_t f ) const
    {
        boundsCheck_( f );
        return fileMeta_[ f ].powersInvRadius;
    }

    int getMaxDegreeSH( std::size_t f ) const
    {
        boundsCheck_( f );
        return fileMeta_[ f ].maxDegreeSH;
    }

    // Column access methods
    Eigen::VectorXd columnForNM( std::size_t f, int n, int m ) const
    {
        auto result = findColumn_( f, n, m );
        bool ok = result.first;
        int col = result.second;
        if(!ok) throw std::out_of_range( "columnForNM: (n,m) not found" );
        return polyCoefficients_[ f ].col( col );
    }

    double value( std::size_t f, Eigen::Index termIndex, int n, int m ) const
    {
        auto result = findColumn_( f, n, m );
        bool ok = result.first;
        int col = result.second;
        if(!ok) throw std::out_of_range( "value: (n,m) not found" );
        if(termIndex < 0 || termIndex >= polyCoefficients_[ f ].rows( ))
            throw std::out_of_range( "value: termIndex out of range" );
        return polyCoefficients_[ f ]( termIndex, col );
    }

    void clear( )
    {
        numPolyCoefFiles_ = 0;
        polyCoefficients_.clear( );
        SHDegreeAndOrderIndices_.clear( );
        fileMeta_.clear( );
        nmToColCache_.clear( );
    }

protected:
    // Allow friend classes to populate data
    friend class ComaPolyDatasetReader;

    void setData( std::size_t numFiles,
                  std::vector< Eigen::MatrixXd > polyCoeffs,
                  std::vector< Eigen::ArrayXXi > shIndices,
                  std::vector< FileMeta > meta )
    {
        numPolyCoefFiles_ = numFiles;
        polyCoefficients_ = std::move( polyCoeffs );
        SHDegreeAndOrderIndices_ = std::move( shIndices );
        fileMeta_ = std::move( meta );

        // Build caches
        nmToColCache_.resize( numFiles );
        for(std::size_t f = 0; f < numFiles; ++f)
            buildNmMap_( f );
    }

private:
    struct PairHash
    {
        std::size_t operator()( const std::pair< int, int >& p ) const noexcept
        {
            return ( static_cast< std::size_t >(p.first) << 32 ) ^
                    static_cast< std::size_t >(p.second);
        }
    };

    void boundsCheck_( std::size_t f ) const
    {
        if(f >= numPolyCoefFiles_)
            throw std::out_of_range( "file index out of range" );
    }

    void buildNmMap_( std::size_t f )
    {
        nmToColCache_[ f ].clear( );
        const auto& sh = SHDegreeAndOrderIndices_[ f ];
        for(Eigen::Index c = 0; c < sh.cols( ); ++c)
        {
            int n = sh( 0, c );
            int m = sh( 1, c );
            nmToColCache_[ f ][ { n, m } ] = static_cast< int >(c);
        }
    }

    std::pair< bool, int > findColumn_( std::size_t f, int n, int m ) const
    {
        boundsCheck_( f );
        const auto& map = nmToColCache_[ f ];
        auto it = map.find( { n, m } );
        if(it == map.end( )) return { false, -1 };
        return { true, it->second };
    }

    std::size_t numPolyCoefFiles_{ 0 };
    std::vector< Eigen::MatrixXd > polyCoefficients_;
    std::vector< Eigen::ArrayXXi > SHDegreeAndOrderIndices_;
    std::vector< FileMeta > fileMeta_;
    std::vector< std::unordered_map< std::pair< int, int >, int, PairHash > > nmToColCache_;
};

// ============= I/O Components (Separate from data) =============

// ---- Reader for Poly Coefficients ----
class ComaPolyDatasetReader
{
public:
    static ComaPolyDataset readFromFiles( const std::vector< std::string >& filePaths )
    {


        // Check that all files exist before attempting to read
        std::vector< std::string > missingFiles;
        for(const auto& path : filePaths)
        {
            std::ifstream testFile( path );
            const bool canOpen = testFile.is_open( );
            testFile.close( );
            if(!canOpen)
                missingFiles.push_back( path );
        }

        if(!missingFiles.empty( ))
        {
            std::string errorMsg = "ComaPolyDatasetReader: the following file(s) do not exist or cannot be opened:\n";
            for(const auto& path : missingFiles)
                errorMsg += "  - " + path + "\n";
            throw std::runtime_error( errorMsg );
        }

        if(filePaths.empty( ))
                    throw std::invalid_argument( "ComaPolyDatasetReader: empty file list" );
        const std::size_t n = filePaths.size( );
        std::vector< Eigen::MatrixXd > polyCoefficients( n );
        std::vector< Eigen::ArrayXXi > SHDegreeAndOrderIndices( n );
        std::vector< ComaPolyDataset::FileMeta > fileMeta( n );

        for(std::size_t fileIdx = 0; fileIdx < n; ++fileIdx)
        {
            readSingleFile( filePaths[ fileIdx ],
                            fileIdx,
                            polyCoefficients,
                            SHDegreeAndOrderIndices,
                            fileMeta );
        }

        ComaPolyDataset dataset;
        dataset.setData( n,
                         std::move( polyCoefficients ),
                         std::move( SHDegreeAndOrderIndices ),
                         std::move( fileMeta ) );
        return dataset;
    }

private:
    static void readSingleFile( const std::string& filePath,
                                std::size_t fileIdx,
                                std::vector< Eigen::MatrixXd >& polyCoefficients,
                                std::vector< Eigen::ArrayXXi >& SHDegreeAndOrderIndices,
                                std::vector< ComaPolyDataset::FileMeta >& fileMeta )
    {
        // ===== Implementation adapted from old ComaPolyDataset::readInputFiles =====
        fileMeta[ fileIdx ].sourcePath = filePath;

        std::ifstream file( filePath );
        if(!file.is_open( ))
        {
            throw std::runtime_error( "ComaPolyDatasetReader: could not open file '" + filePath + "'" );
        }

        std::string line;
        std::vector< std::string > tokens;

        int maxDegreeSH = 0;
        Eigen::Index numTerms = 0, numCoefs = 0, numRadialTerms = 0, numIntervals = 0;
        Eigen::ArrayXd powers;

        // ----- Parse header -----
        while(std::getline( file, line ))
        {
            if(line.empty( )) continue;
            if(line[ 0 ] != '#') break;

            std::string headerLine = line.substr( 1 );
            boost::trim( headerLine );
            boost::split( tokens, headerLine, boost::is_any_of( ", \t" ), boost::token_compress_on );
            if(tokens.empty( )) continue;

            const std::string& key = tokens[ 0 ];

            if(boost::iequals( key, "N(SH)" ))
            {
                maxDegreeSH = std::stoi( line.substr( line.find_last_of( " \t" ) + 1 ) );
                numCoefs = ( maxDegreeSH + 1 ) * ( maxDegreeSH + 1 );
            }
            else if(boost::icontains( key, "PWRS" ))
            {
                std::string tail = line.substr( line.find( "PWRS" ) );
                boost::trim( tail );
                std::vector< std::string > pwrtok;
                boost::split( pwrtok, tail, boost::is_any_of( ", \t" ), boost::token_compress_on );
                std::size_t start = ( !pwrtok.empty( ) &&
                            !std::all_of( pwrtok[ 0 ].begin( ), pwrtok[ 0 ].end( ), ::isdigit ) )
                        ? 1
                        : 0;
                auto count = static_cast< Eigen::Index >(pwrtok.size( ) - start);
                powers.resize( count );
                for(Eigen::Index j = 0; j < count; ++j)
                    powers[ j ] = std::stod( pwrtok[ start + j ] );
            }
            else if(boost::iequals( key, "R" ))
            {
                double R_km = std::stod( line.substr( line.find_last_of( " \t" ) + 1 ) );
                fileMeta[ fileIdx ].referenceRadius = R_km * 1000.0; // Convert km to meters
            }
            else if(line.find( "N(r)" ) != std::string::npos && line.find( "N(T)" ) != std::string::npos)
            {
                std::string content = line.substr( 1 ); // strip '#'
                boost::trim( content );
                boost::split( tokens, content, boost::is_any_of( ", \t" ), boost::token_compress_on );

                if(tokens.size( ) >= 2)
                {
                    int a = std::stoi( tokens[ tokens.size( ) - 2 ] );
                    int b = std::stoi( tokens[ tokens.size( ) - 1 ] );

                    if(line.find( "N(r)" ) < line.find( "N(T)" ))
                    {
                        numRadialTerms = a;
                        numIntervals = b;
                    }
                    else
                    {
                        numIntervals = a;
                        numRadialTerms = b;
                    }
                    numTerms = numRadialTerms * numIntervals;
                }
                else
                {
                    std::cerr << "[ERROR] Failed to parse N(T)/N(r) line: " << line << std::endl;
                    std::exit( EXIT_FAILURE );
                }
            }
            else if(line.find( "Start Date:" ) != std::string::npos && line.find( "End Date:" ) != std::string::npos)
            {
                // Parse start and end dates from line like: "# Start Date: 2015/07/21    End Date: 2015/08/21"
                std::string content = line.substr( 1 ); // strip '#'
                boost::trim( content );

                std::size_t startDatePos = content.find( "Start Date:" );
                std::size_t endDatePos = content.find( "End Date:" );

                if(startDatePos != std::string::npos && endDatePos != std::string::npos)
                {
                    // Extract start date
                    std::string startDateStr = content.substr( startDatePos + 11 ); // Skip "Start Date:"
                    std::size_t endPosStart = startDateStr.find( "End Date:" );
                    if(endPosStart != std::string::npos)
                    {
                        startDateStr = startDateStr.substr( 0, endPosStart );
                    }
                    boost::trim( startDateStr );

                    // Extract end date
                    std::string endDateStr = content.substr( endDatePos + 9 ); // Skip "End Date:"
                    boost::trim( endDateStr );

                    // Store the date strings in metadata for reference
                    fileMeta[ fileIdx ].startDateString = startDateStr;
                    fileMeta[ fileIdx ].endDateString = endDateStr;
                }
            }
            else if(line.find( "Start J2000:" ) != std::string::npos && line.find( "End J2000:" ) != std::string::npos)
            {
                // Parse start and end J2000 times from line like: "# Start J2000: 5679.5      End J2000: 5710.5"
                std::string content = line.substr( 1 ); // strip '#'
                boost::trim( content );

                std::size_t startJ2000Pos = content.find( "Start J2000:" );
                std::size_t endJ2000Pos = content.find( "End J2000:" );

                if(startJ2000Pos != std::string::npos && endJ2000Pos != std::string::npos)
                {
                    // Extract start J2000
                    std::string startJ2000Str = content.substr( startJ2000Pos + 12 ); // Skip "Start J2000:"
                    std::size_t endPosStart = startJ2000Str.find( "End J2000:" );
                    if(endPosStart != std::string::npos)
                    {
                        startJ2000Str = startJ2000Str.substr( 0, endPosStart );
                    }
                    boost::trim( startJ2000Str );

                    // Extract end J2000
                    std::string endJ2000Str = content.substr( endJ2000Pos + 10 ); // Skip "End J2000:"
                    boost::trim( endJ2000Str );

                    try
                    {
                        double startJ2000 = std::stod( startJ2000Str );
                        double endJ2000 = std::stod( endJ2000Str );

                        // Add time period to metadata
                        fileMeta[ fileIdx ].timePeriods.emplace_back( startJ2000, endJ2000 );
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << "[WARNING] Failed to parse J2000 times in file: " << filePath
                                << " Error: " << e.what( ) << std::endl;
                    }
                }
            }
        }

        // ----- Validation -----
        if(numTerms <= 0 || numCoefs <= 0 || powers.size( ) == 0)
        {
            std::cerr << "[ERROR] Header parsing failed in file: " << filePath << std::endl;
            std::cerr << "  numTerms = " << numTerms
                    << "\n  numCoefs = " << numCoefs
                    << "\n  powersInvRadius.size() = " << powers.size( ) << std::endl;
            std::exit( EXIT_FAILURE );
        }

        fileMeta[ fileIdx ].maxDegreeSH = maxDegreeSH;
        fileMeta[ fileIdx ].numRadialTerms = numRadialTerms;
        fileMeta[ fileIdx ].numIntervals = numIntervals;
        fileMeta[ fileIdx ].powersInvRadius = powers;

        // --- Allocate containers ---
        Eigen::MatrixXd& currentPolyCoefficients = polyCoefficients[ fileIdx ];
        Eigen::ArrayXXi& currentShDegreeAndOrder = SHDegreeAndOrderIndices[ fileIdx ];

        currentPolyCoefficients.resize( numTerms, numCoefs );
        currentShDegreeAndOrder.resize( 2, numCoefs );

        // ----- Read poly coefficients data block -----
        Eigen::Index coefIndex = -1;
        do
        {
            boost::trim( line );
            if(line.empty( ) || line[ 0 ] == '#') continue;

            boost::split( tokens, line, boost::is_any_of( ", \t" ), boost::token_compress_on );
            if(static_cast< Eigen::Index >(tokens.size( )) == numTerms + 2)
            {
                ++coefIndex;
                currentShDegreeAndOrder( 0, coefIndex ) = std::stoi( tokens[ 0 ] ); // n
                currentShDegreeAndOrder( 1, coefIndex ) = std::stoi( tokens[ 1 ] ); // m
                for(Eigen::Index j = 0; j < numTerms; ++j)
                    currentPolyCoefficients( j, coefIndex ) = std::stod( tokens[ j + 2 ] );
            }
        } while(std::getline( file, line ));

        file.close( );
        // The (n,m)->col map is built by ComaPolyDataset::setData()
    }
};

// ---- Writer for Stokes Coefficients ----
class ComaStokesDatasetWriter
{
public:
    static void writeCsvForFile( const ComaStokesDataset& dataset,
                                 std::size_t f,
                                 const std::string& outputPath )
    {
        if(f >= dataset.nFiles( ))
            throw std::out_of_range( "writeCsvForFile: file index out of range" );

        auto csvEscape = []( std::ostream& os, const std::string& s ) {
            bool needs = false;
            for(char c: s)
            {
                if(c == ',' || c == '"' || c == '\n' || c == '\r')
                {
                    needs = true;
                    break;
                }
            }
            if(!needs)
            {
                os << s;
                return;
            }
            os << '"';
            for(char c: s) os << ( c == '"' ? "\"\"" : std::string( 1, c ) );
            os << '"';
        };

        boost::filesystem::path p( outputPath );
        if(p.has_parent_path( ))
            boost::filesystem::create_directories( p.parent_path( ) );

        std::ofstream os( outputPath, std::ios::binary );
        if(!os)
            throw std::runtime_error( "writeCsvForFile: cannot open " + outputPath );
        os.imbue( std::locale::classic( ) );

        const auto& fm = dataset.files( )[ f ];

        auto set_sci = [&] {
            os.setf( std::ios::scientific, std::ios::floatfield );
            os << std::setprecision( 17 );
        };
        auto set_def = [&] {
            os.setf( std::ios::fmtflags( 0 ), std::ios::floatfield );
            os << std::setprecision( 17 );
        };

        // Row 1: metadata
        os << "meta";
        set_sci( );
        os << ",start_epoch=" << fm.start_epoch;
        os << ",end_epoch=" << fm.end_epoch;
        os << ",max_degree=" << dataset.nmax( );
        os << ",max_order=" << dataset.nmax( );
        os << ",n_radii=" << dataset.nRadii( );
        os << ",n_lons=" << dataset.nLongitudes( );
        os << ",n_coeffs=" << dataset.nCoeffs( );
        os << ",has_reduced_coeffs=" << ( dataset.hasReducedCoeffs( ) ? "true" : "false" );
        os << ",source=";
        csvEscape( os, fm.source_tag );
        os << '\n';

        // Row 2: radii
        os << "radii [meter]";
        for(double r: dataset.radii( ))
        {
            set_def( );
            os << ',' << r;
        }
        os << '\n';

        // Row 3: longitudes
        os << "longitudes [degree]";
        for(double L: dataset.lons( ))
        {
            set_def( );
            os << ',' << L;
        }
        os << '\n';

        // REDUCED COEFFICIENTS SECTION (for radius > reference radius)
        // Only write if reduced coefficients are present (not needed for wind models)
        if( dataset.hasReducedCoeffs( ) )
        {
            os << "\n# REDUCED COEFFICIENTS SECTION (for radius > ref_radius)\n";
            os << "# These coefficients are computed using reducedToTemporalIFFT\n";
            os << "# No radius dimension - coefficients are independent of radius\n";

            for(std::size_t li = 0; li < dataset.nLongitudes( ); ++li)
            {
                os << "ID," << li << ',';
                set_def( );
                os << "l_0=" << dataset.lons( )[ li ] << '\n';

                os << "n,m,C,S\n";

                const auto blk = dataset.reducedBlock( f, li );
                for(std::size_t k = 0; k < dataset.nCoeffs( ); ++k)
                {
                    const auto nm = index_to_nm_deg_major( k );
                    const double C = blk( static_cast< Eigen::Index >(k), 0 );
                    const double S = blk( static_cast< Eigen::Index >(k), 1 );
                os << nm.first << ',' << nm.second << ',';
                set_sci( );
                os << C << ',' << S << '\n';
            }
        }
        } // End of if( dataset.hasReducedCoeffs( ) )

        // REGULAR COEFFICIENTS SECTION (for radius <= ref_radius)
        os << "\n# REGULAR COEFFICIENTS SECTION (for radius <= ref_radius)\n";
        os << "# These coefficients are computed using radialPolyvalAndTemporalIFFT\n";

        for(std::size_t ri = 0; ri < dataset.nRadii( ); ++ri)
        {
            for(std::size_t li = 0; li < dataset.nLongitudes( ); ++li)
            {
                const std::size_t block_id = ri * dataset.nLongitudes( ) + li;

                os << "ID," << block_id << ',';
                set_def( );
                os << "r_0=" << dataset.radii( )[ ri ] << ','
                        << "l_0=" << dataset.lons( )[ li ] << '\n';

                os << "n,m,C,S\n";

                const auto blk = dataset.block( f, ri, li );
                for(std::size_t k = 0; k < dataset.nCoeffs( ); ++k)
                {
                    const auto nm = index_to_nm_deg_major( k );
                    const double C = blk( static_cast< Eigen::Index >(k), 0 );
                    const double S = blk( static_cast< Eigen::Index >(k), 1 );
                    os << nm.first << ',' << nm.second << ',';
                    set_sci( );
                    os << C << ',' << S << '\n';
                }
            }
        }

        os.flush( );
        if(!os) throw std::runtime_error( "writeCsvForFile: write failed" );
    }

    static void writeCsvAll( const ComaStokesDataset& dataset,
                             const std::string& outputDir,
                             const std::string& prefix = "stokes" )
    {
        boost::filesystem::create_directories( outputDir );
        for(std::size_t f = 0; f < dataset.nFiles( ); ++f)
        {
            boost::filesystem::path path = boost::filesystem::path( outputDir ) /
                    ( prefix + "_file" + std::to_string( f ) + ".csv" );
            writeCsvForFile( dataset, f, path.string( ) );
        }
    }
};

// ---- Reader for Stokes Coefficients (future) ----
class ComaStokesDatasetReader
{
public:
    static ComaStokesDataset readFromCsv( const std::string& csvPath )
    {
        std::ifstream ifs( csvPath );
        if(!ifs)
            throw std::runtime_error( "readFromCsv: cannot open " + csvPath );

        // Parse metadata line
        std::string line;
        std::getline( ifs, line );
        if(line.empty( ) || line.substr( 0, 4 ) != "meta")
            throw std::runtime_error( "readFromCsv: invalid metadata line" );

        auto parseMeta = []( const std::string& metaLine ) -> std::tuple< double, double, int, int, int, int, int, bool, double, std::string > {
            std::istringstream ss( metaLine );
            std::string token;
            double start_epoch = 0, end_epoch = 0, ref_radius = 0.0;
            int max_degree = 0, max_order = 0, n_radii = 0, n_lons = 0, n_coeffs = 0;
            bool has_reduced_coeffs = false;
            std::string source;

            while(std::getline( ss, token, ',' ))
            {
                if(token.find( "start_epoch=" ) == 0)
                    start_epoch = std::stod( token.substr( 12 ) );
                else if(token.find( "end_epoch=" ) == 0)
                    end_epoch = std::stod( token.substr( 10 ) );
                else if(token.find( "max_degree=" ) == 0)
                    max_degree = std::stoi( token.substr( 11 ) );
                else if(token.find( "max_order=" ) == 0)
                    max_order = std::stoi( token.substr( 10 ) );
                else if(token.find( "n_radii=" ) == 0)
                    n_radii = std::stoi( token.substr( 8 ) );
                else if(token.find( "n_lons=" ) == 0)
                    n_lons = std::stoi( token.substr( 7 ) );
                else if(token.find( "n_coeffs=" ) == 0)
                    n_coeffs = std::stoi( token.substr( 9 ) );
                else if(token.find( "has_reduced_coeffs=" ) == 0)
                    has_reduced_coeffs = (token.substr( 19 ) == "true");
                else if(token.find( "ref_radius=" ) == 0)
                    ref_radius = std::stod( token.substr( 11 ) );
                else if(token.find( "source=" ) == 0)
                    source = token.substr( 7 );
            }
            return { start_epoch, end_epoch, max_degree, max_order, n_radii, n_lons, n_coeffs, has_reduced_coeffs, ref_radius, source };
        };

        auto meta = parseMeta( line );
        double start_epoch = std::get<0>( meta );
        double end_epoch = std::get<1>( meta );
        int max_degree = std::get<2>( meta );
        int max_order = std::get<3>( meta );
        int n_radii = std::get<4>( meta );
        int n_lons = std::get<5>( meta );
        int n_coeffs = std::get<6>( meta );
        bool has_reduced_coeffs = std::get<7>( meta );
        double ref_radius = std::get<8>( meta );
        std::string source = std::get<9>( meta );

        // Parse radii line
        std::getline( ifs, line );
        std::istringstream radiiStream( line );
        std::string token;
        std::getline( radiiStream, token, ',' ); // skip "radii [meter]"
        std::vector< double > radii;
        while(std::getline( radiiStream, token, ',' ))
            radii.push_back( std::stod( token ) );

        // Parse longitudes line
        std::getline( ifs, line );
        std::istringstream lonsStream( line );
        std::getline( lonsStream, token, ',' ); // skip "longitudes [degree]"
        std::vector< double > lons;
        while(std::getline( lonsStream, token, ',' ))
            lons.push_back( std::stod( token ) );

        // Create dataset with single file
        std::vector< ComaStokesDataset::FileMeta > files( 1 );
        files[ 0 ].start_epoch = start_epoch;
        files[ 0 ].end_epoch = end_epoch;
        files[ 0 ].source_tag = source;
        files[ 0 ].referenceRadius = ref_radius;

        ComaStokesDataset dataset = ComaStokesDataset::create( std::move( files ), radii, lons, max_degree );

        // Parse coefficient blocks (handles both reduced and regular coefficients in any order)
        while(std::getline( ifs, line ))
        {
            // Skip comment lines and empty lines
            if(line.empty() || line[0] == '#')
                continue;

            if(line.find( "ID," ) == 0)
            {
                // Parse ID line: check if it's reduced (only l_0) or regular (r_0 and l_0)
                bool isReduced = (line.find( "r_0=" ) == std::string::npos);

                std::istringstream idStream( line );
                std::string idToken;
                std::getline( idStream, idToken, ',' ); // skip "ID"
                std::getline( idStream, idToken, ',' ); // block_id
                int block_id = std::stoi( idToken );

                // Skip header line "n,m,C,S"
                std::getline( ifs, line );

                if(isReduced)
                {
                    // Reduced coefficients: ID is just the longitude index
                    std::size_t li = block_id;

                    for(int k = 0; k < n_coeffs; ++k)
                    {
                        std::getline( ifs, line );
                        std::istringstream coeffStream( line );
                        std::string nStr, mStr, cStr, sStr;
                        std::getline( coeffStream, nStr, ',' );
                        std::getline( coeffStream, mStr, ',' );
                        std::getline( coeffStream, cStr, ',' );
                        std::getline( coeffStream, sStr, ',' );

                        int n = std::stoi( nStr );
                        int m = std::stoi( mStr );
                        double C = std::stod( cStr );
                        double S = std::stod( sStr );

                        dataset.setReducedCoeff( 0, li, n, m, C, S );
                    }
                }
                else
                {
                    // Regular coefficients: ID is ri * n_lons + li
                    std::size_t ri = block_id / n_lons;
                    std::size_t li = block_id % n_lons;

                    for(int k = 0; k < n_coeffs; ++k)
                    {
                        std::getline( ifs, line );
                        std::istringstream coeffStream( line );
                        std::string nStr, mStr, cStr, sStr;
                        std::getline( coeffStream, nStr, ',' );
                        std::getline( coeffStream, mStr, ',' );
                        std::getline( coeffStream, cStr, ',' );
                        std::getline( coeffStream, sStr, ',' );

                        int n = std::stoi( nStr );
                        int m = std::stoi( mStr );
                        double C = std::stod( cStr );
                        double S = std::stod( sStr );

                        dataset.setCoeff( 0, ri, li, n, m, C, S );
                    }
                }
            }
        }

        return dataset;
    }

    static ComaStokesDataset readFromCsvFolder( const std::string& dir,
                                                const std::string& prefix = "stokes" )
    {
        namespace bf = boost::filesystem;

        if(!bf::exists( dir ) || !bf::is_directory( dir ))
            throw std::runtime_error( "readFromCsvFolder: directory does not exist: " + dir );

        // Find all CSV files with the given prefix
        std::vector< std::string > csvFiles;
        for(bf::directory_iterator it( dir ); it != bf::directory_iterator( ); ++it)
        {
            const std::string filename = it->path( ).filename( ).string( );
            if(filename.find( prefix + "_file" ) == 0 && filename.substr( filename.length( ) - 4 ) == ".csv")
                csvFiles.push_back( it->path( ).string( ) );
        }

        if(csvFiles.empty( ))
            throw std::runtime_error( "readFromCsvFolder: no CSV files found with prefix " + prefix );

        // Sort files by number
        std::sort( csvFiles.begin( ),
                   csvFiles.end( ),
                   [&prefix]( const std::string& a, const std::string& b ) {
                       auto extractNum = [&prefix]( const std::string& path ) -> int {
                           bf::path p( path );
                           std::string name = p.stem( ).string( );
                           std::string numStr = name.substr( prefix.length( ) + 5 ); // +5 for "_file"
                           return std::stoi( numStr );
                       };
                       return extractNum( a ) < extractNum( b );
                   } );

        // Read first file to get structure
        ComaStokesDataset firstDataset = readFromCsv( csvFiles[ 0 ] );

        if(csvFiles.size( ) == 1)
            return firstDataset;

        // Create multi-file dataset
        std::vector< ComaStokesDataset::FileMeta > allFiles;
        for(const auto& csvFile: csvFiles)
        {
            ComaStokesDataset singleDataset = readFromCsv( csvFile );
            allFiles.push_back( singleDataset.files( )[ 0 ] );
        }

        ComaStokesDataset multiDataset = ComaStokesDataset::create(
                std::move( allFiles ),
                firstDataset.radii( ),
                firstDataset.lons( ),
                firstDataset.nmax( ) );

        // Copy data from all files
        for(std::size_t f = 0; f < csvFiles.size( ); ++f)
        {
            ComaStokesDataset singleDataset = readFromCsv( csvFiles[ f ] );
            for(std::size_t ri = 0; ri < multiDataset.nRadii( ); ++ri)
            {
                for(std::size_t li = 0; li < multiDataset.nLongitudes( ); ++li)
                {
                    auto srcBlock = singleDataset.block( 0, ri, li );
                    auto destBlock = multiDataset.block( f, ri, li );
                    destBlock = srcBlock;
                }
            }
        }

        return multiDataset;
    }
};

// ============= Processing/Transformation Components =============

class StokesCoefficientsEvaluator
{
public:
private:
    // Cached frequency index arrays - computed once per unique numFrequencyTerms
    struct FrequencyCache
    {
        Eigen::ArrayXd cosineFrequencies;
        Eigen::ArrayXd sineFrequencies;
        int numCosTerms;
        int numSinTerms;
    };

    // Thread-safe cache using static local variable (initialized on first use)
    static const FrequencyCache& getFrequencyCache( const int numFrequencyTerms )
    {
        // Static map to cache frequency arrays for different numFrequencyTerms values
        static std::unordered_map< int, FrequencyCache > cache;

        // Check if already cached
        auto it = cache.find( numFrequencyTerms );
        if( it != cache.end( ) )
        {
            return it->second;
        }

        // Not cached - compute and store
        FrequencyCache newCache;

        if( numFrequencyTerms % 2 == 0 )
        {
            // Even: cos frequencies 0, 1, ..., N/2
            newCache.numCosTerms = numFrequencyTerms / 2 + 1;
            newCache.cosineFrequencies = Eigen::ArrayXd::LinSpaced(
                newCache.numCosTerms, 0.0, double( numFrequencyTerms / 2 ) );

            // Even: sin frequencies 1, 2, ..., N/2-1
            newCache.numSinTerms = numFrequencyTerms / 2 - 1;
            newCache.sineFrequencies = Eigen::ArrayXd::LinSpaced(
                newCache.numSinTerms, 1.0, double( numFrequencyTerms / 2 - 1 ) );
        }
        else
        {
            // Odd: cos frequencies 0, 1, ..., N/2
            newCache.numCosTerms = numFrequencyTerms / 2 + 1;
            newCache.cosineFrequencies = Eigen::ArrayXd::LinSpaced(
                newCache.numCosTerms, 0.0, double( numFrequencyTerms / 2 ) );

            // Odd: sin frequencies 1, 2, ..., N/2
            newCache.numSinTerms = numFrequencyTerms / 2;
            newCache.sineFrequencies = Eigen::ArrayXd::LinSpaced(
                newCache.numSinTerms, 1.0, double( numFrequencyTerms / 2 ) );
        }

        cache[ numFrequencyTerms ] = newCache;
        return cache[ numFrequencyTerms ];
    }

    // Pre-compute IFFT basis (trigonometric terms) for given solar longitude and number of frequency terms
    static Eigen::RowVectorXd computeIFFTBasis( const int numFrequencyTerms, const double solarLongitude )
    {
        using namespace Eigen;

        // Get cached frequency indices
        const FrequencyCache& freqCache = getFrequencyCache( numFrequencyTerms );

        // Compute scaled frequencies for vectorized trigonometric operations
        const ArrayXd cosFreqsScaled = freqCache.cosineFrequencies * solarLongitude;
        const ArrayXd sinFreqsScaled = freqCache.sineFrequencies * solarLongitude;

        // Compute trigonometric terms using Eigen's vectorized operations
        // Modern compilers will use SIMD instructions (SSE/AVX) for these operations
        ArrayXd cosTerms( freqCache.numCosTerms );
        ArrayXd sinTerms( freqCache.numSinTerms );

        // Vectorized trig computation
        cosTerms = cosFreqsScaled.cos( );
        sinTerms = sinFreqsScaled.sin( );

        // Pre-allocate and fill basis vector using direct assignment (avoid .transpose())
        RowVectorXd basis( numFrequencyTerms );
        basis.head( freqCache.numCosTerms ) = cosTerms;
        basis.tail( freqCache.numSinTerms ) = sinTerms;

        return basis / double( numFrequencyTerms );
    }

    static double radialPolyvalAndTemporalIFFT( const Eigen::Ref<const Eigen::RowVectorXd>& ifftBasis,
                                                const Eigen::Ref<const Eigen::MatrixXd>& polynomialMatrix,
                                                const Eigen::Ref<const Eigen::ArrayXd>& radialPowers )
    {
        return ifftBasis.dot( polynomialMatrix * radialPowers.matrix() );
    }

    static double reducedToTemporalIFFT( const Eigen::Ref<const Eigen::RowVectorXd>& ifftBasis,
                                        const Eigen::Ref<const Eigen::VectorXd>& polynomialVector )
    {
        return ( ifftBasis * polynomialVector ).value( );
    }

public:
    static void evaluate2D(
            const double radius_m,
            // meter
            const double solarLongitude,
            // radians
            const Eigen::MatrixXd& polyCoefficients,
            const Eigen::ArrayXXi& atDegreeAndOrder,
            const Eigen::ArrayXd& atPowersInvRadius,
            double refRadius_m,
            // meter
            Eigen::MatrixXd& cosineCoefficients,
            Eigen::MatrixXd& sineCoefficients,
            int maxDegree,
            int maxOrder,
            bool isLog2Data = true )
    {
        // --- Unit conversion ---
        const double radius_km = radius_m / 1000.0; // Conversion from m to km
        const double refRadius = refRadius_m / 1000.0; // Conversion from m to km

        const int maxDegreeAvailable = atDegreeAndOrder.row( 0 ).maxCoeff( );
        const int maxOrderAvailable = atDegreeAndOrder.row( 1 ).abs( ).maxCoeff( );

        if(maxDegree < 0) maxDegree = maxDegreeAvailable;
        if(maxOrder < 0) maxOrder = maxOrderAvailable;

        if(maxDegree > maxDegreeAvailable || maxOrder > maxOrderAvailable)
        {
            std::ostringstream err;
            err << "[FATAL] Requested maxDegree=" << maxDegree
                    << ", maxOrder=" << maxOrder
                    << " exceeds available (degree=" << maxDegreeAvailable
                    << ", order=" << maxOrderAvailable << ")";
            throw std::runtime_error( err.str( ) );
        }

        const int numRadialTerms = atPowersInvRadius.size( );
        const int numIntervals = polyCoefficients.rows( ) / numRadialTerms;

        cosineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );
        sineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );

        double inverseReferenceRadius = ( refRadius < 1.0e-10 ) ? 0.0 : 1.0 / refRadius;

        if(radius_km <= refRadius || refRadius < 1.0e-10)
        {
            // Pre-compute IFFT basis (trigonometric terms) - computed once per (radius, solarLongitude) pair
            const Eigen::RowVectorXd ifftBasis = computeIFFTBasis( numIntervals, solarLongitude );

            // Pre-compute radial power terms - computed once per radius
            const double radialDistance = 1.0 / radius_km - inverseReferenceRadius;
            const Eigen::ArrayXd radialPowers = pow( radialDistance, atPowersInvRadius );

            for(int coefficientIndex = 0; coefficientIndex < polyCoefficients.cols( ); ++coefficientIndex)
            {
                const int degree = atDegreeAndOrder( 0, coefficientIndex );
                const int order = atDegreeAndOrder( 1, coefficientIndex );
                const int absoluteOrder = std::abs( order );

                if(degree > maxDegree || absoluteOrder > maxOrder)
                    continue;

                // Use Eigen::Map to avoid copying/reshaping polynomial coefficients
                // The column is stored as [T0R0, T0R1, T0R2, T0R3, T1R0, T1R1, ...] (col-major)
                // Map interprets this as a (numIntervals × numRadialTerms) matrix without allocation
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                    polyCoefs( polyCoefficients.col( coefficientIndex ).data( ), numIntervals, numRadialTerms );

                double value = radialPolyvalAndTemporalIFFT( ifftBasis, polyCoefs, radialPowers );

                if(order >= 0)
                    cosineCoefficients( degree, absoluteOrder ) = value;
                else
                    sineCoefficients( degree, absoluteOrder ) = value;
            }
        }
        else
        {
            // Pre-compute IFFT basis (trigonometric terms) - computed once per solarLongitude
            const Eigen::RowVectorXd ifftBasis = computeIFFTBasis( numIntervals, solarLongitude );

            for(int coefficientIndex = 0; coefficientIndex < polyCoefficients.cols( ); ++coefficientIndex)
            {
                const int degree = atDegreeAndOrder( 0, coefficientIndex );
                const int order = atDegreeAndOrder( 1, coefficientIndex );
                const int absoluteOrder = std::abs( order );

                if(degree > maxDegree || absoluteOrder > maxOrder)
                    continue;

                const Eigen::VectorXd polyCoefs =
                        polyCoefficients.block( 0, coefficientIndex, numIntervals, 1 );

                // Use pre-computed basis
                double value = reducedToTemporalIFFT( ifftBasis, polyCoefs );

                if(order >= 0)
                    cosineCoefficients( degree, absoluteOrder ) = value;
                else
                    sineCoefficients( degree, absoluteOrder ) = value;
            }

            // Apply 1/r² decay term
            applyDecayTerm( cosineCoefficients, radius_m, refRadius_m, isLog2Data );
        }
    }

    /// Evaluate polynomial coefficients using only reduced temporal method (no radial component)
    /// This method is used to pre-compute Stokes coefficients for radius > reference radius
    /// The decay term must be applied separately during runtime
    /// @param solarLongitude Solar longitude in radians
    /// @param polyCoefficients Polynomial coefficient array
    /// @param atDegreeAndOrder Degree and order indices
    /// @param atPowersInvRadius Powers of inverse radius array
    /// @param cosineCoefficients Output cosine coefficients matrix
    /// @param sineCoefficients Output sine coefficients matrix
    /// @param maxDegree Maximum degree (-1 for auto)
    /// @param maxOrder Maximum order (-1 for auto)
    static void evaluate2DReduced(
            const double solarLongitude,
            // radians
            const Eigen::MatrixXd& polyCoefficients,
            const Eigen::ArrayXXi& atDegreeAndOrder,
            const Eigen::ArrayXd& atPowersInvRadius,
            Eigen::MatrixXd& cosineCoefficients,
            Eigen::MatrixXd& sineCoefficients,
            int maxDegree,
            int maxOrder )
    {
        const int maxDegreeAvailable = atDegreeAndOrder.row( 0 ).maxCoeff( );
        const int maxOrderAvailable = atDegreeAndOrder.row( 1 ).abs( ).maxCoeff( );

        if(maxDegree < 0) maxDegree = maxDegreeAvailable;
        if(maxOrder < 0) maxOrder = maxOrderAvailable;

        if(maxDegree > maxDegreeAvailable || maxOrder > maxOrderAvailable)
        {
            std::ostringstream err;
            err << "[FATAL] Requested maxDegree=" << maxDegree
                    << ", maxOrder=" << maxOrder
                    << " exceeds available (degree=" << maxDegreeAvailable
                    << ", order=" << maxOrderAvailable << ")";
            throw std::runtime_error( err.str( ) );
        }

        const int numRadialTerms = atPowersInvRadius.size( );
        const int numIntervals = polyCoefficients.rows( ) / numRadialTerms;

        cosineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );
        sineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );

        // Pre-compute IFFT basis (trigonometric terms) - computed once per solarLongitude
        const Eigen::RowVectorXd ifftBasis = computeIFFTBasis( numIntervals, solarLongitude );

        // Process only the reduced temporal coefficients (first row of each interval)
        for(int coefficientIndex = 0; coefficientIndex < polyCoefficients.cols( ); ++coefficientIndex)
        {
            const int degree = atDegreeAndOrder( 0, coefficientIndex );
            const int order = atDegreeAndOrder( 1, coefficientIndex );
            const int absoluteOrder = std::abs( order );

            if(degree > maxDegree || absoluteOrder > maxOrder)
                continue;

            // Extract the reduced polynomial coefficients (first row of each interval)
            const Eigen::VectorXd polyCoefs =
                    polyCoefficients.block( 0, coefficientIndex, numIntervals, 1 );

            // Use pre-computed basis for evaluation (no radius component)
            const double value = reducedToTemporalIFFT( ifftBasis, polyCoefs );

            if(order >= 0)
                cosineCoefficients( degree, absoluteOrder ) = value;
            else
                sineCoefficients( degree, absoluteOrder ) = value;
        }
    }

    /// Apply 1/r² decay term when radius exceeds reference radius.
    /// For log2-transformed data: additive in log2 space, C(0,0) += 2*log2(R_ref/r).
    /// For linear data: multiplicative, C(0,0) *= (R_ref/r)².
    /// @param cosineCoefficients Matrix of cosine coefficients to modify
    /// @param radius_m Current radius in meters
    /// @param referenceRadius_m Reference radius in meters
    /// @param isLog2Data Whether coefficients represent log2-transformed data (default true)
    static void applyDecayTerm( Eigen::MatrixXd& cosineCoefficients, double radius_m, double referenceRadius_m,
                                bool isLog2Data = true )
    {
        if ( isLog2Data )
        {
            // In log2 space: log2(n * (R_ref/r)²) = log2(n) + 2*log2(R_ref/r)
            const double radius_km = radius_m / 1000.0;
            const double referenceRadius_km = referenceRadius_m / 1000.0;
            cosineCoefficients( 0, 0 ) += 2.0 *
                    std::log2( ( referenceRadius_km < 1.0e-10 ) ? 1.0 / radius_km : referenceRadius_km / radius_km );
        }
        else
        {
            // In linear space: n * (R_ref/r)²
            const double ratio = ( referenceRadius_m < 1.0e-7 ) ? 0.0 : referenceRadius_m / radius_m;
            cosineCoefficients( 0, 0 ) *= ratio * ratio;
        }
    }
};

class ComaDatasetTransformer
{
public:
    static ComaStokesDataset transformPolyToStokes(
            const ComaPolyDataset& polyDataset,
            const std::vector< double >& radii_m,
            const std::vector< double >& solLongitudes_deg,
            const int requestedMaxDegree = -1,
            const int requestedMaxOrder = -1,
            const bool computeReducedCoeffs = true,
            const bool isLog2Data = true )
    {
        // Validate inputs
        validateTransformInputs( polyDataset,
                                 radii_m,
                                 solLongitudes_deg,
                                 requestedMaxDegree,
                                 requestedMaxOrder );

        // Determine effective maxima
        auto maxima = determineEffectiveMaxima(
                polyDataset,
                requestedMaxDegree,
                requestedMaxOrder );
        int effectiveMaxDegree = maxima.first;
        int effectiveMaxOrder = maxima.second;

        // Build file metadata
        auto files = buildFileMeta( polyDataset );

        // Convert solar longitudes from degrees to radians
        // The dataset will store radians internally for use with interpolators
        std::vector<double> solLongitudes_rad(solLongitudes_deg.size());
        for (std::size_t i = 0; i < solLongitudes_deg.size(); ++i)
        {
            solLongitudes_rad[i] = solLongitudes_deg[i] * mathematical_constants::PI / 180.0;
        }

        // Create empty Stokes dataset (with longitudes in radians)
        ComaStokesDataset stokesDataset = ComaStokesDataset::create(
                std::move( files ),
                radii_m,
                solLongitudes_rad,
                effectiveMaxDegree,
                computeReducedCoeffs );

        // Fill dataset with computed coefficients
        fillStokesDataset( polyDataset, stokesDataset, effectiveMaxDegree, effectiveMaxOrder, computeReducedCoeffs, isLog2Data );

        return stokesDataset;
    }

private:
    static void validateTransformInputs(
            const ComaPolyDataset& polyDataset,
            const std::vector< double >& radii_m,
            const std::vector< double >& solLongitudes_deg,
            const int requestedMaxDegree,
            const int requestedMaxOrder )
    {
        if(radii_m.empty( ) || solLongitudes_deg.empty( ))
            throw std::invalid_argument( "transformPolyToStokes: radii and longitudes must be non-empty." );
        if(radii_m.size( ) < 2)
            throw std::invalid_argument( "transformPolyToStokes: at least two radii are necessary for proper interpolation." );
        if(solLongitudes_deg.size( ) < 2)
            throw std::invalid_argument( "transformPolyToStokes: at least two solar longitudes are necessary for proper interpolation." );
        if(requestedMaxDegree < -1 || requestedMaxOrder < -1)
            throw std::invalid_argument( "transformPolyToStokes: requested maxima must be >= -1." );
        if(polyDataset.getNumFiles( ) == 0)
            throw std::invalid_argument( "transformPolyToStokes: polyDataset has no files." );
    }

    static std::pair< int, int > determineEffectiveMaxima(
            const ComaPolyDataset& polyDataset,
            const int requestedMaxDegree,
            const int requestedMaxOrder )
    {
        int globalMaxDegree = 0;
        int globalMaxOrder = 0;

        const std::size_t numFiles = polyDataset.getNumFiles( );
        std::vector< int > perFileMaxDegree( numFiles, 0 );
        std::vector< int > perFileMaxOrder( numFiles, 0 );

        for(std::size_t fileIndex = 0; fileIndex < numFiles; ++fileIndex)
        {
            const int fileMaxDegree = polyDataset.getMaxDegreeSH( fileIndex );
            const int fileMaxOrder = polyDataset.getSHDegreeAndOrderIndices( fileIndex ).row( 1 ).abs( ).maxCoeff( );

            perFileMaxDegree[ fileIndex ] = fileMaxDegree;
            perFileMaxOrder[ fileIndex ] = fileMaxOrder;

            globalMaxDegree = std::max( globalMaxDegree, fileMaxDegree );
            globalMaxOrder = std::max( globalMaxOrder, fileMaxOrder );
        }

        const int effectiveMaxDegree = requestedMaxDegree < 0 ? globalMaxDegree : requestedMaxDegree;
        const int effectiveMaxOrder = requestedMaxOrder < 0 ? globalMaxOrder : requestedMaxOrder;

        if(effectiveMaxDegree < 0 || effectiveMaxOrder < 0)
            throw std::invalid_argument( "determineEffectiveMaxima: resolved negative maxima." );
        // Ensure every file can support requested maxima
        for(std::size_t fileIndex = 0; fileIndex < numFiles; ++fileIndex)
        {
            if(effectiveMaxDegree > perFileMaxDegree[ fileIndex ] || effectiveMaxOrder > perFileMaxOrder[ fileIndex ])
            {
                std::ostringstream oss;
                oss << "Requested (nmax=" << effectiveMaxDegree << ", mmax=" << effectiveMaxOrder
                        << ") exceeds availability in file #" << fileIndex
                        << " [" << polyDataset.getFileMeta( fileIndex ).sourcePath << "]: "
                        << "(maxDegree=" << perFileMaxDegree[ fileIndex ]
                        << ", maxOrder=" << perFileMaxOrder[ fileIndex ] << ").";
                throw std::invalid_argument( oss.str( ) );
            }
        }

        return { effectiveMaxDegree, effectiveMaxOrder };
    }

    static std::vector< ComaStokesDataset::FileMeta > buildFileMeta(
            const ComaPolyDataset& polyDataset )
    {
        std::vector< ComaStokesDataset::FileMeta > files;
        files.reserve( polyDataset.getNumFiles( ) );
        for(std::size_t fileIndex = 0; fileIndex < polyDataset.getNumFiles( ); ++fileIndex)
        {
            ComaStokesDataset::FileMeta fileMeta{};
            fileMeta.source_tag = polyDataset.getFileMeta( fileIndex ).sourcePath;
            fileMeta.referenceRadius = polyDataset.getFileMeta( fileIndex ).referenceRadius;
            if(!polyDataset.getFileMeta( fileIndex ).timePeriods.empty( ))
            {
                fileMeta.start_epoch = polyDataset.getFileMeta( fileIndex ).timePeriods.front( ).first;
                fileMeta.end_epoch = polyDataset.getFileMeta( fileIndex ).timePeriods.front( ).second;
            }
            else
            {
                fileMeta.start_epoch = 0.0;
                fileMeta.end_epoch = 0.0;
            }
            files.push_back( std::move( fileMeta ) );
        }
        return files;
    }

    static void fillStokesDataset(
            const ComaPolyDataset& polyDataset,
            ComaStokesDataset& stokesDataset,
            const int effectiveMaxDegree,
            const int effectiveMaxOrder,
            const bool computeReducedCoeffs,
            const bool isLog2Data = true )
    {
        Eigen::MatrixXd cosineCoefficients( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        Eigen::MatrixXd sineCoefficients( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );

        // Only allocate reduced coefficient matrices if needed
        Eigen::MatrixXd reducedCosineCoefficients;
        Eigen::MatrixXd reducedSineCoefficients;
        if( computeReducedCoeffs )
        {
            reducedCosineCoefficients.resize( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
            reducedSineCoefficients.resize( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        }

        for(std::size_t fileIndex = 0; fileIndex < stokesDataset.nFiles( ); ++fileIndex)
        {
            const auto& polynomialCoefficients = polyDataset.getPolyCoefficients( fileIndex );
            const auto& degreeAndOrderIndices = polyDataset.getSHDegreeAndOrderIndices( fileIndex );
            const auto& inversePowers = polyDataset.getPowersInvRadius( fileIndex );
            const double referenceRadius = polyDataset.getReferenceRadius( fileIndex ); // (unit as in files; matches old code)

            // Fill regular coefficients (for radius <= reference radius)
            for(std::size_t radiusIndex = 0; radiusIndex < stokesDataset.nRadii( ); ++radiusIndex)
            {
                const double radius_m = stokesDataset.radii( )[ radiusIndex ]; // radius in meters

                for(std::size_t longitudeIndex = 0; longitudeIndex < stokesDataset.nLongitudes( ); ++longitudeIndex)
                {
                    const double solarLongitudeRadians = stokesDataset.lons( )[ longitudeIndex ]; // already in radians

                    // Evaluate regular coefficients using existing method
                    StokesCoefficientsEvaluator::evaluate2D(
                            radius_m,
                            solarLongitudeRadians,
                            polynomialCoefficients,
                            degreeAndOrderIndices,
                            inversePowers,
                            referenceRadius,
                            cosineCoefficients,
                            sineCoefficients,
                            effectiveMaxDegree,
                            effectiveMaxOrder,
                            isLog2Data );

                    // Store regular coefficients using batch write via block reference
                    auto coeffBlock = stokesDataset.block( fileIndex, radiusIndex, longitudeIndex );
                    for(int degree = 0; degree <= effectiveMaxDegree; ++degree)
                    {
                        const int maxOrderForDegree = std::min( degree, effectiveMaxOrder );
                        for(int order = 0; order <= maxOrderForDegree; ++order)
                        {
                            if(degree >= cosineCoefficients.rows( ) || order >= cosineCoefficients.cols( ))
                                throw std::runtime_error( "Evaluator returned C/S smaller than requested (nmax,mmax)." );
                            const std::size_t coefficientIndex = nm_to_index_deg_major( degree, order );
                            coeffBlock( static_cast< Eigen::Index >(coefficientIndex), 0 ) = cosineCoefficients( degree, order );
                            coeffBlock( static_cast< Eigen::Index >(coefficientIndex), 1 ) = sineCoefficients( degree, order );
                        }
                    }
                }
            }

            // Fill reduced coefficients (computed once per solar longitude, no radius dependency)
            // Only compute and store if needed (normal coma model requires them, wind model does not)
            if( computeReducedCoeffs )
            {
                for(std::size_t longitudeIndex = 0; longitudeIndex < stokesDataset.nLongitudes( ); ++longitudeIndex)
                {
                    const double solarLongitudeRadians = stokesDataset.lons( )[ longitudeIndex ]; // already in radians


                    reducedCosineCoefficients.setZero();
                    reducedSineCoefficients.setZero();

                    StokesCoefficientsEvaluator::evaluate2DReduced(
                        solarLongitudeRadians,
                        polynomialCoefficients,
                        degreeAndOrderIndices,
                        inversePowers,
                        reducedCosineCoefficients,
                        reducedSineCoefficients,
                        effectiveMaxDegree,
                        effectiveMaxOrder );

                    // Store reduced coefficients using batch write via block reference
                    auto reducedCoeffBlock = stokesDataset.reducedBlock( fileIndex, longitudeIndex );
                    for(int degree = 0; degree <= effectiveMaxDegree; ++degree)
                    {
                        const int maxOrderForDegree = std::min( degree, effectiveMaxOrder );
                        for(int order = 0; order <= maxOrderForDegree; ++order)
                        {
                            if(degree >= reducedCosineCoefficients.rows( ) || order >= reducedCosineCoefficients.cols( ))
                                throw std::runtime_error( "Evaluator returned reducedC/S smaller than requested (nmax,mmax)." );
                            const std::size_t coefficientIndex = nm_to_index_deg_major( degree, order );
                            reducedCoeffBlock( static_cast< Eigen::Index >(coefficientIndex), 0 ) = reducedCosineCoefficients( degree, order );
                            reducedCoeffBlock( static_cast< Eigen::Index >(coefficientIndex), 1 ) = reducedSineCoefficients( degree, order );
                        }
                    }
                }
            }
        }
    }
};

// ============= Wind Model Dataset Collection =============

/*!
 * \brief Collection of three coma datasets for wind model in modified vertical frame.
 *
 * This class holds three datasets (one for each component of the wind velocity in the modified vertical frame)
 * that are used together to construct a ComaWindModel. All three datasets must be of
 * the same type (either all polynomial or all Stokes coefficients).
 *
 * The wind velocity components in these input datasets are defined in a modified vertical frame:
 *   - X-component: Meridional direction (North, in meridian plane)
 *   - Y-component: Zonal direction (West, completing the right-handed frame)
 *   - Z-component: Radial direction pointing OUTWARD from the nucleus (opposite to standard vertical frame)
 */
class ComaWindDatasetCollection
{
public:
    //! Type alias for the dataset variant
    using DataVariant = boost::variant<ComaPolyDataset, ComaStokesDataset>;

    /*!
     * \brief Factory method to create collection from polynomial datasets.
     * \param xPolyDataset Polynomial dataset for x-component wind
     * \param yPolyDataset Polynomial dataset for y-component wind
     * \param zPolyDataset Polynomial dataset for z-component wind
     * \return ComaWindDatasetCollection containing the three datasets
     */
    static ComaWindDatasetCollection createFromPoly(
        ComaPolyDataset xPolyDataset,
        ComaPolyDataset yPolyDataset,
        ComaPolyDataset zPolyDataset)
    {
        ComaWindDatasetCollection collection;
        collection.xData_ = std::move(xPolyDataset);
        collection.yData_ = std::move(yPolyDataset);
        collection.zData_ = std::move(zPolyDataset);
        collection.dataType_ = DataType::Polynomial;
        return collection;
    }

    /*!
     * \brief Factory method to create collection from Stokes datasets.
     * \param xStokesDataset Stokes dataset for x-component wind
     * \param yStokesDataset Stokes dataset for y-component wind
     * \param zStokesDataset Stokes dataset for z-component wind
     * \return ComaWindDatasetCollection containing the three datasets
     */
    static ComaWindDatasetCollection createFromStokes(
        ComaStokesDataset xStokesDataset,
        ComaStokesDataset yStokesDataset,
        ComaStokesDataset zStokesDataset)
    {
        ComaWindDatasetCollection collection;
        collection.xData_ = std::move(xStokesDataset);
        collection.yData_ = std::move(yStokesDataset);
        collection.zData_ = std::move(zStokesDataset);
        collection.dataType_ = DataType::Stokes;
        return collection;
    }

    /*!
     * \brief Check if the collection contains polynomial datasets.
     * \return True if datasets are polynomial type
     */
    bool isPoly() const
    {
        return dataType_ == DataType::Polynomial;
    }

    /*!
     * \brief Check if the collection contains Stokes datasets.
     * \return True if datasets are Stokes type
     */
    bool isStokes() const
    {
        return dataType_ == DataType::Stokes;
    }

    /*!
     * \brief Get x-component polynomial dataset.
     * \return Reference to x-component polynomial dataset
     * \throws std::runtime_error if datasets are not polynomial type
     */
    const ComaPolyDataset& getXPolyDataset() const
    {
        if (auto* p = boost::get<ComaPolyDataset>(&xData_))
            return *p;
        throw std::runtime_error("X-component data does not contain polynomial dataset");
    }

    /*!
     * \brief Get y-component polynomial dataset.
     * \return Reference to y-component polynomial dataset
     * \throws std::runtime_error if datasets are not polynomial type
     */
    const ComaPolyDataset& getYPolyDataset() const
    {
        if (auto* p = boost::get<ComaPolyDataset>(&yData_))
            return *p;
        throw std::runtime_error("Y-component data does not contain polynomial dataset");
    }

    /*!
     * \brief Get z-component polynomial dataset.
     * \return Reference to z-component polynomial dataset
     * \throws std::runtime_error if datasets are not polynomial type
     */
    const ComaPolyDataset& getZPolyDataset() const
    {
        if (auto* p = boost::get<ComaPolyDataset>(&zData_))
            return *p;
        throw std::runtime_error("Z-component data does not contain polynomial dataset");
    }

    /*!
     * \brief Get x-component Stokes dataset.
     * \return Reference to x-component Stokes dataset
     * \throws std::runtime_error if datasets are not Stokes type
     */
    const ComaStokesDataset& getXStokesDataset() const
    {
        if (auto* p = boost::get<ComaStokesDataset>(&xData_))
            return *p;
        throw std::runtime_error("X-component data does not contain Stokes dataset");
    }

    /*!
     * \brief Get y-component Stokes dataset.
     * \return Reference to y-component Stokes dataset
     * \throws std::runtime_error if datasets are not Stokes type
     */
    const ComaStokesDataset& getYStokesDataset() const
    {
        if (auto* p = boost::get<ComaStokesDataset>(&yData_))
            return *p;
        throw std::runtime_error("Y-component data does not contain Stokes dataset");
    }

    /*!
     * \brief Get z-component Stokes dataset.
     * \return Reference to z-component Stokes dataset
     * \throws std::runtime_error if datasets are not Stokes type
     */
    const ComaStokesDataset& getZStokesDataset() const
    {
        if (auto* p = boost::get<ComaStokesDataset>(&zData_))
            return *p;
        throw std::runtime_error("Z-component data does not contain Stokes dataset");
    }

private:
    //! Enumeration for dataset type
    enum class DataType
    {
        Polynomial,
        Stokes
    };

    //! Type of datasets stored in this collection
    DataType dataType_;

    //! X-component dataset (either polynomial or Stokes)
    DataVariant xData_;

    //! Y-component dataset (either polynomial or Stokes)
    DataVariant yData_;

    //! Z-component dataset (either polynomial or Stokes)
    DataVariant zData_;
};

// ============= High-Level Processing Interface =============

class ComaModelFileProcessor
{
public:
    enum class FileType
    {
        PolyCoefficients,
        StokesCoefficients
    };

    // Constructor for polynomial coefficient files
    explicit ComaModelFileProcessor( std::vector< std::string > filePaths ) :
        filePaths_( std::move( filePaths ) ), fileType_( FileType::PolyCoefficients )
    {
        if(filePaths_.empty( ))
        {
            throw std::invalid_argument( "ComaModelFileProcessor: empty file list provided to constructor" );
        }
    }

    // Constructor for Stokes coefficient files (SH files)
    explicit ComaModelFileProcessor( const std::string& inputDir, const std::string& prefix = "stokes" ) :
        fileType_( FileType::StokesCoefficients ), shInputDir_( inputDir ), shPrefix_( prefix ),
        preloadedSHDataset_( ComaStokesDatasetReader::readFromCsvFolder( inputDir, prefix ) )
    {
    }

    // Create poly dataset from files (only available for poly coefficient files)
    ComaPolyDataset createPolyCoefDataset( ) const
    {
        if(fileType_ != FileType::PolyCoefficients)
        {
            throw std::runtime_error( "createPolyCoefDataset: not available when processor is constructed from SH files. "
                    "Use a processor constructed from polynomial coefficient files instead." );
        }
        if(filePaths_.empty( ))
        {
            throw std::invalid_argument( "ComaPolyDataset createPolyCoefDataset( ): empty file list. "
                    "This should never happen as the constructor validates the file list." );
        }
        return ComaPolyDatasetReader::readFromFiles( filePaths_ );
    }

    // Create Stokes dataset (parameterless version for Stokes coefficient files)
    ComaStokesDataset createSHDataset() const
    {
        if(fileType_ != FileType::StokesCoefficients)
        {
            throw std::runtime_error(
                "createSHDataset(): parameterless version only available when processor is constructed from "
                "Stokes coefficient CSV files. Use createSHDataset(radii_m, sol_longitudes_deg, ...) for "
                "polynomial coefficient files." );
        }
        return preloadedSHDataset_;
    }

    // Create Stokes dataset (with parameters for polynomial coefficient files)
    ComaStokesDataset createSHDataset(
            const std::vector< double >& radii_m,
            const std::vector< double >& solLongitudes_deg,
            const int requestedMaxDegree = -1,
            const int requestedMaxOrder = -1,
            const bool computeReducedCoeffs = true,
            const bool isLog2Data = true ) const
    {
        if(fileType_ != FileType::PolyCoefficients)
        {
            throw std::runtime_error(
                "createSHDataset(radii_m, sol_longitudes_deg, ...): only available when processor is constructed "
                "from polynomial coefficient files. Use parameterless createSHDataset() for Stokes coefficient CSV files." );
        }

        // Transform poly coefficients to Stokes coefficients
        const ComaPolyDataset polyDataset = createPolyCoefDataset(  );
        return ComaDatasetTransformer::transformPolyToStokes(
                polyDataset,
                radii_m,
                solLongitudes_deg,
                requestedMaxDegree,
                requestedMaxOrder,
                computeReducedCoeffs,
                isLog2Data );
    }

    // Create SH files (combines dataset creation and writing)
    void createSHFiles(
            const std::string& outputDir,
            const std::vector< double >& radii_m,
            const std::vector< double >& solLongitudes_deg,
            const int requestedMaxDegree = -1,
            const int requestedMaxOrder = -1,
            const bool computeReducedCoeffs = true,
            const bool isLog2Data = true ) const
    {
        const ComaStokesDataset stokesDataset = createSHDataset(
                radii_m,
                solLongitudes_deg,
                requestedMaxDegree,
                requestedMaxOrder,
                computeReducedCoeffs,
                isLog2Data );

        ComaStokesDatasetWriter::writeCsvAll( stokesDataset, outputDir );
    }

    // Get the file type of this processor
    FileType getFileType( ) const
    {
        return fileType_;
    }

private:
    std::vector< std::string > filePaths_;
    FileType fileType_;

    // For SH file processor
    std::string shInputDir_;
    std::string shPrefix_;
    ComaStokesDataset preloadedSHDataset_;
};

// Stream operator for FileType enum (needed for Boost Test)
inline std::ostream& operator<<( std::ostream& os, const ComaModelFileProcessor::FileType& type )
{
    switch(type)
    {
        case ComaModelFileProcessor::FileType::PolyCoefficients:
            return os << "PolyCoefficients";
        case ComaModelFileProcessor::FileType::StokesCoefficients:
            return os << "StokesCoefficients";
        default:
            return os << "Unknown";
    }
}

/*!
 * \brief Processor for creating wind model datasets from three component file sources.
 *
 * This class manages the creation of ComaWindDatasetCollection from three sets of files
 * (one for each component in the modified vertical frame: X=meridional/North, Y=zonal/West, Z=radial outward).
 * It internally uses three ComaModelFileProcessor instances and provides a simplified interface for wind model setup.
 *
 * \note The Z-component represents radial wind pointing OUTWARD from the nucleus, opposite to standard vertical frame.
 */
class ComaWindModelFileProcessor
{
public:
    /*!
     * \brief Constructor for polynomial coefficient files.
     * \param xFilePaths Vector of file paths for X-component (meridional/North) polynomial coefficients
     * \param yFilePaths Vector of file paths for Y-component (zonal/West) polynomial coefficients
     * \param zFilePaths Vector of file paths for Z-component (radial outward) polynomial coefficients
     * \throws std::invalid_argument if any file list is empty
     */
    ComaWindModelFileProcessor(
        std::vector<std::string> xFilePaths,
        std::vector<std::string> yFilePaths,
        std::vector<std::string> zFilePaths)
        : xProcessor_(std::move(xFilePaths)),
          yProcessor_(std::move(yFilePaths)),
          zProcessor_(std::move(zFilePaths)),
          isPolyType_(true)
    {
        if (xProcessor_.getFileType() != ComaModelFileProcessor::FileType::PolyCoefficients ||
            yProcessor_.getFileType() != ComaModelFileProcessor::FileType::PolyCoefficients ||
            zProcessor_.getFileType() != ComaModelFileProcessor::FileType::PolyCoefficients)
        {
            throw std::runtime_error("ComaWindModelFileProcessor: All processors must be polynomial type for this constructor");
        }
    }

    /*!
     * \brief Constructor for Stokes coefficient directories.
     * \param xInputDir Directory containing X-component (meridional/North) Stokes CSV files
     * \param yInputDir Directory containing Y-component (zonal/West) Stokes CSV files
     * \param zInputDir Directory containing Z-component (radial outward) Stokes CSV files
     * \param prefix File prefix for the CSV files (default: "stokes")
     */
    ComaWindModelFileProcessor(
        const std::string& xInputDir,
        const std::string& yInputDir,
        const std::string& zInputDir,
        const std::string& prefix = "stokes")
        : xProcessor_(xInputDir, prefix),
          yProcessor_(yInputDir, prefix),
          zProcessor_(zInputDir, prefix),
          isPolyType_(false)
    {
        if (xProcessor_.getFileType() != ComaModelFileProcessor::FileType::StokesCoefficients ||
            yProcessor_.getFileType() != ComaModelFileProcessor::FileType::StokesCoefficients ||
            zProcessor_.getFileType() != ComaModelFileProcessor::FileType::StokesCoefficients)
        {
            throw std::runtime_error("ComaWindModelFileProcessor: All processors must be Stokes type for this constructor");
        }
    }

    /*!
     * \brief Create polynomial coefficient dataset collection for all three components.
     * \return ComaWindDatasetCollection containing x, y, z polynomial datasets
     * \throws std::runtime_error if processor was constructed from Stokes files
     */
    ComaWindDatasetCollection createPolyCoefDatasets() const
    {
        if (!isPolyType_)
        {
            throw std::runtime_error("createPolyCoefDatasets: not available when processor is constructed from SH files. "
                                   "Use a processor constructed from polynomial coefficient files instead.");
        }

        ComaPolyDataset xData = xProcessor_.createPolyCoefDataset();
        ComaPolyDataset yData = yProcessor_.createPolyCoefDataset();
        ComaPolyDataset zData = zProcessor_.createPolyCoefDataset();

        return ComaWindDatasetCollection::createFromPoly(
            std::move(xData),
            std::move(yData),
            std::move(zData));
    }

    /*!
     * \brief Create Stokes coefficient dataset collection (parameterless version for Stokes CSV files).
     * \return ComaWindDatasetCollection containing x, y, z Stokes datasets from preloaded CSV files
     * \throws std::runtime_error if processor was constructed from polynomial files
     */
    ComaWindDatasetCollection createSHDatasets() const
    {
        if (!isStokesType())
        {
            throw std::runtime_error(
                "createSHDatasets(): parameterless version only available when processor is constructed from "
                "Stokes coefficient CSV files. Use createSHDatasets(radii_m, sol_longitudes_deg, ...) for "
                "polynomial coefficient files.");
        }

        ComaStokesDataset xData = xProcessor_.createSHDataset();
        ComaStokesDataset yData = yProcessor_.createSHDataset();
        ComaStokesDataset zData = zProcessor_.createSHDataset();

        return ComaWindDatasetCollection::createFromStokes(
            std::move(xData),
            std::move(yData),
            std::move(zData));
    }

    /*!
     * \brief Create Stokes coefficient dataset collection (with parameters for polynomial files).
     * \param radii_m Vector of radii at which to evaluate Stokes coefficients [m]
     * \param solLongitudes_deg Vector of solar longitudes [degrees]
     * \param requestedMaxDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedMaxOrder Maximum spherical harmonic order (-1 for auto)
     * \return ComaWindDatasetCollection containing x, y, z Stokes datasets
     * \throws std::runtime_error if processor was constructed from Stokes CSV files
     */
    ComaWindDatasetCollection createSHDatasets(
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree = -1,
        const int requestedMaxOrder = -1) const
    {
        if (!isPolyType())
        {
            throw std::runtime_error(
                "createSHDatasets(radii_m, sol_longitudes_deg, ...): only available when processor is constructed "
                "from polynomial coefficient files. Use parameterless createSHDatasets() for Stokes coefficient CSV files.");
        }

        ComaStokesDataset xData = xProcessor_.createSHDataset(radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder);
        ComaStokesDataset yData = yProcessor_.createSHDataset(radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder);
        ComaStokesDataset zData = zProcessor_.createSHDataset(radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder);

        return ComaWindDatasetCollection::createFromStokes(
            std::move(xData),
            std::move(yData),
            std::move(zData));
    }

    /*!
     * \brief Create Stokes coefficient CSV files for all three components.
     * \param xOutputDir Output directory for x-component files
     * \param yOutputDir Output directory for y-component files
     * \param zOutputDir Output directory for z-component files
     * \param radii_m Vector of radii at which to evaluate Stokes coefficients [m]
     * \param solLongitudes_deg Vector of solar longitudes [degrees]
     * \param requestedMaxDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedMaxOrder Maximum spherical harmonic order (-1 for auto)
     * \note Wind models do not require reduced Stokes coefficients, so they are not computed/saved
     */
    void createSHFiles(
        const std::string& xOutputDir,
        const std::string& yOutputDir,
        const std::string& zOutputDir,
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree = -1,
        const int requestedMaxOrder = -1) const
    {
        // Wind models do not need reduced coefficients - pass false to skip computation and save storage
        xProcessor_.createSHFiles(xOutputDir, radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder, false);
        yProcessor_.createSHFiles(yOutputDir, radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder, false);
        zProcessor_.createSHFiles(zOutputDir, radii_m, solLongitudes_deg, requestedMaxDegree, requestedMaxOrder, false);
    }

    /*!
     * \brief Check if this processor works with polynomial coefficient files.
     * \return True if constructed from polynomial files
     */
    bool isPolyType() const { return isPolyType_; }

    /*!
     * \brief Check if this processor works with Stokes coefficient files.
     * \return True if constructed from Stokes files
     */
    bool isStokesType() const { return !isPolyType_; }

private:
    //! Processor for x-component files
    ComaModelFileProcessor xProcessor_;

    //! Processor for y-component files
    ComaModelFileProcessor yProcessor_;

    //! Processor for z-component files
    ComaModelFileProcessor zProcessor_;

    //! Flag indicating if processors work with polynomial (true) or Stokes (false) files
    bool isPolyType_;
};

// Factory functions for ComaWindModelFileProcessor
inline std::shared_ptr<ComaWindModelFileProcessor> comaWindModelFileProcessorFromPolyFiles(
    const std::vector<std::string>& xFilePaths,
    const std::vector<std::string>& yFilePaths,
    const std::vector<std::string>& zFilePaths)
{
    return std::make_shared<ComaWindModelFileProcessor>(xFilePaths, yFilePaths, zFilePaths);
}

inline std::shared_ptr<ComaWindModelFileProcessor> comaWindModelFileProcessorFromSHFiles(
    const std::string& xInputDir,
    const std::string& yInputDir,
    const std::string& zInputDir,
    const std::string& prefix = "stokes")
{
    return std::make_shared<ComaWindModelFileProcessor>(xInputDir, yInputDir, zInputDir, prefix);
}

// === Coma processing: factory-style helpers (header-inline) ===
inline std::shared_ptr< ComaModelFileProcessor > comaModelFileProcessorFromPolyFiles(
        const std::vector< std::string >& filePaths )
{
    return std::make_shared< ComaModelFileProcessor >( filePaths );
}

inline std::shared_ptr< ComaModelFileProcessor > comaModelFileProcessorFromSHFiles(
        const std::string& inputDir,
        const std::string& prefix = "stokes" )
{
    return std::make_shared< ComaModelFileProcessor >( inputDir, prefix );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif // TUDAT_COMA_MODEL_INPUT_OUTPUT_H
