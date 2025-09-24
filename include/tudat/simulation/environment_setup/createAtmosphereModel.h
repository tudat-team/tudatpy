/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEATMOSPHEREMODEL_H
#define TUDAT_CREATEATMOSPHEREMODEL_H

#include <string>
#include <map>

#include <memory>
#include <boost/date_time/posix_time/time_period.hpp>

#include "body.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/basics/identityElements.h"

namespace tudat
{
namespace simulation_setup
{
using namespace aerodynamics;


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


/**
 * \class ComaStokesDataset
 * \brief Data handler for storing and managing sets of Stokes coefficients for the coma model.
 *
 * This class provides a structured container for spherical harmonic Stokes coefficients
 * (Cnm, Snm) evaluated across multiple input files, radii, and solar longitudes. It
 * maintains contiguous Eigen-based storage for efficiency while exposing convenient
 * accessors for individual blocks and coefficients.
 *
 * Key features:
 *  - Metadata storage for each file (validity epochs, descriptive tag).
 *  - Vectors of all radii and solar longitudes at which coefficients are evaluated.
 *  - Degree-major ordering of coefficients with closed-form (n,m) <-> index mapping.
 *  - Access to whole blocks of coefficients for (file, radius, longitude) combinations.
 *  - Access to individual coefficients by (n,m).
 *  - CSV export functionality: each file can be written to a human-readable CSV with
 *    metadata, radii, longitudes, and coefficient blocks.
 *
 * Planned extensions:
 *  - CSV reader: reconstruct a StokesDataset from one or more previously written CSVs.
 *  - Optional support for HDF5 or other binary formats for large expansions.
 *
 * Example usage:
 * \code{.cpp}
 * std::vector<FileMeta> files = {
 *     { 2.015e9, 2.0150864e9, "simulation_run_1" }
 * };
 * std::vector<double> radii = { 1000.0, 2000.0 };
 * std::vector<double> longitudes = { 0.0, 180.0 };
 *
 * auto dataset = StokesDataset::create(files, radii, longitudes, 4);
 * dataset.setCoeff(0, 0, 0, 2, 0, -0.484e-03, 0.0); // C20 term
 * dataset.writeCsvAll("output_dir/");
 * \endcode
 *
 * \see FileMeta
 */
class ComaStokesDataset
{
public:
    // -------- Construction --------
    static ComaStokesDataset create( std::vector< FileMeta > files,
                                 std::vector< double > radii,
                                 std::vector< double > lons,
                                 int nmax )
    {
        if(files.empty( ) || radii.empty( ) || lons.empty( ) || nmax < 0)
            throw std::runtime_error( "StokesDataset: invalid metadata." );

        ComaStokesDataset g;
        g.files_ = std::move( files );
        g.radii_ = std::move( radii );
        g.lons_ = std::move( lons );
        g.nmax_ = nmax;

        g.n_files_ = g.files_.size( );
        g.n_radii_ = g.radii_.size( );
        g.n_lons_ = g.lons_.size( );
        g.n_coeffs_ = static_cast< std::size_t >(( nmax + 1 ) * ( nmax + 2 ) / 2);

        const std::size_t totalRows = g.n_files_ * g.n_radii_ * g.n_lons_ * g.n_coeffs_;
        g.data_.setZero( static_cast< Eigen::Index >(totalRows), 2 );
        return g;
    }

    // -------- Sizes / metadata --------
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

    // -------- Block access (nCoeffs x 2 view) --------
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

    // -------- Single coefficient access via (n,m) --------
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

    // -------- CSV writing (one CSV per file) --------
    void writeCsvForFile( std::size_t f, const std::string& out_path ) const
    {
        if(f >= n_files_) throw std::out_of_range( "writeCsvForFile: file index OOR" );

        boost::filesystem::path p( out_path );
        if(p.has_parent_path( )) boost::filesystem::create_directories( p.parent_path( ) );

        std::ofstream os( out_path, std::ios::binary );
        if(!os) throw std::runtime_error( "writeCsvForFile: cannot open " + out_path );
        os.imbue( std::locale::classic( ) );

        const auto& fm = files_[ f ];

        // Row 1: metadata (key=value pairs)
        os << "meta";
        os << ",start_epoch=" << std::setprecision( 17 ) << std::scientific << fm.start_epoch;
        os << ",end_epoch=" << std::setprecision( 17 ) << std::scientific << fm.end_epoch;
        os << ",max_degree=" << nmax_;
        os << ",max_order=" << nmax_;
        os << ",n_radii=" << n_radii_;
        os << ",n_lons=" << n_lons_;
        os << ",n_coeffs=" << n_coeffs_;
        os << ",source=";
        csvEscape_( os, fm.source_tag );
        os << '\n';

        // Row 2: radii
        os << "radii";
        for(double r: radii_) os << ',' << std::setprecision( 17 ) << std::scientific << r;
        os << '\n';

        // Row 3: longitudes
        os << "longitudes";
        for(double L: lons_) os << ',' << std::setprecision( 17 ) << std::scientific << L;
        os << '\n';

        // Blocks: for each (ri, li)
        for(std::size_t ri = 0; ri < n_radii_; ++ri)
        {
            for(std::size_t li = 0; li < n_lons_; ++li)
            {
                const std::size_t block_id = ri * n_lons_ + li;

                // Block header
                os << "ID," << block_id << ','
                        << std::setprecision( 17 ) << std::scientific << radii_[ ri ] << ','
                        << std::setprecision( 17 ) << std::scientific << lons_[ li ] << '\n';

                // Column header
                os << "n,m,C,S\n";

                const auto blk = block( f, ri, li ); // n_coeffs_ x 2
                for(std::size_t k = 0; k < n_coeffs_; ++k)
                {
                    const auto nm = index_to_nm_deg_major( k );
                    const double C = blk( static_cast< Eigen::Index >(k), 0 );
                    const double S = blk( static_cast< Eigen::Index >(k), 1 );
                    os << nm.first << ',' << nm.second << ','
                            << std::setprecision( 17 ) << std::scientific << C << ','
                            << std::setprecision( 17 ) << std::scientific << S << '\n';
                }
            }
        }

        os.flush( );
        if(!os) throw std::runtime_error( "writeCsvForFile: write failed" );
    }

    void writeCsvAll( const std::string& out_dir, const std::string& prefix = "stokes" ) const
    {
        boost::filesystem::create_directories( out_dir );
        for(std::size_t f = 0; f < n_files_; ++f)
        {
            boost::filesystem::path path = boost::filesystem::path( out_dir )
                    / ( prefix + "_file" + std::to_string( f ) + ".csv" );
            writeCsvForFile( f, path.string( ) );
        }
    }

    // -------- (Optional) raw data access --------
    const StokesBlock& data( ) const
    {
        return data_;
    }

    StokesBlock& data( )
    {
        return data_;
    }

    // -------- Reader hooks (to implement later) --------
    // static StokesDataset readCsvFiles(const std::vector<std::string>& paths);
    // static StokesDataset readCsvFolder(const std::string& dir, const std::string& prefix);

private:
    // Layout: ((((f * n_radii_) + r) * n_lons_) + l) * n_coeffs_ + k
    std::size_t startRow_( std::size_t f, std::size_t r, std::size_t l ) const
    {
        if(f >= n_files_ || r >= n_radii_ || l >= n_lons_)
            throw std::out_of_range( "StokesDataset: index out of range." );
        const std::size_t cell = ( ( f * n_radii_ ) + r ) * n_lons_ + l;
        return cell * n_coeffs_;
    }

    static void csvEscape_( std::ofstream& os, const std::string& s )
    {
        bool needs = false;
        for(char c: s) if(c == ',' || c == '"' || c == '\n' || c == '\r')
        {
            needs = true;
            break;
        }
        if(!needs)
        {
            os << s;
            return;
        }
        os << '"';
        for(char c: s) os << ( c == '"' ? "\"\"" : std::string( 1, c ) );
        os << '"';
    }

    // ---- Metadata
    std::vector< FileMeta > files_;
    std::vector< double > radii_;
    std::vector< double > lons_;
    int nmax_{};

    // ---- Dimensions
    std::size_t n_files_{}, n_radii_{}, n_lons_{}, n_coeffs_{};

    // ---- Storage
    StokesBlock data_;
};


/**
 * \class ComaPolyDataset
 * \brief Data handler for polynomial coefficients used to build Stokes coefficients.
 *
 * Reads a list of "PolyCoefficient" files and stores, per file:
 *  - polyCoefficients  : (numTerms x numCoefs) matrix
 *  - SHDegreeAndOrder  : (2 x numCoefs) array; row 0 = degree n, row 1 = order m
 *  - referenceRadius   : scalar
 *  - powersInvRadius   : vector<double> of inverse-radius powers (from "PWRS" header)
 *  - timePeriods       : vector of (start_epoch, end_epoch) pairs [placeholder]
 *
 * Convenience:
 *  - Fast lookup from (n,m) -> column index via a cached map per file.
 *  - Accessors to get a full coefficient *column* (all terms for one (n,m)),
 *    or a specific entry (termIndex, (n,m)).
 */
class ComaPolyDataset
{
public:
    // Construct empty
    ComaPolyDataset( ) = default;

    // Construct directly from a list of file paths
    explicit ComaPolyDataset( const std::vector< std::string >& filePathList )
    {
        readInputFiles( filePathList );
    }

    struct FileMeta
    {
        double referenceRadius{}; //!< From header "R"
        Eigen::VectorXd powersInvRadius; //!< From "PWRS"
        std::vector< std::pair< double, double > > timePeriods; //!< TODO: parsed if available
        int maxDegreeSH{}; //!< From "N(SH)"
        Eigen::Index numRadialTerms{}; //!< From "N(r)"
        Eigen::Index numIntervals{}; //!< From "N(T)"
        std::string sourcePath; //!< File path
    };

    // -------- Reading API --------
    void readInputFiles( const std::vector< std::string >& filePathList )
    {
        clear( );

        const std::size_t n = filePathList.size( );
        numPolyCoefFiles_ = n;

        polyCoefficients_.resize( n );
        SHDegreeAndOrderIndices_.resize( n );
        fileMeta_.resize( n );

        // internal maps sized; built lazily
        nmToColCache_.clear( );
        nmToColCache_.resize( n );

        for(std::size_t fileIdx = 0; fileIdx < n; ++fileIdx)
        {
            const std::string& currentFile = filePathList[ fileIdx ];
            fileMeta_[ fileIdx ].sourcePath = currentFile;

            std::ifstream file( currentFile );
            if(!file.is_open( ))
            {
                std::cerr << "[ERROR] Could not open file '" << currentFile << "'.\n";
                std::exit( EXIT_FAILURE );
            }

            std::string line;
            std::vector< std::string > tokens;

            int maxDegreeSH = 0;
            Eigen::Index numTerms = 0, numCoefs = 0, numRadialTerms = 0, numIntervals = 0;

            Eigen::VectorXd powers; // temp; moved to meta

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
                    double R = std::stod( line.substr( line.find_last_of( " \t" ) + 1 ) );
                    fileMeta_[ fileIdx ].referenceRadius = R;
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
                // OPTIONAL: parse time period lines here if present in your format;
                // set fileMeta_[fileIdx].timePeriods accordingly.
            }

            // ----- Validation -----
            if(numTerms <= 0 || numCoefs <= 0 || powers.size( ) == 0)
            {
                std::cerr << "[ERROR] Header parsing failed in file: " << currentFile << std::endl;
                std::cerr << "  numTerms = " << numTerms
                        << "\n  numCoefs = " << numCoefs
                        << "\n  powersInvRadius.size() = " << powers.size( ) << std::endl;
                std::exit( EXIT_FAILURE );
            }

            fileMeta_[ fileIdx ].maxDegreeSH = maxDegreeSH;
            fileMeta_[ fileIdx ].numRadialTerms = numRadialTerms;
            fileMeta_[ fileIdx ].numIntervals = numIntervals;
            fileMeta_[ fileIdx ].powersInvRadius = powers;

            // --- Allocate containers ---
            Eigen::MatrixXd& currentPolyCoefficients = polyCoefficients_[ fileIdx ];
            Eigen::ArrayXXi& currentShDegreeAndOrder = SHDegreeAndOrderIndices_[ fileIdx ];

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

            // Build (n,m)->col map for this file
            buildNmMap_( fileIdx );
        }
    }

    // -------- Accessors --------
    std::size_t getNumFiles( ) const
    {
        return numPolyCoefFiles_;
    }

    const FileMeta& getFileMeta( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ];
    }

    //! Get the polynomial coefficient matrix of a specific input file.
    /*!
     * Matrix has dimensions (number of radial terms × number of SH coefficients).
     */
    const Eigen::MatrixXd& getPolyCoefficients( std::size_t f ) const
    {
        boundsFile_( f );
        return polyCoefficients_[ f ];
    }

    //! Get the spherical harmonics degree and order indices of a specified input file.
    /*!
     * Each inner vector holds a flattened (2 × N) list of degree and order pairs, where N is the number of coefficients.
     * Index `2*k` contains degree, `2*k+1` contains order.
     */
    const Eigen::ArrayXXi& getSHDegreeAndOrderIndices( std::size_t f ) const
    {
        boundsFile_( f );
        return SHDegreeAndOrderIndices_[ f ];
    }

    //! Get the reference radius for a specific input file.
    /*!
     * Reference radius specifies the radius until which DSMC model is valid.
     */
    double getReferenceRadius( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ].referenceRadius;
    }

    //! Get the inverse-radius powers.
    /*!
     * These powers define the radial dependency of the polynomial basis.
     */
    const Eigen::VectorXd& getPowersInvRadius( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ].powersInvRadius;
    }

    //! Get maximum degree and order of SH expansion for which poly coefficients are valid
    int getMaxDegreeSH( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ].maxDegreeSH;
    }

    Eigen::Index getNumRadialTerms( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ].numRadialTerms;
    }

    Eigen::Index getNumIntervals( std::size_t f ) const
    {
        boundsFile_( f );
        return fileMeta_[ f ].numIntervals;
    }


    /// Get the (numTerms x 1) column vector for a given (n,m).
    Eigen::VectorXd columnForNM( std::size_t f, int n, int m ) const
    {
        auto [ ok, col ] = findColumn_( f, n, m );
        if(!ok) throw std::out_of_range( "columnForNM: (n,m) not found in file" );
        return polyCoefficients_[ f ].col( col );
    }

    /// Get a single value: termIndex in [0, N(r)*N(T)), for given (n,m).
    double value( std::size_t f, Eigen::Index termIndex, int n, int m ) const
    {
        auto [ ok, col ] = findColumn_( f, n, m );
        if(!ok) throw std::out_of_range( "value: (n,m) not found in file" );
        if(termIndex < 0 || termIndex >= polyCoefficients_[ f ].rows( ))
            throw std::out_of_range( "value: termIndex out of range" );
        return polyCoefficients_[ f ]( termIndex, col );
    }

    /// Convenience: (radialIdx, intervalIdx) -> termIndex
    static inline Eigen::Index termIndex( Eigen::Index radialIdx,
                                          Eigen::Index intervalIdx,
                                          Eigen::Index numIntervals )
    {
        return radialIdx * numIntervals + intervalIdx;
    }

    void clear( )
    {
        numPolyCoefFiles_ = 0;
        polyCoefficients_.clear( );
        SHDegreeAndOrderIndices_.clear( );
        fileMeta_.clear( );
        nmToColCache_.clear( );
    }

private:
    // Map key for (n,m)
    struct PairHash
    {
        std::size_t operator()( const std::pair< int, int >& p ) const noexcept
        {
            // Simple hash combine for (n,m)
            return ( static_cast< std::size_t >(p.first) << 32 ) ^ static_cast< std::size_t >(p.second);
        }
    };

    void boundsFile_( std::size_t f ) const
    {
        if(f >= numPolyCoefFiles_) throw std::out_of_range( "file index out of range" );
    }

    void buildNmMap_( std::size_t f )
    {
        boundsFile_( f );
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
        boundsFile_( f );
        const auto& map = nmToColCache_[ f ];
        auto it = map.find( { n, m } );
        if(it == map.end( )) return { false, -1 };
        return { true, it->second };
    }

    // -------- Data --------
    std::size_t numPolyCoefFiles_{ 0 };

    // Per-file containers
    std::vector< Eigen::MatrixXd > polyCoefficients_; // (numTerms x numCoefs)
    std::vector< Eigen::ArrayXXi > SHDegreeAndOrderIndices_; // 2 x numCoefs
    std::vector< FileMeta > fileMeta_;

    // Per-file (n,m)->column cache
    std::vector< std::unordered_map< std::pair< int, int >, int, PairHash > > nmToColCache_;
};



class PolyCoefFileProcessing
{
public:
    /// Construct from a list of PolyCoefficient file paths.
    explicit PolyCoefFileProcessing(std::vector<std::string> filePaths)
        : filePaths_(std::move(filePaths))
    {
        if (filePaths_.empty()) {
            throw std::invalid_argument("PolycoefFileProcessing: file path list is empty.");
        }
    }

    /// Access the raw file list used to build datasets.
    const std::vector<std::string>& getFilePaths() const { return filePaths_; }

    /**
     * \brief Read all PolyCoefficient files and return a populated ComaPolyDataset.
     * \details This simply forwards the stored file list to ComaPolyDataset’s reader/ctor.
     */
    ComaPolyDataset createPolyCoefDataset() const
    {
        return ComaPolyDataset(filePaths_);
    }

    /**
     * \brief (Stub) Create SH/Stokes CSV files for given settings.
     * \param nmax           Maximum requested SH degree
     * \param mmax           Maximum requested SH order
     * \param radii          Vector of radii where SH will be evaluated
     * \param solLongitudes  Vector of solar longitudes where SH will be evaluated
     * \param outputDir      Directory to write output files
     *
     * TODO: Implement using ComaPolyDataset → ComaStokesDataset pipeline + CSV writer.
     */
    void createSHFiles(int /*nmax*/,
                       int /*mmax*/,
                       const std::vector<double>& /*radii*/,
                       const std::vector<double>& /*solLongitudes*/,
                       const std::string& /*outputDir*/) const
    {
        throw std::logic_error("PolycoefFileProcessing::createSHFiles: not implemented yet.");
    }

    /**
     * \brief (Stub) Compute and return a populated ComaStokesDataset.
     * \param nmax           Maximum requested SH degree
     * \param mmax           Maximum requested SH order
     * \param radii          Vector of radii where SH will be evaluated
     * \param solLongitudes  Vector of solar longitudes where SH will be evaluated
     *
     * TODO: Implement by transforming ComaPolyDataset → ComaStokesDataset and filling coefficients.
     */
    ComaStokesDataset createSHDataset(int /*nmax*/,
                                      int /*mmax*/,
                                      const std::vector<double>& /*radii*/,
                                      const std::vector<double>& /*solLongitudes*/) const
    {
        throw std::logic_error("PolycoefFileProcessing::createSHDataset: not implemented yet.");
    }

private:
    std::vector<std::string> filePaths_;
};



//  List of wind models available in simulations
/*
 *  List of wind models available in simulations. Wind models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum WindModelTypes { constant_wind_model, custom_wind_model, coma_wind_model };

//  Class for providing settings for wind model.
/*
 *  Class for providing settings for automatic wind model creation. This class is a
 *  functional (base) class for settings of wind models that require no information in
 *  addition to their type. Wind model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(WindModelSettings.__docstring__)
class WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windModelType Type of wind model that is to be created
     */
    WindModelSettings( const WindModelTypes windModelType,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        windModelType_( windModelType ), associatedFrame_( associatedFrame )
    {
    }

    //  Destructor
    virtual ~WindModelSettings( )
    {
    }

    //  Function to retrieve type of wind model that is to be created
    /*
     * Function to retrieve type of wind model that is to be created
     * \return Type of wind model that is to be created
     */
    WindModelTypes getWindModelType( )
    {
        return windModelType_;
    }

    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

protected:
    //  Type of wind model that is to be created
    WindModelTypes windModelType_;

    reference_frames::AerodynamicsReferenceFrames associatedFrame_;
};

class ConstantWindModelSettings : public WindModelSettings
{
public:
    ConstantWindModelSettings( const Eigen::Vector3d constantWindVelocity,
                               const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModelSettings( constant_wind_model, associatedFrame ), constantWindVelocity_( constantWindVelocity )
    {
    }

    Eigen::Vector3d getConstantWindVelocity( )
    {
        return constantWindVelocity_;
    }

private:
    Eigen::Vector3d constantWindVelocity_;
};

//  Class to define settings for a custom, user-defined, wind model
class CustomWindModelSettings : public WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     */
    CustomWindModelSettings( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                             const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModelSettings( custom_wind_model, associatedFrame ), windFunction_( windFunction )
    {
    }

    //  Destructor
    ~CustomWindModelSettings( )
    {
    }

    //  Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
     * \return Function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > getWindFunction( )
    {
        return windFunction_;
    }

    //  Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
     * \param windFunction New function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    void setWindFunction( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction )
    {
        windFunction_ = windFunction;
    }

protected:
    //  Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

//  List of atmosphere models available in simulations
/*
 *  List of atmosphere models available in simulations. Atmosphere models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum AtmosphereTypes
{
    exponential_atmosphere,
    custom_constant_temperature_atmosphere,
    tabulated_atmosphere,
    nrlmsise00,
    mars_dtm_atmosphere,
    scaled_atmosphere,
    coma_model
};

//  Class for providing settings for atmosphere model.
/*
 *  Class for providing settings for automatic atmosphere model creation. This class is a
 *  functional (base) class for settings of atmosphere models that require no information in
 *  addition to their type. Atmosphere model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(AtmosphereSettings.__docstring__)
class AtmosphereSettings
{
public:
    //  Constructor, sets type of atmosphere model.
    /*
     *  Constructor, sets type of atmosphere model. Settings for atmosphere models requiring
     *  additional information should be defined in a derived class.
     *  \param atmosphereType Type of atmosphere model that is to be created.
     */
    AtmosphereSettings( const AtmosphereTypes atmosphereType ):
        atmosphereType_( atmosphereType )
    {
    }

    //  Destructor
    virtual ~AtmosphereSettings( )
    {
    }

    //  Function to return type of atmosphere model that is to be created.
    /*
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( )
    {
        return atmosphereType_;
    }

    //  Function to return settings for the atmosphere's wind model.
    /*
     *  Function to return settings for the atmosphere's wind model.
     *  \return Settings for the atmosphere's wind model.
     */
    std::shared_ptr< WindModelSettings > getWindSettings( )
    {
        return windSettings_;
    }

    //  Function to (re)set settings for the atmosphere's wind model.
    /*
     *  Function to (re)set settings for the atmosphere's wind model.
     *  \param windSettings Settings for the atmosphere's wind model.
     */
    void setWindSettings( const std::shared_ptr< WindModelSettings > windSettings )
    {
        windSettings_ = windSettings;
    }

private:
    //   Type of atmosphere model that is to be created.
    AtmosphereTypes atmosphereType_;

    //  Settings for the atmosphere's wind model.
    std::shared_ptr< WindModelSettings > windSettings_;
};

//  AtmosphereSettings for defining an exponential atmosphere.
//! @get_docstring(ExponentialAtmosphereSettings.__docstring__)
class ExponentialAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor, taking full atmosphere model parameters.
     *  \param densityScaleHeight Scale height for density profile of atmosphere.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at ground level.
     *  \param specificGasConstant Specific gas constant for (constant) atmospheric chemical
     *  composition.
     *  \param ratioOfSpecificHeats Ratio of specific heats for (constant) atmospheric chemical
     *  composition.
     */
    ExponentialAtmosphereSettings( const double densityScaleHeight,
                                   const double constantTemperature,
                                   const double densityAtZeroAltitude,
                                   const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                   const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( exponential_atmosphere ), densityScaleHeight_( densityScaleHeight ),
        constantTemperature_( constantTemperature ), densityAtZeroAltitude_( densityAtZeroAltitude ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        bodyWithPredefinedExponentialAtmosphere_( undefined_body )
    {
    }

    //  Default constructor.
    /*
     *  Default constructor, taking only the name of the body for which to load the predefined
     *  exponential atmosphere model parameters.
     *  \param bodyWithPredefinedExponentialAtmosphere Enumeration denoting the name of the body for which the
     *  predefined atmosphere model is to be loaded.
     */
    ExponentialAtmosphereSettings( const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere ):
        AtmosphereSettings( exponential_atmosphere ), bodyWithPredefinedExponentialAtmosphere_( bodyWithPredefinedExponentialAtmosphere )
    {
        // Check that the body name inserted is available
        switch(bodyWithPredefinedExponentialAtmosphere)
        {
            case earth:
            case mars:
                // all is good
                break;
            default:
                throw std::runtime_error(
                        "Error while creating exponential atmosphere. The body name provided "
                        "does not match any predefined atmosphere model. Available models for: "
                        "Earth, Mars." );
        }
    }

    //  Function to return scale heigh for density profile of atmosphere.
    /*
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( )
    {
        return densityScaleHeight_;
    }

    //  Function to return constant atmospheric temperature.
    /*
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Function to return atmospheric density at ground level.
    /*
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( )
    {
        return densityAtZeroAltitude_;
    }

    //  Function to return specific gas constant for (constant) atmospheric chemical
    /*
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats for (constant) atmospheric chemical
    /*
     *  Function to return ratio of specific heats for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return the name of the body for which to load the predefined
    //  atmosphere model parameters.
    BodiesWithPredefinedExponentialAtmospheres getBodyName( )
    {
        return bodyWithPredefinedExponentialAtmosphere_;
    }

private:
    //  Scale heigh for density profile of atmosphere.
    double densityScaleHeight_;

    //  Constant atmospheric temperature.
    double constantTemperature_;

    //  Atmospheric density at ground level.
    double densityAtZeroAltitude_;

    //  Specific gas constant for (constant) atmospheric chemical
    double specificGasConstant_;

    //  Ratio of specific heats for (constant) atmospheric chemical
    double ratioOfSpecificHeats_;

    //  Enumeration denoting the name of the body for which to load the predefined
    //  atmosphere model.
    BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere_;
};

//  AtmosphereSettings for defining custom constant temperature atmosphere.
class CustomConstantTemperatureAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Typedef for density function.
    typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;

    //  Default constructor.
    /*
     *  Default constructor setting all parameters manually.
     *  \param densityFunction Function to retireve the density at the current altitude,
     *      longitude, latitude and time.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    CustomConstantTemperatureAtmosphereSettings( const DensityFunction& densityFunction,
                                                 const double constantTemperature,
                                                 const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                                 const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunction_( densityFunction ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    {
    }

    //  Constructor.
    /*
     *  Constructor which uses one of the built-in density functions as input.
     *  \param densityFunctionType Enumeration denoting which density function to implement.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param modelSpecificParameters Vector of parameters to be used to set up the density
     *      function. Both meaning and number of parameters depends on the model.
     */
    CustomConstantTemperatureAtmosphereSettings( const AvailableConstantTemperatureAtmosphereModels densityFunctionType,
                                                 const double constantTemperature,
                                                 const double specificGasConstant,
                                                 const double ratioOfSpecificHeats,
                                                 const std::vector< double >& modelSpecificParameters ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunctionType_( densityFunctionType ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ), modelSpecificParameters_( modelSpecificParameters )
    {
    }

    //  Get the function to compute the density at the current conditions.
    /*
     *  Get the function to compute the density at the current conditions.
     *  \return Function to compute the density at the current conditions.
     */
    DensityFunction getDensityFunction( )
    {
        return densityFunction_;
    }

    //  Get the type of function to compute the density at the current conditions.
    /*
     *  Get the type of function to compute the density at the current conditions.
     *  \return Type of function to compute the density at the current conditions.
     */
    AvailableConstantTemperatureAtmosphereModels getDensityFunctionType( )
    {
        return densityFunctionType_;
    }

    //  Get constant temperature.
    /*
     *  Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     *  Kelvin.
     *  \return constantTemperature Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Get specific gas constant.
    /*
     *  Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Get ratio of specific heats.
    /*
     *  Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Get model specific parameters.
    /*
     *  Get model specific parameters.
     *  \return Vector of parameters to be used to set up the density function.
     */
    std::vector< double > getModelSpecificParameters( )
    {
        return modelSpecificParameters_;
    }

private:
    //  Function to compute the density at the current conditions.
    /*
     *  Function to compute the density at the current conditions. Note that the independent
     *  variables need to be altitude, longitude, latitude and time, in this precise order.
     */
    DensityFunction densityFunction_;

    //  Enumeration denoting which density function to implement.
    /*
     *  Enumeration denoting which density function to implement
     */
    AvailableConstantTemperatureAtmosphereModels densityFunctionType_;

    //  Constant temperature.
    /*
     *  The atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    const double constantTemperature_;

    //  Specific gas constant.
    /*
     *  Specific gas constant of the atmosphere, its value is assumed constant, due to the
     *  assumption of constant atmospheric composition.
     */
    const double specificGasConstant_;

    //  Ratio of specific heats at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    const double ratioOfSpecificHeats_;

    //  Vector of parameters to be used to set up the density function.
    /*
     *  Vector of parameters to be used to set up the density function. Both meaning and number of parameters depends on the model.
     */
    std::vector< double > modelSpecificParameters_;
};

//  AtmosphereSettings for defining an NRLMSISE00 atmosphere reading space weather data from a text file.
class NRLMSISE00AtmosphereSettings : public AtmosphereSettings
{
public:
    //  Constructor.
    /*
     *  Constructor.
     *  \param spaceWeatherFile File containing space weather data, as in
     *  https://celestrak.com/SpaceData/sw19571001.txt
     */
    NRLMSISE00AtmosphereSettings( const std::string& spaceWeatherFile,
                                  const bool useStormConditions = false,
                                  const bool useAnomalousOxygen = true ):
        AtmosphereSettings( nrlmsise00 ), spaceWeatherFile_( spaceWeatherFile ), useStormConditions_( useStormConditions ),
        useAnomalousOxygen_( useAnomalousOxygen )
    {
    }

    //  Function to return file containing space weather data.
    /*
     *  Function to return file containing space weather data.
     *  \return Filename containing space weather data.
     */
    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

    //  Function to return geomagnetic activity setting.
    /*
     *  Function to return geomagnetic activity setting.
     *  \return Geomagnetic activity value.
     */
    bool getUseStormConditions( )
    {
        return useStormConditions_;
    }

    bool getUseAnomalousOxygen( )
    {
        return useAnomalousOxygen_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;

    //  Boolean denoting how to deal with geomagnetic activity setting.
    /*
     *  Controls the geomagnetic activity behavior.
     *  If true, it uses full vector of Ap values (switch for geomagnetic activity in nrlmsise to -1).
     *  If false, it only uses daily value of Ap (switch for geomagnetic activity in nrlmsise to 0).
     */
    bool useStormConditions_;

    bool useAnomalousOxygen_;
};

class MarsDtmAtmosphereSettings : public AtmosphereSettings
{
public:
    MarsDtmAtmosphereSettings( const std::string& spaceWeatherFile = "" ):
        AtmosphereSettings( mars_dtm_atmosphere ), spaceWeatherFile_( spaceWeatherFile )
    {
    }

    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;
};


//  AtmosphereSettings for defining an atmosphere with tabulated data from file.
// //! @get_docstring(TabulatedAtmosphereSettings.__docstring__)
class TabulatedAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::map< int, std::string >& atmosphereTableFile,
            const std::vector< AtmosphereIndependentVariables >& independentVariablesNames = { altitude_dependent_atmosphere },
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = {},
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ), boundaryHandling_( boundaryHandling ),
        defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with single boundary handling parameters.
    /*
     *  Constructor with single boundary handling parameters. The specifier is assumed to be the same for
     *  each (in)dependent variable.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const double specificGasConstant,
                                 const double ratioOfSpecificHeats,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                specificGasConstant,
                ratioOfSpecificHeats,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor compatible with old version.
    /*
     *  Constructor compatible with old version.
     *  \param atmosphereTableFile File containing atmospheric properties. The file name of the atmosphere table. The file
     *      should contain four columns of data, containing altitude (first column), and the associated density, pressure and
     *      density values in the second, third and fourth columns.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::string& atmosphereTableFile,
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
            const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings( { { 0, atmosphereTableFile } },
                                     { altitude_dependent_atmosphere },
                                     dependentVariablesNames,
                                     specificGasConstant,
                                     ratioOfSpecificHeats,
                                     { boundaryHandling },
                                     std::vector< std::vector< std::pair< double, double > > >(
                                             dependentVariablesNames.size( ),
                                             std::vector< std::pair< double, double > >(
                                                     1,
                                                     std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector).
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< double >& defaultExtrapolationValue ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling )
    {
        // Assign default values
        defaultExtrapolationValue_.resize( dependentVariablesNames.size( ) );
        for(unsigned int i = 0; i < dependentVariablesNames.size( ); i++)
        {
            for(unsigned int j = 0; j < independentVariablesNames.size( ); j++)
            {
                if(boundaryHandling_.at( j ) == interpolators::use_default_value ||
                    boundaryHandling_.at( j ) == interpolators::use_default_value_with_warning)
                {
                    defaultExtrapolationValue_.at( i ).push_back(
                            std::make_pair( defaultExtrapolationValue.at( i ), defaultExtrapolationValue.at( i ) ) );
                }
                else
                {
                    defaultExtrapolationValue_.at( i ).push_back( std::make_pair( IdentityElement::getAdditionIdentity< double >( ),
                                                                                  IdentityElement::getAdditionIdentity< double >( ) ) );
                }
            }
        }
    }

    //  Constructor with no specific gas constant nor ratio of specific heats, and with
    //  single boundary handling parameters.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector). Only one boundary handling parameter is specified, which is then repeated for
     *  dimension.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                1.4,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Map of filenames containing atmospheric properties.
     */
    std::map< int, std::string > getAtmosphereFile( )
    {
        return atmosphereFile_;
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Filename containing atmospheric properties.
     */
    std::string getAtmosphereFile( const unsigned int fileIndex )
    {
        return atmosphereFile_.at( fileIndex );
    }

    //  Function to return independent variables names.
    /*
     *  Function to return independent variables names.
     *  \return Independent variables.
     */
    std::vector< AtmosphereIndependentVariables > getIndependentVariables( )
    {
        return independentVariables_;
    }

    //  Function to return dependent variables names.
    /*
     *  Function to return dependent variables names.
     *  \return Dependent variables.
     */
    std::vector< AtmosphereDependentVariables > getDependentVariables( )
    {
        return dependentVariables_;
    }

    //  Function to return specific gas constant of the atmosphere.
    /*
     *  Function to return specific gas constant of the atmosphere.
     *  \return Specific gas constant of the atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats of the atmosphere.
    /*
     *  Function to return ratio of specific heats of the atmosphere.
     *  \return Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return boundary handling method.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< interpolators::BoundaryInterpolationType > getBoundaryHandling( )
    {
        return boundaryHandling_;
    }

    //  Function to return default extrapolation value.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< std::vector< std::pair< double, double > > > getDefaultExtrapolationValue( )
    {
        return defaultExtrapolationValue_;
    }

private:
    //  File containing atmospheric properties.
    /*
     *  File containing atmospheric properties, file should contain
     *  columns of atmospheric data with at least density, pressure and temperature,
     *  (whose order is specified in dependentVariables), and with at least one
     *  indendent variables.
     */
    std::map< int, std::string > atmosphereFile_;

    //  A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //  A vector of strings containing the names of the variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //  Specific gas constant of the atmosphere.
    /*
     * Specific gas constant of the atmosphere.
     */
    double specificGasConstant_;

    //  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double ratioOfSpecificHeats_;

    //  Behavior of interpolator when independent variable is outside range.
    /*
     *  Behavior of interpolator when independent variable is outside range.
     */
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling_;

    //  Default value to be used for extrapolation.
    /*
     *  Default value to be used for extrapolation.
     */
    std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValue_;
};

class ScaledAtmosphereSettings : public AtmosphereSettings
{
public:
    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                              const double scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( [ = ]( const double ) {
            return scaling;
        } ),
        isScalingAbsolute_( isScalingAbsolute )
    {
    }

    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                              const std::function< double( const double ) > scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( scaling ), isScalingAbsolute_( isScalingAbsolute )
    {
    }

    std::shared_ptr< AtmosphereSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< double( const double ) > getScaling( )
    {
        return scaling_;
    }

    bool getIsScalingAbsolute( )
    {
        return isScalingAbsolute_;
    }

protected:
    std::shared_ptr< AtmosphereSettings > baseSettings_;

    std::function< double( const double ) > scaling_;

    bool isScalingAbsolute_;
};


class ComaSettings final : public AtmosphereSettings
{
public:
    // Constructor taking PolyCoefficient file list as input.
    ComaSettings( std::vector< std::string > polyCoefFileList, int requestedDegree, int requestedOrder ):
        AtmosphereSettings( coma_model ),
        polyCoefData_( std::move( polyCoefFileList ) ),
        requestedDegree_( requestedDegree ),
        requestedOrder_( requestedOrder )
    {
    }

    //! Get the poly coefficient data handler object
    ComaPolyDataset getPolyDataset( ) const
    {
        return polyCoefData_;
    }


    //! Get the list of polynomial coefficient matrices, one per input file.
    /*!
     * Each matrix has dimensions (number of radial terms × number of SH coefficients) for one file.
     */
    const std::vector< Eigen::MatrixXd >& getPolyCoefficients( ) const
    {
        return polyCoefficients_;
    }

    //! Get the list of spherical harmonics degree and order indices for each input file.
    /*!
     * Each inner vector holds a flattened (2 × N) list of degree and order pairs, where N is the number of coefficients.
     * Index `2*k` contains degree, `2*k+1` contains order.
     */
    const std::vector< Eigen::ArrayXXi >& getSHDegreeAndOrder( ) const
    {
        return SHDegreeAndOrderIndices_;
    }

    //! Get the list of reference radii used in each input file.
    /*!
     * Each inner vector corresponds to one input file. Typically contains one element per file.
     */
    const std::vector< double >& getReferenceRadius( ) const
    {
        return referenceRadius_;
    }

    //! Get the list of inverse-radius powers used in each input file.
    /*!
     * These powers define the radial dependency of the polynomial basis. Each inner vector corresponds to one file.
     */
    const std::vector< Eigen::VectorXd >& getPowersInvRadius( ) const
    {
        return powersInvRadius_;
    }

    const std::vector< std::vector< double > >& getTimePeriods( ) const
    {
        return TimePeriods_;
    }

    // Get requested Degree
    const int& getRequestedDegree( ) const
    {
        return requestedDegree_;
    }

    // Get requested Order
    const int& getRequestedOrder( ) const
    {
        return requestedOrder_;
    }

private:
    // Spherical harmonics model.
    // SphericalHarmonicsModel sphericalHarmonicsModel_ = customModel;

    // poly coefficient data handler
    ComaPolyDataset polyCoefData_;

    // Path of loaded coma input files.
    std::vector< std::string > polyCoefFileList_;

    // Maximum degree used for computation
    int requestedDegree_;

    // requested order used for computation
    int requestedOrder_;

    // number of Poly Coefficient input files used
    int numPolyCoefFiles_;

    // PolyCoefficients of input files
    std::vector< Eigen::MatrixXd > polyCoefficients_;

    // Spherical Harmonics Degree and Order indices of input tables
    std::vector< Eigen::ArrayXXi > SHDegreeAndOrderIndices_;

    // Reference radius of input table
    std::vector< double > referenceRadius_;

    // Power Inverse Radius of input tables
    std::vector< Eigen::VectorXd > powersInvRadius_;

    // Time periods where input tables are valid
    std::vector< std::vector< double > > TimePeriods_;


    //  A function reading PolyCoefficient files
    /*
     * A function reading and processing a PolyCoefficient file list. Extracts:
     * - polyCoefficients, used to compute stokes coefficients
     * - SH degree and order indices
     * - reference radius, designating the radius from which a 1/r^2 law will be used to compute the density
     * - power inverse radius, used to compute stokes coefficients
     * - time periods, for which each file is valid
     * */
    void readInputFiles( const std::vector< std::string >& filePathList_ )
    {
        const std::size_t n = filePathList_.size( );

        // Pre-allocate
        numPolyCoefFiles_ = n;
        polyCoefficients_.resize( n );
        SHDegreeAndOrderIndices_.resize( n );
        referenceRadius_.resize( n );
        powersInvRadius_.resize( n );
        TimePeriods_.resize( n );

        for(std::size_t fileIdx = 0; fileIdx < n; ++fileIdx)
        {
            const std::string& currentFile = filePathList_[ fileIdx ];
            Eigen::MatrixXd& currentPolyCoefficients = polyCoefficients_[ fileIdx ];
            Eigen::ArrayXXi& currentShDegreeAndOrder = SHDegreeAndOrderIndices_[ fileIdx ];
            double& currentReferenceRadius = referenceRadius_[ fileIdx ];
            Eigen::VectorXd& currentPowersInvRadius = powersInvRadius_[ fileIdx ];
            std::vector< double >& currentTimePeriod = TimePeriods_[ fileIdx ]; // TODO: read current time period

            std::ifstream file( currentFile );
            if(!file.is_open( ))
            {
                std::cerr << "[ERROR] Could not open file '" << currentFile << "'." << std::endl;
                std::exit( EXIT_FAILURE );
            }

            std::string line;
            std::vector< std::string > tokens;
            int maxDegreeSH = 0;
            Eigen::Index numTerms = 0, numCoefs = 0, numRadialTerms = 0, numIntervals = 0;

            // Parse header
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
                    boost::split( tokens, tail, boost::is_any_of( ", \t" ), boost::token_compress_on );
                    std::size_t start = ( !tokens.empty( ) && !std::all_of( tokens[ 0 ].begin( ), tokens[ 0 ].end( ), ::isdigit ) ) ? 1 : 0;
                    auto count = static_cast< Eigen::Index >(tokens.size( ) - start);
                    currentPowersInvRadius.resize( count );
                    for(Eigen::Index j = 0; j < count; ++j)
                        currentPowersInvRadius[ j ] = std::stod( tokens[ start + j ] );
                }
                else if(boost::iequals( key, "R" ))
                {
                    double R = std::stod( line.substr( line.find_last_of( " \t" ) + 1 ) );
                    currentReferenceRadius = R;
                }
                else if(line.find( "N(r)" ) != std::string::npos && line.find( "N(T)" ) != std::string::npos)
                {
                    std::string content = line.substr( 1 ); // strip '#'
                    boost::trim( content );
                    boost::split( tokens, content, boost::is_any_of( ", \t" ), boost::token_compress_on );

                    // Find last two tokens (assumed to be the numbers)
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
            }

            // Validation
            if(numTerms <= 0 || numCoefs <= 0 || currentPowersInvRadius.size( ) == 0)
            {
                std::cerr << "[ERROR] Header parsing failed in file: " << currentFile << std::endl;
                std::cerr << "  numTerms = " << numTerms << "\n  numCoefs = " << numCoefs
                        << "\n  powersInvRadius.size() = " << currentPowersInvRadius.size( ) << std::endl;
                std::exit( EXIT_FAILURE );
            }

            currentPolyCoefficients.resize( numTerms, numCoefs );
            currentShDegreeAndOrder.resize( 2, numCoefs );

            // Read data block
            Eigen::Index coefIndex = -1;
            do
            {
                boost::trim( line );
                if(line.empty( ) || line[ 0 ] == '#') continue;

                boost::split( tokens, line, boost::is_any_of( ", \t" ), boost::token_compress_on );
                if(static_cast< Eigen::Index >(tokens.size( )) == numTerms + 2)
                {
                    ++coefIndex;
                    currentShDegreeAndOrder( 0, coefIndex ) = std::stoi( tokens[ 0 ] );
                    currentShDegreeAndOrder( 1, coefIndex ) = std::stoi( tokens[ 1 ] );
                    for(Eigen::Index j = 0; j < numTerms; ++j)
                        currentPolyCoefficients( j, coefIndex ) = std::stod( tokens[ j + 2 ] );
                }
            } while(std::getline( file, line ));

            file.close( );
        }
    }
};


//! @get_docstring(exponentialAtmosphereSettings,2)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const double densityScaleHeight,
        const double densityAtZeroAltitude,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< ExponentialAtmosphereSettings >(
            densityScaleHeight,
            constantTemperature,
            densityAtZeroAltitude,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(exponentialAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const double densityScaleHeight,
                                                                            const double densityAtZeroAltitude )
{
    return std::make_shared< ExponentialAtmosphereSettings >( densityScaleHeight, TUDAT_NAN, densityAtZeroAltitude, TUDAT_NAN, TUDAT_NAN );
}

//! @get_docstring(exponentialAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const std::string& bodyName )
{
    BodiesWithPredefinedExponentialAtmospheres bodyId;
    if(bodyName == "Earth")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::earth;
    }
    else if(bodyName == "Mars")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::mars;
    }
    else
    {
        throw std::runtime_error(
                "Error while creating exponential atmosphere. The body name provided "
                "does not match any predefined atmosphere model. Available models for: "
                "Earth, Mars." );
    }
    return std::make_shared< ExponentialAtmosphereSettings >( bodyId );
}

//! @get_docstring(nrlmsise00AtmosphereSettings)
inline std::shared_ptr< AtmosphereSettings > nrlmsise00AtmosphereSettings( const std::string dataFile = paths::getSpaceWeatherDataPath( ) +
                                                                                   "/sw19571001.txt",
                                                                           const bool useStormConditions = true,
                                                                           const bool useAnomalousOxygen = true )
{
    return std::make_shared< NRLMSISE00AtmosphereSettings >( dataFile, useStormConditions, useAnomalousOxygen );
}

inline std::shared_ptr< AtmosphereSettings > marsDtmAtmosphereSettings( )
{
    return std::make_shared< MarsDtmAtmosphereSettings >( );
}


typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;
//! @get_docstring(customConstantTemperatureAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const std::function< double( const double ) > densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    DensityFunction fullDensityFunction = [ = ]( const double altitude, const double, const double, const double ) {
        return densityFunction( altitude );
    };
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            fullDensityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(customConstantTemperatureAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const DensityFunction densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            densityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(scaledAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const std::function< double( const double ) > scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const double scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

inline std::shared_ptr< AtmosphereSettings > tabulatedAtmosphereSettings(
        const std::string& atmosphereTableFile,
        const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                       pressure_dependent_atmosphere,
                                                                                       temperature_dependent_atmosphere },
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile,
                                                            dependentVariablesNames,
                                                            specificGasConstant,
                                                            ratioOfSpecificHeats,
                                                            interpolators::throw_exception_at_boundary );
}


inline std::shared_ptr< AtmosphereSettings > comaSettings(
        std::vector< std::string > polyList,
        const int requestedDegree,
        const int requestedOrder )
{
    return std::make_shared< ComaSettings >( polyList, requestedDegree, requestedOrder );
}


//! @get_docstring(customWindModelSettings)
inline std::shared_ptr< WindModelSettings > customWindModelSettings(
        const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< CustomWindModelSettings >( windFunction, associatedFrame );
}

//! @get_docstring(constantWindModelSettings)
inline std::shared_ptr< WindModelSettings > constantWindModelSettings(
        const Eigen::Vector3d constantWindVelocity,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< ConstantWindModelSettings >( constantWindVelocity, associatedFrame );
}

//  Function to create a wind model.
/*
 *  Function to create a wind model based on model-specific settings for the wind model.
 *  \param windSettings Settings for the wind model that is to be created, defined
 *  a pointer to an object of class (derived from) WindModelSettings.
 *  \param body Name of the body for which the wind model is to be created.
 *  \return Wind model created according to settings in windSettings.
 */
std::shared_ptr< aerodynamics::WindModel > createWindModel( const std::shared_ptr< WindModelSettings > windSettings,
                                                            const std::string& body );

//  Function to create an atmosphere model.
/*
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel( const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
                                                                        const std::string& body,
                                                                        const SystemOfBodies& bodies );
} // namespace simulation_setup
} // namespace tudat

#endif  // TUDAT_CREATEATMOSPHEREMODEL_H