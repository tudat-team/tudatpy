#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/io/basicInputOutput.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace tudat
{
namespace unit_tests
{
using namespace simulation_setup;

// ==================== Test Fixtures ====================

struct PolyDatasetConfig {
    std::string label;
    boost::filesystem::path polyFile;
    boost::filesystem::path stokesFileLon0;
    boost::filesystem::path stokesFileLon30;
    int maxDegree;
    int maxOrder;
};

struct TestDataPaths {
    boost::filesystem::path testFile;
    boost::filesystem::path rotatedTestFile;
    boost::filesystem::path dataDir;
    boost::filesystem::path outputDir;
    boost::filesystem::path stokesFileLegacyLon0;
    boost::filesystem::path stokesFileLegacyLon30;
    boost::filesystem::path stokesFileRotatedLon0;
    boost::filesystem::path stokesFileRotatedLon30;

    TestDataPaths( )
    {
        dataDir = boost::filesystem::path( tudat::paths::getTudatTestDataPath( ) ) / "coma";
        const boost::filesystem::path polyDir = dataDir / "density" / "polynomial";
        const boost::filesystem::path stokesDir = dataDir / "density" / "stokes";

        testFile = polyDir / "input_poly_coef_test_file.txt";
        rotatedTestFile = polyDir / "input_poly_coef_test_file_rotated.txt";
        outputDir = dataDir / "output" / "density";

        stokesFileLegacyLon0 = stokesDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-000_04km.txt";
        stokesFileLegacyLon30 = stokesDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-030_10km.txt";
        stokesFileRotatedLon0 = stokesDir / "SH-d30-rp3-fft12__r_cometFixed_ep10-000_04km_rotated.txt";
        stokesFileRotatedLon30 = stokesDir / "SH-d30-rp3-fft12__r_cometFixed_ep10-030_10km_rotated.txt";

        // Ensure required files exist
        const std::vector< boost::filesystem::path > requiredFiles = {
            testFile, rotatedTestFile, stokesFileLegacyLon0, stokesFileLegacyLon30, stokesFileRotatedLon0, stokesFileRotatedLon30
        };

        for( const auto& path : requiredFiles )
        {
            if( !boost::filesystem::exists( path ) )
            {
                throw std::runtime_error( "Test data file not found: " + path.string( ) );
            }
        }

        // Clean and create output directory
        if( boost::filesystem::exists( outputDir ) )
        {
            boost::filesystem::remove_all( outputDir );
        }
        boost::filesystem::create_directories( outputDir );
    }

    std::vector< PolyDatasetConfig > getPolyDatasetConfigs( ) const
    {
        return { { "legacy", testFile, stokesFileLegacyLon0, stokesFileLegacyLon30, 10, 10 },
                 { "rotated", rotatedTestFile, stokesFileRotatedLon0, stokesFileRotatedLon30, 30, 30 } };
    }

    ~TestDataPaths( )
    {
        // Optional: cleanup output dir after tests
        // boost::filesystem::remove_all(outputDir);
    }
};

// Expected values from your original tests
struct ExpectedPolyValues {
    static constexpr int numTerms = 48;
    static constexpr int numCoeffs = 121;
    static constexpr double refRadius = 10000.0;
    static constexpr int maxDegree = 10;
    static constexpr int maxOrder = 10;

    // Sample poly coefficient values
    static constexpr double polyCoef_0_0 = 6.262302500423528E+02;
    static constexpr double polyCoef_3_1 = -4.459813047951577E-02;
    static constexpr double polyCoef_47_120 = -1.049515320967208E+02;
    static constexpr double polyCoef_10_22 = 1.287417812579956E-01;
};

// C++14 requires definitions for static constexpr members that are ODR-used
constexpr int ExpectedPolyValues::numTerms;
constexpr int ExpectedPolyValues::numCoeffs;
constexpr double ExpectedPolyValues::refRadius;
constexpr int ExpectedPolyValues::maxDegree;
constexpr int ExpectedPolyValues::maxOrder;
constexpr double ExpectedPolyValues::polyCoef_0_0;
constexpr double ExpectedPolyValues::polyCoef_3_1;
constexpr double ExpectedPolyValues::polyCoef_47_120;
constexpr double ExpectedPolyValues::polyCoef_10_22;

// Helper function to read Stokes coefficient data from test files
struct StokesTestData {
    double solarLongitude;         // degrees
    double radius;                 // meters
    Eigen::MatrixXd cosineCoeffs;  // (degree+1) x (order+1)
    Eigen::MatrixXd sineCoeffs;    // (degree+1) x (order+1)

    static StokesTestData readFromFile( const std::string& filepath, int maxDegree = 10 )
    {
        StokesTestData data;
        std::ifstream file( filepath );
        if( !file.is_open( ) )
        {
            throw std::runtime_error( "Cannot open test data file: " + filepath );
        }

        std::string line;

        // Read first line: solar longitude and radius
        std::getline( file, line );
        std::istringstream iss( line );
        std::string dummy;
        iss >> dummy >> dummy >> dummy;  // skip "# sol r"
        iss >> data.solarLongitude >> data.radius;
        data.radius *= 1000.0;  // Convert km to meters

        // Read second line: header (GM R values)
        std::getline( file, line );

        // Initialize coefficient matrices
        data.cosineCoeffs = Eigen::MatrixXd::Zero( maxDegree + 1, maxDegree + 1 );
        data.sineCoeffs = Eigen::MatrixXd::Zero( maxDegree + 1, maxDegree + 1 );

        // Read coefficient data
        int degree, order;
        double cosine, sine;
        while( file >> degree >> order >> cosine >> sine )
        {
            if( degree <= maxDegree && order <= degree )
            {
                data.cosineCoeffs( degree, order ) = cosine;
                data.sineCoeffs( degree, order ) = sine;
            }
        }

        file.close( );
        return data;
    }
};

// ==================== Data Model Tests ====================

/**
 * @brief Test suite for ComaStokesDataset data structure operations.
 *
 * Tests the core data model functionality including dataset creation,
 * coefficient storage/retrieval, and bounds checking. These tests verify
 * that the underlying data structure correctly manages spherical harmonicple
 * coefficient data.
 */
BOOST_AUTO_TEST_SUITE( test_data_models )

/**
 * @brief Verifies ComaStokesDataset can be created with proper dimensions.
 *
 * Input/Setup:
 * - Creates a dataset with 2 file metadata entries, 3 radii (+1 ref), 4 longitudes, and nmax=10
 *
 * Expected Behavior:
 * - Dataset reports correct dimensions for files, radii, longitudes, and coefficients
 * - Coefficient get/set operations work correctly and return exact values
 * - Block access returns properly sized coefficient blocks
 * - Coefficient matrices have correct dimensions and store values accurately
 */
BOOST_AUTO_TEST_CASE( test_stokes_dataset_creation )
{
    // Test pure data model creation
    const std::vector< ComaStokesDataset::FileMeta > files = { { 2.015e9, 2.0150864e9, "test_file_1", 10000.0 },
                                                               { 2.016e9, 2.0160864e9, "test_file_2", 10000.0 } };
    const std::vector< double > radii = { 1000.0, 2000.0, 3000.0 };
    const std::vector< double > lons = { 0.0, 30.0, 60.0, 90.0 };
    constexpr int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create( files, radii, lons, nmax );

    // Verify metadata
    BOOST_CHECK_EQUAL( dataset.nFiles( ), 2 );
    BOOST_CHECK_EQUAL( dataset.nRadii( ), 4 );  // 3 from radii vector + added ref radius
    BOOST_CHECK_EQUAL( dataset.nLongitudes( ), 4 );
    BOOST_CHECK_EQUAL( dataset.nmax( ), nmax );

    // Expected number of coefficients for degree 10
    constexpr std::size_t expectedCoeffs = ( nmax + 1 ) * ( nmax + 2 ) / 2;
    BOOST_CHECK_EQUAL( dataset.nCoeffs( ), expectedCoeffs );

    // Test coefficient setting and getting
    dataset.setCoeff( 0, 0, 0, 2, 1, 0.5, -0.3 );
    auto coeff_pair = dataset.getCoeff( 0, 0, 0, 2, 1 );
    double C = coeff_pair.first;
    double S = coeff_pair.second;
    BOOST_CHECK_CLOSE( C, 0.5, 1e-12 );
    BOOST_CHECK_CLOSE( S, -0.3, 1e-12 );

    // Test block access
    const auto block = dataset.block( 0, 0, 0 );
    BOOST_CHECK_EQUAL( block.rows( ), expectedCoeffs );
    BOOST_CHECK_EQUAL( block.cols( ), 2 );

    // Test coefficient matrices
    dataset.setCoeff( 0, 1, 1, 3, 2, 1.5, 2.5 );
    auto matrices = dataset.getCoefficientMatrices( 0, 1, 1 );
    auto cosineMatrix = matrices.first;
    auto sineMatrix = matrices.second;
    BOOST_CHECK_EQUAL( cosineMatrix.rows( ), nmax + 1 );
    BOOST_CHECK_EQUAL( cosineMatrix.cols( ), nmax + 1 );
    BOOST_CHECK_CLOSE( cosineMatrix( 3, 2 ), 1.5, 1e-12 );
    BOOST_CHECK_CLOSE( sineMatrix( 3, 2 ), 2.5, 1e-12 );
}

/**
 * @brief Verifies that out-of-bounds access throws std::out_of_range exceptions.
 *
 * Input/Setup:
 * - Creates a minimal dataset with 1 file, 1 radius (+1 ref), 1 longitude, and nmax=5
 *
 * Expected Behavior:
 * - Accessing file index >= nFiles throws std::out_of_range
 * - Accessing radius index >= nRadii throws std::out_of_range
 * - Accessing longitude index >= nLongitudes throws std::out_of_range
 * - Accessing degree n > nmax throws std::out_of_range
 * - Accessing order m > n throws std::out_of_range
 */
BOOST_AUTO_TEST_CASE( test_stokes_dataset_bounds_checking )
{
    const std::vector< ComaStokesDataset::FileMeta > files = { { 0, 0, "test", 10000.0 } };
    const std::vector< double > radii = { 1000.0 };  // + ref. Radius
    const std::vector< double > lons = { 0.0 };
    constexpr int nmax = 5;

    ComaStokesDataset dataset = ComaStokesDataset::create( files, radii, lons, nmax );

    // Test out of bounds access
    BOOST_CHECK_THROW( dataset.setCoeff( 1, 0, 0, 0, 0, 0, 0 ), std::out_of_range );  // file OOR
    BOOST_CHECK_THROW( dataset.setCoeff( 0, 2, 0, 0, 0, 0, 0 ), std::out_of_range );  // radius OOR
    BOOST_CHECK_THROW( dataset.setCoeff( 0, 0, 1, 0, 0, 0, 0 ), std::out_of_range );  // longitude OOR
    BOOST_CHECK_THROW( dataset.setCoeff( 0, 0, 0, 6, 0, 0, 0 ), std::out_of_range );  // n > nmax
    BOOST_CHECK_THROW( dataset.setCoeff( 0, 0, 0, 5, 6, 0, 0 ), std::out_of_range );  // m > n
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== I/O Component Tests ====================

/**
 * @brief Test suite for file I/O reader and writer components.
 *
 * Tests the ComaPolyDatasetReader, ComaStokesDatasetReader, and ComaStokesDatasetWriter
 * classes. Verifies correct reading of polynomial coefficient files, writing of Stokes
 * datasets to CSV format, and round-trip consistency for CSV read/write operations.
 */
BOOST_AUTO_TEST_SUITE( test_io_components )

/**
 * @brief Verifies ComaPolyDatasetReader correctly reads polynomial coefficient files.
 *
 * Input/Setup:
 * - Reads the test polynomial coefficient file using ComaPolyDatasetReader::readFromFiles
 *
 * Expected Behavior:
 * - Dataset reports 1 file loaded
 * - Polynomial coefficient matrix has correct dimensions (48 rows x 121 columns)
 * - Specific coefficient values match expected reference values at positions (0,0), (3,1), (47,120), (10,22)
 * - SH degree/order indices have correct dimensions and values
 * - Reference radius and max degree are correctly extracted
 * - Powers of inverse radius are correctly read
 * - Column access by (n,m) returns correct values
 */
BOOST_FIXTURE_TEST_CASE( test_poly_dataset_reader, TestDataPaths )
{
    const std::vector< std::string > files = { testFile.string( ) };

    // Test reader functionality
    const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles( files );

    BOOST_CHECK_EQUAL( dataset.getNumFiles( ), 1 );

    // Check dimensions
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients( 0 );
    BOOST_CHECK_EQUAL( polyCoefs.rows( ), ExpectedPolyValues::numTerms );
    BOOST_CHECK_EQUAL( polyCoefs.cols( ), ExpectedPolyValues::numCoeffs );

    // Check specific coefficient values
    BOOST_CHECK_CLOSE( polyCoefs( 0, 0 ), ExpectedPolyValues::polyCoef_0_0, 1e-12 );
    BOOST_CHECK_CLOSE( polyCoefs( 3, 1 ), ExpectedPolyValues::polyCoef_3_1, 1e-12 );
    BOOST_CHECK_CLOSE( polyCoefs( 47, 120 ), ExpectedPolyValues::polyCoef_47_120, 1e-12 );
    BOOST_CHECK_CLOSE( polyCoefs( 10, 22 ), ExpectedPolyValues::polyCoef_10_22, 1e-12 );

    // Check SH degree and order indices
    const Eigen::ArrayXXi& indices = dataset.getSHDegreeAndOrderIndices( 0 );
    BOOST_CHECK_EQUAL( indices.rows( ), 2 );
    BOOST_CHECK_EQUAL( indices.cols( ), ExpectedPolyValues::numCoeffs );
    BOOST_CHECK_EQUAL( indices( 0, 120 ), 10 );   // degree
    BOOST_CHECK_EQUAL( indices( 1, 120 ), -10 );  // order

    // Check metadata
    BOOST_CHECK_EQUAL( dataset.getReferenceRadius( 0 ), ExpectedPolyValues::refRadius );
    BOOST_CHECK_EQUAL( dataset.getMaxDegreeSH( 0 ), ExpectedPolyValues::maxDegree );

    // Check powers
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius( 0 );
    BOOST_CHECK_EQUAL( powers.rows( ), 4 );
    BOOST_CHECK_EQUAL( powers[ 0 ], 0 );
    BOOST_CHECK_EQUAL( powers[ 3 ], 3 );

    // Test column access by (n,m)
    const Eigen::VectorXd col = dataset.columnForNM( 0, 3, 1 );
    BOOST_CHECK_EQUAL( col.rows( ), ExpectedPolyValues::numTerms );

    // Test value access
    const double val = dataset.value( 0, 10, 3, 1 );
    BOOST_CHECK_CLOSE( val, col[ 10 ], 1e-12 );
}

/**
 * @brief Verifies ComaStokesDatasetWriter correctly writes Stokes dataset to CSV format.
 *
 * Input/Setup:
 * - Creates a small Stokes dataset with 1 file, 2 radii, 2 longitudes, and nmax=10
 * - Sets specific coefficient values at various indices
 *
 * Expected Behavior:
 * - CSV file is created at the specified output path
 * - Meta line contains correct metadata (epochs, max_degree, n_radii, n_lons)
 * - Radii line contains the correct radius values in meters
 * - Longitudes line contains the correct longitude values in degrees
 */
BOOST_FIXTURE_TEST_CASE( test_stokes_dataset_writer, TestDataPaths )
{
    // Create a small dataset for testing
    std::vector< ComaStokesDataset::FileMeta > files = { { 0.0, 1.0, "test_source", 10000 } };
    std::vector< double > radii = { 6000.0, 10000.0 };
    std::vector< double > lons = { 0.0, 30.0 };
    int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create( files, radii, lons, nmax );

    // Set some test values
    dataset.setCoeff( 0, 0, 0, 0, 0, 1.0, 0.0 );
    dataset.setCoeff( 0, 0, 1, 2, 1, 0.5, -0.5 );
    dataset.setCoeff( 0, 1, 0, 3, 3, 0.25, 0.75 );

    // Write to CSV
    boost::filesystem::path csvPath = outputDir / "test_stokes.csv";
    ComaStokesDatasetWriter::writeCsvForFile( dataset, 0, csvPath.string( ) );

    BOOST_CHECK( boost::filesystem::exists( csvPath ) );

    // Verify CSV content
    std::ifstream ifs( csvPath.string( ) );
    BOOST_REQUIRE( ifs.is_open( ) );

    std::string line;

    // Check meta line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "meta" ) != std::string::npos );
    BOOST_CHECK( line.find( "start_epoch=0" ) != std::string::npos );
    BOOST_CHECK( line.find( "end_epoch=1" ) != std::string::npos );
    BOOST_CHECK( line.find( "max_degree=10" ) != std::string::npos );
    BOOST_CHECK( line.find( "n_radii=2" ) != std::string::npos );
    BOOST_CHECK( line.find( "n_lons=2" ) != std::string::npos );

    // Check radii line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "radii [meter]" ) != std::string::npos );
    BOOST_CHECK( line.find( "6000" ) != std::string::npos );
    BOOST_CHECK( line.find( "10000" ) != std::string::npos );

    // Check longitudes line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "longitudes [degree]" ) != std::string::npos );
    BOOST_CHECK( line.find( "0" ) != std::string::npos );
    BOOST_CHECK( line.find( "30" ) != std::string::npos );

    ifs.close( );
}

/**
 * @brief Verifies round-trip: write a Stokes dataset to CSV, read it back, verify values match.
 *
 * Input/Setup:
 * - Creates a Stokes dataset with 1 file, 2 radii, 2 longitudes, and nmax=8
 * - Sets specific coefficient values at various (file, radius, longitude, n, m) indices
 * - Writes the dataset to CSV using ComaStokesDatasetWriter
 *
 * Expected Behavior:
 * - ComaStokesDatasetReader::readFromCsv correctly reads the written file
 * - Read dataset has same structure (nFiles, nRadii, nLongitudes, nmax)
 * - File metadata (epochs, source_tag) is preserved
 * - Radii and longitude values are preserved with high precision
 * - All coefficient values match the original values within tolerance
 */
BOOST_FIXTURE_TEST_CASE( test_stokes_dataset_reader_from_csv, TestDataPaths )
{
    // First, create and write a test dataset
    std::vector< ComaStokesDataset::FileMeta > files = { { 2.015e9, 2.0150864e9, "test_source", 10000.0 } };
    std::vector< double > radii = { 6000.0, 10000.0 };
    std::vector< double > lons = { 0.0, 30.0 };
    int nmax = 8;

    ComaStokesDataset originalDataset = ComaStokesDataset::create( files, radii, lons, nmax );

    // Set some known test values
    originalDataset.setCoeff( 0, 0, 0, 0, 0, 54.0, 0.0 );
    originalDataset.setCoeff( 0, 0, 0, 2, 1, -0.232, 0.138 );
    originalDataset.setCoeff( 0, 0, 1, 0, 0, 53.93, 0.0 );
    originalDataset.setCoeff( 0, 1, 0, 1, 1, -1.75, 0.407 );
    originalDataset.setCoeff( 0, 1, 1, 3, 2, -0.084, -0.026 );

    // Write to CSV
    boost::filesystem::path csvPath = outputDir / "test_reader.csv";
    ComaStokesDatasetWriter::writeCsvForFile( originalDataset, 0, csvPath.string( ) );

    // Now test reading it back
    ComaStokesDataset readDataset = ComaStokesDatasetReader::readFromCsv( csvPath.string( ) );

    // Verify structure
    BOOST_CHECK_EQUAL( readDataset.nFiles( ), 1 );
    BOOST_CHECK_EQUAL( readDataset.nRadii( ), 2 );
    BOOST_CHECK_EQUAL( readDataset.nLongitudes( ), 2 );
    BOOST_CHECK_EQUAL( readDataset.nmax( ), nmax );

    // Verify metadata
    const auto& filesMeta = readDataset.files( );
    BOOST_CHECK_CLOSE( filesMeta[ 0 ].startEpoch, 2.015e9, 1e-6 );
    BOOST_CHECK_CLOSE( filesMeta[ 0 ].endEpoch, 2.0150864e9, 1e-6 );
    BOOST_CHECK_EQUAL( filesMeta[ 0 ].sourceTag, "test_source" );

    // Verify radii and longitudes
    const auto& readRadii = readDataset.radii( );
    const auto& readLons = readDataset.lons( );
    BOOST_CHECK_CLOSE( readRadii[ 0 ], 6000.0, 1e-10 );
    BOOST_CHECK_CLOSE( readRadii[ 1 ], 10000.0, 1e-10 );
    BOOST_CHECK_CLOSE( readLons[ 0 ], 0.0, 1e-10 );
    BOOST_CHECK_CLOSE( readLons[ 1 ], 30.0, 1e-10 );

    // Verify coefficient values
    auto coeff_0_0_r0_l0 = readDataset.getCoeff( 0, 0, 0, 0, 0 );
    double C_0_0_r0_l0 = coeff_0_0_r0_l0.first;
    double S_0_0_r0_l0 = coeff_0_0_r0_l0.second;
    BOOST_CHECK_CLOSE( C_0_0_r0_l0, 54.0, 1e-10 );
    BOOST_CHECK_CLOSE( S_0_0_r0_l0, 0.0, 1e-10 );

    auto coeff_2_1_r0_l0 = readDataset.getCoeff( 0, 0, 0, 2, 1 );
    double C_2_1_r0_l0 = coeff_2_1_r0_l0.first;
    double S_2_1_r0_l0 = coeff_2_1_r0_l0.second;
    BOOST_CHECK_CLOSE( C_2_1_r0_l0, -0.232, 1e-10 );
    BOOST_CHECK_CLOSE( S_2_1_r0_l0, 0.138, 1e-10 );

    auto coeff_0_0_r0_l1 = readDataset.getCoeff( 0, 0, 1, 0, 0 );
    double C_0_0_r0_l1 = coeff_0_0_r0_l1.first;
    double S_0_0_r0_l1 = coeff_0_0_r0_l1.second;
    BOOST_CHECK_CLOSE( C_0_0_r0_l1, 53.93, 1e-10 );

    auto coeff_1_1_r1_l0 = readDataset.getCoeff( 0, 1, 0, 1, 1 );
    double C_1_1_r1_l0 = coeff_1_1_r1_l0.first;
    double S_1_1_r1_l0 = coeff_1_1_r1_l0.second;
    BOOST_CHECK_CLOSE( C_1_1_r1_l0, -1.75, 1e-10 );
    BOOST_CHECK_CLOSE( S_1_1_r1_l0, 0.407, 1e-10 );

    auto coeff_3_2_r1_l1 = readDataset.getCoeff( 0, 1, 1, 3, 2 );
    double C_3_2_r1_l1 = coeff_3_2_r1_l1.first;
    double S_3_2_r1_l1 = coeff_3_2_r1_l1.second;
    BOOST_CHECK_CLOSE( C_3_2_r1_l1, -0.084, 1e-10 );
    BOOST_CHECK_CLOSE( S_3_2_r1_l1, -0.026, 1e-10 );
}

/**
 * @brief Verifies reading multiple Stokes CSV files from a folder and combining into a single dataset.
 *
 * Input/Setup:
 * - Creates two separate Stokes datasets representing different time epochs
 * - Each dataset has different coefficient values to verify correct file association
 * - Writes both datasets to separate CSV files in a common folder
 *
 * Expected Behavior:
 * - ComaStokesDatasetReader::readFromCsvFolder correctly reads all matching files
 * - Combined dataset has nFiles=2 with shared structure (nRadii, nLongitudes, nmax)
 * - File metadata for each file is correctly preserved (epochs, source_tag)
 * - Coefficient values from each file are correctly associated with file index
 */
BOOST_FIXTURE_TEST_CASE( test_stokes_dataset_reader_from_csv_folder, TestDataPaths )
{
    // Create multiple test datasets (simulating multiple time epochs)
    std::vector< ComaStokesDataset::FileMeta > files1 = { { 0.0, 1.0, "test_file_1", 10000.0 } };
    std::vector< ComaStokesDataset::FileMeta > files2 = { { 1.0, 2.0, "test_file_2", 10000.0 } };

    std::vector< double > radii = { 6000.0, 10000.0 };
    std::vector< double > lons = { 0.0, 45.0 };
    int nmax = 5;

    ComaStokesDataset dataset1 = ComaStokesDataset::create( files1, radii, lons, nmax );
    ComaStokesDataset dataset2 = ComaStokesDataset::create( files2, radii, lons, nmax );

    // Set different coefficients for each dataset
    dataset1.setCoeff( 0, 0, 0, 0, 0, 10.0, 0.0 );
    dataset1.setCoeff( 0, 0, 1, 1, 1, 1.5, -0.5 );
    dataset1.setCoeff( 0, 1, 0, 2, 0, 0.25, 0.0 );

    dataset2.setCoeff( 0, 0, 0, 0, 0, 20.0, 0.0 );
    dataset2.setCoeff( 0, 0, 1, 1, 1, 2.5, -1.0 );
    dataset2.setCoeff( 0, 1, 0, 2, 0, 0.5, 0.0 );

    // Create folder for multi-file test
    boost::filesystem::path folderPath = outputDir / "multi_file_test";
    boost::filesystem::create_directories( folderPath );

    // Write both datasets to separate CSV files
    ComaStokesDatasetWriter::writeCsvForFile( dataset1, 0, ( folderPath / "test_file0.csv" ).string( ) );
    ComaStokesDatasetWriter::writeCsvForFile( dataset2, 0, ( folderPath / "test_file1.csv" ).string( ) );

    // Test reading from folder
    ComaStokesDataset multiDataset = ComaStokesDatasetReader::readFromCsvFolder( folderPath.string( ), "test" );

    // Verify structure
    BOOST_CHECK_EQUAL( multiDataset.nFiles( ), 2 );
    BOOST_CHECK_EQUAL( multiDataset.nRadii( ), 2 );
    BOOST_CHECK_EQUAL( multiDataset.nLongitudes( ), 2 );
    BOOST_CHECK_EQUAL( multiDataset.nmax( ), nmax );

    // Verify file metadata
    const auto& filesMeta = multiDataset.files( );
    BOOST_CHECK_CLOSE( filesMeta[ 0 ].startEpoch, 0.0, 1e-10 );
    BOOST_CHECK_CLOSE( filesMeta[ 0 ].endEpoch, 1.0, 1e-10 );
    BOOST_CHECK_EQUAL( filesMeta[ 0 ].sourceTag, "test_file_1" );
    BOOST_CHECK_CLOSE( filesMeta[ 1 ].startEpoch, 1.0, 1e-10 );
    BOOST_CHECK_CLOSE( filesMeta[ 1 ].endEpoch, 2.0, 1e-10 );
    BOOST_CHECK_EQUAL( filesMeta[ 1 ].sourceTag, "test_file_2" );

    // Verify coefficient values from first file
    auto coeff_1_0_0 = multiDataset.getCoeff( 0, 0, 0, 0, 0 );
    double C1_0_0 = coeff_1_0_0.first;
    double S1_0_0 = coeff_1_0_0.second;
    BOOST_CHECK_CLOSE( C1_0_0, 10.0, 1e-10 );

    auto coeff_1_1_1 = multiDataset.getCoeff( 0, 0, 1, 1, 1 );
    double C1_1_1 = coeff_1_1_1.first;
    double S1_1_1 = coeff_1_1_1.second;
    BOOST_CHECK_CLOSE( C1_1_1, 1.5, 1e-10 );
    BOOST_CHECK_CLOSE( S1_1_1, -0.5, 1e-10 );

    // Verify coefficient values from second file
    auto coeff_2_0_0 = multiDataset.getCoeff( 1, 0, 0, 0, 0 );
    double C2_0_0 = coeff_2_0_0.first;
    double S2_0_0 = coeff_2_0_0.second;
    BOOST_CHECK_CLOSE( C2_0_0, 20.0, 1e-10 );

    auto coeff_2_1_1 = multiDataset.getCoeff( 1, 0, 1, 1, 1 );
    double C2_1_1 = coeff_2_1_1.first;
    double S2_1_1 = coeff_2_1_1.second;
    BOOST_CHECK_CLOSE( C2_1_1, 2.5, 1e-10 );
    BOOST_CHECK_CLOSE( S2_1_1, -1.0, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== Transformation/Processing Tests ====================

/**
 * @brief Test suite for coefficient evaluation and dataset transformation.
 *
 * Tests the StokesCoefficientsEvaluator and ComaDatasetTransformer classes.
 * Verifies correct computation of Stokes coefficients from polynomial data,
 * including truncation handling and error conditions for invalid parameters.
 */
BOOST_AUTO_TEST_SUITE( test_processing_components )

/**
 * @brief Verifies StokesCoefficientsEvaluator::evaluate2D computes correct Stokes coefficients.
 *
 * Input/Setup:
 * - Reads polynomial coefficient data from test files (legacy and rotated datasets)
 * - Loads reference Stokes coefficient data from pre-computed test files
 *
 * Expected Behavior:
 * - Test Case 1 (0° solar longitude, 4km radius): All computed cosine and sine coefficients
 *   match reference values within 1e-10 relative tolerance
 * - Test Case 2 (30° solar longitude, 10km radius): All computed coefficients match reference
 * - Truncated evaluation (maxDegree=5, maxOrder=3): Returns correctly sized matrices with
 *   matching coefficient values for the truncated range
 * - Invalid degree/order requests: Throws std::runtime_error when requested maxima exceed
 *   available data
 */
BOOST_FIXTURE_TEST_CASE( test_stokes_coefficients_evaluator, TestDataPaths )
{
    const auto configs = getPolyDatasetConfigs( );

    for( const auto& config : configs )
    {
        BOOST_TEST_CONTEXT( "dataset=" << config.label )
        {
            const std::vector< std::string > files = { config.polyFile.string( ) };
            const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles( files );

            const int maxDegree = config.maxDegree;
            const int maxOrder = config.maxOrder;

            const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients( 0 );
            const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices( 0 );
            const Eigen::VectorXd& powers = dataset.getPowersInvRadius( 0 );
            const double refRadius = dataset.getReferenceRadius( 0 );

            const Eigen::ArrayXXd polyCoefficients = polyCoefs.array( );

            StokesTestData expectedData1 = StokesTestData::readFromFile( config.stokesFileLon0.string( ), maxDegree );
            StokesTestData expectedData2 = StokesTestData::readFromFile( config.stokesFileLon30.string( ), maxDegree );

            // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
            {
                const double radius = expectedData1.radius;
                const double solarLongitude = expectedData1.solarLongitude * mathematical_constants::PI / 180.0;

                Eigen::MatrixXd cosineCoefficients, sineCoefficients;

                simulation_setup::StokesCoefficientsEvaluator::evaluate2D( radius,
                                                                           solarLongitude,
                                                                           polyCoefficients,
                                                                           shIndices,
                                                                           powers,
                                                                           refRadius,
                                                                           cosineCoefficients,
                                                                           sineCoefficients,
                                                                           maxDegree,
                                                                           maxOrder );

                BOOST_CHECK_EQUAL( cosineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( cosineCoefficients.cols( ), maxOrder + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.cols( ), maxOrder + 1 );

                int numCoeffsChecked = 0;
                for( int n = 0; n <= maxDegree; ++n )
                {
                    for( int m = 0; m <= n; ++m )
                    {
                        BOOST_CHECK_MESSAGE( std::abs( cosineCoefficients( n, m ) - expectedData1.cosineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData1.cosineCoeffs( n, m ) ), 1e-10 ) <
                                                     1e-10,
                                             "Cosine coefficient C(" << n << "," << m
                                                                     << ") mismatch: computed = " << cosineCoefficients( n, m )
                                                                     << ", expected = " << expectedData1.cosineCoeffs( n, m ) );

                        BOOST_CHECK_MESSAGE( std::abs( sineCoefficients( n, m ) - expectedData1.sineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData1.sineCoeffs( n, m ) ), 1e-10 ) <
                                                     1e-10,
                                             "Sine coefficient S(" << n << "," << m << ") mismatch: computed = " << sineCoefficients( n, m )
                                                                   << ", expected = " << expectedData1.sineCoeffs( n, m ) );

                        numCoeffsChecked++;
                    }
                }
                const int expectedNumCoeffs = ( maxDegree + 1 ) * ( maxDegree + 2 ) / 2;
                BOOST_CHECK_EQUAL( numCoeffsChecked, expectedNumCoeffs );
                BOOST_TEST_MESSAGE( "Test Case 1: Verified all " << numCoeffsChecked << " cosine and " << numCoeffsChecked
                                                                 << " sine coefficients (total: " << 2 * numCoeffsChecked << ")" );
            }

            // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
            {
                const double radius = expectedData2.radius;
                const double solarLongitude = expectedData2.solarLongitude * mathematical_constants::PI / 180.0;

                Eigen::MatrixXd cosineCoefficients, sineCoefficients;

                simulation_setup::StokesCoefficientsEvaluator::evaluate2D( radius,
                                                                           solarLongitude,
                                                                           polyCoefficients,
                                                                           shIndices,
                                                                           powers,
                                                                           refRadius,
                                                                           cosineCoefficients,
                                                                           sineCoefficients,
                                                                           maxDegree,
                                                                           maxOrder );

                BOOST_CHECK_EQUAL( cosineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( cosineCoefficients.cols( ), maxOrder + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.cols( ), maxOrder + 1 );

                int numCoeffsChecked = 0;
                for( int n = 0; n <= maxDegree; ++n )
                {
                    for( int m = 0; m <= n; ++m )
                    {
                        BOOST_CHECK_MESSAGE( std::abs( cosineCoefficients( n, m ) - expectedData2.cosineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData2.cosineCoeffs( n, m ) ), 1e-10 ) <
                                                     1e-10,
                                             "Cosine coefficient C(" << n << "," << m
                                                                     << ") mismatch: computed = " << cosineCoefficients( n, m )
                                                                     << ", expected = " << expectedData2.cosineCoeffs( n, m ) );

                        BOOST_CHECK_MESSAGE( std::abs( sineCoefficients( n, m ) - expectedData2.sineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData2.sineCoeffs( n, m ) ), 1e-10 ) <
                                                     9.9e-9,
                                             "Sine coefficient S(" << n << "," << m << ") mismatch: computed = " << sineCoefficients( n, m )
                                                                   << ", expected = " << expectedData2.sineCoeffs( n, m ) );

                        numCoeffsChecked++;
                    }
                }
                const int expectedNumCoeffs = ( maxDegree + 1 ) * ( maxDegree + 2 ) / 2;
                BOOST_CHECK_EQUAL( numCoeffsChecked, expectedNumCoeffs );
                BOOST_TEST_MESSAGE( "Test Case 2: Verified all " << numCoeffsChecked << " cosine and " << numCoeffsChecked
                                                                 << " sine coefficients (total: " << 2 * numCoeffsChecked << ")" );
            }

            // ========== Test truncated degree/order ==========
            {
                Eigen::MatrixXd truncatedCosine, truncatedSine;
                const int truncatedMaxDegree = 5;
                const int truncatedMaxOrder = 3;
                const double radius = expectedData1.radius;
                const double solarLongitude = expectedData1.solarLongitude * mathematical_constants::PI / 180.0;

                simulation_setup::StokesCoefficientsEvaluator::evaluate2D( radius,
                                                                           solarLongitude,
                                                                           polyCoefficients,
                                                                           shIndices,
                                                                           powers,
                                                                           refRadius,
                                                                           truncatedCosine,
                                                                           truncatedSine,
                                                                           truncatedMaxDegree,
                                                                           truncatedMaxOrder );

                BOOST_CHECK_EQUAL( truncatedCosine.rows( ), truncatedMaxDegree + 1 );
                BOOST_CHECK_EQUAL( truncatedCosine.cols( ), truncatedMaxOrder + 1 );
                BOOST_CHECK_EQUAL( truncatedSine.rows( ), truncatedMaxDegree + 1 );
                BOOST_CHECK_EQUAL( truncatedSine.cols( ), truncatedMaxOrder + 1 );

                for( int n = 0; n <= truncatedMaxDegree; ++n )
                {
                    for( int m = 0; m <= std::min( n, truncatedMaxOrder ); ++m )
                    {
                        BOOST_CHECK_CLOSE( truncatedCosine( n, m ), expectedData1.cosineCoeffs( n, m ), 1e-10 );
                        BOOST_CHECK_CLOSE( truncatedSine( n, m ), expectedData1.sineCoeffs( n, m ), 1e-10 );
                    }
                }
            }

            // ========== Test exceeding available maxima ==========
            {
                Eigen::MatrixXd dummyCosine, dummySine;
                const double radius = expectedData1.radius;
                const double solarLongitude = expectedData1.solarLongitude * mathematical_constants::PI / 180.0;
                const int invalidMaxDegree = config.maxDegree + 5;
                const int invalidMaxOrder = config.maxOrder + 5;

                BOOST_CHECK_THROW( simulation_setup::StokesCoefficientsEvaluator::evaluate2D( radius,
                                                                                              solarLongitude,
                                                                                              polyCoefficients,
                                                                                              shIndices,
                                                                                              powers,
                                                                                              refRadius,
                                                                                              dummyCosine,
                                                                                              dummySine,
                                                                                              invalidMaxDegree,
                                                                                              config.maxOrder ),
                                   std::runtime_error );

                BOOST_CHECK_THROW( simulation_setup::StokesCoefficientsEvaluator::evaluate2D( radius,
                                                                                              solarLongitude,
                                                                                              polyCoefficients,
                                                                                              shIndices,
                                                                                              powers,
                                                                                              refRadius,
                                                                                              dummyCosine,
                                                                                              dummySine,
                                                                                              config.maxDegree,
                                                                                              invalidMaxOrder ),
                                   std::runtime_error );
            }
        }
    }
}

/**
 * @brief Verifies ComaDatasetTransformer correctly transforms polynomial dataset to Stokes dataset.
 *
 * Input/Setup:
 * - Reads polynomial coefficient data from test files (legacy and rotated datasets)
 * - Uses radii = {4000, 10000} meters and longitudes = {0, 30} degrees
 * - Loads reference Stokes data for comparison
 *
 * Expected Behavior:
 * - Transformed dataset has correct structure (nFiles, nRadii, nLongitudes, nmax)
 * - Coefficient matrices at each (radius, longitude) point match reference values
 * - Truncated transformation (maxDegree=6, maxOrder=4) produces correctly sized output
 * - Invalid maxDegree/maxOrder throws std::invalid_argument
 */
BOOST_FIXTURE_TEST_CASE( test_dataset_transformer, TestDataPaths )
{
    const auto configs = getPolyDatasetConfigs( );
    const std::vector< double > radii_m = { 4000.0, 10000.0 };
    const std::vector< double > longitudes_deg = { 0.0, 30.0 };

    for( const auto& config : configs )
    {
        BOOST_TEST_CONTEXT( "dataset=" << config.label )
        {
            std::vector< std::string > files = { config.polyFile.string( ) };
            ComaPolyDataset polyDataset = ComaPolyDatasetReader::readFromFiles( files );

            StokesTestData expectedData1 = StokesTestData::readFromFile( config.stokesFileLon0.string( ), config.maxDegree );
            StokesTestData expectedData2 = StokesTestData::readFromFile( config.stokesFileLon30.string( ), config.maxDegree );

            ComaStokesDataset stokesDataset = ComaDatasetTransformer::transformPolyToStokes( polyDataset, radii_m, longitudes_deg );

            BOOST_CHECK_EQUAL( stokesDataset.nFiles( ), 1 );
            BOOST_CHECK_EQUAL( stokesDataset.nRadii( ), radii_m.size( ) );
            BOOST_CHECK_EQUAL( stokesDataset.nLongitudes( ), longitudes_deg.size( ) );
            BOOST_CHECK_EQUAL( stokesDataset.nmax( ), config.maxDegree );

            auto verifyCoefficients = [ & ]( int radiusIndex, int longitudeIndex, const StokesTestData& expectedData ) {
                auto matrices = stokesDataset.getCoefficientMatrices( 0, radiusIndex, longitudeIndex );
                const auto cosineMatrix = matrices.first;
                const auto sineMatrix = matrices.second;

                BOOST_CHECK_EQUAL( cosineMatrix.rows( ), config.maxDegree + 1 );
                BOOST_CHECK_EQUAL( cosineMatrix.cols( ), config.maxOrder + 1 );
                BOOST_CHECK_EQUAL( sineMatrix.rows( ), config.maxDegree + 1 );
                BOOST_CHECK_EQUAL( sineMatrix.cols( ), config.maxOrder + 1 );

                for( int n = 0; n <= config.maxDegree; ++n )
                {
                    for( int m = 0; m <= n; ++m )
                    {
                        BOOST_CHECK_MESSAGE( std::abs( cosineMatrix( n, m ) - expectedData.cosineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData.cosineCoeffs( n, m ) ), 1e-10 ) <
                                                     1e-10,
                                             "Cosine coefficient mismatch at (" << n << "," << m << ")" );
                        BOOST_CHECK_MESSAGE( std::abs( sineMatrix( n, m ) - expectedData.sineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData.sineCoeffs( n, m ) ), 1e-10 ) <
                                                     9.9e-9,
                                             "Sine coefficient mismatch at (" << n << "," << m << ")" );
                    }
                }
            };

            verifyCoefficients( 0, 0, expectedData1 );
            verifyCoefficients( 1, 1, expectedData2 );

            ComaStokesDataset truncatedDataset =
                    ComaDatasetTransformer::transformPolyToStokes( polyDataset, radii_m, longitudes_deg, 6, 4 );

            BOOST_CHECK_EQUAL( truncatedDataset.nmax( ), 6 );

            BOOST_CHECK_THROW( ComaDatasetTransformer::transformPolyToStokes(
                                       polyDataset, radii_m, longitudes_deg, config.maxDegree + 5, config.maxOrder ),
                               std::invalid_argument );

            BOOST_CHECK_THROW( ComaDatasetTransformer::transformPolyToStokes(
                                       polyDataset, radii_m, longitudes_deg, config.maxDegree, config.maxOrder + 5 ),
                               std::invalid_argument );
        }
    }
}

/**
 * @brief Verifies ComaModelFileProcessor::createSHDataset produces correct Stokes coefficients.
 *
 * Input/Setup:
 * - Creates ComaModelFileProcessor from test polynomial coefficient files
 * - Uses radii = {4000, 10000} meters and longitudes = {0, 30} degrees
 * - Loads reference Stokes coefficient data from pre-computed test files
 *
 * Expected Behavior:
 * - Created dataset has correct structure (nFiles=1, nRadii=2, nLongitudes=2, nmax as specified)
 * - Coefficient matrices at (ri=0, li=0) match 0°/4km reference values
 * - Coefficient matrices at (ri=1, li=1) match 30°/10km reference values
 * - All cosine and sine coefficients verified against reference data
 */
BOOST_FIXTURE_TEST_CASE( test_stokes_dataset_creation_via_processor, TestDataPaths )
{
    const auto configs = getPolyDatasetConfigs( );
    const std::vector< double > radii_m = { 4000.0, 10000.0 };
    const std::vector< double > longitudes_deg = { 0.0, 30.0 };

    for( const auto& config : configs )
    {
        BOOST_TEST_CONTEXT( "dataset=" << config.label )
        {
            const std::vector< std::string > files = { config.polyFile.string( ) };
            const ComaModelFileProcessor processor( files );

            const int maxDegree = config.maxDegree;
            const int maxOrder = config.maxOrder;

            const ComaStokesDataset stokesDataset = processor.createSHDataset( radii_m, longitudes_deg, maxDegree, maxOrder );

            BOOST_CHECK_EQUAL( stokesDataset.nFiles( ), 1 );
            BOOST_CHECK_EQUAL( stokesDataset.nRadii( ), radii_m.size( ) );
            BOOST_CHECK_EQUAL( stokesDataset.nLongitudes( ), longitudes_deg.size( ) );
            BOOST_CHECK_EQUAL( stokesDataset.nmax( ), maxDegree );

            StokesTestData expectedData1 = StokesTestData::readFromFile( config.stokesFileLon0.string( ), maxDegree );
            StokesTestData expectedData2 = StokesTestData::readFromFile( config.stokesFileLon30.string( ), maxDegree );

            auto verifyCoefficients = [ & ]( int radiusIndex, int longitudeIndex, const StokesTestData& expectedData ) {
                auto matrices = stokesDataset.getCoefficientMatrices( 0, radiusIndex, longitudeIndex );
                const auto cosineCoefficients = matrices.first;
                const auto sineCoefficients = matrices.second;

                BOOST_CHECK_EQUAL( cosineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( cosineCoefficients.cols( ), maxOrder + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.rows( ), maxDegree + 1 );
                BOOST_CHECK_EQUAL( sineCoefficients.cols( ), maxOrder + 1 );

                int numCoeffsChecked = 0;
                for( int n = 0; n <= maxDegree; ++n )
                {
                    for( int m = 0; m <= n; ++m )
                    {
                        BOOST_CHECK_MESSAGE( std::abs( cosineCoefficients( n, m ) - expectedData.cosineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData.cosineCoeffs( n, m ) ), 1e-10 ) <
                                                     1e-10,
                                             "Cosine coefficient C(" << n << "," << m
                                                                     << ") mismatch: computed = " << cosineCoefficients( n, m )
                                                                     << ", expected = " << expectedData.cosineCoeffs( n, m ) );

                        BOOST_CHECK_MESSAGE( std::abs( sineCoefficients( n, m ) - expectedData.sineCoeffs( n, m ) ) /
                                                             std::max( std::abs( expectedData.sineCoeffs( n, m ) ), 1e-10 ) <
                                                     9.9e-10,
                                             "Sine coefficient S(" << n << "," << m << ") mismatch: computed = " << sineCoefficients( n, m )
                                                                   << ", expected = " << expectedData.sineCoeffs( n, m ) );

                        numCoeffsChecked++;
                    }
                }

                const int expectedNumCoeffs = ( maxDegree + 1 ) * ( maxDegree + 2 ) / 2;
                BOOST_CHECK_EQUAL( numCoeffsChecked, expectedNumCoeffs );
                BOOST_TEST_MESSAGE( "Verified " << numCoeffsChecked << " cosine and sine coefficients for indices (" << radiusIndex << ","
                                                << longitudeIndex << ")" );
            };

            verifyCoefficients( 0, 0, expectedData1 );
            verifyCoefficients( 1, 1, expectedData2 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== ComaModel Tests ====================

/**
 * @brief Test suite for ComaModel density calculation functionality.
 *
 * Tests the ComaModel class which computes gas number density and mass density
 * in a cometary coma using spherical harmonic expansions. Validates against
 * reference residual files and Python-computed reference values.
 */
BOOST_AUTO_TEST_SUITE( test_coma_model )

BOOST_AUTO_TEST_CASE( test_solar_longitude_dependent_variable_update_requirements )
{
    const std::string cometName = "Churyumov_Gerasimenko";
    const std::string sunName = "Sun";

    SystemOfBodies bodies( "SSB", "ECLIPJ2000" );
    bodies.createEmptyBody( cometName );
    bodies.createEmptyBody( sunName );

    bodies.at( cometName )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) ) );
    bodies.at( cometName )
            ->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
                    Eigen::Quaterniond::Identity( ), "ECLIPJ2000", "Churyumov_Gerasimenko_Fixed" ) );

    Eigen::Vector6d sunState = Eigen::Vector6d::Zero( );
    sunState( 0 ) = 1.0E11;
    bodies.at( sunName )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( sunState ) );

    const std::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSettings =
            propagators::solarLongitudeDependentVariable( cometName );

    const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > updateSettings =
            propagators::createEnvironmentUpdaterSettingsForDependentVariables( dependentVariableSettings, bodies );

    BOOST_CHECK_NO_THROW( propagators::checkValidityOfRequiredEnvironmentUpdates( updateSettings, bodies ) );
    BOOST_CHECK_EQUAL( updateSettings.count( propagators::vehicle_flight_conditions_update ), 0 );

    const auto containsBody = []( const std::vector< std::string >& bodyList, const std::string& bodyName ) {
        return std::find( bodyList.begin( ), bodyList.end( ), bodyName ) != bodyList.end( );
    };

    BOOST_REQUIRE_EQUAL( updateSettings.count( propagators::body_translational_state_update ), 1 );
    const std::vector< std::string >& translationalUpdates = updateSettings.at( propagators::body_translational_state_update );
    BOOST_CHECK_EQUAL( translationalUpdates.size( ), 2 );
    BOOST_CHECK( containsBody( translationalUpdates, cometName ) );
    BOOST_CHECK( containsBody( translationalUpdates, sunName ) );
    BOOST_CHECK( !containsBody( translationalUpdates, "" ) );

    BOOST_REQUIRE_EQUAL( updateSettings.count( propagators::body_rotational_state_update ), 1 );
    const std::vector< std::string >& rotationalUpdates = updateSettings.at( propagators::body_rotational_state_update );
    BOOST_CHECK_EQUAL( rotationalUpdates.size( ), 1 );
    BOOST_CHECK( containsBody( rotationalUpdates, cometName ) );
    BOOST_CHECK( !containsBody( rotationalUpdates, sunName ) );
    BOOST_CHECK( !containsBody( rotationalUpdates, "" ) );
}

/**
 * @brief Verifies ComaModel::getTotalNumberDensity returns correct log2(density) values.
 *
 * Input/Setup:
 * - Creates ComaModel from polynomial coefficients and pre-computed Stokes dataset
 * - Tests with both polynomial and Stokes dataset input types
 * - Uses reference residual files containing expected log2(number_density) values
 *
 * Expected Behavior:
 * - Test Case 1 (0° solar longitude, 4km radius): 100 randomly selected points
 *   match reference log2(density) values within 1e-10 tolerance
 * - Test Case 2 (30° solar longitude, 10km radius): 50 randomly selected points
 *   match reference values
 * - Both polynomial and Stokes dataset inputs produce identical results
 */
BOOST_FIXTURE_TEST_CASE( test_coma_model_number_density, TestDataPaths )
{
    // Load polynomial coefficients from test data file
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaPolyDataset polyDataset = ComaPolyDatasetReader::readFromFiles( files );

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Molecular weight (arbitrary for number density test, but needed for constructor)
    const double molecularWeight = 0.018;  // kg/mol (water vapor)

    // Create Stokes dataset from polynomial dataset
    // Use the test points from both test files as the grid
    const std::vector< double > radii_m = { 4000.0, 10000.0 };   // 4 km, 10 km
    const std::vector< double > longitudes_deg = { 0.0, 30.0 };  // 0°, 30°
    const ComaModelFileProcessor processor( files );
    const ComaStokesDataset stokesDataset = processor.createSHDataset( radii_m, longitudes_deg, maxDegree, maxOrder );

    // Test with both polynomial and Stokes datasets
    for( int datasetType = 0; datasetType < 2; ++datasetType )
    {
        const std::string datasetTypeName = ( datasetType == 0 ) ? "Polynomial" : "Stokes";
        BOOST_TEST_MESSAGE( "========== Testing with " << datasetTypeName << " dataset ==========" );

        // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
        {
            const boost::filesystem::path residualsFile1 = dataDir / "density" / "residual" / "residual_r_cometFixed_ep10-000_04km.txt";
            const double testRadius = 4000.0;       // 4 km in meters
            const double testSolarLongitude = 0.0;  // 0 degrees in radians
            const double testTime = 490708800;      // s since J2000

            // Define state functions for solar longitude = 0°
            // Solar longitude is calculated as atan2(y, x) of the Sun direction in body-fixed frame
            // For solar longitude = 0°, we need Sun along the +X axis in body-fixed frame

            // Sun state function - returns position of Sun along +X axis
            auto sunStateFunction = []( ) -> Eigen::Vector6d {
                Eigen::Vector6d state = Eigen::Vector6d::Zero( );
                state.segment( 0, 3 ) = Eigen::Vector3d( 1.0e11, 0.0, 0.0 );  // Sun along +X axis
                return state;
            };

            // Comet state function - returns position of comet at origin
            auto cometStateFunction = []( ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); };

            // Comet rotation function - returns identity matrix
            auto cometRotationFunction = []( ) -> Eigen::Matrix3d { return Eigen::Matrix3d::Identity( ); };

            // Create ComaModel based on dataset type
            std::unique_ptr< ComaModel > comaModel;
            if( datasetType == 0 )
            {
                // Create with polynomial coefficients
                comaModel = std::make_unique< ComaModel >(
                        polyDataset, molecularWeight, sunStateFunction, cometStateFunction, cometRotationFunction, maxDegree, maxOrder );
            }
            else
            {
                // Create with Stokes coefficients
                comaModel = std::make_unique< ComaModel >(
                        stokesDataset, molecularWeight, sunStateFunction, cometStateFunction, cometRotationFunction, maxDegree, maxOrder );
            }

            // Read all data points from residuals file
            std::ifstream file1( residualsFile1.string( ) );
            BOOST_REQUIRE_MESSAGE( file1.is_open( ), "Cannot open residuals file: " + residualsFile1.string( ) );

            struct DataPoint {
                double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
            };
            std::vector< DataPoint > allPoints;

            std::string line;
            while( std::getline( file1, line ) )
            {
                std::istringstream iss( line );
                DataPoint point;
                if( iss >> point.longitude_deg >> point.latitude_deg >> point.originalData >> point.shEvaluation >> point.difference )
                {
                    allPoints.push_back( point );
                }
            }
            file1.close( );

            // Randomly select 100 points for testing
            const int numTestPoints = 100;
            std::vector< int > selectedIndices;
            std::srand( 12345 );  // Fixed seed for reproducibility

            if( allPoints.size( ) <= numTestPoints )
            {
                for( size_t i = 0; i < allPoints.size( ); ++i ) selectedIndices.push_back( i );
            }
            else
            {
                std::vector< int > allIndices;
                for( size_t i = 0; i < allPoints.size( ); ++i ) allIndices.push_back( i );

                for( int i = 0; i < numTestPoints; ++i )
                {
                    int randomIndex = std::rand( ) % allIndices.size( );
                    selectedIndices.push_back( allIndices[ randomIndex ] );
                    allIndices.erase( allIndices.begin( ) + randomIndex );
                }
            }

            // Test the selected points
            int failCount = 0;
            for( int idx : selectedIndices )
            {
                const DataPoint& point = allPoints[ idx ];

                // Convert to radians
                const double longitude_rad = point.longitude_deg * mathematical_constants::PI / 180.0;
                const double latitude_rad = point.latitude_deg * mathematical_constants::PI / 180.0;

                // Call getTotalNumberDensity from ComaModel (returns actual number density)
                const double computedNumberDensity = comaModel->getTotalNumberDensity( testRadius, longitude_rad, latitude_rad, testTime );

                // Convert to log2 for comparison since test file contains log2(number_density)
                const double computedNumberDensityLog2 = std::log2( computedNumberDensity );
                const double expectedNumberDensityLog2 = point.shEvaluation;

                const double tolerance = 1e-10;
                const double diff = std::abs( computedNumberDensityLog2 - expectedNumberDensityLog2 );

                if( diff > tolerance )
                {
                    failCount++;
                    if( failCount <= 10 )  // Only report first 10 failures
                    {
                        BOOST_TEST_MESSAGE( "Case 1 - Point " << idx << " mismatch: "
                                                              << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                                              << "computed_log2=" << computedNumberDensityLog2
                                                              << ", expected_log2=" << expectedNumberDensityLog2 << ", diff=" << diff );
                    }
                }

                BOOST_CHECK_MESSAGE( diff <= tolerance,
                                     "Number density mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°" );
            }

            BOOST_TEST_MESSAGE( "Test Case 1 (" << datasetTypeName << ", 0°, 4km): Verified " << selectedIndices.size( )
                                                << " randomly selected points out of " << allPoints.size( ) << " total, " << failCount
                                                << " failures" );
            BOOST_CHECK_EQUAL( failCount, 0 );
        }

        // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
        {
            const boost::filesystem::path residualsFile2 = dataDir / "density" / "residual" / "residual_r_cometFixed_ep10-030_10km.txt";
            const double testRadius = 10000.0;                                            // 10 km in meters
            const double testSolarLongitude = 30.0 * mathematical_constants::PI / 180.0;  // 30 degrees in radians
            const double testTime = 490708800;                                            // s since J2000

            // Define state functions for solar longitude = 30°
            // Solar longitude is calculated as atan2(y, x) of the Sun direction in body-fixed frame
            // For solar longitude = 30°, we need Sun at angle 30° from +X axis in XY plane
            // This means: x = cos(30°), y = sin(30°), z = 0

            const double solarDistance = 1.0e11;  // Distance to Sun in meters
            const double sunX = solarDistance * std::cos( testSolarLongitude );
            const double sunY = solarDistance * std::sin( testSolarLongitude );

            // Sun state function - returns position of Sun at 30° angle
            auto sunStateFunction = [ sunX, sunY ]( ) -> Eigen::Vector6d {
                Eigen::Vector6d state = Eigen::Vector6d::Zero( );
                state.segment( 0, 3 ) = Eigen::Vector3d( sunX, sunY, 0.0 );
                return state;
            };

            // Comet state function - returns position of comet at origin
            auto cometStateFunction = []( ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); };

            // Comet rotation function - returns identity matrix
            auto cometRotationFunction = []( ) -> Eigen::Matrix3d { return Eigen::Matrix3d::Identity( ); };

            // Create ComaModel based on dataset type
            std::unique_ptr< ComaModel > comaModel;
            if( datasetType == 0 )
            {
                // Create with polynomial coefficients
                comaModel = std::make_unique< ComaModel >(
                        polyDataset, molecularWeight, sunStateFunction, cometStateFunction, cometRotationFunction, maxDegree, maxOrder );
            }
            else
            {
                // Create with Stokes coefficients
                comaModel = std::make_unique< ComaModel >(
                        stokesDataset, molecularWeight, sunStateFunction, cometStateFunction, cometRotationFunction, maxDegree, maxOrder );
            }

            // Read all data points from residuals file
            std::ifstream file2( residualsFile2.string( ) );
            BOOST_REQUIRE_MESSAGE( file2.is_open( ), "Cannot open residuals file: " + residualsFile2.string( ) );

            struct DataPoint {
                double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
            };
            std::vector< DataPoint > allPoints;

            std::string line;
            while( std::getline( file2, line ) )
            {
                std::istringstream iss( line );
                DataPoint point;
                if( iss >> point.longitude_deg >> point.latitude_deg >> point.originalData >> point.shEvaluation >> point.difference )
                {
                    allPoints.push_back( point );
                }
            }
            file2.close( );

            // Randomly select 50 points for testing
            const int numTestPoints = 50;
            std::vector< int > selectedIndices;
            std::srand( 67890 );  // Different seed from case 1

            if( allPoints.size( ) <= numTestPoints )
            {
                for( size_t i = 0; i < allPoints.size( ); ++i ) selectedIndices.push_back( i );
            }
            else
            {
                std::vector< int > allIndices;
                for( size_t i = 0; i < allPoints.size( ); ++i ) allIndices.push_back( i );

                for( int i = 0; i < numTestPoints; ++i )
                {
                    int randomIndex = std::rand( ) % allIndices.size( );
                    selectedIndices.push_back( allIndices[ randomIndex ] );
                    allIndices.erase( allIndices.begin( ) + randomIndex );
                }
            }

            // Test the selected points
            int failCount = 0;
            for( int idx : selectedIndices )
            {
                const DataPoint& point = allPoints[ idx ];

                // Convert to radians
                const double longitude_rad = point.longitude_deg * mathematical_constants::PI / 180.0;
                const double latitude_rad = point.latitude_deg * mathematical_constants::PI / 180.0;

                // Call getTotalNumberDensity from ComaModel (returns actual number density)
                const double computedNumberDensity = comaModel->getTotalNumberDensity( testRadius, longitude_rad, latitude_rad, testTime );

                // Convert to log2 for comparison since test file contains log2(number_density)
                const double computedNumberDensityLog2 = std::log2( computedNumberDensity );
                const double expectedNumberDensityLog2 = point.shEvaluation;

                const double tolerance = 1e-10;
                const double diff = std::abs( computedNumberDensityLog2 - expectedNumberDensityLog2 );

                if( diff > tolerance )
                {
                    failCount++;
                    if( failCount <= 10 )  // Only report first 10 failures
                    {
                        BOOST_TEST_MESSAGE( "Case 2 - Point " << idx << " mismatch: "
                                                              << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                                              << "computed_log2=" << computedNumberDensityLog2
                                                              << ", expected_log2=" << expectedNumberDensityLog2 << ", diff=" << diff );
                    }
                }

                BOOST_CHECK_MESSAGE( diff <= tolerance,
                                     "Number density mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°" );
            }

            BOOST_TEST_MESSAGE( "Test Case 2 (" << datasetTypeName << ", 30°, 10km): Verified " << selectedIndices.size( )
                                                << " randomly selected points out of " << allPoints.size( ) << " total, " << failCount
                                                << " failures" );
            BOOST_CHECK_EQUAL( failCount, 0 );
        }
    }
}

/**
 * @brief Validates the complete density pipeline against reference values computed from Python interface.
 *
 * Input/Setup:
 * - Creates Stokes dataset with dense grid (100 radii from 4-10km, 37 longitudes from 0-360°)
 * - Reads reference values from file containing: time, radial distance, latitude, longitude,
 *   solar longitude, and expected density
 *
 * Expected Behavior:
 * - ComaModel::getDensity returns mass density values (kg/m³) matching Python reference
 *   within 1e-8 relative tolerance
 * - 100 randomly selected test points all pass validation
 * - Verifies consistency between C++ and Python implementations of the coma model
 */
BOOST_FIXTURE_TEST_CASE( test_coma_model_density_validation_from_python, TestDataPaths )
{
    // This test validates the entire pipeline by using reference values computed from the Python interface.
    // The reference_values.txt file contains: time, radial distance, latitude, longitude, solar longitude,
    // density, and wind velocity components. This test only uses the density column.
    // We use these values as input to calculate density with the verified code and validate against the reference density.

    // Load polynomial coefficients from test data file
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaModelFileProcessor processor( files );

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Molecular weight (kg/mol, water vapor)
    const double molecularWeight = 0.018;

    // Create Stokes dataset from polynomial dataset
    // Use the same grid as in Python: radii from 4 km to 10 km (100 points), sol_longitude 0 to 360° (37 points)
    std::vector< double > radii_m;
    for( int i = 0; i < 100; ++i )
    {
        radii_m.push_back( 4000.0 + i * ( 10000.0 - 4000.0 ) / 99.0 );
    }

    std::vector< double > longitudes_deg;
    for( int i = 0; i < 37; ++i )
    {
        longitudes_deg.push_back( i * 360.0 / 36.0 );
    }

    const ComaStokesDataset stokesDataset = processor.createSHDataset( radii_m, longitudes_deg, maxDegree, maxOrder );

    const boost::filesystem::path referenceFile = dataDir / "reference_values.txt";

    // Read reference values file
    std::ifstream file( referenceFile.string( ) );
    BOOST_REQUIRE_MESSAGE( file.is_open( ), "Cannot open reference values file: " + referenceFile.string( ) );

    struct ReferencePoint {
        double time;
        double radialDistance;
        double latitude;
        double longitude;
        double solarLongitude;
        double density;
    };
    std::vector< ReferencePoint > allPoints;

    std::string line;
    // Skip header line
    std::getline( file, line );

    while( std::getline( file, line ) )
    {
        std::istringstream iss( line );
        ReferencePoint point;
        if( iss >> point.time >> point.radialDistance >> point.latitude >> point.longitude >> point.solarLongitude >> point.density )
        {
            allPoints.push_back( point );
        }
    }
    file.close( );

    BOOST_REQUIRE_MESSAGE( allPoints.size( ) > 0, "No data points found in reference values file" );
    BOOST_TEST_MESSAGE( "Loaded " << allPoints.size( ) << " reference data points" );

    // Randomly select 100 points for testing (or all if less than 100)
    const int numTestPoints = 100;
    std::vector< int > selectedIndices;
    std::srand( 54321 );  // Fixed seed for reproducibility

    if( allPoints.size( ) <= numTestPoints )
    {
        for( size_t i = 0; i < allPoints.size( ); ++i ) selectedIndices.push_back( i );
    }
    else
    {
        std::vector< int > allIndices;
        for( size_t i = 0; i < allPoints.size( ); ++i ) allIndices.push_back( i );

        for( int i = 0; i < numTestPoints; ++i )
        {
            int randomIndex = std::rand( ) % allIndices.size( );
            selectedIndices.push_back( allIndices[ randomIndex ] );
            allIndices.erase( allIndices.begin( ) + randomIndex );
        }
    }

    // Test with Stokes dataset
    int failCount = 0;
    for( int idx : selectedIndices )
    {
        const ReferencePoint& point = allPoints[ idx ];

        // Define state functions based on the reference solar longitude
        const double solarDistance = 1.0e11;
        const double sunX = solarDistance * std::cos( point.solarLongitude );
        const double sunY = solarDistance * std::sin( point.solarLongitude );

        auto sunStateFunction = [ sunX, sunY ]( ) -> Eigen::Vector6d {
            Eigen::Vector6d state = Eigen::Vector6d::Zero( );
            state.segment( 0, 3 ) = Eigen::Vector3d( sunX, sunY, 0.0 );
            return state;
        };

        auto cometStateFunction = []( ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); };

        auto cometRotationFunction = []( ) -> Eigen::Matrix3d { return Eigen::Matrix3d::Identity( ); };

        // Create ComaModel with Stokes coefficients
        ComaModel comaModel(
                stokesDataset, molecularWeight, sunStateFunction, cometStateFunction, cometRotationFunction, maxDegree, maxOrder );

        // Calculate density using the ComaModel
        // getDensity() returns mass density in kg/m³ directly
        const double computedDensity = comaModel.getDensity( point.radialDistance, point.longitude, point.latitude, point.time );

        // Calculate relative error
        const double relativeTolerance = 1e-8;  // 1e-8 relative error
        const double relativeError = std::abs( computedDensity - point.density ) / std::max( std::abs( point.density ), 1e-30 );

        if( relativeError > relativeTolerance )
        {
            failCount++;
            if( failCount <= 10 )  // Only report first 10 failures
            {
                BOOST_TEST_MESSAGE( "Point " << idx << " mismatch: "
                                             << "r=" << point.radialDistance << " m, "
                                             << "lat=" << point.latitude * 180.0 / mathematical_constants::PI << "°, "
                                             << "lon=" << point.longitude * 180.0 / mathematical_constants::PI << "°, "
                                             << "sol=" << point.solarLongitude * 180.0 / mathematical_constants::PI << "°, "
                                             << "computed=" << computedDensity << " kg/m³, "
                                             << "expected=" << point.density << " kg/m³, "
                                             << "rel_error=" << relativeError );
            }
        }

        BOOST_CHECK_MESSAGE( relativeError <= relativeTolerance,
                             "Density mismatch at r=" << point.radialDistance << " m, "
                                                      << "lat=" << point.latitude * 180.0 / mathematical_constants::PI << "°, "
                                                      << "lon=" << point.longitude * 180.0 / mathematical_constants::PI << "°" );
    }

    BOOST_TEST_MESSAGE( "Python reference validation: Verified " << selectedIndices.size( ) << " randomly selected points out of "
                                                                 << allPoints.size( ) << " total, " << failCount << " failures" );
    BOOST_CHECK_EQUAL( failCount, 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== High-Level Interface Tests ====================

/**
 * @brief Test suite for ComaModelFileProcessor API tests.
 *
 * Tests the high-level ComaModelFileProcessor class which provides a unified interface
 * for creating polynomial and Stokes datasets from files. Includes validation tests
 * for constructor parameters, file type detection, and error handling.
 */
BOOST_AUTO_TEST_SUITE( test_high_level_interface )

/**
 * @brief Verifies ComaModelFileProcessor correctly creates polynomial dataset from input files.
 *
 * Input/Setup:
 * - Creates processor from single polynomial coefficient test file
 *
 * Expected Behavior:
 * - createPolyCoefDataset() returns dataset with 1 file
 * - File metadata contains correct source path
 * - Reference radius and max degree match expected values
 */
BOOST_FIXTURE_TEST_CASE( test_poly_coef_processor_create_poly_dataset, TestDataPaths )
{
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaModelFileProcessor processor( files );

    const ComaPolyDataset polyDataset = processor.createPolyCoefDataset( );

    BOOST_CHECK_EQUAL( polyDataset.getNumFiles( ), 1 );
    BOOST_CHECK_EQUAL( polyDataset.getFileMeta( 0 ).sourcePath, files[ 0 ] );
    BOOST_CHECK_EQUAL( polyDataset.getReferenceRadius( 0 ), ExpectedPolyValues::refRadius );
    BOOST_CHECK_EQUAL( polyDataset.getMaxDegreeSH( 0 ), ExpectedPolyValues::maxDegree );
}

/**
 * @brief Verifies ComaModelFileProcessor::createSHDataset produces correct coefficient values.
 *
 * Input/Setup:
 * - Creates processor from polynomial coefficient test file
 * - Uses radii = {4000, 10000} meters and longitudes = {0, 30} degrees
 * - Loads reference Stokes data from pre-computed test files
 *
 * Expected Behavior:
 * - Dataset structure is correct (nFiles=1, nRadii=2, nLongitudes=2, nmax=10)
 * - Coefficient values at (ri=0, li=0) match 4km/0° reference data
 * - Coefficient values at (ri=1, li=1) match 10km/30° reference data
 */
BOOST_FIXTURE_TEST_CASE( test_poly_coef_processor_create_sh_dataset, TestDataPaths )
{
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaModelFileProcessor processor( files );

    // Use the same test points as in test_stokes_coefficients_evaluator
    const std::vector< double > radii_m = { 4000.0, 10000.0 };
    const std::vector< double > longitudes_deg = { 0.0, 30.0 };

    const boost::filesystem::path testFile1 = stokesFileLegacyLon0;
    const boost::filesystem::path testFile2 = stokesFileLegacyLon30;
    StokesTestData expectedData1 = StokesTestData::readFromFile( testFile1.string( ), 10 );
    StokesTestData expectedData2 = StokesTestData::readFromFile( testFile2.string( ), 10 );

    // Test with auto-selected maxima
    const ComaStokesDataset dataset = processor.createSHDataset( radii_m, longitudes_deg );

    BOOST_CHECK_EQUAL( dataset.nFiles( ), 1 );
    BOOST_CHECK_EQUAL( dataset.nRadii( ), 2 );
    BOOST_CHECK_EQUAL( dataset.nLongitudes( ), 2 );
    BOOST_CHECK_EQUAL( dataset.nmax( ), 10 );

    // Verify computed coefficients match expected values at (ri=0, li=0) -> 4000m, 0deg
    auto matrix_pair_0 = dataset.getCoefficientMatrices( 0, 0, 0 );
    auto cosineMatrix0 = matrix_pair_0.first;
    auto sineMatrix0 = matrix_pair_0.second;
    BOOST_CHECK_CLOSE( cosineMatrix0( 0, 0 ), expectedData1.cosineCoeffs( 0, 0 ), 1e-10 );
    BOOST_CHECK_CLOSE( cosineMatrix0( 3, 1 ), expectedData1.cosineCoeffs( 3, 1 ), 1e-10 );
    BOOST_CHECK_CLOSE( sineMatrix0( 6, 2 ), expectedData1.sineCoeffs( 6, 2 ), 1e-10 );
    BOOST_CHECK_CLOSE( sineMatrix0( 7, 5 ), expectedData1.sineCoeffs( 7, 5 ), 1e-10 );

    // Verify computed coefficients at (ri=1, li=1) -> 10000m, 30deg
    auto matrix_pair_1 = dataset.getCoefficientMatrices( 0, 1, 1 );
    auto cosineMatrix1 = matrix_pair_1.first;
    auto sineMatrix1 = matrix_pair_1.second;
    BOOST_CHECK_CLOSE( cosineMatrix1( 0, 0 ), expectedData2.cosineCoeffs( 0, 0 ), 1e-10 );
    BOOST_CHECK_CLOSE( cosineMatrix1( 3, 1 ), expectedData2.cosineCoeffs( 3, 1 ), 1e-10 );
}

/**
 * @brief Verifies ComaModelFileProcessor::createSHFiles creates correctly formatted CSV files.
 *
 * Input/Setup:
 * - Creates processor from polynomial coefficient test file
 * - Calls createSHFiles with radii = {6000, 10000} and longitudes = {0, 30}
 *
 * Expected Behavior:
 * - CSV file is created at expected path (outputDir/stokes_file0.csv)
 * - Meta line contains correct max_degree and max_order values
 * - Radii line contains the specified radius values
 * - Longitudes line header is present
 */
BOOST_FIXTURE_TEST_CASE( test_poly_coef_processor_create_sh_files, TestDataPaths )
{
    std::vector< std::string > files = { testFile.string( ) };
    ComaModelFileProcessor processor( files );

    std::vector< double > radii_m = { 6000.0, 10000.0 };
    std::vector< double > longitudes_deg = { 0.0, 30.0 };

    // Clean output directory
    if( boost::filesystem::exists( outputDir ) )
    {
        boost::filesystem::remove_all( outputDir );
    }

    // Generate CSV files
    processor.createSHFiles( outputDir.string( ), radii_m, longitudes_deg );

    // Check that file was created
    boost::filesystem::path expectedFile = outputDir / "stokes_file0.csv";
    BOOST_CHECK( boost::filesystem::exists( expectedFile ) );

    // Verify file content structure
    std::ifstream ifs( expectedFile.string( ) );
    BOOST_REQUIRE( ifs.is_open( ) );

    std::string line;

    // Meta line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "meta" ) != std::string::npos );
    BOOST_CHECK( line.find( "max_degree=10" ) != std::string::npos );
    BOOST_CHECK( line.find( "max_order=10" ) != std::string::npos );

    // Radii line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "radii [meter]" ) != std::string::npos );
    BOOST_CHECK( line.find( "6000" ) != std::string::npos );
    BOOST_CHECK( line.find( "10000" ) != std::string::npos );

    // Longitudes line
    BOOST_REQUIRE( std::getline( ifs, line ) );
    BOOST_CHECK( line.find( "longitudes [degree]" ) != std::string::npos );

    ifs.close( );
}

/**
 * @brief Verifies round-trip: create SH files from poly data, read back via new processor, verify values.
 *
 * Input/Setup:
 * - Creates original Stokes dataset from polynomial coefficients
 * - Writes dataset to CSV files using createSHFiles
 * - Creates new processor from the SH files directory
 *
 * Expected Behavior:
 * - New processor correctly identifies file type as StokesCoefficients
 * - createSHDataset() (parameterless) reads files and returns dataset with matching structure
 * - Radii and longitude values are preserved
 * - Coefficient values match original within tolerance
 * - Custom file prefix works correctly
 * - Calling createPolyCoefDataset() on SH processor throws std::runtime_error
 */
BOOST_FIXTURE_TEST_CASE( test_sh_processor_from_existing_files, TestDataPaths )
{
    std::vector< std::string > files = { testFile.string( ) };
    ComaModelFileProcessor processor( files );

    std::vector< double > radii_m = { 6000.0, 10000.0 };
    std::vector< double > longitudes_deg = { 0.0, 30.0 };

    // First, create the original SH dataset for comparison
    ComaStokesDataset originalDataset = processor.createSHDataset( radii_m, longitudes_deg, 8, 6 );

    // Create SH files directory
    boost::filesystem::path shFilesDir = outputDir / "sh_files_test";
    boost::filesystem::create_directories( shFilesDir );

    // Generate CSV files using createSHFiles
    processor.createSHFiles( shFilesDir.string( ), radii_m, longitudes_deg, 8, 6 );

    // Verify files were created
    boost::filesystem::path expectedFile = shFilesDir / "stokes_file0.csv";
    BOOST_CHECK( boost::filesystem::exists( expectedFile ) );

    // Now test creating a new processor from SH files and reading the dataset
    ComaModelFileProcessor shProcessor( shFilesDir.string( ) );
    ComaStokesDataset readDataset = shProcessor.createSHDataset( );

    // Verify structure matches original
    BOOST_CHECK_EQUAL( readDataset.nFiles( ), originalDataset.nFiles( ) );
    BOOST_CHECK_EQUAL( readDataset.nRadii( ), originalDataset.nRadii( ) );
    BOOST_CHECK_EQUAL( readDataset.nLongitudes( ), originalDataset.nLongitudes( ) );
    BOOST_CHECK_EQUAL( readDataset.nmax( ), originalDataset.nmax( ) );

    // Verify radii and longitudes match
    const auto& readRadii = readDataset.radii( );
    const auto& readLons = readDataset.lons( );
    const auto& origRadii = originalDataset.radii( );
    const auto& origLons = originalDataset.lons( );

    for( std::size_t i = 0; i < readRadii.size( ); ++i )
    {
        BOOST_CHECK_CLOSE( readRadii[ i ], origRadii[ i ], 1e-10 );
    }
    for( std::size_t i = 0; i < readLons.size( ); ++i )
    {
        BOOST_CHECK_CLOSE( readLons[ i ], origLons[ i ], 1e-10 );
    }

    // Verify selected coefficient values match between original and read datasets
    // Test a few specific coefficients across different radii and longitudes
    auto orig_coeff_0_0_r0_l0 = originalDataset.getCoeff( 0, 0, 0, 0, 0 );
    double orig_C_0_0_r0_l0 = orig_coeff_0_0_r0_l0.first;
    double orig_S_0_0_r0_l0 = orig_coeff_0_0_r0_l0.second;
    auto read_coeff_0_0_r0_l0 = readDataset.getCoeff( 0, 0, 0, 0, 0 );
    double read_C_0_0_r0_l0 = read_coeff_0_0_r0_l0.first;
    double read_S_0_0_r0_l0 = read_coeff_0_0_r0_l0.second;
    BOOST_CHECK_CLOSE( read_C_0_0_r0_l0, orig_C_0_0_r0_l0, 1e-10 );
    BOOST_CHECK_CLOSE( read_S_0_0_r0_l0, orig_S_0_0_r0_l0, 1e-10 );

    auto orig_coeff_2_1_r0_l1 = originalDataset.getCoeff( 0, 0, 1, 2, 1 );
    double orig_C_2_1_r0_l1 = orig_coeff_2_1_r0_l1.first;
    double orig_S_2_1_r0_l1 = orig_coeff_2_1_r0_l1.second;
    auto read_coeff_2_1_r0_l1 = readDataset.getCoeff( 0, 0, 1, 2, 1 );
    double read_C_2_1_r0_l1 = read_coeff_2_1_r0_l1.first;
    double read_S_2_1_r0_l1 = read_coeff_2_1_r0_l1.second;
    BOOST_CHECK_CLOSE( read_C_2_1_r0_l1, orig_C_2_1_r0_l1, 1e-10 );
    BOOST_CHECK_CLOSE( read_S_2_1_r0_l1, orig_S_2_1_r0_l1, 1e-10 );

    auto orig_coeff_3_2_r1_l0 = originalDataset.getCoeff( 0, 1, 0, 3, 2 );
    double orig_C_3_2_r1_l0 = orig_coeff_3_2_r1_l0.first;
    double orig_S_3_2_r1_l0 = orig_coeff_3_2_r1_l0.second;
    auto read_coeff_3_2_r1_l0 = readDataset.getCoeff( 0, 1, 0, 3, 2 );
    double read_C_3_2_r1_l0 = read_coeff_3_2_r1_l0.first;
    double read_S_3_2_r1_l0 = read_coeff_3_2_r1_l0.second;
    BOOST_CHECK_CLOSE( read_C_3_2_r1_l0, orig_C_3_2_r1_l0, 1e-10 );
    BOOST_CHECK_CLOSE( read_S_3_2_r1_l0, orig_S_3_2_r1_l0, 1e-10 );

    auto orig_coeff_5_4_r1_l1 = originalDataset.getCoeff( 0, 1, 1, 5, 4 );
    double orig_C_5_4_r1_l1 = orig_coeff_5_4_r1_l1.first;
    double orig_S_5_4_r1_l1 = orig_coeff_5_4_r1_l1.second;
    auto read_coeff_5_4_r1_l1 = readDataset.getCoeff( 0, 1, 1, 5, 4 );
    double read_C_5_4_r1_l1 = read_coeff_5_4_r1_l1.first;
    double read_S_5_4_r1_l1 = read_coeff_5_4_r1_l1.second;
    BOOST_CHECK_CLOSE( read_C_5_4_r1_l1, orig_C_5_4_r1_l1, 1e-10 );
    BOOST_CHECK_CLOSE( read_S_5_4_r1_l1, orig_S_5_4_r1_l1, 1e-10 );

    // Test with custom prefix
    boost::filesystem::path customPrefixDir = outputDir / "custom_prefix_test";
    boost::filesystem::create_directories( customPrefixDir );

    // Write CSV files with custom prefix
    ComaStokesDatasetWriter::writeCsvAll( originalDataset, customPrefixDir.string( ), "custom" );

    // Read back with custom prefix using new processor
    ComaModelFileProcessor customShProcessor( customPrefixDir.string( ), "custom" );
    ComaStokesDataset customReadDataset = customShProcessor.createSHDataset( );

    // Verify it matches the original
    BOOST_CHECK_EQUAL( customReadDataset.nFiles( ), originalDataset.nFiles( ) );
    BOOST_CHECK_EQUAL( customReadDataset.nmax( ), originalDataset.nmax( ) );

    // Verify one coefficient value
    auto custom_coeff_0_0 = customReadDataset.getCoeff( 0, 0, 0, 0, 0 );
    double custom_C_0_0 = custom_coeff_0_0.first;
    double custom_S_0_0 = custom_coeff_0_0.second;
    BOOST_CHECK_CLOSE( custom_C_0_0, orig_C_0_0_r0_l0, 1e-10 );

    // Test error handling: calling createPolyCoefDataset on SH processor should throw
    BOOST_CHECK_THROW( shProcessor.createPolyCoefDataset( ), std::runtime_error );
}

/**
 * @brief Verifies createSHDataset throws std::invalid_argument for invalid degree/order or empty inputs.
 *
 * Input/Setup:
 * - Creates processor from valid polynomial coefficient test file
 *
 * Expected Behavior:
 * - Requesting maxDegree > available data throws std::invalid_argument
 * - Requesting maxOrder > available data throws std::invalid_argument
 * - Empty radii vector throws std::invalid_argument
 * - Empty longitudes vector throws std::invalid_argument
 */
BOOST_FIXTURE_TEST_CASE( test_poly_coef_processor_validation, TestDataPaths )
{
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaModelFileProcessor processor( files );

    const std::vector< double > radii_m = { 6000.0 };
    const std::vector< double > longitudes_deg = { 30.0 };

    // Test with invalid degree/order requests
    BOOST_CHECK_THROW( processor.createSHDataset( radii_m, longitudes_deg, 15, 10 ), std::invalid_argument );

    BOOST_CHECK_THROW( processor.createSHDataset( radii_m, longitudes_deg, 10, 15 ), std::invalid_argument );

    // Test with empty inputs
    const std::vector< double > emptyVec;
    BOOST_CHECK_THROW( processor.createSHDataset( emptyVec, longitudes_deg ), std::invalid_argument );

    BOOST_CHECK_THROW( processor.createSHDataset( radii_m, emptyVec ), std::invalid_argument );
}

/**
 * @brief Verifies constructor throws std::invalid_argument for empty file list.
 *
 * Input/Setup:
 * - Attempts to create processor with empty std::vector<std::string>
 *
 * Expected Behavior:
 * - Constructor throws std::invalid_argument
 */
BOOST_AUTO_TEST_CASE( test_poly_coef_processor_constructor_validation )
{
    // Test with empty file list
    const std::vector< std::string > emptyFiles;
    BOOST_CHECK_THROW( ComaModelFileProcessor processor( emptyFiles ), std::invalid_argument );
}

/**
 * @brief Verifies constructor throws std::runtime_error for non-existent or empty directory.
 *
 * Input/Setup:
 * - Attempts to create processor with non-existent directory path
 * - Attempts to create processor with empty directory (no SH files)
 *
 * Expected Behavior:
 * - Non-existent directory throws std::runtime_error
 * - Empty directory (no matching SH CSV files) throws std::runtime_error
 */
BOOST_FIXTURE_TEST_CASE( test_sh_processor_constructor_validation, TestDataPaths )
{
    // Test with non-existent directory
    const std::string nonExistentDir = "/path/that/does/not/exist";
    BOOST_CHECK_THROW( ComaModelFileProcessor processor( nonExistentDir ), std::runtime_error );

    // Test with empty directory (no SH files)
    boost::filesystem::path emptyDir = outputDir / "empty_dir";
    boost::filesystem::create_directories( emptyDir );
    BOOST_CHECK_THROW( ComaModelFileProcessor processor( emptyDir.string( ) ), std::runtime_error );
}

/**
 * @brief Verifies file type detection and method availability based on processor type.
 *
 * Input/Setup:
 * - Creates polynomial processor from test file
 * - Creates SH processor from directory containing SH CSV files
 *
 * Expected Behavior:
 * - Poly processor: getFileType() returns PolyCoefficients
 * - Poly processor: createPolyCoefDataset() succeeds
 * - Poly processor: createSHDataset(params) succeeds
 * - Poly processor: createSHDataset() (parameterless) throws std::runtime_error
 * - SH processor: getFileType() returns StokesCoefficients
 * - SH processor: createPolyCoefDataset() throws std::runtime_error
 * - SH processor: createSHDataset() (parameterless) succeeds and returns non-empty dataset
 * - SH processor: createSHDataset(params) throws std::runtime_error
 */
BOOST_FIXTURE_TEST_CASE( test_processor_file_type_behavior, TestDataPaths )
{
    // Create test SH files first
    const std::vector< std::string > files = { testFile.string( ) };
    ComaModelFileProcessor polyProcessor( files );
    std::vector< double > radii_m = { 6000.0, 10000.0 };   // Need at least 2 radii for interpolation
    std::vector< double > longitudes_deg = { 0.0, 30.0 };  // Need at least 2 longitudes for interpolation

    boost::filesystem::path shTestDir = outputDir / "file_type_test";
    boost::filesystem::create_directories( shTestDir );
    polyProcessor.createSHFiles( shTestDir.string( ), radii_m, longitudes_deg, 5, 5 );

    // Test poly processor behavior
    ComaModelFileProcessor polyProc( files );
    BOOST_CHECK_EQUAL( polyProc.getFileType( ), ComaModelFileProcessor::FileType::PolyCoefficients );
    BOOST_CHECK_NO_THROW( polyProc.createPolyCoefDataset( ) );
    BOOST_CHECK_NO_THROW( polyProc.createSHDataset( radii_m, longitudes_deg ) );

    // Test SH processor behavior
    ComaModelFileProcessor shProc( shTestDir.string( ) );
    BOOST_CHECK_EQUAL( shProc.getFileType( ), ComaModelFileProcessor::FileType::StokesCoefficients );
    BOOST_CHECK_THROW( shProc.createPolyCoefDataset( ), std::runtime_error );

    // Test that parameterless version works on SH processor
    BOOST_CHECK_NO_THROW( shProc.createSHDataset( ) );
    ComaStokesDataset shDataset = shProc.createSHDataset( );
    BOOST_CHECK_GT( shDataset.nRadii( ), 0 );

    // Test that parameterized version throws error on SH processor
    BOOST_CHECK_THROW( shProc.createSHDataset( { 1000.0 }, { 45.0 } ), std::runtime_error );
    BOOST_CHECK_THROW( shProc.createSHDataset( radii_m, longitudes_deg ), std::runtime_error );

    // Test that parameterless version throws error on Poly processor
    BOOST_CHECK_THROW( polyProc.createSHDataset( ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== Integration Tests ====================

/**
 * @brief Test suite for end-to-end pipeline integration tests.
 *
 * Tests the complete workflow from polynomial coefficient files through
 * Stokes dataset transformation, CSV file output, and reading back the
 * processed data. Validates full round-trip data integrity.
 */
BOOST_AUTO_TEST_SUITE( test_integration )

/**
 * @brief End-to-end test: poly files → poly dataset → Stokes dataset → CSV files → read back → verify.
 *
 * Input/Setup:
 * - Starts with polynomial coefficient test file
 * - Uses radii = {6000, 10000} meters and longitudes = {0, 30} degrees
 * - Creates output directory for CSV files
 *
 * Expected Behavior:
 * - Step 1: createPolyCoefDataset returns dataset with 1 file
 * - Step 2: createSHDataset returns Stokes dataset with nmax=8 and correct coefficient count
 * - Step 3: writeCsvAll creates output CSV file
 * - Step 4: Output file exists at expected location
 * - Step 5: readFromCsv returns dataset with matching structure and coefficient values
 * - Step 6: Creating new SH processor from folder returns matching dataset
 */
BOOST_FIXTURE_TEST_CASE( test_full_pipeline, TestDataPaths )
{
    // Test the complete pipeline from files to CSV output
    std::vector< std::string > files = { testFile.string( ) };
    ComaModelFileProcessor processor( files );

    // Step 1: Create poly dataset
    ComaPolyDataset polyDataset = processor.createPolyCoefDataset( );
    BOOST_CHECK_EQUAL( polyDataset.getNumFiles( ), 1 );

    // Step 2: Transform to Stokes dataset
    std::vector< double > radii_m = { 6000.0, 10000.0 };
    std::vector< double > longitudes_deg = { 0.0, 30.0 };
    ComaStokesDataset stokesDataset = processor.createSHDataset( radii_m, longitudes_deg, 8, 6 );

    BOOST_CHECK_EQUAL( stokesDataset.nmax( ), 8 );
    BOOST_CHECK_LE( stokesDataset.nCoeffs( ), ( 8 + 1 ) * ( 8 + 2 ) / 2 );

    // Step 3: Write to files
    boost::filesystem::path integratedOutput = outputDir / "integrated";
    boost::filesystem::create_directories( integratedOutput );

    ComaStokesDatasetWriter::writeCsvAll( stokesDataset, integratedOutput.string( ), "integrated" );

    // Step 4: Verify output
    boost::filesystem::path outputFile = integratedOutput / "integrated_file0.csv";
    BOOST_CHECK( boost::filesystem::exists( outputFile ) );

    // Step 5: Test reading back the written data
    ComaStokesDataset readDataset = ComaStokesDatasetReader::readFromCsv( outputFile.string( ) );
    BOOST_CHECK_EQUAL( readDataset.nmax( ), stokesDataset.nmax( ) );
    BOOST_CHECK_EQUAL( readDataset.nFiles( ), stokesDataset.nFiles( ) );
    BOOST_CHECK_EQUAL( readDataset.nRadii( ), stokesDataset.nRadii( ) );
    BOOST_CHECK_EQUAL( readDataset.nLongitudes( ), stokesDataset.nLongitudes( ) );

    // Verify some coefficient values match
    auto orig_pair = stokesDataset.getCoeff( 0, 0, 0, 0, 0 );
    double orig_coeff = orig_pair.first;
    double orig_sine = orig_pair.second;
    auto read_pair = readDataset.getCoeff( 0, 0, 0, 0, 0 );
    double read_coeff = read_pair.first;
    double read_sine = read_pair.second;
    BOOST_CHECK_CLOSE( read_coeff, orig_coeff, 1e-10 );
    BOOST_CHECK_CLOSE( read_sine, orig_sine, 1e-10 );

    // Step 6: Test reading from folder using new SH processor
    ComaModelFileProcessor integratedShProcessor( integratedOutput.string( ), "integrated" );
    ComaStokesDataset folderReadDataset = integratedShProcessor.createSHDataset( );
    BOOST_CHECK_EQUAL( folderReadDataset.nmax( ), stokesDataset.nmax( ) );

    // Verify coefficient values from folder read match original
    auto folder_pair = folderReadDataset.getCoeff( 0, 0, 0, 0, 0 );
    double folder_coeff = folder_pair.first;
    double folder_sine = folder_pair.second;
    BOOST_CHECK_CLOSE( folder_coeff, orig_coeff, 1e-10 );
    BOOST_CHECK_CLOSE( folder_sine, orig_sine, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

// ==================== SurfaceSphericalHarmonicsCalculator Tests ====================

/**
 * @brief Test suite for SurfaceSphericalHarmonicsCalculator surface evaluation.
 *
 * Tests the SurfaceSphericalHarmonicsCalculator class which computes spherical harmonic
 * expansions at specific latitude/longitude coordinates. Validates against
 * reference residual data from pre-computed test files.
 */
BOOST_AUTO_TEST_SUITE( test_spherical_harmonics_calculator )

/**
 * @brief Verifies SurfaceSphericalHarmonicsCalculator::calculateSurfaceSphericalHarmonics computes correct values.
 *
 * Input/Setup:
 * - Reads polynomial coefficients and computes Stokes coefficients
 * - Uses reference residual files containing expected SH evaluation values at various lat/lon points
 *
 * Expected Behavior:
 * - Test Case 1 (0° solar longitude, 4km radius): 100 randomly selected points match
 *   reference SH evaluation values within 1e-10 absolute tolerance
 * - Test Case 2 (30° solar longitude, 10km radius): 100 randomly selected points match
 *   reference values
 * - All tested points report zero failures
 */
BOOST_FIXTURE_TEST_CASE( test_calculate_surface_spherical_harmonics, TestDataPaths )
{
    // Step 1: Get polynomial coefficients from test data file
    const std::vector< std::string > files = { testFile.string( ) };
    const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles( files );

    // Step 2: Get polynomial data from dataset
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients( 0 );
    const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices( 0 );
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius( 0 );
    const double refRadius = dataset.getReferenceRadius( 0 );
    const Eigen::ArrayXXd polyCoefficients = polyCoefs.array( );

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Create SurfaceSphericalHarmonicsCalculator instance
    SurfaceSphericalHarmonicsCalculator calculator;

    // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
    {
        const boost::filesystem::path residualsFile1 = dataDir / "density" / "residual" / "residual_r_cometFixed_ep10-000_04km.txt";
        const double testRadius = 4000.0;                                            // 4 km in meters
        const double testSolarLongitude = 0.0 * mathematical_constants::PI / 180.0;  // 0 degrees in radians

        // Calculate Stokes coefficients for this case
        Eigen::MatrixXd cosineCoefficients, sineCoefficients;
        simulation_setup::StokesCoefficientsEvaluator::evaluate2D( testRadius,
                                                                   testSolarLongitude,
                                                                   polyCoefficients,
                                                                   shIndices,
                                                                   powers,
                                                                   refRadius,
                                                                   cosineCoefficients,
                                                                   sineCoefficients,
                                                                   maxDegree,
                                                                   maxOrder );

        // Read all data points from file
        std::ifstream file1( residualsFile1.string( ) );
        BOOST_REQUIRE_MESSAGE( file1.is_open( ), "Cannot open residuals file: " + residualsFile1.string( ) );

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector< DataPoint > allPoints;

        std::string line;
        while( std::getline( file1, line ) )
        {
            std::istringstream iss( line );
            DataPoint point;
            if( iss >> point.longitude_deg >> point.latitude_deg >> point.originalData >> point.shEvaluation >> point.difference )
            {
                allPoints.push_back( point );
            }
        }
        file1.close( );

        // Randomly select 100 points
        const int numTestPoints = 100;
        std::vector< int > selectedIndices;
        std::srand( 12345 );  // Fixed seed for reproducibility

        if( allPoints.size( ) <= numTestPoints )
        {
            // Use all points if we have fewer than requested
            for( size_t i = 0; i < allPoints.size( ); ++i ) selectedIndices.push_back( i );
        }
        else
        {
            // Randomly select numTestPoints indices without replacement
            std::vector< int > allIndices;
            for( size_t i = 0; i < allPoints.size( ); ++i ) allIndices.push_back( i );

            for( int i = 0; i < numTestPoints; ++i )
            {
                int randomIndex = std::rand( ) % allIndices.size( );
                selectedIndices.push_back( allIndices[ randomIndex ] );
                allIndices.erase( allIndices.begin( ) + randomIndex );
            }
        }

        // Test the selected points
        int failCount = 0;
        for( int idx : selectedIndices )
        {
            const DataPoint& point = allPoints[ idx ];

            // Convert to radians (latitude is already geodetic latitude in the file)
            const double longitude_rad = point.longitude_deg * mathematical_constants::PI / 180.0;
            const double latitude_rad = point.latitude_deg * mathematical_constants::PI / 180.0;

            // Calculate using our implementation
            const double computedResult = calculator.calculateSurfaceSphericalHarmonics(
                    sineCoefficients, cosineCoefficients, latitude_rad, longitude_rad, maxDegree, maxOrder );

            // Compare with SH evaluation from file (column 4)
            const double tolerance = 1e-10;
            const double diff = std::abs( computedResult - point.shEvaluation );

            if( diff > tolerance )
            {
                failCount++;
                if( failCount <= 10 )  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE( "Case 1 - Point " << idx << " mismatch: "
                                                          << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                                          << "computed=" << computedResult << ", expected=" << point.shEvaluation
                                                          << ", diff=" << diff );
                }
            }

            BOOST_CHECK_MESSAGE( diff <= tolerance,
                                 "SH evaluation mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°" );
        }

        BOOST_TEST_MESSAGE( "Test Case 1 (0°, 4km): Verified " << selectedIndices.size( ) << " randomly selected points out of "
                                                               << allPoints.size( ) << " total, " << failCount << " failures" );
        BOOST_CHECK_EQUAL( failCount, 0 );
    }

    // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
    {
        const boost::filesystem::path residualsFile2 = dataDir / "density" / "residual" / "residual_r_cometFixed_ep10-030_10km.txt";
        const double testRadius = 10000.0;                                            // 10 km in meters
        const double testSolarLongitude = 30.0 * mathematical_constants::PI / 180.0;  // 30 degrees in radians

        // Calculate Stokes coefficients for this case
        Eigen::MatrixXd cosineCoefficients, sineCoefficients;
        simulation_setup::StokesCoefficientsEvaluator::evaluate2D( testRadius,
                                                                   testSolarLongitude,
                                                                   polyCoefficients,
                                                                   shIndices,
                                                                   powers,
                                                                   refRadius,
                                                                   cosineCoefficients,
                                                                   sineCoefficients,
                                                                   maxDegree,
                                                                   maxOrder );

        // Read all data points from file
        std::ifstream file2( residualsFile2.string( ) );
        BOOST_REQUIRE_MESSAGE( file2.is_open( ), "Cannot open residuals file: " + residualsFile2.string( ) );

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector< DataPoint > allPoints;

        std::string line;
        while( std::getline( file2, line ) )
        {
            std::istringstream iss( line );
            DataPoint point;
            if( iss >> point.longitude_deg >> point.latitude_deg >> point.originalData >> point.shEvaluation >> point.difference )
            {
                allPoints.push_back( point );
            }
        }
        file2.close( );

        // Randomly select 100 points
        const int numTestPoints = 100;
        std::vector< int > selectedIndices;
        std::srand( 67890 );  // Fixed seed for reproducibility (different from case 1)

        if( allPoints.size( ) <= numTestPoints )
        {
            // Use all points if we have fewer than requested
            for( size_t i = 0; i < allPoints.size( ); ++i ) selectedIndices.push_back( i );
        }
        else
        {
            // Randomly select numTestPoints indices without replacement
            std::vector< int > allIndices;
            for( size_t i = 0; i < allPoints.size( ); ++i ) allIndices.push_back( i );

            for( int i = 0; i < numTestPoints; ++i )
            {
                int randomIndex = std::rand( ) % allIndices.size( );
                selectedIndices.push_back( allIndices[ randomIndex ] );
                allIndices.erase( allIndices.begin( ) + randomIndex );
            }
        }

        // Test the selected points
        int failCount = 0;
        for( int idx : selectedIndices )
        {
            const DataPoint& point = allPoints[ idx ];

            // Convert to radians (latitude is already geodetic latitude in the file)
            const double longitude_rad = point.longitude_deg * mathematical_constants::PI / 180.0;
            const double latitude_rad = point.latitude_deg * mathematical_constants::PI / 180.0;

            // Calculate using our implementation
            const double computedResult = calculator.calculateSurfaceSphericalHarmonics(
                    sineCoefficients, cosineCoefficients, latitude_rad, longitude_rad, maxDegree, maxOrder );

            // Compare with SH evaluation from file (column 4)
            const double tolerance = 1e-10;
            const double diff = std::abs( computedResult - point.shEvaluation );

            if( diff > tolerance )
            {
                failCount++;
                if( failCount <= 10 )  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE( "Case 2 - Point " << idx << " mismatch: "
                                                          << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                                          << "computed=" << computedResult << ", expected=" << point.shEvaluation
                                                          << ", diff=" << diff );
                }
            }

            BOOST_CHECK_MESSAGE( diff <= tolerance,
                                 "SH evaluation mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°" );
        }

        BOOST_TEST_MESSAGE( "Test Case 2 (30°, 10km): Verified " << selectedIndices.size( ) << " randomly selected points out of "
                                                                 << allPoints.size( ) << " total, " << failCount << " failures" );
        BOOST_CHECK_EQUAL( failCount, 0 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
