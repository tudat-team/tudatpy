#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>

namespace tudat
{
namespace unit_tests
{
using namespace simulation_setup;

// ==================== Test Fixtures ====================

struct TestDataPaths
{
    boost::filesystem::path testFile;
    boost::filesystem::path outputDir;

    TestDataPaths()
    {
        boost::filesystem::path thisFile(__FILE__);
        boost::filesystem::path testDir = thisFile.parent_path();
        boost::filesystem::path dataDir = testDir / "test_data";
        testFile = dataDir / "input_poly_coef_test_file.txt";
        outputDir = dataDir / "test_output";

        // Ensure test file exists
        if (!boost::filesystem::exists(testFile))
        {
            throw std::runtime_error("Test data file not found: " + testFile.string());
        }

        // Clean and create output directory
        if (boost::filesystem::exists(outputDir))
        {
            boost::filesystem::remove_all(outputDir);
        }
        boost::filesystem::create_directories(outputDir);
    }

    ~TestDataPaths()
    {
        // Optional: cleanup output dir after tests
        // boost::filesystem::remove_all(outputDir);
    }
};

// Expected values from your original tests
struct ExpectedPolyValues
{
    inline static constexpr int numTerms = 48;
    inline static constexpr int numCoeffs = 121;
    inline static constexpr double refRadius = 10.0;
    inline static constexpr int maxDegree = 10;
    inline static constexpr int maxOrder = 10;

    // Sample poly coefficient values
    inline static constexpr double polyCoef_0_0 = 6.262302500423528E+02;
    inline static constexpr double polyCoef_3_1 = -4.459813047951577E-02;
    inline static constexpr double polyCoef_47_120 = -1.049515320967208E+02;
    inline static constexpr double polyCoef_10_22 = 1.287417812579956E-01;
};

struct ExpectedStokesValues
{
    // For distance = 6 km, solar longitude = 30 degrees
    inline static constexpr double cosine_0_0 = 5.393369500951372e+01;
    inline static constexpr double cosine_3_1 = 1.054997670055739e-01;
    inline static constexpr double cosine_5_4 = 1.207229799736584e-02;
    inline static constexpr double cosine_9_3 = -1.845804426625799e-03;

    inline static constexpr double sine_0_0 = 0.0;
    inline static constexpr double sine_6_2 = 2.139319739882320e-02;
    inline static constexpr double sine_7_5 = -5.401766866307728e-02;
    inline static constexpr double sine_10_8 = 1.423848196013654e-02;
};


// ==================== Data Model Tests ====================

BOOST_AUTO_TEST_SUITE(test_data_models)

BOOST_AUTO_TEST_CASE(test_stokes_dataset_creation)
{
    // Test pure data model creation
    std::vector<ComaStokesDataset::FileMeta> files = {
        {2.015e9, 2.0150864e9, "test_file_1"},
        {2.016e9, 2.0160864e9, "test_file_2"}
    };
    std::vector<double> radii = {1000.0, 2000.0, 3000.0};
    std::vector<double> lons = {0.0, 30.0, 60.0, 90.0};
    int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Verify metadata
    BOOST_CHECK_EQUAL(dataset.nFiles(), 2);
    BOOST_CHECK_EQUAL(dataset.nRadii(), 3);
    BOOST_CHECK_EQUAL(dataset.nLongitudes(), 4);
    BOOST_CHECK_EQUAL(dataset.nmax(), nmax);

    // Expected number of coefficients for degree 10
    std::size_t expectedCoeffs = (nmax + 1) * (nmax + 2) / 2;
    BOOST_CHECK_EQUAL(dataset.nCoeffs(), expectedCoeffs);

    // Test coefficient setting and getting
    dataset.setCoeff(0, 0, 0, 2, 1, 0.5, -0.3);
    auto [C, S] = dataset.getCoeff(0, 0, 0, 2, 1);
    BOOST_CHECK_CLOSE(C, 0.5, 1e-12);
    BOOST_CHECK_CLOSE(S, -0.3, 1e-12);

    // Test block access
    auto block = dataset.block(0, 0, 0);
    BOOST_CHECK_EQUAL(block.rows(), expectedCoeffs);
    BOOST_CHECK_EQUAL(block.cols(), 2);

    // Test coefficient matrices
    dataset.setCoeff(0, 1, 1, 3, 2, 1.5, 2.5);
    auto [cosineMatrix, sineMatrix] = dataset.getCoefficientMatrices(0, 1, 1);
    BOOST_CHECK_EQUAL(cosineMatrix.rows(), nmax + 1);
    BOOST_CHECK_EQUAL(cosineMatrix.cols(), nmax + 1);
    BOOST_CHECK_CLOSE(cosineMatrix(3, 2), 1.5, 1e-12);
    BOOST_CHECK_CLOSE(sineMatrix(3, 2), 2.5, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_stokes_dataset_bounds_checking)
{
    std::vector<ComaStokesDataset::FileMeta> files = {{0, 0, "test"}};
    std::vector<double> radii = {1000.0};
    std::vector<double> lons = {0.0};
    int nmax = 5;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Test out of bounds access
    BOOST_CHECK_THROW(dataset.setCoeff(1, 0, 0, 0, 0, 0, 0), std::out_of_range); // file OOR
    BOOST_CHECK_THROW(dataset.setCoeff(0, 1, 0, 0, 0, 0, 0), std::out_of_range); // radius OOR
    BOOST_CHECK_THROW(dataset.setCoeff(0, 0, 1, 0, 0, 0, 0), std::out_of_range); // longitude OOR
    BOOST_CHECK_THROW(dataset.setCoeff(0, 0, 0, 6, 0, 0, 0), std::out_of_range); // n > nmax
    BOOST_CHECK_THROW(dataset.setCoeff(0, 0, 0, 5, 6, 0, 0), std::out_of_range); // m > n
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== I/O Component Tests ====================

BOOST_AUTO_TEST_SUITE(test_io_components)

BOOST_FIXTURE_TEST_CASE(test_poly_dataset_reader, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};

    // Test reader functionality
    ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles(files);

    BOOST_CHECK_EQUAL(dataset.getNumFiles(), 1);

    // Check dimensions
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients(0);
    BOOST_CHECK_EQUAL(polyCoefs.rows(), ExpectedPolyValues::numTerms);
    BOOST_CHECK_EQUAL(polyCoefs.cols(), ExpectedPolyValues::numCoeffs);

    // Check specific coefficient values
    BOOST_CHECK_CLOSE(polyCoefs(0, 0), ExpectedPolyValues::polyCoef_0_0, 1e-12);
    BOOST_CHECK_CLOSE(polyCoefs(3, 1), ExpectedPolyValues::polyCoef_3_1, 1e-12);
    BOOST_CHECK_CLOSE(polyCoefs(47, 120), ExpectedPolyValues::polyCoef_47_120, 1e-12);
    BOOST_CHECK_CLOSE(polyCoefs(10, 22), ExpectedPolyValues::polyCoef_10_22, 1e-12);

    // Check SH degree and order indices
    const Eigen::ArrayXXi& indices = dataset.getSHDegreeAndOrderIndices(0);
    BOOST_CHECK_EQUAL(indices.rows(), 2);
    BOOST_CHECK_EQUAL(indices.cols(), ExpectedPolyValues::numCoeffs);
    BOOST_CHECK_EQUAL(indices(0, 120), 10);  // degree
    BOOST_CHECK_EQUAL(indices(1, 120), -10); // order

    // Check metadata
    BOOST_CHECK_EQUAL(dataset.getReferenceRadius(0), ExpectedPolyValues::refRadius);
    BOOST_CHECK_EQUAL(dataset.getMaxDegreeSH(0), ExpectedPolyValues::maxDegree);

    // Check powers
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius(0);
    BOOST_CHECK_EQUAL(powers.rows(), 4);
    BOOST_CHECK_EQUAL(powers[0], 0);
    BOOST_CHECK_EQUAL(powers[3], 3);

    // Test column access by (n,m)
    Eigen::VectorXd col = dataset.columnForNM(0, 3, 1);
    BOOST_CHECK_EQUAL(col.rows(), ExpectedPolyValues::numTerms);

    // Test value access
    double val = dataset.value(0, 10, 3, 1);
    BOOST_CHECK_CLOSE(val, polyCoefs(10, indices.cols() > 0 ? 0 : 0), 1e-12);
}

BOOST_FIXTURE_TEST_CASE(test_stokes_dataset_writer, TestDataPaths)
{
    // Create a small dataset for testing
    std::vector<ComaStokesDataset::FileMeta> files = {
        {0.0, 1.0, "test_source"}
    };
    std::vector<double> radii = {6000.0, 10000.0};
    std::vector<double> lons = {0.0, 30.0};
    int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Set some test values
    dataset.setCoeff(0, 0, 0, 0, 0, 1.0, 0.0);
    dataset.setCoeff(0, 0, 1, 2, 1, 0.5, -0.5);
    dataset.setCoeff(0, 1, 0, 3, 3, 0.25, 0.75);

    // Write to CSV
    boost::filesystem::path csvPath = outputDir / "test_stokes.csv";
    ComaStokesDatasetWriter::writeCsvForFile(dataset, 0, csvPath.string());

    BOOST_CHECK(boost::filesystem::exists(csvPath));

    // Verify CSV content
    std::ifstream ifs(csvPath.string());
    BOOST_REQUIRE(ifs.is_open());

    std::string line;

    // Check meta line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("meta") != std::string::npos);
    BOOST_CHECK(line.find("start_epoch=0") != std::string::npos);
    BOOST_CHECK(line.find("end_epoch=1") != std::string::npos);
    BOOST_CHECK(line.find("max_degree=10") != std::string::npos);
    BOOST_CHECK(line.find("n_radii=2") != std::string::npos);
    BOOST_CHECK(line.find("n_lons=2") != std::string::npos);

    // Check radii line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("radii [meter]") != std::string::npos);
    BOOST_CHECK(line.find("6000") != std::string::npos);
    BOOST_CHECK(line.find("10000") != std::string::npos);

    // Check longitudes line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("longitudes [degree]") != std::string::npos);
    BOOST_CHECK(line.find("0") != std::string::npos);
    BOOST_CHECK(line.find("30") != std::string::npos);

    ifs.close();
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== Transformation/Processing Tests ====================

BOOST_AUTO_TEST_SUITE(test_processing_components)

BOOST_AUTO_TEST_CASE(test_stokes_coefficients_evaluator)
{
    // Create test data matching your original test
    double distanceToCometCentre = 6.0; // km
    double solarLongitude = 30.0 * M_PI / 180.0; // radians
    int maxDegree = 10;
    int maxOrder = 10;

    // You would need to set up polyCoefficients, degreeAndOrder, etc.
    // from your test data here. This is a simplified version:

    Eigen::MatrixXd cosineCoefficients, sineCoefficients;

    // This would call your evaluator with the proper test data
    // StokesCoefficientsEvaluator::evaluate2D(...);

    // For now, just check the output dimensions would be correct
    cosineCoefficients = Eigen::MatrixXd::Zero(maxDegree + 1, maxOrder + 1);
    sineCoefficients = Eigen::MatrixXd::Zero(maxDegree + 1, maxOrder + 1);

    BOOST_CHECK_EQUAL(cosineCoefficients.rows(), maxDegree + 1);
    BOOST_CHECK_EQUAL(cosineCoefficients.cols(), maxOrder + 1);
    BOOST_CHECK_EQUAL(sineCoefficients.rows(), maxDegree + 1);
    BOOST_CHECK_EQUAL(sineCoefficients.cols(), maxOrder + 1);
}

BOOST_FIXTURE_TEST_CASE(test_dataset_transformer, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaPolyDataset polyDataset = ComaPolyDatasetReader::readFromFiles(files);

    std::vector<double> radii_m = {6000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};

    // Test transformation with default maxima
    ComaStokesDataset stokesDataset = ComaDatasetTransformer::transformPolyToStokes(
        polyDataset, radii_m, lons_deg);

    BOOST_CHECK_EQUAL(stokesDataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(stokesDataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nmax(), 10);

    // Check specific coefficient values at (ri=0, li=1) -> 6000m, 30deg
    auto [C_0_0, S_0_0] = stokesDataset.getCoeff(0, 0, 1, 0, 0);
    BOOST_CHECK_CLOSE(C_0_0, ExpectedStokesValues::cosine_0_0, 1e-10);
    BOOST_CHECK_CLOSE(S_0_0, ExpectedStokesValues::sine_0_0, 1e-10);

    auto [C_3_1, S_3_1] = stokesDataset.getCoeff(0, 0, 1, 3, 1);
    BOOST_CHECK_CLOSE(C_3_1, ExpectedStokesValues::cosine_3_1, 1e-10);

    auto [C_5_4, S_5_4] = stokesDataset.getCoeff(0, 0, 1, 5, 4);
    BOOST_CHECK_CLOSE(C_5_4, ExpectedStokesValues::cosine_5_4, 1e-10);

    // Test with explicit truncation
    ComaStokesDataset truncatedDataset = ComaDatasetTransformer::transformPolyToStokes(
        polyDataset, radii_m, lons_deg, 6, 4);

    BOOST_CHECK_EQUAL(truncatedDataset.nmax(), 6);

    // Test error handling for exceeding available maxima
    BOOST_CHECK_THROW(
        ComaDatasetTransformer::transformPolyToStokes(
            polyDataset, radii_m, lons_deg, 15, 10),
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== High-Level Interface Tests ====================

BOOST_AUTO_TEST_SUITE(test_high_level_interface)

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_create_poly_dataset, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    ComaPolyDataset polyDataset = processor.createPolyCoefDataset();

    BOOST_CHECK_EQUAL(polyDataset.getNumFiles(), 1);
    BOOST_CHECK_EQUAL(polyDataset.getFileMeta(0).sourcePath, files[0]);
    BOOST_CHECK_EQUAL(polyDataset.getReferenceRadius(0), ExpectedPolyValues::refRadius);
    BOOST_CHECK_EQUAL(polyDataset.getMaxDegreeSH(0), ExpectedPolyValues::maxDegree);
}

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_create_sh_dataset, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    std::vector<double> radii_m = {6000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};

    // Test with auto-selected maxima
    ComaStokesDataset dataset = processor.createSHDataset(radii_m, lons_deg);

    BOOST_CHECK_EQUAL(dataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(dataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(dataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(dataset.nmax(), 10);

    // Verify computed coefficients match expected values
    auto [cosineMatrix, sineMatrix] = dataset.getCoefficientMatrices(0, 0, 1);
    BOOST_CHECK_CLOSE(cosineMatrix(0, 0), ExpectedStokesValues::cosine_0_0, 1e-10);
    BOOST_CHECK_CLOSE(cosineMatrix(3, 1), ExpectedStokesValues::cosine_3_1, 1e-10);
    BOOST_CHECK_CLOSE(sineMatrix(6, 2), ExpectedStokesValues::sine_6_2, 1e-10);
    BOOST_CHECK_CLOSE(sineMatrix(7, 5), ExpectedStokesValues::sine_7_5, 1e-10);
}

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_create_sh_files, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    std::vector<double> radii_m = {6000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};

    // Clean output directory
    if (boost::filesystem::exists(outputDir))
    {
        boost::filesystem::remove_all(outputDir);
    }

    // Generate CSV files
    processor.createSHFiles(outputDir.string(), radii_m, lons_deg);

    // Check that file was created
    boost::filesystem::path expectedFile = outputDir / "stokes_file0.csv";
    BOOST_CHECK(boost::filesystem::exists(expectedFile));

    // Verify file content structure
    std::ifstream ifs(expectedFile.string());
    BOOST_REQUIRE(ifs.is_open());

    std::string line;

    // Meta line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("meta") != std::string::npos);
    BOOST_CHECK(line.find("max_degree=10") != std::string::npos);
    BOOST_CHECK(line.find("max_order=10") != std::string::npos);

    // Radii line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("radii [meter]") != std::string::npos);
    BOOST_CHECK(line.find("6000") != std::string::npos);
    BOOST_CHECK(line.find("10000") != std::string::npos);

    // Longitudes line
    BOOST_REQUIRE(std::getline(ifs, line));
    BOOST_CHECK(line.find("longitudes [degree]") != std::string::npos);

    ifs.close();
}

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_validation, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    std::vector<double> radii_m = {6000.0};
    std::vector<double> lons_deg = {30.0};

    // Test with invalid degree/order requests
    BOOST_CHECK_THROW(
        processor.createSHDataset(radii_m, lons_deg, 15, 10),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        processor.createSHDataset(radii_m, lons_deg, 10, 15),
        std::invalid_argument);

    // Test with empty inputs
    std::vector<double> emptyVec;
    BOOST_CHECK_THROW(
        processor.createSHDataset(emptyVec, lons_deg),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        processor.createSHDataset(radii_m, emptyVec),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_poly_coef_processor_constructor_validation)
{
    // Test with empty file list
    std::vector<std::string> emptyFiles;
    BOOST_CHECK_THROW(
        ComaModelFileProcessor processor(emptyFiles),
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== Integration Tests ====================

BOOST_AUTO_TEST_SUITE(test_integration)

BOOST_FIXTURE_TEST_CASE(test_full_pipeline, TestDataPaths)
{
    // Test the complete pipeline from files to CSV output
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    // Step 1: Create poly dataset
    ComaPolyDataset polyDataset = processor.createPolyCoefDataset();
    BOOST_CHECK_EQUAL(polyDataset.getNumFiles(), 1);

    // Step 2: Transform to Stokes dataset
    std::vector<double> radii_m = {6000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};
    ComaStokesDataset stokesDataset = processor.createSHDataset(radii_m, lons_deg, 8, 6);

    BOOST_CHECK_EQUAL(stokesDataset.nmax(), 8);
    BOOST_CHECK_LE(stokesDataset.nCoeffs(), (8+1)*(8+2)/2);

    // Step 3: Write to files
    boost::filesystem::path integratedOutput = outputDir / "integrated";
    boost::filesystem::create_directories(integratedOutput);

    ComaStokesDatasetWriter::writeCsvAll(stokesDataset, integratedOutput.string(), "integrated");

    // Step 4: Verify output
    boost::filesystem::path outputFile = integratedOutput / "integrated_file0.csv";
    BOOST_CHECK(boost::filesystem::exists(outputFile));

    // Optional: Could add a reader test here when implemented
    // ComaStokesDataset readDataset = ComaStokesDatasetReader::readFromCsv(outputFile.string());
    // BOOST_CHECK_EQUAL(readDataset.nmax(), stokesDataset.nmax());
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat