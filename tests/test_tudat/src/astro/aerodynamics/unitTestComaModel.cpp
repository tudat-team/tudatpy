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
        const boost::filesystem::path thisFile(__FILE__);
        const boost::filesystem::path testDir = thisFile.parent_path();
        const boost::filesystem::path dataDir = testDir / "test_data";
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
    const std::vector<ComaStokesDataset::FileMeta> files = {
        {2.015e9, 2.0150864e9, "test_file_1"},
        {2.016e9, 2.0160864e9, "test_file_2"}
    };
    const std::vector<double> radii = {1000.0, 2000.0, 3000.0};
    const std::vector<double> lons = {0.0, 30.0, 60.0, 90.0};
    constexpr int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Verify metadata
    BOOST_CHECK_EQUAL(dataset.nFiles(), 2);
    BOOST_CHECK_EQUAL(dataset.nRadii(), 3);
    BOOST_CHECK_EQUAL(dataset.nLongitudes(), 4);
    BOOST_CHECK_EQUAL(dataset.nmax(), nmax);

    // Expected number of coefficients for degree 10
    constexpr std::size_t expectedCoeffs = (nmax + 1) * (nmax + 2) / 2;
    BOOST_CHECK_EQUAL(dataset.nCoeffs(), expectedCoeffs);

    // Test coefficient setting and getting
    dataset.setCoeff(0, 0, 0, 2, 1, 0.5, -0.3);
    auto [C, S] = dataset.getCoeff(0, 0, 0, 2, 1);
    BOOST_CHECK_CLOSE(C, 0.5, 1e-12);
    BOOST_CHECK_CLOSE(S, -0.3, 1e-12);

    // Test block access
    const auto block = dataset.block(0, 0, 0);
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
    const std::vector<ComaStokesDataset::FileMeta> files = {{0, 0, "test"}};
    const std::vector<double> radii = {1000.0};
    const std::vector<double> lons = {0.0};
    constexpr int nmax = 5;

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
    const std::vector<std::string> files = {testFile.string()};

    // Test reader functionality
    const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles(files);

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
    const Eigen::VectorXd col = dataset.columnForNM(0, 3, 1);
    BOOST_CHECK_EQUAL(col.rows(), ExpectedPolyValues::numTerms);

    // Test value access
    const double val = dataset.value(0, 10, 3, 1);
    BOOST_CHECK_CLOSE(val, col[10], 1e-12);
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

BOOST_FIXTURE_TEST_CASE(test_stokes_dataset_reader_from_csv, TestDataPaths)
{
    // First, create and write a test dataset
    std::vector<ComaStokesDataset::FileMeta> files = {
        {2.015e9, 2.0150864e9, "test_source"}
    };
    std::vector<double> radii = {6000.0, 10000.0};
    std::vector<double> lons = {0.0, 30.0};
    int nmax = 8;

    ComaStokesDataset originalDataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Set some known test values
    originalDataset.setCoeff(0, 0, 0, 0, 0, 54.0, 0.0);
    originalDataset.setCoeff(0, 0, 0, 2, 1, -0.232, 0.138);
    originalDataset.setCoeff(0, 0, 1, 0, 0, 53.93, 0.0);
    originalDataset.setCoeff(0, 1, 0, 1, 1, -1.75, 0.407);
    originalDataset.setCoeff(0, 1, 1, 3, 2, -0.084, -0.026);

    // Write to CSV
    boost::filesystem::path csvPath = outputDir / "test_reader.csv";
    ComaStokesDatasetWriter::writeCsvForFile(originalDataset, 0, csvPath.string());

    // Now test reading it back
    ComaStokesDataset readDataset = ComaStokesDatasetReader::readFromCsv(csvPath.string());

    // Verify structure
    BOOST_CHECK_EQUAL(readDataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(readDataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(readDataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(readDataset.nmax(), nmax);

    // Verify metadata
    const auto& filesMeta = readDataset.files();
    BOOST_CHECK_CLOSE(filesMeta[0].start_epoch, 2.015e9, 1e-6);
    BOOST_CHECK_CLOSE(filesMeta[0].end_epoch, 2.0150864e9, 1e-6);
    BOOST_CHECK_EQUAL(filesMeta[0].source_tag, "test_source");

    // Verify radii and longitudes
    const auto& readRadii = readDataset.radii();
    const auto& readLons = readDataset.lons();
    BOOST_CHECK_CLOSE(readRadii[0], 6000.0, 1e-10);
    BOOST_CHECK_CLOSE(readRadii[1], 10000.0, 1e-10);
    BOOST_CHECK_CLOSE(readLons[0], 0.0, 1e-10);
    BOOST_CHECK_CLOSE(readLons[1], 30.0, 1e-10);

    // Verify coefficient values
    auto [C_0_0_r0_l0, S_0_0_r0_l0] = readDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(C_0_0_r0_l0, 54.0, 1e-10);
    BOOST_CHECK_CLOSE(S_0_0_r0_l0, 0.0, 1e-10);

    auto [C_2_1_r0_l0, S_2_1_r0_l0] = readDataset.getCoeff(0, 0, 0, 2, 1);
    BOOST_CHECK_CLOSE(C_2_1_r0_l0, -0.232, 1e-10);
    BOOST_CHECK_CLOSE(S_2_1_r0_l0, 0.138, 1e-10);

    auto [C_0_0_r0_l1, S_0_0_r0_l1] = readDataset.getCoeff(0, 0, 1, 0, 0);
    BOOST_CHECK_CLOSE(C_0_0_r0_l1, 53.93, 1e-10);

    auto [C_1_1_r1_l0, S_1_1_r1_l0] = readDataset.getCoeff(0, 1, 0, 1, 1);
    BOOST_CHECK_CLOSE(C_1_1_r1_l0, -1.75, 1e-10);
    BOOST_CHECK_CLOSE(S_1_1_r1_l0, 0.407, 1e-10);

    auto [C_3_2_r1_l1, S_3_2_r1_l1] = readDataset.getCoeff(0, 1, 1, 3, 2);
    BOOST_CHECK_CLOSE(C_3_2_r1_l1, -0.084, 1e-10);
    BOOST_CHECK_CLOSE(S_3_2_r1_l1, -0.026, 1e-10);
}

BOOST_FIXTURE_TEST_CASE(test_stokes_dataset_reader_from_csv_folder, TestDataPaths)
{
    // Create multiple test datasets (simulating multiple time epochs)
    std::vector<ComaStokesDataset::FileMeta> files1 = {
        {0.0, 1.0, "test_file_1"}
    };
    std::vector<ComaStokesDataset::FileMeta> files2 = {
        {1.0, 2.0, "test_file_2"}
    };

    std::vector<double> radii = {6000.0, 8000.0};
    std::vector<double> lons = {0.0, 45.0};
    int nmax = 5;

    ComaStokesDataset dataset1 = ComaStokesDataset::create(files1, radii, lons, nmax);
    ComaStokesDataset dataset2 = ComaStokesDataset::create(files2, radii, lons, nmax);

    // Set different coefficients for each dataset
    dataset1.setCoeff(0, 0, 0, 0, 0, 10.0, 0.0);
    dataset1.setCoeff(0, 0, 1, 1, 1, 1.5, -0.5);
    dataset1.setCoeff(0, 1, 0, 2, 0, 0.25, 0.0);

    dataset2.setCoeff(0, 0, 0, 0, 0, 20.0, 0.0);
    dataset2.setCoeff(0, 0, 1, 1, 1, 2.5, -1.0);
    dataset2.setCoeff(0, 1, 0, 2, 0, 0.5, 0.0);

    // Create folder for multi-file test
    boost::filesystem::path folderPath = outputDir / "multi_file_test";
    boost::filesystem::create_directories(folderPath);

    // Write both datasets to separate CSV files
    ComaStokesDatasetWriter::writeCsvForFile(dataset1, 0, (folderPath / "test_file0.csv").string());
    ComaStokesDatasetWriter::writeCsvForFile(dataset2, 0, (folderPath / "test_file1.csv").string());

    // Test reading from folder
    ComaStokesDataset multiDataset = ComaStokesDatasetReader::readFromCsvFolder(folderPath.string(), "test");

    // Verify structure
    BOOST_CHECK_EQUAL(multiDataset.nFiles(), 2);
    BOOST_CHECK_EQUAL(multiDataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(multiDataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(multiDataset.nmax(), nmax);

    // Verify file metadata
    const auto& filesMeta = multiDataset.files();
    BOOST_CHECK_CLOSE(filesMeta[0].start_epoch, 0.0, 1e-10);
    BOOST_CHECK_CLOSE(filesMeta[0].end_epoch, 1.0, 1e-10);
    BOOST_CHECK_EQUAL(filesMeta[0].source_tag, "test_file_1");
    BOOST_CHECK_CLOSE(filesMeta[1].start_epoch, 1.0, 1e-10);
    BOOST_CHECK_CLOSE(filesMeta[1].end_epoch, 2.0, 1e-10);
    BOOST_CHECK_EQUAL(filesMeta[1].source_tag, "test_file_2");

    // Verify coefficient values from first file
    auto [C1_0_0, S1_0_0] = multiDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(C1_0_0, 10.0, 1e-10);

    auto [C1_1_1, S1_1_1] = multiDataset.getCoeff(0, 0, 1, 1, 1);
    BOOST_CHECK_CLOSE(C1_1_1, 1.5, 1e-10);
    BOOST_CHECK_CLOSE(S1_1_1, -0.5, 1e-10);

    // Verify coefficient values from second file
    auto [C2_0_0, S2_0_0] = multiDataset.getCoeff(1, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(C2_0_0, 20.0, 1e-10);

    auto [C2_1_1, S2_1_1] = multiDataset.getCoeff(1, 0, 1, 1, 1);
    BOOST_CHECK_CLOSE(C2_1_1, 2.5, 1e-10);
    BOOST_CHECK_CLOSE(S2_1_1, -1.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== Transformation/Processing Tests ====================

BOOST_AUTO_TEST_SUITE(test_processing_components)

BOOST_FIXTURE_TEST_CASE(test_stokes_coefficients_evaluator, TestDataPaths)
{
    // Load test data to get the polynomial coefficients and metadata
    const std::vector<std::string> files = {testFile.string()};
    const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles(files);

    // Test parameters matching ExpectedStokesValues
    const double radius = 6000.0; // m (distance to comet center)
    const double solarLongitude = 30.0 * M_PI / 180.0; // 30 degrees in radians
    const int maxDegree = 10;
    const int maxOrder = 10;

    // Get the polynomial data from the dataset
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients(0);
    const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices(0);
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius(0);
    const double refRadius = dataset.getReferenceRadius(0);

    // Convert to ArrayXXd as required by the evaluator
    const Eigen::ArrayXXd polyCoefficients = polyCoefs.array();

    // Output matrices
    Eigen::MatrixXd cosineCoefficients, sineCoefficients;

    // Call the evaluator
    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius,
        solarLongitude,
        polyCoefficients,
        shIndices,
        powers,
        refRadius,
        cosineCoefficients,
        sineCoefficients,
        maxDegree,
        maxOrder
    );

    // Check output dimensions
    BOOST_CHECK_EQUAL(cosineCoefficients.rows(), maxDegree + 1);
    BOOST_CHECK_EQUAL(cosineCoefficients.cols(), maxOrder + 1);
    BOOST_CHECK_EQUAL(sineCoefficients.rows(), maxDegree + 1);
    BOOST_CHECK_EQUAL(sineCoefficients.cols(), maxOrder + 1);

    // Verify specific coefficient values against expected values
    BOOST_CHECK_CLOSE(cosineCoefficients(0, 0), ExpectedStokesValues::cosine_0_0, 1e-10);
    BOOST_CHECK_CLOSE(cosineCoefficients(3, 1), ExpectedStokesValues::cosine_3_1, 1e-10);
    BOOST_CHECK_CLOSE(cosineCoefficients(5, 4), ExpectedStokesValues::cosine_5_4, 1e-10);
    BOOST_CHECK_CLOSE(cosineCoefficients(9, 3), ExpectedStokesValues::cosine_9_3, 1e-10);

    // Check sine coefficients
    BOOST_CHECK_CLOSE(sineCoefficients(0, 0), ExpectedStokesValues::sine_0_0, 1e-10);
    BOOST_CHECK_CLOSE(sineCoefficients(6, 2), ExpectedStokesValues::sine_6_2, 1e-10);
    BOOST_CHECK_CLOSE(sineCoefficients(7, 5), ExpectedStokesValues::sine_7_5, 1e-10);
    BOOST_CHECK_CLOSE(sineCoefficients(10, 8), ExpectedStokesValues::sine_10_8, 1e-10);

    // Test with truncated degree/order
    Eigen::MatrixXd truncatedCosine, truncatedSine;
    const int truncatedMaxDegree = 5;
    const int truncatedMaxOrder = 3;

    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius,
        solarLongitude,
        polyCoefficients,
        shIndices,
        powers,
        refRadius,
        truncatedCosine,
        truncatedSine,
        truncatedMaxDegree,
        truncatedMaxOrder
    );

    // Check truncated dimensions
    BOOST_CHECK_EQUAL(truncatedCosine.rows(), truncatedMaxDegree + 1);
    BOOST_CHECK_EQUAL(truncatedCosine.cols(), truncatedMaxOrder + 1);
    BOOST_CHECK_EQUAL(truncatedSine.rows(), truncatedMaxDegree + 1);
    BOOST_CHECK_EQUAL(truncatedSine.cols(), truncatedMaxOrder + 1);

    // Verify that truncated results match the corresponding elements from full evaluation
    BOOST_CHECK_CLOSE(truncatedCosine(0, 0), cosineCoefficients(0, 0), 1e-12);
    BOOST_CHECK_CLOSE(truncatedCosine(3, 1), cosineCoefficients(3, 1), 1e-12);
    BOOST_CHECK_CLOSE(truncatedSine(0, 0), sineCoefficients(0, 0), 1e-12);

    // Test error handling for excessive degree/order requests
    Eigen::MatrixXd dummyCosine, dummySine;
    BOOST_CHECK_THROW(
        simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
            radius, solarLongitude, polyCoefficients, shIndices, powers, refRadius,
            dummyCosine, dummySine, 15, 10), // Exceeds available maxDegree (10)
        std::runtime_error
    );

    BOOST_CHECK_THROW(
        simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
            radius, solarLongitude, polyCoefficients, shIndices, powers, refRadius,
            dummyCosine, dummySine, 10, 15), // Exceeds available maxOrder (10)
        std::runtime_error
    );
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
    const std::vector<std::string> files = {testFile.string()};
    const ComaModelFileProcessor processor(files);

    const ComaPolyDataset polyDataset = processor.createPolyCoefDataset();

    BOOST_CHECK_EQUAL(polyDataset.getNumFiles(), 1);
    BOOST_CHECK_EQUAL(polyDataset.getFileMeta(0).sourcePath, files[0]);
    BOOST_CHECK_EQUAL(polyDataset.getReferenceRadius(0), ExpectedPolyValues::refRadius);
    BOOST_CHECK_EQUAL(polyDataset.getMaxDegreeSH(0), ExpectedPolyValues::maxDegree);
}

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_create_sh_dataset, TestDataPaths)
{
    const std::vector<std::string> files = {testFile.string()};
    const ComaModelFileProcessor processor(files);

    const std::vector<double> radii_m = {6000.0, 10000.0};
    const std::vector<double> lons_deg = {0.0, 30.0};

    // Test with auto-selected maxima
    const ComaStokesDataset dataset = processor.createSHDataset(radii_m, lons_deg);

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

BOOST_FIXTURE_TEST_CASE(test_sh_processor_from_existing_files, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor processor(files);

    std::vector<double> radii_m = {6000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};

    // First, create the original SH dataset for comparison
    ComaStokesDataset originalDataset = processor.createSHDataset(radii_m, lons_deg, 8, 6);

    // Create SH files directory
    boost::filesystem::path shFilesDir = outputDir / "sh_files_test";
    boost::filesystem::create_directories(shFilesDir);

    // Generate CSV files using createSHFiles
    processor.createSHFiles(shFilesDir.string(), radii_m, lons_deg, 8, 6);

    // Verify files were created
    boost::filesystem::path expectedFile = shFilesDir / "stokes_file0.csv";
    BOOST_CHECK(boost::filesystem::exists(expectedFile));

    // Now test creating a new processor from SH files and reading the dataset
    ComaModelFileProcessor shProcessor(shFilesDir.string());
    ComaStokesDataset readDataset = shProcessor.createSHDataset({}, {});

    // Verify structure matches original
    BOOST_CHECK_EQUAL(readDataset.nFiles(), originalDataset.nFiles());
    BOOST_CHECK_EQUAL(readDataset.nRadii(), originalDataset.nRadii());
    BOOST_CHECK_EQUAL(readDataset.nLongitudes(), originalDataset.nLongitudes());
    BOOST_CHECK_EQUAL(readDataset.nmax(), originalDataset.nmax());

    // Verify radii and longitudes match
    const auto& readRadii = readDataset.radii();
    const auto& readLons = readDataset.lons();
    const auto& origRadii = originalDataset.radii();
    const auto& origLons = originalDataset.lons();

    for (std::size_t i = 0; i < readRadii.size(); ++i)
    {
        BOOST_CHECK_CLOSE(readRadii[i], origRadii[i], 1e-10);
    }
    for (std::size_t i = 0; i < readLons.size(); ++i)
    {
        BOOST_CHECK_CLOSE(readLons[i], origLons[i], 1e-10);
    }

    // Verify selected coefficient values match between original and read datasets
    // Test a few specific coefficients across different radii and longitudes
    auto [orig_C_0_0_r0_l0, orig_S_0_0_r0_l0] = originalDataset.getCoeff(0, 0, 0, 0, 0);
    auto [read_C_0_0_r0_l0, read_S_0_0_r0_l0] = readDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(read_C_0_0_r0_l0, orig_C_0_0_r0_l0, 1e-10);
    BOOST_CHECK_CLOSE(read_S_0_0_r0_l0, orig_S_0_0_r0_l0, 1e-10);

    auto [orig_C_2_1_r0_l1, orig_S_2_1_r0_l1] = originalDataset.getCoeff(0, 0, 1, 2, 1);
    auto [read_C_2_1_r0_l1, read_S_2_1_r0_l1] = readDataset.getCoeff(0, 0, 1, 2, 1);
    BOOST_CHECK_CLOSE(read_C_2_1_r0_l1, orig_C_2_1_r0_l1, 1e-10);
    BOOST_CHECK_CLOSE(read_S_2_1_r0_l1, orig_S_2_1_r0_l1, 1e-10);

    auto [orig_C_3_2_r1_l0, orig_S_3_2_r1_l0] = originalDataset.getCoeff(0, 1, 0, 3, 2);
    auto [read_C_3_2_r1_l0, read_S_3_2_r1_l0] = readDataset.getCoeff(0, 1, 0, 3, 2);
    BOOST_CHECK_CLOSE(read_C_3_2_r1_l0, orig_C_3_2_r1_l0, 1e-10);
    BOOST_CHECK_CLOSE(read_S_3_2_r1_l0, orig_S_3_2_r1_l0, 1e-10);

    auto [orig_C_5_4_r1_l1, orig_S_5_4_r1_l1] = originalDataset.getCoeff(0, 1, 1, 5, 4);
    auto [read_C_5_4_r1_l1, read_S_5_4_r1_l1] = readDataset.getCoeff(0, 1, 1, 5, 4);
    BOOST_CHECK_CLOSE(read_C_5_4_r1_l1, orig_C_5_4_r1_l1, 1e-10);
    BOOST_CHECK_CLOSE(read_S_5_4_r1_l1, orig_S_5_4_r1_l1, 1e-10);

    // Test with custom prefix
    boost::filesystem::path customPrefixDir = outputDir / "custom_prefix_test";
    boost::filesystem::create_directories(customPrefixDir);

    // Write CSV files with custom prefix
    ComaStokesDatasetWriter::writeCsvAll(originalDataset, customPrefixDir.string(), "custom");

    // Read back with custom prefix using new processor
    ComaModelFileProcessor customShProcessor(customPrefixDir.string(), "custom");
    ComaStokesDataset customReadDataset = customShProcessor.createSHDataset({}, {});

    // Verify it matches the original
    BOOST_CHECK_EQUAL(customReadDataset.nFiles(), originalDataset.nFiles());
    BOOST_CHECK_EQUAL(customReadDataset.nmax(), originalDataset.nmax());

    // Verify one coefficient value
    auto [custom_C_0_0, custom_S_0_0] = customReadDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(custom_C_0_0, orig_C_0_0_r0_l0, 1e-10);

    // Test error handling: calling createPolyCoefDataset on SH processor should throw
    BOOST_CHECK_THROW(
        shProcessor.createPolyCoefDataset(),
        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_poly_coef_processor_validation, TestDataPaths)
{
    const std::vector<std::string> files = {testFile.string()};
    const ComaModelFileProcessor processor(files);

    const std::vector<double> radii_m = {6000.0};
    const std::vector<double> lons_deg = {30.0};

    // Test with invalid degree/order requests
    BOOST_CHECK_THROW(
        processor.createSHDataset(radii_m, lons_deg, 15, 10),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        processor.createSHDataset(radii_m, lons_deg, 10, 15),
        std::invalid_argument);

    // Test with empty inputs
    const std::vector<double> emptyVec;
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
    const std::vector<std::string> emptyFiles;
    BOOST_CHECK_THROW(
        ComaModelFileProcessor processor(emptyFiles),
        std::invalid_argument);
}

BOOST_FIXTURE_TEST_CASE(test_sh_processor_constructor_validation, TestDataPaths)
{
    // Test with non-existent directory
    const std::string nonExistentDir = "/path/that/does/not/exist";
    BOOST_CHECK_THROW(
        ComaModelFileProcessor processor(nonExistentDir),
        std::runtime_error);

    // Test with empty directory (no SH files)
    boost::filesystem::path emptyDir = outputDir / "empty_dir";
    boost::filesystem::create_directories(emptyDir);
    BOOST_CHECK_THROW(
        ComaModelFileProcessor processor(emptyDir.string()),
        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(test_processor_file_type_behavior, TestDataPaths)
{
    // Create test SH files first
    const std::vector<std::string> files = {testFile.string()};
    ComaModelFileProcessor polyProcessor(files);
    std::vector<double> radii_m = {6000.0};
    std::vector<double> lons_deg = {0.0};

    boost::filesystem::path shTestDir = outputDir / "file_type_test";
    boost::filesystem::create_directories(shTestDir);
    polyProcessor.createSHFiles(shTestDir.string(), radii_m, lons_deg, 5, 5);

    // Test poly processor behavior
    ComaModelFileProcessor polyProc(files);
    BOOST_CHECK_EQUAL(polyProc.getFileType(), ComaModelFileProcessor::FileType::PolyCoefficients);
    BOOST_CHECK_NO_THROW(polyProc.createPolyCoefDataset());
    BOOST_CHECK_NO_THROW(polyProc.createSHDataset(radii_m, lons_deg));

    // Test SH processor behavior
    ComaModelFileProcessor shProc(shTestDir.string());
    BOOST_CHECK_EQUAL(shProc.getFileType(), ComaModelFileProcessor::FileType::StokesCoefficients);
    BOOST_CHECK_THROW(shProc.createPolyCoefDataset(), std::runtime_error);
    BOOST_CHECK_NO_THROW(shProc.createSHDataset({}, {})); // Parameters ignored

    // Verify that SH processor ignores createSHDataset parameters
    ComaStokesDataset shDataset1 = shProc.createSHDataset({1000.0}, {45.0});
    ComaStokesDataset shDataset2 = shProc.createSHDataset({2000.0, 3000.0}, {90.0, 180.0});

    // Both should return the same preloaded dataset
    BOOST_CHECK_EQUAL(shDataset1.nRadii(), shDataset2.nRadii());
    BOOST_CHECK_EQUAL(shDataset1.nLongitudes(), shDataset2.nLongitudes());
    BOOST_CHECK_EQUAL(shDataset1.nmax(), shDataset2.nmax());
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

    // Step 5: Test reading back the written data
    ComaStokesDataset readDataset = ComaStokesDatasetReader::readFromCsv(outputFile.string());
    BOOST_CHECK_EQUAL(readDataset.nmax(), stokesDataset.nmax());
    BOOST_CHECK_EQUAL(readDataset.nFiles(), stokesDataset.nFiles());
    BOOST_CHECK_EQUAL(readDataset.nRadii(), stokesDataset.nRadii());
    BOOST_CHECK_EQUAL(readDataset.nLongitudes(), stokesDataset.nLongitudes());

    // Verify some coefficient values match
    auto [orig_coeff, orig_sine] = stokesDataset.getCoeff(0, 0, 0, 0, 0);
    auto [read_coeff, read_sine] = readDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(read_coeff, orig_coeff, 1e-10);
    BOOST_CHECK_CLOSE(read_sine, orig_sine, 1e-10);

    // Step 6: Test reading from folder using new SH processor
    ComaModelFileProcessor integratedShProcessor(integratedOutput.string(), "integrated");
    ComaStokesDataset folderReadDataset = integratedShProcessor.createSHDataset({}, {});
    BOOST_CHECK_EQUAL(folderReadDataset.nmax(), stokesDataset.nmax());

    // Verify coefficient values from folder read match original
    auto [folder_coeff, folder_sine] = folderReadDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(folder_coeff, orig_coeff, 1e-10);
    BOOST_CHECK_CLOSE(folder_sine, orig_sine, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== SphericalHarmonicsCalculator Tests ====================

BOOST_AUTO_TEST_SUITE(test_spherical_harmonics_calculator)

BOOST_FIXTURE_TEST_CASE(test_calculate_surface_spherical_harmonics, TestDataPaths)
{
    // Step 1: Get polynomial coefficients from test data file
    const std::vector<std::string> files = {testFile.string()};
    const ComaPolyDataset dataset = ComaPolyDatasetReader::readFromFiles(files);

    // Step 2: Set up test parameters
    const double testRadius = 4000.0; //
    const double testSolarLongitude = 0.0 * M_PI / 180.0; // radians
    const int maxDegree = 10;
    const int maxOrder = 10;

    const double testLatitude = (90.0 - 11.5) * M_PI / 180.0; // converted from co-latitude to latitude
    const double testLongitude = -179.2 * M_PI / 180.0; // radians
    const double expectedResult = 1.54525E+16;

    const double testLatitude2 = (90.0 - 78.3) * M_PI / 180.0; // converted from co-latitude to latitude
    const double testLongitude2 = 110.9 * M_PI / 180.0; // radians
    const double expectedResult2 = 2.57521E+17;

    // Step 3: Get polynomial data from dataset
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients(0);
    const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices(0);
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius(0);
    const double refRadius = dataset.getReferenceRadius(0);

    // Step 4: Use evaluate2D to calculate Stokes coefficients from polynomial data
    const Eigen::ArrayXXd polyCoefficients = polyCoefs.array();
    Eigen::MatrixXd cosineCoefficients, sineCoefficients;

    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        testRadius,
        testSolarLongitude,
        polyCoefficients,
        shIndices,
        powers,
        refRadius,
        cosineCoefficients,
        sineCoefficients,
        maxDegree,
        maxOrder
    );

    // Step 5: Create SphericalHarmonicsCalculator instance
    SphericalHarmonicsCalculator calculator;

    // Step 6: Test calculateSurfaceSphericalHarmonics
    const double actualResult = calculator.calculateSurfaceSphericalHarmonics(
        sineCoefficients,
        cosineCoefficients,
        testLatitude,
        testLongitude,
        maxDegree,
        maxOrder
    );

    // Step 7: Validate result against expected value
    BOOST_CHECK_CLOSE( actualResult, std::log2(expectedResult), 1 );
    BOOST_CHECK_CLOSE( pow(2, actualResult), expectedResult, 10 );


    const double actualResult2 = calculator.calculateSurfaceSphericalHarmonics(
        sineCoefficients,
        cosineCoefficients,
        testLatitude2,
        testLongitude2,
        maxDegree,
        maxOrder
    );

    BOOST_CHECK_CLOSE( actualResult2, std::log2(expectedResult2), 1 );
    BOOST_CHECK_CLOSE( pow(2, actualResult2), expectedResult2, 10 );

}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat