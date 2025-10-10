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
    inline static constexpr double refRadius = 10000.0;
    inline static constexpr int maxDegree = 10;
    inline static constexpr int maxOrder = 10;

    // Sample poly coefficient values
    inline static constexpr double polyCoef_0_0 = 6.262302500423528E+02;
    inline static constexpr double polyCoef_3_1 = -4.459813047951577E-02;
    inline static constexpr double polyCoef_47_120 = -1.049515320967208E+02;
    inline static constexpr double polyCoef_10_22 = 1.287417812579956E-01;
};

// Helper function to read Stokes coefficient data from test files
struct StokesTestData
{
    double solarLongitude; // degrees
    double radius; // meters
    Eigen::MatrixXd cosineCoeffs; // (degree+1) x (order+1)
    Eigen::MatrixXd sineCoeffs;   // (degree+1) x (order+1)

    static StokesTestData readFromFile(const std::string& filepath, int maxDegree = 10)
    {
        StokesTestData data;
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open test data file: " + filepath);
        }

        std::string line;

        // Read first line: solar longitude and radius
        std::getline(file, line);
        std::istringstream iss(line);
        std::string dummy;
        iss >> dummy >> dummy >> dummy; // skip "# sol r"
        iss >> data.solarLongitude >> data.radius;
        data.radius *= 1000.0; // Convert km to meters

        // Read second line: header (GM R values)
        std::getline(file, line);

        // Initialize coefficient matrices
        data.cosineCoeffs = Eigen::MatrixXd::Zero(maxDegree + 1, maxDegree + 1);
        data.sineCoeffs = Eigen::MatrixXd::Zero(maxDegree + 1, maxDegree + 1);

        // Read coefficient data
        int degree, order;
        double cosine, sine;
        while (file >> degree >> order >> cosine >> sine)
        {
            if (degree <= maxDegree && order <= degree)
            {
                data.cosineCoeffs(degree, order) = cosine;
                data.sineCoeffs(degree, order) = sine;
            }
        }

        file.close();
        return data;
    }
};


// ==================== Data Model Tests ====================

BOOST_AUTO_TEST_SUITE(test_data_models)

BOOST_AUTO_TEST_CASE(test_stokes_dataset_creation)
{
    // Test pure data model creation
    const std::vector<ComaStokesDataset::FileMeta> files = {
        {2.015e9, 2.0150864e9, "test_file_1", 10000.0},
        {2.016e9, 2.0160864e9, "test_file_2", 10000.0}
    };
    const std::vector<double> radii = {1000.0, 2000.0, 3000.0};
    const std::vector<double> lons = {0.0, 30.0, 60.0, 90.0};
    constexpr int nmax = 10;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Verify metadata
    BOOST_CHECK_EQUAL(dataset.nFiles(), 2);
    BOOST_CHECK_EQUAL(dataset.nRadii(), 4); // 3 from radii vector + added ref radius
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
    const std::vector<ComaStokesDataset::FileMeta> files = {{0, 0, "test", 10000.0}};
    const std::vector<double> radii = {1000.0}; // + ref. Radius
    const std::vector<double> lons = {0.0};
    constexpr int nmax = 5;

    ComaStokesDataset dataset = ComaStokesDataset::create(files, radii, lons, nmax);

    // Test out of bounds access
    BOOST_CHECK_THROW(dataset.setCoeff(1, 0, 0, 0, 0, 0, 0), std::out_of_range); // file OOR
    BOOST_CHECK_THROW(dataset.setCoeff(0, 2, 0, 0, 0, 0, 0), std::out_of_range); // radius OOR
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
        {0.0, 1.0, "test_source", 10000}
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
        {2.015e9, 2.0150864e9, "test_source", 10000.0}
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
        {0.0, 1.0, "test_file_1", 10000.0}
    };
    std::vector<ComaStokesDataset::FileMeta> files2 = {
        {1.0, 2.0, "test_file_2", 10000.0}
    };

    std::vector<double> radii = {6000.0, 10000.0};
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

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Get the polynomial data from the dataset
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients(0);
    const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices(0);
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius(0);
    const double refRadius = dataset.getReferenceRadius(0);

    // Convert to ArrayXXd as required by the evaluator
    const Eigen::ArrayXXd polyCoefficients = polyCoefs.array();

    // Construct paths to test data files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";
    const boost::filesystem::path testFile1 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-000_04km.txt";
    const boost::filesystem::path testFile2 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-030_10km.txt";

    // Load expected values from both test files
    StokesTestData expectedData1 = StokesTestData::readFromFile(testFile1.string(), maxDegree);
    StokesTestData expectedData2 = StokesTestData::readFromFile(testFile2.string(), maxDegree);

    // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
    {
        const double radius = expectedData1.radius; // 4000 m
        const double solarLongitude = expectedData1.solarLongitude * M_PI / 180.0; // 0° in radians

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

        // Compare ALL coefficients against expected values from file
        int numCoeffsChecked = 0;
        for (int n = 0; n <= maxDegree; ++n)
        {
            for (int m = 0; m <= n; ++m)
            {
                // Check cosine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(cosineCoefficients(n, m) - expectedData1.cosineCoeffs(n, m)) /
                    std::max(std::abs(expectedData1.cosineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Cosine coefficient C(" << n << "," << m << ") mismatch: computed = "
                    << cosineCoefficients(n, m) << ", expected = " << expectedData1.cosineCoeffs(n, m)
                );

                // Check sine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(sineCoefficients(n, m) - expectedData1.sineCoeffs(n, m)) /
                    std::max(std::abs(expectedData1.sineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Sine coefficient S(" << n << "," << m << ") mismatch: computed = "
                    << sineCoefficients(n, m) << ", expected = " << expectedData1.sineCoeffs(n, m)
                );

                numCoeffsChecked++;
            }
        }
        // Verify we checked the expected number of coefficients: sum from m=0 to n for n=0 to maxDegree
        const int expectedNumCoeffs = (maxDegree + 1) * (maxDegree + 2) / 2;
        BOOST_CHECK_EQUAL(numCoeffsChecked, expectedNumCoeffs);
        BOOST_TEST_MESSAGE("Test Case 1: Verified all " << numCoeffsChecked << " cosine and "
                          << numCoeffsChecked << " sine coefficients (total: " << 2*numCoeffsChecked << ")");
    }

    // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
    {
        const double radius = expectedData2.radius; // 10000 m
        const double solarLongitude = expectedData2.solarLongitude * M_PI / 180.0; // 30° in radians

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

        // Compare ALL coefficients against expected values from file
        int numCoeffsChecked = 0;
        for (int n = 0; n <= maxDegree; ++n)
        {
            for (int m = 0; m <= n; ++m)
            {
                // Check cosine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(cosineCoefficients(n, m) - expectedData2.cosineCoeffs(n, m)) /
                    std::max(std::abs(expectedData2.cosineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Cosine coefficient C(" << n << "," << m << ") mismatch: computed = "
                    << cosineCoefficients(n, m) << ", expected = " << expectedData2.cosineCoeffs(n, m)
                );

                // Check sine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(sineCoefficients(n, m) - expectedData2.sineCoeffs(n, m)) /
                    std::max(std::abs(expectedData2.sineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Sine coefficient S(" << n << "," << m << ") mismatch: computed = "
                    << sineCoefficients(n, m) << ", expected = " << expectedData2.sineCoeffs(n, m)
                );

                numCoeffsChecked++;
            }
        }
        // Verify we checked the expected number of coefficients: sum from m=0 to n for n=0 to maxDegree
        const int expectedNumCoeffs = (maxDegree + 1) * (maxDegree + 2) / 2;
        BOOST_CHECK_EQUAL(numCoeffsChecked, expectedNumCoeffs);
        BOOST_TEST_MESSAGE("Test Case 2: Verified all " << numCoeffsChecked << " cosine and "
                          << numCoeffsChecked << " sine coefficients (total: " << 2*numCoeffsChecked << ")");
    }

    // ========== Test truncated degree/order ==========
    {
        Eigen::MatrixXd truncatedCosine, truncatedSine;
        const int truncatedMaxDegree = 5;
        const int truncatedMaxOrder = 3;
        const double radius = expectedData1.radius;
        const double solarLongitude = expectedData1.solarLongitude * M_PI / 180.0;

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

        // Verify that truncated results match the corresponding elements from expected data
        for (int n = 0; n <= truncatedMaxDegree; ++n)
        {
            for (int m = 0; m <= std::min(n, truncatedMaxOrder); ++m)
            {
                BOOST_CHECK_CLOSE(truncatedCosine(n, m), expectedData1.cosineCoeffs(n, m), 1e-12);
                BOOST_CHECK_CLOSE(truncatedSine(n, m), expectedData1.sineCoeffs(n, m), 1e-12);
            }
        }
    }

    // ========== Test error handling ==========
    {
        Eigen::MatrixXd dummyCosine, dummySine;
        const double radius = expectedData1.radius;
        const double solarLongitude = expectedData1.solarLongitude * M_PI / 180.0;

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
}

BOOST_FIXTURE_TEST_CASE(test_dataset_transformer, TestDataPaths)
{
    std::vector<std::string> files = {testFile.string()};
    ComaPolyDataset polyDataset = ComaPolyDatasetReader::readFromFiles(files);

    // Use the same test points as in test_stokes_coefficients_evaluator
    std::vector<double> radii_m = {4000.0, 10000.0};
    std::vector<double> lons_deg = {0.0, 30.0};

    // Load expected values from test files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";
    const boost::filesystem::path testFile1 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-000_04km.txt";
    const boost::filesystem::path testFile2 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-030_10km.txt";
    StokesTestData expectedData1 = StokesTestData::readFromFile(testFile1.string(), 10);
    StokesTestData expectedData2 = StokesTestData::readFromFile(testFile2.string(), 10);

    // Test transformation with default maxima
    ComaStokesDataset stokesDataset = ComaDatasetTransformer::transformPolyToStokes(
        polyDataset, radii_m, lons_deg);

    BOOST_CHECK_EQUAL(stokesDataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(stokesDataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nmax(), 10);

    // Check specific coefficient values at (ri=0, li=0) -> 4000m, 0deg
    auto [C_0_0_r0l0, S_0_0_r0l0] = stokesDataset.getCoeff(0, 0, 0, 0, 0);
    BOOST_CHECK_CLOSE(C_0_0_r0l0, expectedData1.cosineCoeffs(0, 0), 1e-10);
    BOOST_CHECK_CLOSE(S_0_0_r0l0, expectedData1.sineCoeffs(0, 0), 1e-10);

    auto [C_3_1_r0l0, S_3_1_r0l0] = stokesDataset.getCoeff(0, 0, 0, 3, 1);
    BOOST_CHECK_CLOSE(C_3_1_r0l0, expectedData1.cosineCoeffs(3, 1), 1e-10);

    auto [C_5_4_r0l0, S_5_4_r0l0] = stokesDataset.getCoeff(0, 0, 0, 5, 4);
    BOOST_CHECK_CLOSE(C_5_4_r0l0, expectedData1.cosineCoeffs(5, 4), 1e-10);

    // Check coefficient values at (ri=1, li=1) -> 10000m, 30deg
    auto [C_0_0_r1l1, S_0_0_r1l1] = stokesDataset.getCoeff(0, 1, 1, 0, 0);
    BOOST_CHECK_CLOSE(C_0_0_r1l1, expectedData2.cosineCoeffs(0, 0), 1e-10);
    BOOST_CHECK_CLOSE(S_0_0_r1l1, expectedData2.sineCoeffs(0, 0), 1e-10);

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

BOOST_FIXTURE_TEST_CASE(test_stokes_dataset_creation_via_processor, TestDataPaths)
{
    // This test is similar to test_stokes_coefficients_evaluator, but instead of directly
    // using the evaluate function, it uses the ComaModelFileProcessor class to create
    // an SH dataset, then extracts and compares the Stokes coefficients.

    const std::vector<std::string> files = {testFile.string()};
    const ComaModelFileProcessor processor(files);

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Create radii and longitudes vectors - need at least 2 of each
    // Include the test points: (0°, 4km) and (30°, 10km)
    const std::vector<double> radii_m = {4000.0, 10000.0};  // 4 km, 10 km
    const std::vector<double> lons_deg = {0.0, 30.0};       // 0°, 30°

    // Create SH dataset using the processor
    const ComaStokesDataset stokesDataset = processor.createSHDataset(radii_m, lons_deg, maxDegree, maxOrder);

    // Verify dataset structure
    BOOST_CHECK_EQUAL(stokesDataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(stokesDataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(stokesDataset.nmax(), maxDegree);

    // Construct paths to test data files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";
    const boost::filesystem::path testFile1 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-000_04km.txt";
    const boost::filesystem::path testFile2 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-030_10km.txt";

    // Load expected values from both test files
    StokesTestData expectedData1 = StokesTestData::readFromFile(testFile1.string(), maxDegree);
    StokesTestData expectedData2 = StokesTestData::readFromFile(testFile2.string(), maxDegree);

    // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
    // This corresponds to dataset indices: radius_index=0, longitude_index=0
    {
        // Get coefficient matrices from the dataset
        auto [cosineCoefficients, sineCoefficients] = stokesDataset.getCoefficientMatrices(0, 0, 0);

        // Check output dimensions
        BOOST_CHECK_EQUAL(cosineCoefficients.rows(), maxDegree + 1);
        BOOST_CHECK_EQUAL(cosineCoefficients.cols(), maxOrder + 1);
        BOOST_CHECK_EQUAL(sineCoefficients.rows(), maxDegree + 1);
        BOOST_CHECK_EQUAL(sineCoefficients.cols(), maxOrder + 1);

        // Compare ALL coefficients against expected values from file
        int numCoeffsChecked = 0;
        for (int n = 0; n <= maxDegree; ++n)
        {
            for (int m = 0; m <= n; ++m)
            {
                // Check cosine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(cosineCoefficients(n, m) - expectedData1.cosineCoeffs(n, m)) /
                    std::max(std::abs(expectedData1.cosineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Cosine coefficient C(" << n << "," << m << ") mismatch: computed = "
                    << cosineCoefficients(n, m) << ", expected = " << expectedData1.cosineCoeffs(n, m)
                );

                // Check sine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(sineCoefficients(n, m) - expectedData1.sineCoeffs(n, m)) /
                    std::max(std::abs(expectedData1.sineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Sine coefficient S(" << n << "," << m << ") mismatch: computed = "
                    << sineCoefficients(n, m) << ", expected = " << expectedData1.sineCoeffs(n, m)
                );

                numCoeffsChecked++;
            }
        }
        // Verify we checked the expected number of coefficients: sum from m=0 to n for n=0 to maxDegree
        const int expectedNumCoeffs = (maxDegree + 1) * (maxDegree + 2) / 2;
        BOOST_CHECK_EQUAL(numCoeffsChecked, expectedNumCoeffs);
        BOOST_TEST_MESSAGE("Test Case 1 (via processor): Verified all " << numCoeffsChecked << " cosine and "
                          << numCoeffsChecked << " sine coefficients (total: " << 2*numCoeffsChecked << ")");
    }

    // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
    // This corresponds to dataset indices: radius_index=1, longitude_index=1
    {
        // Get coefficient matrices from the dataset
        auto [cosineCoefficients, sineCoefficients] = stokesDataset.getCoefficientMatrices(0, 1, 1);

        // Check output dimensions
        BOOST_CHECK_EQUAL(cosineCoefficients.rows(), maxDegree + 1);
        BOOST_CHECK_EQUAL(cosineCoefficients.cols(), maxOrder + 1);
        BOOST_CHECK_EQUAL(sineCoefficients.rows(), maxDegree + 1);
        BOOST_CHECK_EQUAL(sineCoefficients.cols(), maxOrder + 1);

        // Compare ALL coefficients against expected values from file
        int numCoeffsChecked = 0;
        for (int n = 0; n <= maxDegree; ++n)
        {
            for (int m = 0; m <= n; ++m)
            {
                // Check cosine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(cosineCoefficients(n, m) - expectedData2.cosineCoeffs(n, m)) /
                    std::max(std::abs(expectedData2.cosineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Cosine coefficient C(" << n << "," << m << ") mismatch: computed = "
                    << cosineCoefficients(n, m) << ", expected = " << expectedData2.cosineCoeffs(n, m)
                );

                // Check sine coefficient
                BOOST_CHECK_MESSAGE(
                    std::abs(sineCoefficients(n, m) - expectedData2.sineCoeffs(n, m)) /
                    std::max(std::abs(expectedData2.sineCoeffs(n, m)), 1e-10) < 1e-10,
                    "Sine coefficient S(" << n << "," << m << ") mismatch: computed = "
                    << sineCoefficients(n, m) << ", expected = " << expectedData2.sineCoeffs(n, m)
                );

                numCoeffsChecked++;
            }
        }
        // Verify we checked the expected number of coefficients: sum from m=0 to n for n=0 to maxDegree
        const int expectedNumCoeffs = (maxDegree + 1) * (maxDegree + 2) / 2;
        BOOST_CHECK_EQUAL(numCoeffsChecked, expectedNumCoeffs);
        BOOST_TEST_MESSAGE("Test Case 2 (via processor): Verified all " << numCoeffsChecked << " cosine and "
                          << numCoeffsChecked << " sine coefficients (total: " << 2*numCoeffsChecked << ")");
    }
}

BOOST_AUTO_TEST_SUITE_END()

// ==================== ComaModel Tests ====================

BOOST_AUTO_TEST_SUITE(test_coma_model)

BOOST_FIXTURE_TEST_CASE(test_coma_model_number_density, TestDataPaths)
{
    // Load polynomial coefficients from test data file
    const std::vector<std::string> files = {testFile.string()};
    const ComaPolyDataset polyDataset = ComaPolyDatasetReader::readFromFiles(files);

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Molecular weight (arbitrary for number density test, but needed for constructor)
    const double molecularWeight = 0.018;  // kg/mol (water vapor)

    // Create Stokes dataset from polynomial dataset
    // Use the test points from both test files as the grid
    const std::vector<double> radii_m = {4000.0, 10000.0};  // 4 km, 10 km
    const std::vector<double> lons_deg = {0.0, 30.0};       // 0°, 30°
    const ComaModelFileProcessor processor(files);
    const ComaStokesDataset stokesDataset = processor.createSHDataset(radii_m, lons_deg, maxDegree, maxOrder);

    // Construct paths to residuals test data files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";

    // Test with both polynomial and Stokes datasets
    for (int datasetType = 0; datasetType < 2; ++datasetType)
    {
        const std::string datasetTypeName = (datasetType == 0) ? "Polynomial" : "Stokes";
        BOOST_TEST_MESSAGE("========== Testing with " << datasetTypeName << " dataset ==========");

        // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
        {
            const boost::filesystem::path residualsFile1 = dataDir / "residual_r_cometFixed_ep10-000_04km.txt";
            const double testRadius = 4000.0;  // 4 km in meters
            const double testSolarLongitude = 0.0;  // 0 degrees in radians
            const double testTime = 490708800;   // s since J2000

            // Define state functions for solar longitude = 0°
            // Solar longitude is calculated as atan2(y, x) of the Sun direction in body-fixed frame
            // For solar longitude = 0°, we need Sun along the +X axis in body-fixed frame

            // Sun state function - returns position of Sun along +X axis
            auto sunStateFunction = []() -> Eigen::Vector6d {
                Eigen::Vector6d state = Eigen::Vector6d::Zero();
                state.segment(0, 3) = Eigen::Vector3d(1.0e11, 0.0, 0.0);  // Sun along +X axis
                return state;
            };

            // Comet state function - returns position of comet at origin
            auto cometStateFunction = []() -> Eigen::Vector6d {
                return Eigen::Vector6d::Zero();
            };

            // Comet rotation function - returns identity matrix
            auto cometRotationFunction = []() -> Eigen::Matrix3d {
                return Eigen::Matrix3d::Identity();
            };

            // Create ComaModel based on dataset type
            std::unique_ptr<ComaModel> comaModel;
            if (datasetType == 0)
            {
                // Create with polynomial coefficients
                comaModel = std::make_unique<ComaModel>(
                    polyDataset,
                    molecularWeight,
                    sunStateFunction,
                    cometStateFunction,
                    cometRotationFunction,
                    maxDegree,
                    maxOrder
                );
            }
            else
            {
                // Create with Stokes coefficients
                comaModel = std::make_unique<ComaModel>(
                    stokesDataset,
                    molecularWeight,
                    sunStateFunction,
                    cometStateFunction,
                    cometRotationFunction,
                    maxDegree,
                    maxOrder
                );
            }

        // Read all data points from residuals file
        std::ifstream file1(residualsFile1.string());
        BOOST_REQUIRE_MESSAGE(file1.is_open(), "Cannot open residuals file: " + residualsFile1.string());

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector<DataPoint> allPoints;

        std::string line;
        while (std::getline(file1, line))
        {
            std::istringstream iss(line);
            DataPoint point;
            if (iss >> point.longitude_deg >> point.latitude_deg >> point.originalData
                    >> point.shEvaluation >> point.difference)
            {
                allPoints.push_back(point);
            }
        }
        file1.close();

        // Randomly select 100 points for testing
        const int numTestPoints = 100;
        std::vector<int> selectedIndices;
        std::srand(12345);  // Fixed seed for reproducibility

        if (allPoints.size() <= numTestPoints)
        {
            for (size_t i = 0; i < allPoints.size(); ++i)
                selectedIndices.push_back(i);
        }
        else
        {
            std::vector<int> allIndices;
            for (size_t i = 0; i < allPoints.size(); ++i)
                allIndices.push_back(i);

            for (int i = 0; i < numTestPoints; ++i)
            {
                int randomIndex = std::rand() % allIndices.size();
                selectedIndices.push_back(allIndices[randomIndex]);
                allIndices.erase(allIndices.begin() + randomIndex);
            }
        }

        // Test the selected points
        int failCount = 0;
        for (int idx : selectedIndices)
        {
            const DataPoint& point = allPoints[idx];

            // Convert to radians
            const double longitude_rad = point.longitude_deg * M_PI / 180.0;
            const double latitude_rad = point.latitude_deg * M_PI / 180.0;

            // Call getNumberDensity from ComaModel
            const double computedNumberDensity = comaModel->getNumberDensity(
                testRadius,
                longitude_rad,
                latitude_rad,
                testTime
            );

            // The expected value is in the 4th column (shEvaluation)
            // Note: The test file contains log2(number_density), so we need to compare directly
            const double expectedNumberDensity = point.shEvaluation;

            const double tolerance = 1e-10;
            const double diff = std::abs(computedNumberDensity - expectedNumberDensity);

            if (diff > tolerance)
            {
                failCount++;
                if (failCount <= 10)  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE("Case 1 - Point " << idx << " mismatch: "
                                      << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                      << "computed=" << computedNumberDensity << ", expected=" << expectedNumberDensity
                                      << ", diff=" << diff);
                }
            }

            BOOST_CHECK_MESSAGE(
                diff <= tolerance,
                "Number density mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°"
            );
        }

            BOOST_TEST_MESSAGE("Test Case 1 (" << datasetTypeName << ", 0°, 4km): Verified " << selectedIndices.size()
                              << " randomly selected points out of " << allPoints.size()
                              << " total, " << failCount << " failures");
            BOOST_CHECK_EQUAL(failCount, 0);
        }

        // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
        {
            const boost::filesystem::path residualsFile2 = dataDir / "residual_r_cometFixed_ep10-030_10km.txt";
            const double testRadius = 10000.0;  // 10 km in meters
            const double testSolarLongitude = 30.0 * M_PI / 180.0;  // 30 degrees in radians
            const double testTime = 490708800;   // s since J2000

            // Define state functions for solar longitude = 30°
            // Solar longitude is calculated as atan2(y, x) of the Sun direction in body-fixed frame
            // For solar longitude = 30°, we need Sun at angle 30° from +X axis in XY plane
            // This means: x = cos(30°), y = sin(30°), z = 0

            const double solarDistance = 1.0e11;  // Distance to Sun in meters
            const double sunX = solarDistance * std::cos(testSolarLongitude);
            const double sunY = solarDistance * std::sin(testSolarLongitude);

            // Sun state function - returns position of Sun at 30° angle
            auto sunStateFunction = [sunX, sunY]() -> Eigen::Vector6d {
                Eigen::Vector6d state = Eigen::Vector6d::Zero();
                state.segment(0, 3) = Eigen::Vector3d(sunX, sunY, 0.0);
                return state;
            };

            // Comet state function - returns position of comet at origin
            auto cometStateFunction = []() -> Eigen::Vector6d {
                return Eigen::Vector6d::Zero();
            };

            // Comet rotation function - returns identity matrix
            auto cometRotationFunction = []() -> Eigen::Matrix3d {
                return Eigen::Matrix3d::Identity();
            };

            // Create ComaModel based on dataset type
            std::unique_ptr<ComaModel> comaModel;
            if (datasetType == 0)
            {
                // Create with polynomial coefficients
                comaModel = std::make_unique<ComaModel>(
                    polyDataset,
                    molecularWeight,
                    sunStateFunction,
                    cometStateFunction,
                    cometRotationFunction,
                    maxDegree,
                    maxOrder
                );
            }
            else
            {
                // Create with Stokes coefficients
                comaModel = std::make_unique<ComaModel>(
                    stokesDataset,
                    molecularWeight,
                    sunStateFunction,
                    cometStateFunction,
                    cometRotationFunction,
                    maxDegree,
                    maxOrder
                );
            }

        // Read all data points from residuals file
        std::ifstream file2(residualsFile2.string());
        BOOST_REQUIRE_MESSAGE(file2.is_open(), "Cannot open residuals file: " + residualsFile2.string());

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector<DataPoint> allPoints;

        std::string line;
        while (std::getline(file2, line))
        {
            std::istringstream iss(line);
            DataPoint point;
            if (iss >> point.longitude_deg >> point.latitude_deg >> point.originalData
                    >> point.shEvaluation >> point.difference)
            {
                allPoints.push_back(point);
            }
        }
        file2.close();

        // Randomly select 50 points for testing
        const int numTestPoints = 50;
        std::vector<int> selectedIndices;
        std::srand(67890);  // Different seed from case 1

        if (allPoints.size() <= numTestPoints)
        {
            for (size_t i = 0; i < allPoints.size(); ++i)
                selectedIndices.push_back(i);
        }
        else
        {
            std::vector<int> allIndices;
            for (size_t i = 0; i < allPoints.size(); ++i)
                allIndices.push_back(i);

            for (int i = 0; i < numTestPoints; ++i)
            {
                int randomIndex = std::rand() % allIndices.size();
                selectedIndices.push_back(allIndices[randomIndex]);
                allIndices.erase(allIndices.begin() + randomIndex);
            }
        }

        // Test the selected points
        int failCount = 0;
        for (int idx : selectedIndices)
        {
            const DataPoint& point = allPoints[idx];

            // Convert to radians
            const double longitude_rad = point.longitude_deg * M_PI / 180.0;
            const double latitude_rad = point.latitude_deg * M_PI / 180.0;

            // Call getNumberDensity from ComaModel
            const double computedNumberDensity = comaModel->getNumberDensity(
                testRadius,
                longitude_rad,
                latitude_rad,
                testTime
            );

            // The expected value is in the 4th column (shEvaluation)
            const double expectedNumberDensity = point.shEvaluation;

            const double tolerance = 1e-10;
            const double diff = std::abs(computedNumberDensity - expectedNumberDensity);

            if (diff > tolerance)
            {
                failCount++;
                if (failCount <= 10)  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE("Case 2 - Point " << idx << " mismatch: "
                                      << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                      << "computed=" << computedNumberDensity << ", expected=" << expectedNumberDensity
                                      << ", diff=" << diff);
                }
            }

            BOOST_CHECK_MESSAGE(
                diff <= tolerance,
                "Number density mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°"
            );
        }

            BOOST_TEST_MESSAGE("Test Case 2 (" << datasetTypeName << ", 30°, 10km): Verified " << selectedIndices.size()
                              << " randomly selected points out of " << allPoints.size()
                              << " total, " << failCount << " failures");
            BOOST_CHECK_EQUAL(failCount, 0);
        }
    }
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

    // Use the same test points as in test_stokes_coefficients_evaluator
    const std::vector<double> radii_m = {4000.0, 10000.0};
    const std::vector<double> lons_deg = {0.0, 30.0};

    // Load expected values from test files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";
    const boost::filesystem::path testFile1 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-000_04km.txt";
    const boost::filesystem::path testFile2 = dataDir / "SH-d10-rp3-fft12__r_cometFixed_ep10-030_10km.txt";
    StokesTestData expectedData1 = StokesTestData::readFromFile(testFile1.string(), 10);
    StokesTestData expectedData2 = StokesTestData::readFromFile(testFile2.string(), 10);

    // Test with auto-selected maxima
    const ComaStokesDataset dataset = processor.createSHDataset(radii_m, lons_deg);

    BOOST_CHECK_EQUAL(dataset.nFiles(), 1);
    BOOST_CHECK_EQUAL(dataset.nRadii(), 2);
    BOOST_CHECK_EQUAL(dataset.nLongitudes(), 2);
    BOOST_CHECK_EQUAL(dataset.nmax(), 10);

    // Verify computed coefficients match expected values at (ri=0, li=0) -> 4000m, 0deg
    auto [cosineMatrix0, sineMatrix0] = dataset.getCoefficientMatrices(0, 0, 0);
    BOOST_CHECK_CLOSE(cosineMatrix0(0, 0), expectedData1.cosineCoeffs(0, 0), 1e-10);
    BOOST_CHECK_CLOSE(cosineMatrix0(3, 1), expectedData1.cosineCoeffs(3, 1), 1e-10);
    BOOST_CHECK_CLOSE(sineMatrix0(6, 2), expectedData1.sineCoeffs(6, 2), 1e-10);
    BOOST_CHECK_CLOSE(sineMatrix0(7, 5), expectedData1.sineCoeffs(7, 5), 1e-10);

    // Verify computed coefficients at (ri=1, li=1) -> 10000m, 30deg
    auto [cosineMatrix1, sineMatrix1] = dataset.getCoefficientMatrices(0, 1, 1);
    BOOST_CHECK_CLOSE(cosineMatrix1(0, 0), expectedData2.cosineCoeffs(0, 0), 1e-10);
    BOOST_CHECK_CLOSE(cosineMatrix1(3, 1), expectedData2.cosineCoeffs(3, 1), 1e-10);
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
    std::vector<double> radii_m = {6000.0, 10000.0};  // Need at least 2 radii for interpolation
    std::vector<double> lons_deg = {0.0, 30.0};       // Need at least 2 longitudes for interpolation

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

    // Step 2: Get polynomial data from dataset
    const Eigen::MatrixXd& polyCoefs = dataset.getPolyCoefficients(0);
    const Eigen::ArrayXXi& shIndices = dataset.getSHDegreeAndOrderIndices(0);
    const Eigen::VectorXd& powers = dataset.getPowersInvRadius(0);
    const double refRadius = dataset.getReferenceRadius(0);
    const Eigen::ArrayXXd polyCoefficients = polyCoefs.array();

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Create SphericalHarmonicsCalculator instance
    SphericalHarmonicsCalculator calculator;

    // Construct paths to residuals test data files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";

    // ========== Test Case 1: solar longitude = 0°, radius = 4 km ==========
    {
        const boost::filesystem::path residualsFile1 = dataDir / "residual_r_cometFixed_ep10-000_04km.txt";
        const double testRadius = 4000.0;  // 4 km in meters
        const double testSolarLongitude = 0.0 * M_PI / 180.0;  // 0 degrees in radians

        // Calculate Stokes coefficients for this case
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

        // Read all data points from file
        std::ifstream file1(residualsFile1.string());
        BOOST_REQUIRE_MESSAGE(file1.is_open(), "Cannot open residuals file: " + residualsFile1.string());

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector<DataPoint> allPoints;

        std::string line;
        while (std::getline(file1, line))
        {
            std::istringstream iss(line);
            DataPoint point;
            if (iss >> point.longitude_deg >> point.latitude_deg >> point.originalData
                    >> point.shEvaluation >> point.difference)
            {
                allPoints.push_back(point);
            }
        }
        file1.close();

        // Randomly select 100 points
        const int numTestPoints = 100;
        std::vector<int> selectedIndices;
        std::srand(12345);  // Fixed seed for reproducibility

        if (allPoints.size() <= numTestPoints)
        {
            // Use all points if we have fewer than requested
            for (size_t i = 0; i < allPoints.size(); ++i)
                selectedIndices.push_back(i);
        }
        else
        {
            // Randomly select numTestPoints indices without replacement
            std::vector<int> allIndices;
            for (size_t i = 0; i < allPoints.size(); ++i)
                allIndices.push_back(i);

            for (int i = 0; i < numTestPoints; ++i)
            {
                int randomIndex = std::rand() % allIndices.size();
                selectedIndices.push_back(allIndices[randomIndex]);
                allIndices.erase(allIndices.begin() + randomIndex);
            }
        }

        // Test the selected points
        int failCount = 0;
        for (int idx : selectedIndices)
        {
            const DataPoint& point = allPoints[idx];

            // Convert to radians (latitude is already geodetic latitude in the file)
            const double longitude_rad = point.longitude_deg * M_PI / 180.0;
            const double latitude_rad = point.latitude_deg * M_PI / 180.0;

            // Calculate using our implementation
            const double computedResult = calculator.calculateSurfaceSphericalHarmonics(
                sineCoefficients,
                cosineCoefficients,
                latitude_rad,
                longitude_rad,
                maxDegree,
                maxOrder
            );

            // Compare with SH evaluation from file (column 4)
            const double tolerance = 1e-10;
            const double diff = std::abs(computedResult - point.shEvaluation);

            if (diff > tolerance)
            {
                failCount++;
                if (failCount <= 10)  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE("Case 1 - Point " << idx << " mismatch: "
                                      << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                      << "computed=" << computedResult << ", expected=" << point.shEvaluation
                                      << ", diff=" << diff);
                }
            }

            BOOST_CHECK_MESSAGE(
                diff <= tolerance,
                "SH evaluation mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°"
            );
        }

        BOOST_TEST_MESSAGE("Test Case 1 (0°, 4km): Verified " << selectedIndices.size()
                          << " randomly selected points out of " << allPoints.size()
                          << " total, " << failCount << " failures");
        BOOST_CHECK_EQUAL(failCount, 0);
    }

    // ========== Test Case 2: solar longitude = 30°, radius = 10 km ==========
    {
        const boost::filesystem::path residualsFile2 = dataDir / "residual_r_cometFixed_ep10-030_10km.txt";
        const double testRadius = 10000.0;  // 10 km in meters
        const double testSolarLongitude = 30.0 * M_PI / 180.0;  // 30 degrees in radians

        // Calculate Stokes coefficients for this case
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

        // Read all data points from file
        std::ifstream file2(residualsFile2.string());
        BOOST_REQUIRE_MESSAGE(file2.is_open(), "Cannot open residuals file: " + residualsFile2.string());

        struct DataPoint {
            double longitude_deg, latitude_deg, originalData, shEvaluation, difference;
        };
        std::vector<DataPoint> allPoints;

        std::string line;
        while (std::getline(file2, line))
        {
            std::istringstream iss(line);
            DataPoint point;
            if (iss >> point.longitude_deg >> point.latitude_deg >> point.originalData
                    >> point.shEvaluation >> point.difference)
            {
                allPoints.push_back(point);
            }
        }
        file2.close();

        // Randomly select 100 points
        const int numTestPoints = 100;
        std::vector<int> selectedIndices;
        std::srand(67890);  // Fixed seed for reproducibility (different from case 1)

        if (allPoints.size() <= numTestPoints)
        {
            // Use all points if we have fewer than requested
            for (size_t i = 0; i < allPoints.size(); ++i)
                selectedIndices.push_back(i);
        }
        else
        {
            // Randomly select numTestPoints indices without replacement
            std::vector<int> allIndices;
            for (size_t i = 0; i < allPoints.size(); ++i)
                allIndices.push_back(i);

            for (int i = 0; i < numTestPoints; ++i)
            {
                int randomIndex = std::rand() % allIndices.size();
                selectedIndices.push_back(allIndices[randomIndex]);
                allIndices.erase(allIndices.begin() + randomIndex);
            }
        }

        // Test the selected points
        int failCount = 0;
        for (int idx : selectedIndices)
        {
            const DataPoint& point = allPoints[idx];

            // Convert to radians (latitude is already geodetic latitude in the file)
            const double longitude_rad = point.longitude_deg * M_PI / 180.0;
            const double latitude_rad = point.latitude_deg * M_PI / 180.0;

            // Calculate using our implementation
            const double computedResult = calculator.calculateSurfaceSphericalHarmonics(
                sineCoefficients,
                cosineCoefficients,
                latitude_rad,
                longitude_rad,
                maxDegree,
                maxOrder
            );

            // Compare with SH evaluation from file (column 4)
            const double tolerance = 1e-10;
            const double diff = std::abs(computedResult - point.shEvaluation);

            if (diff > tolerance)
            {
                failCount++;
                if (failCount <= 10)  // Only report first 10 failures
                {
                    BOOST_TEST_MESSAGE("Case 2 - Point " << idx << " mismatch: "
                                      << "lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°, "
                                      << "computed=" << computedResult << ", expected=" << point.shEvaluation
                                      << ", diff=" << diff);
                }
            }

            BOOST_CHECK_MESSAGE(
                diff <= tolerance,
                "SH evaluation mismatch at lon=" << point.longitude_deg << "°, lat=" << point.latitude_deg << "°"
            );
        }

        BOOST_TEST_MESSAGE("Test Case 2 (30°, 10km): Verified " << selectedIndices.size()
                          << " randomly selected points out of " << allPoints.size()
                          << " total, " << failCount << " failures");
        BOOST_CHECK_EQUAL(failCount, 0);
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat