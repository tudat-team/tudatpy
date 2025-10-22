#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

namespace tudat
{
namespace unit_tests
{
using namespace simulation_setup;
using namespace aerodynamics;

// ==================== Test Fixtures ====================

struct WindTestDataPaths
{
    boost::filesystem::path windFileX;
    boost::filesystem::path windFileY;
    boost::filesystem::path windFileZ;
    boost::filesystem::path outputDir;

    WindTestDataPaths()
    {
        const boost::filesystem::path thisFile(__FILE__);
        const boost::filesystem::path testDir = thisFile.parent_path();
        const boost::filesystem::path dataDir = testDir / "test_data";

        windFileX = dataDir / "input_poly_wind_x.txt";
        windFileY = dataDir / "input_poly_wind_y.txt";
        windFileZ = dataDir / "input_poly_wind_z.txt";
        outputDir = dataDir / "test_output_wind";

        // Ensure test files exist
        if (!boost::filesystem::exists(windFileX))
        {
            throw std::runtime_error("Wind test data file not found: " + windFileX.string());
        }
        if (!boost::filesystem::exists(windFileY))
        {
            throw std::runtime_error("Wind test data file not found: " + windFileY.string());
        }
        if (!boost::filesystem::exists(windFileZ))
        {
            throw std::runtime_error("Wind test data file not found: " + windFileZ.string());
        }

        // Clean and create output directory
        if (boost::filesystem::exists(outputDir))
        {
            boost::filesystem::remove_all(outputDir);
        }
        boost::filesystem::create_directories(outputDir);
    }

    ~WindTestDataPaths()
    {
        // Optional: cleanup output dir after tests
        // boost::filesystem::remove_all(outputDir);
    }
};

// Expected values from wind input files
struct ExpectedWindValues
{
    static constexpr int maxDegree = 10;
    static constexpr int maxOrder = 10;
    static constexpr double refRadius = 10.0; // in meters (converted from km in file)
    static constexpr int numTimeTerms = 12;
    static constexpr int numRadiiTerms = 2;
};

// C++14 requires definitions for static constexpr members that are ODR-used
constexpr int ExpectedWindValues::maxDegree;
constexpr int ExpectedWindValues::maxOrder;
constexpr double ExpectedWindValues::refRadius;
constexpr int ExpectedWindValues::numTimeTerms;
constexpr int ExpectedWindValues::numRadiiTerms;

// ==================== ComaWindModelFileProcessor Tests ====================

BOOST_AUTO_TEST_SUITE(test_coma_wind_model_file_processor)

BOOST_FIXTURE_TEST_CASE(test_wind_processor_construction, WindTestDataPaths)
{
    // Test construction with polynomial coefficient files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    // Should construct without throwing
    BOOST_CHECK_NO_THROW(
        ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles)
    );

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Verify it's recognized as polynomial type
    BOOST_CHECK(processor.isPolyType());
    BOOST_CHECK(!processor.isStokesType());
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_factory_function, WindTestDataPaths)
{
    // Test the factory function for polynomial files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    auto processor = comaWindModelFileProcessorFromPolyFiles(xFiles, yFiles, zFiles);

    BOOST_CHECK(processor != nullptr);
    BOOST_CHECK(processor->isPolyType());
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_create_poly_datasets, WindTestDataPaths)
{
    // Test creating polynomial coefficient datasets
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Create polynomial datasets
    ComaWindDatasetCollection windDatasets = processor.createPolyCoefDatasets();

    // Verify the datasets were created
    BOOST_CHECK(windDatasets.isPoly());
    BOOST_CHECK(!windDatasets.isStokes());

    // Verify all three components are present
    const ComaPolyDataset& xData = windDatasets.getXPolyDataset();
    const ComaPolyDataset& yData = windDatasets.getYPolyDataset();
    const ComaPolyDataset& zData = windDatasets.getZPolyDataset();

    // Check that each component has data from one file
    BOOST_CHECK_EQUAL(xData.getNumFiles(), 1);
    BOOST_CHECK_EQUAL(yData.getNumFiles(), 1);
    BOOST_CHECK_EQUAL(zData.getNumFiles(), 1);

    // Verify reference radius (should be 10 km = 10000 m)
    BOOST_CHECK_CLOSE(xData.getReferenceRadius(0), 10000.0, 1e-10);
    BOOST_CHECK_CLOSE(yData.getReferenceRadius(0), 10000.0, 1e-10);
    BOOST_CHECK_CLOSE(zData.getReferenceRadius(0), 10000.0, 1e-10);

    // Verify maximum degree
    BOOST_CHECK_EQUAL(xData.getMaxDegreeSH(0), ExpectedWindValues::maxDegree);
    BOOST_CHECK_EQUAL(yData.getMaxDegreeSH(0), ExpectedWindValues::maxDegree);
    BOOST_CHECK_EQUAL(zData.getMaxDegreeSH(0), ExpectedWindValues::maxDegree);
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_create_sh_datasets_from_poly, WindTestDataPaths)
{
    // Test creating Stokes coefficient datasets from polynomial files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Define radii and solar longitudes for evaluation
    const std::vector<double> radii_m = {4000.0, 9000.0};  // 4 km, 10 km
    const std::vector<double> solLongitudes_deg = {0.0, 30.0};  // 0°, 30°

    // Create Stokes datasets
    ComaWindDatasetCollection windDatasets = processor.createSHDatasets(
        radii_m, solLongitudes_deg,
        ExpectedWindValues::maxDegree,
        ExpectedWindValues::maxOrder
    );

    // Verify the datasets were created
    BOOST_CHECK(!windDatasets.isPoly());
    BOOST_CHECK(windDatasets.isStokes());

    // Get the three components
    const ComaStokesDataset& xStokes = windDatasets.getXStokesDataset();
    const ComaStokesDataset& yStokes = windDatasets.getYStokesDataset();
    const ComaStokesDataset& zStokes = windDatasets.getZStokesDataset();

    // Verify structure of each component
    BOOST_CHECK_EQUAL(xStokes.nFiles(), 1);
    BOOST_CHECK_EQUAL(yStokes.nFiles(), 1);
    BOOST_CHECK_EQUAL(zStokes.nFiles(), 1);

    // Verify radii (should have reference radius + requested radii)
    BOOST_CHECK_EQUAL(xStokes.nRadii(), radii_m.size() + 1); // +1 for reference radius
    BOOST_CHECK_EQUAL(yStokes.nRadii(), radii_m.size() + 1);
    BOOST_CHECK_EQUAL(zStokes.nRadii(), radii_m.size() + 1);

    // Verify longitudes
    BOOST_CHECK_EQUAL(xStokes.nLongitudes(), solLongitudes_deg.size());
    BOOST_CHECK_EQUAL(yStokes.nLongitudes(), solLongitudes_deg.size());
    BOOST_CHECK_EQUAL(zStokes.nLongitudes(), solLongitudes_deg.size());

    // Verify max degree
    BOOST_CHECK_EQUAL(xStokes.nmax(), ExpectedWindValues::maxDegree);
    BOOST_CHECK_EQUAL(yStokes.nmax(), ExpectedWindValues::maxDegree);
    BOOST_CHECK_EQUAL(zStokes.nmax(), ExpectedWindValues::maxDegree);
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_create_sh_files, WindTestDataPaths)
{
    // Test creating Stokes coefficient CSV files from polynomial files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Define radii and solar longitudes
    const std::vector<double> radii_m = {6000.0, 10000.0};
    const std::vector<double> solLongitudes_deg = {0.0, 30.0};

    // Create output directories for each component
    boost::filesystem::path xOutputDir = outputDir / "wind_x_stokes";
    boost::filesystem::path yOutputDir = outputDir / "wind_y_stokes";
    boost::filesystem::path zOutputDir = outputDir / "wind_z_stokes";

    // Create SH files
    processor.createSHFiles(
        xOutputDir.string(),
        yOutputDir.string(),
        zOutputDir.string(),
        radii_m,
        solLongitudes_deg,
        8, 6  // max degree and order
    );

    // Verify directories were created
    BOOST_CHECK(boost::filesystem::exists(xOutputDir));
    BOOST_CHECK(boost::filesystem::exists(yOutputDir));
    BOOST_CHECK(boost::filesystem::exists(zOutputDir));

    // Verify CSV files were created (default prefix is "stokes")
    boost::filesystem::path xCsvFile = xOutputDir / "stokes_file0.csv";
    boost::filesystem::path yCsvFile = yOutputDir / "stokes_file0.csv";
    boost::filesystem::path zCsvFile = zOutputDir / "stokes_file0.csv";

    BOOST_CHECK(boost::filesystem::exists(xCsvFile));
    BOOST_CHECK(boost::filesystem::exists(yCsvFile));
    BOOST_CHECK(boost::filesystem::exists(zCsvFile));
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_roundtrip_poly_to_stokes, WindTestDataPaths)
{
    // Test the complete roundtrip: poly files -> SH datasets -> verify consistency
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    const std::vector<double> radii_m = {6000.0, 10000.0};
    const std::vector<double> solLongitudes_deg = {0.0, 45.0, 90.0};

    // Create Stokes datasets
    ComaWindDatasetCollection windDatasets = processor.createSHDatasets(
        radii_m, solLongitudes_deg, 8, 6
    );

    const ComaStokesDataset& xStokes = windDatasets.getXStokesDataset();
    const ComaStokesDataset& yStokes = windDatasets.getYStokesDataset();
    const ComaStokesDataset& zStokes = windDatasets.getZStokesDataset();

    // All three components should have the same structure
    BOOST_CHECK_EQUAL(xStokes.nRadii(), yStokes.nRadii());
    BOOST_CHECK_EQUAL(xStokes.nRadii(), zStokes.nRadii());

    BOOST_CHECK_EQUAL(xStokes.nLongitudes(), yStokes.nLongitudes());
    BOOST_CHECK_EQUAL(xStokes.nLongitudes(), zStokes.nLongitudes());

    BOOST_CHECK_EQUAL(xStokes.nmax(), yStokes.nmax());
    BOOST_CHECK_EQUAL(xStokes.nmax(), zStokes.nmax());

    // Verify we can access coefficients without throwing
    BOOST_CHECK_NO_THROW(xStokes.getCoeff(0, 0, 0, 0, 0));
    BOOST_CHECK_NO_THROW(yStokes.getCoeff(0, 0, 0, 0, 0));
    BOOST_CHECK_NO_THROW(zStokes.getCoeff(0, 0, 0, 0, 0));

    // For these test files with zero coefficients, all values should be zero
    auto xCoeff = xStokes.getCoeff(0, 0, 0, 2, 1);
    auto yCoeff = yStokes.getCoeff(0, 0, 0, 2, 1);
    auto zCoeff = zStokes.getCoeff(0, 0, 0, 2, 1);

    BOOST_CHECK_SMALL(xCoeff.first, 1e-15);  // Cosine coefficient
    BOOST_CHECK_SMALL(xCoeff.second, 1e-15); // Sine coefficient
    BOOST_CHECK_SMALL(yCoeff.first, 1e-15);
    BOOST_CHECK_SMALL(yCoeff.second, 1e-15);
    BOOST_CHECK_SMALL(zCoeff.first, 1e-15);
    BOOST_CHECK_SMALL(zCoeff.second, 1e-15);
}

BOOST_FIXTURE_TEST_CASE(test_wind_processor_error_handling, WindTestDataPaths)
{
    // Test error handling for mismatched file types
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Trying to call createSHDatasets() without parameters on poly processor should throw
    BOOST_CHECK_THROW(
        processor.createSHDatasets(),
        std::runtime_error
    );

    // Trying to call createPolyCoefDatasets() on Stokes processor should throw
    // First create some Stokes files
    boost::filesystem::path xOutputDir = outputDir / "wind_x_sh_test";
    boost::filesystem::path yOutputDir = outputDir / "wind_y_sh_test";
    boost::filesystem::path zOutputDir = outputDir / "wind_z_sh_test";

    processor.createSHFiles(
        xOutputDir.string(),
        yOutputDir.string(),
        zOutputDir.string(),
        {6000.0, 10000.0},
        {0.0, 30.0},
        5, 5
    );

    // Create a Stokes-based processor
    ComaWindModelFileProcessor stokesProcessor(
        xOutputDir.string(),
        yOutputDir.string(),
        zOutputDir.string()
    );

    BOOST_CHECK(stokesProcessor.isStokesType());
    BOOST_CHECK(!stokesProcessor.isPolyType());

    // Should throw when trying to create poly datasets
    BOOST_CHECK_THROW(
        stokesProcessor.createPolyCoefDatasets(),
        std::runtime_error
    );

    // Should throw when trying to create SH datasets with parameters
    BOOST_CHECK_THROW(
        stokesProcessor.createSHDatasets({1000.0}, {45.0}),
        std::runtime_error
    );

    // But should work without parameters
    BOOST_CHECK_NO_THROW(
        stokesProcessor.createSHDatasets()
    );
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
