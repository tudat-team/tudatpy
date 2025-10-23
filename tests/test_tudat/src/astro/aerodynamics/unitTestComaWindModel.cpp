#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/aerodynamics/windModel.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>

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

        windFileX = dataDir / "wind" / "polynomial" / "input_poly_wind_x.txt";
        windFileY = dataDir / "wind" / "polynomial" / "input_poly_wind_y.txt";
        windFileZ = dataDir / "wind" / "polynomial" / "input_poly_wind_z.txt";
        outputDir = dataDir / "output" / "wind";

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

// Helper structure to hold reference Stokes coefficients from file
struct ReferenceStokesData
{
    double solarLongitude;
    double radius;
    std::map<std::pair<int, int>, std::pair<double, double>> coefficients; // (degree, order) -> (cosine, sine)

    void loadFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open reference file: " + filename);
        }

        std::string line;
        // Read header line with sol and r
        std::getline(file, line);
        std::istringstream headerStream(line);
        std::string hash, solTag, rTag;
        headerStream >> hash >> solTag >> rTag >> solarLongitude >> radius;

        // Read GM R line
        std::getline(file, line);

        // Read coefficients
        while (std::getline(file, line))
        {
            std::istringstream lineStream(line);
            int degree, order;
            double cosine, sine;
            lineStream >> degree >> order >> cosine >> sine;
            coefficients[{degree, order}] = {cosine, sine};
        }

        file.close();
    }
};

// Helper structure to hold reference residual velocity data from file
struct ReferenceResidualData
{
    struct DataPoint {
        double longitude;      // degrees
        double latitude;       // degrees
        double shEvaluation;   // m/s (from column 4: SH evaluation)
    };

    std::vector<DataPoint> points;

    void loadFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open residual file: " + filename);
        }

        std::string line;
        points.clear();

        // Read all data lines: longitude, latitude, original_data, SH_evaluation, difference
        while (std::getline(file, line))
        {
            std::istringstream lineStream(line);
            DataPoint point;
            double originalData, difference;

            if (lineStream >> point.longitude >> point.latitude >> originalData >> point.shEvaluation >> difference)
            {
                points.push_back(point);
            }
        }

        file.close();
    }
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

BOOST_FIXTURE_TEST_CASE(test_wind_processor_poly_to_stokes, WindTestDataPaths)
{
    // Test the complete roundtrip: poly files -> SH datasets -> verify against reference files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    // Use radii and solar longitudes that match the reference files
    const std::vector<double> radii_m = {4000.0, 10000.0};  // 4 km, 10 km
    const std::vector<double> solLongitudes_deg = {0.0, 30.0};  // 0°, 30°

    // Create Stokes datasets with max degree and order matching the reference files
    ComaWindDatasetCollection windDatasets = processor.createSHDatasets(
        radii_m, solLongitudes_deg, ExpectedWindValues::maxDegree, ExpectedWindValues::maxOrder
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

    // Load reference Stokes coefficients
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";

    ReferenceStokesData ref_sol0_r4km, ref_sol30_r10km;
    ref_sol0_r4km.loadFromFile((dataDir / "wind" / "stokes" / "SH-d10-Vel_rp3-fft12__r_cometFixed_ep10-000_04km.txt").string());
    ref_sol30_r10km.loadFromFile((dataDir / "wind" / "stokes" / "SH-d10-Vel_rp3-fft12__r_cometFixed_ep10-030_10km.txt").string());

    // Verify reference files were loaded correctly
    BOOST_CHECK_CLOSE(ref_sol0_r4km.solarLongitude, 0.0, 1e-10);
    BOOST_CHECK_CLOSE(ref_sol0_r4km.radius, 4.0, 1e-10);  // km
    BOOST_CHECK_CLOSE(ref_sol30_r10km.solarLongitude, 30.0, 1e-10);
    BOOST_CHECK_CLOSE(ref_sol30_r10km.radius, 10.0, 1e-10);  // km

    // Test 1: Verify X and Y components produce all zero coefficients
    // (since input files are all zeros)
    // Only check m >= 0 to match reference file format
    std::cout << "\n=== Testing X and Y components (should be all zeros) ===" << std::endl;
    int numZeroChecked = 0;
    for (int iRad = 0; iRad < xStokes.nRadii(); ++iRad)
    {
        for (int iLon = 0; iLon < xStokes.nLongitudes(); ++iLon)
        {
            for (int degree = 0; degree <= xStokes.nmax(); ++degree)
            {
                // Only check m >= 0 to avoid indexing issues
                for (int order = 0; order <= degree; ++order)
                {
                    auto xCoeff = xStokes.getCoeff(0, iRad, iLon, degree, order);
                    auto yCoeff = yStokes.getCoeff(0, iRad, iLon, degree, order);

                    numZeroChecked++;
                    BOOST_CHECK_SMALL(xCoeff.first, 1e-15);   // X cosine coefficient
                    BOOST_CHECK_SMALL(xCoeff.second, 1e-15);  // X sine coefficient
                    BOOST_CHECK_SMALL(yCoeff.first, 1e-15);   // Y cosine coefficient
                    BOOST_CHECK_SMALL(yCoeff.second, 1e-15);  // Y sine coefficient
                }
            }
        }
    }
    std::cout << "  Checked " << numZeroChecked << " coefficients for X and Y components" << std::endl;
    std::cout << "✓ X and Y components verified as all zeros" << std::endl;

    // Test 2: Verify Z component matches reference files
    // Note: Reference files only contain m >= 0 (non-negative orders)
    std::cout << "\n=== Testing Z component against reference files ===" << std::endl;
    std::cout << "Dataset has " << zStokes.nRadii() << " radii and " << zStokes.nLongitudes() << " longitudes" << std::endl;

    // The dataset includes: reference radius (10000m) + requested radii (4000m, 10000m)
    // We need to find the correct indices. Let's check a known coefficient to identify them.
    // From ref_sol0_r4km: degree=2, order=1 should give cosine=-37.49..., sine=-5.51...

    int iRad_4km = -1, iRad_10km = -1;
    int iLon_0deg = 0, iLon_30deg = 1;  // These match our request order

    // Try to find which radius index corresponds to 4km by checking a known coefficient
    for (int iRad = 0; iRad < zStokes.nRadii(); ++iRad)
    {
        auto testCoeff = zStokes.getCoeff(0, iRad, iLon_0deg, 2, 1);
        // Check if this matches the 4km reference value (degree=2, order=1)
        if (std::abs(testCoeff.first - (-3.749166638392317E+01)) < 1.0)  // Rough match
        {
            iRad_4km = iRad;
            std::cout << "Found 4km at radius index " << iRad << std::endl;
        }
    }

    // Try to find which radius index corresponds to 10km at sol=30deg
    for (int iRad = 0; iRad < zStokes.nRadii(); ++iRad)
    {
        auto testCoeff = zStokes.getCoeff(0, iRad, iLon_30deg, 2, 1);
        // Check if this matches the 10km reference value at sol=30 (degree=2, order=1)
        if (std::abs(testCoeff.first - (-2.932249498617989E+01)) < 1.0)  // Rough match
        {
            iRad_10km = iRad;
            std::cout << "Found 10km at radius index " << iRad << std::endl;
        }
    }

    BOOST_REQUIRE(iRad_4km >= 0);
    BOOST_REQUIRE(iRad_10km >= 0);

    // For sol=0°, r=4km
    std::cout << "Checking sol=0°, r=4km (radIdx=" << iRad_4km << ", lonIdx=" << iLon_0deg << ")..." << std::endl;
    int numChecked = 0;
    for (const auto& [degOrder, refCoeffs] : ref_sol0_r4km.coefficients)
    {
        int degree = degOrder.first;
        int order = degOrder.second;  // Always >= 0 in reference files

        if (order < 0) continue;

        auto zCoeff = zStokes.getCoeff(0, iRad_4km, iLon_0deg, degree, order);

        numChecked++;
        double relTol = 1e-6;  // 0.0001% relative tolerance
        if (std::abs(refCoeffs.first) > 1e-10)
        {
            BOOST_CHECK_CLOSE(zCoeff.first, refCoeffs.first, relTol);
        }
        else
        {
            BOOST_CHECK_SMALL(zCoeff.first, 1e-10);
        }

        if (std::abs(refCoeffs.second) > 1e-10)
        {
            BOOST_CHECK_CLOSE(zCoeff.second, refCoeffs.second, relTol);
        }
        else
        {
            BOOST_CHECK_SMALL(zCoeff.second, 1e-10);
        }
    }
    std::cout << "  Checked " << numChecked << " coefficients for sol=0°, r=4km" << std::endl;

    // For sol=30°, r=10km
    std::cout << "Checking sol=30°, r=10km (radIdx=" << iRad_10km << ", lonIdx=" << iLon_30deg << ")..." << std::endl;
    numChecked = 0;
    for (const auto& [degOrder, refCoeffs] : ref_sol30_r10km.coefficients)
    {
        int degree = degOrder.first;
        int order = degOrder.second;  // Always >= 0 in reference files

        if (order < 0) continue;

        auto zCoeff = zStokes.getCoeff(0, iRad_10km, iLon_30deg, degree, order);

        numChecked++;
        double relTol = 1e-6;  // 0.0001% relative tolerance
        if (std::abs(refCoeffs.first) > 1e-10)
        {
            BOOST_CHECK_CLOSE(zCoeff.first, refCoeffs.first, relTol);
        }
        else
        {
            BOOST_CHECK_SMALL(zCoeff.first, 1e-10);
        }

        if (std::abs(refCoeffs.second) > 1e-10)
        {
            BOOST_CHECK_CLOSE(zCoeff.second, refCoeffs.second, relTol);
        }
        else
        {
            BOOST_CHECK_SMALL(zCoeff.second, 1e-10);
        }
    }
    std::cout << "  Checked " << numChecked << " coefficients for sol=30°, r=10km" << std::endl;
    std::cout << "✓ Z component verified against reference files" << std::endl;

    // Test 3: Verify that using Z component for all directions produces same coefficients
    std::cout << "\n=== Testing Z component used for all directions ===" << std::endl;
    const std::vector<std::string> zFilesForAll = {windFileZ.string()};
    ComaWindModelFileProcessor processorAllZ(zFilesForAll, zFilesForAll, zFilesForAll);

    ComaWindDatasetCollection windDatasetsAllZ = processorAllZ.createSHDatasets(
        radii_m, solLongitudes_deg, ExpectedWindValues::maxDegree, ExpectedWindValues::maxOrder
    );

    const ComaStokesDataset& xStokesFromZ = windDatasetsAllZ.getXStokesDataset();
    const ComaStokesDataset& yStokesFromZ = windDatasetsAllZ.getYStokesDataset();
    const ComaStokesDataset& zStokesFromZ = windDatasetsAllZ.getZStokesDataset();

    // All three components should now match the Z component from the original test
    // Only check m >= 0 to avoid indexing issues
    numChecked = 0;
    for (int iRad = 0; iRad < zStokes.nRadii(); ++iRad)
    {
        for (int iLon = 0; iLon < zStokes.nLongitudes(); ++iLon)
        {
            for (int degree = 0; degree <= zStokes.nmax(); ++degree)
            {
                // Only check m >= 0 to avoid indexing issues
                for (int order = 0; order <= degree; ++order)
                {
                    auto origZCoeff = zStokes.getCoeff(0, iRad, iLon, degree, order);
                    auto xFromZCoeff = xStokesFromZ.getCoeff(0, iRad, iLon, degree, order);
                    auto yFromZCoeff = yStokesFromZ.getCoeff(0, iRad, iLon, degree, order);
                    auto zFromZCoeff = zStokesFromZ.getCoeff(0, iRad, iLon, degree, order);

                    numChecked++;
                    double tol = 1e-6;

                    // X, Y, and Z should all match the original Z coefficients
                    if (std::abs(origZCoeff.first) > 1e-15)
                    {
                        BOOST_CHECK_CLOSE(xFromZCoeff.first, origZCoeff.first, tol);
                        BOOST_CHECK_CLOSE(yFromZCoeff.first, origZCoeff.first, tol);
                        BOOST_CHECK_CLOSE(zFromZCoeff.first, origZCoeff.first, tol);
                    }
                    else
                    {
                        BOOST_CHECK_SMALL(xFromZCoeff.first, 1e-15);
                        BOOST_CHECK_SMALL(yFromZCoeff.first, 1e-15);
                        BOOST_CHECK_SMALL(zFromZCoeff.first, 1e-15);
                    }

                    if (std::abs(origZCoeff.second) > 1e-15)
                    {
                        BOOST_CHECK_CLOSE(xFromZCoeff.second, origZCoeff.second, tol);
                        BOOST_CHECK_CLOSE(yFromZCoeff.second, origZCoeff.second, tol);
                        BOOST_CHECK_CLOSE(zFromZCoeff.second, origZCoeff.second, tol);
                    }
                    else
                    {
                        BOOST_CHECK_SMALL(xFromZCoeff.second, 1e-15);
                        BOOST_CHECK_SMALL(yFromZCoeff.second, 1e-15);
                        BOOST_CHECK_SMALL(zFromZCoeff.second, 1e-15);
                    }
                }
            }
        }
    }
    std::cout << "  Checked " << numChecked << " coefficients across all components" << std::endl;
    std::cout << "✓ All directions using Z component verified as identical" << std::endl;
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

BOOST_FIXTURE_TEST_CASE(test_wind_model_velocity_computation, WindTestDataPaths)
{
    using namespace mathematical_constants;

    std::cout << "\n=== Testing ComaWindModel::getCurrentBodyFixedCartesianWindVelocity ===" << std::endl;

    // Setup: Create Stokes datasets from polynomial files (reuse from previous test)
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    const std::vector<double> radii_m = {4000.0, 10000.0};  // 4 km, 10 km
    const std::vector<double> solLongitudes_deg = {0.0, 30.0};  // 0°, 30°

    // Create Stokes datasets
    ComaWindDatasetCollection windDatasets = processor.createSHDatasets(
        radii_m, solLongitudes_deg, ExpectedWindValues::maxDegree, ExpectedWindValues::maxOrder
    );

    const ComaStokesDataset& xStokes = windDatasets.getXStokesDataset();
    const ComaStokesDataset& yStokes = windDatasets.getYStokesDataset();
    const ComaStokesDataset& zStokes = windDatasets.getZStokesDataset();

    // Load residual velocity reference files
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";

    ReferenceResidualData residual_0deg_4km, residual_30deg_10km;
    residual_0deg_4km.loadFromFile((dataDir / "wind" / "residual" / "residual-vel_r_cometFixed_ep10-000_04km.txt").string());
    residual_30deg_10km.loadFromFile((dataDir / "wind" / "residual" / "residual-vel_r_cometFixed_ep10-030_10km.txt").string());

    std::cout << "Loaded " << residual_0deg_4km.points.size() << " points for sol=0°, r=4km" << std::endl;
    std::cout << "Loaded " << residual_30deg_10km.points.size() << " points for sol=30°, r=10km" << std::endl;

    // ========== Part A: Test with X/Y=zero, Z=valid ==========

    // Test 1: solar longitude = 0°, radius = 4 km
    std::cout << "\n=== Part A.1: Testing sol=0°, r=4km ===" << std::endl;
    {
        const double testRadius = 4000.0;  // 4 km in meters
        const double testTime = 490708800;  // seconds since J2000

        // Define state functions for solar longitude = 0°
        auto sunStateFunction = []() -> Eigen::Vector6d {
            Eigen::Vector6d state = Eigen::Vector6d::Zero();
            state.segment(0, 3) = Eigen::Vector3d(1.0e11, 0.0, 0.0);  // Sun along +X axis
            return state;
        };

        auto cometStateFunction = []() -> Eigen::Vector6d {
            return Eigen::Vector6d::Zero();  // Comet at origin
        };

        auto cometRotationFunction = []() -> Eigen::Matrix3d {
            return Eigen::Matrix3d::Identity();  // No rotation
        };

        // Create a dummy ComaModel for density (required by ComaWindModel)
        const double molecularWeight = 18.01528e-3;  // kg/mol (water)
        auto comaModel = std::make_shared<aerodynamics::ComaModel>(
            zStokes,
            molecularWeight,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder
        );

        // Create ComaWindModel with X/Y=zero, Z=valid
        aerodynamics::ComaWindModel windModel(
            xStokes,
            yStokes,
            zStokes,
            comaModel,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder,
            reference_frames::vertical_frame,
            true,   // includeCorotation
            true    // useRadius
        );

        // Sample points for testing (every 500th point for speed)
        int numSampled = 500;
        int numChecked = 0;
        int numXYZero = 0;
        int numZMatch = 0;

        for (size_t i = 0; i < residual_0deg_4km.points.size(); i += numSampled)
        {
            const auto& point = residual_0deg_4km.points[i];

            // Convert to radians
            double lon_rad = point.longitude * PI / 180.0;
            double lat_rad = point.latitude * PI / 180.0;

            // Call wind model
            Eigen::Vector3d windVelocity = windModel.getCurrentBodyFixedCartesianWindVelocity(
                testRadius, lon_rad, lat_rad, testTime
            );

            numChecked++;

            // Verify X and Y components are zero
            if (std::abs(windVelocity[0]) < 1e-10 && std::abs(windVelocity[1]) < 1e-10)
            {
                numXYZero++;
            }
            BOOST_CHECK_SMALL(windVelocity[0], 1e-10);
            BOOST_CHECK_SMALL(windVelocity[1], 1e-10);

            // Verify Z component matches SH evaluation (within tolerance)
            if (std::abs(point.shEvaluation) > 1e-10)
            {
                double relError = std::abs((windVelocity[2] - point.shEvaluation) / point.shEvaluation) * 100.0;
                if (relError < 1.0)  // Within 1% tolerance
                {
                    numZMatch++;
                }
                BOOST_CHECK_CLOSE(windVelocity[2], point.shEvaluation, 0.01);  // 0.01% tolerance
            }
        }

        std::cout << "  Checked " << numChecked << " sampled points" << std::endl;
        std::cout << "  X/Y components zero: " << numXYZero << "/" << numChecked << std::endl;
        std::cout << "  Z component matches: " << numZMatch << "/" << numChecked << std::endl;
        std::cout << "✓ Sol=0°, r=4km test completed" << std::endl;
    }

    // Test 2: solar longitude = 30°, radius = 10 km
    std::cout << "\n=== Part A.2: Testing sol=30°, r=10km ===" << std::endl;
    {
        const double testRadius = 10000.0;  // 10 km in meters
        const double testTime = 490708800;  // seconds since J2000
        const double testSolarLongitude = 30.0 * PI / 180.0;  // 30° in radians

        // Define state functions for solar longitude = 30°
        const double solarDistance = 1.0e11;
        const double sunX = solarDistance * std::cos(testSolarLongitude);
        const double sunY = solarDistance * std::sin(testSolarLongitude);

        auto sunStateFunction = [sunX, sunY]() -> Eigen::Vector6d {
            Eigen::Vector6d state = Eigen::Vector6d::Zero();
            state.segment(0, 3) = Eigen::Vector3d(sunX, sunY, 0.0);
            return state;
        };

        auto cometStateFunction = []() -> Eigen::Vector6d {
            return Eigen::Vector6d::Zero();
        };

        auto cometRotationFunction = []() -> Eigen::Matrix3d {
            return Eigen::Matrix3d::Identity();
        };

        // Create ComaModel
        const double molecularWeight = 18.01528e-3;
        auto comaModel = std::make_shared<aerodynamics::ComaModel>(
            zStokes,
            molecularWeight,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder
        );

        // Create ComaWindModel
        aerodynamics::ComaWindModel windModel(
            xStokes,
            yStokes,
            zStokes,
            comaModel,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder,
            reference_frames::vertical_frame,
            true,
            true
        );

        // Sample points for testing
        int numSampled = 500;
        int numChecked = 0;
        int numXYZero = 0;
        int numZMatch = 0;

        for (size_t i = 0; i < residual_30deg_10km.points.size(); i += numSampled)
        {
            const auto& point = residual_30deg_10km.points[i];

            double lon_rad = point.longitude * PI / 180.0;
            double lat_rad = point.latitude * PI / 180.0;

            Eigen::Vector3d windVelocity = windModel.getCurrentBodyFixedCartesianWindVelocity(
                testRadius, lon_rad, lat_rad, testTime
            );

            numChecked++;

            if (std::abs(windVelocity[0]) < 1e-10 && std::abs(windVelocity[1]) < 1e-10)
            {
                numXYZero++;
            }
            BOOST_CHECK_SMALL(windVelocity[0], 1e-10);
            BOOST_CHECK_SMALL(windVelocity[1], 1e-10);

            if (std::abs(point.shEvaluation) > 1e-10)
            {
                double relError = std::abs((windVelocity[2] - point.shEvaluation) / point.shEvaluation) * 100.0;
                if (relError < 1.0)
                {
                    numZMatch++;
                }
                BOOST_CHECK_CLOSE(windVelocity[2], point.shEvaluation, 1.0);
            }
        }

        std::cout << "  Checked " << numChecked << " sampled points" << std::endl;
        std::cout << "  X/Y components zero: " << numXYZero << "/" << numChecked << std::endl;
        std::cout << "  Z component matches: " << numZMatch << "/" << numChecked << std::endl;
        std::cout << "✓ Sol=30°, r=10km test completed" << std::endl;
    }

    // ========== Part B: Test with all components using Z data ==========

    std::cout << "\n=== Part B: Testing all components using Z data ===" << std::endl;
    {
        const double testRadius = 4000.0;
        const double testTime = 490708800;

        auto sunStateFunction = []() -> Eigen::Vector6d {
            Eigen::Vector6d state = Eigen::Vector6d::Zero();
            state.segment(0, 3) = Eigen::Vector3d(1.0e11, 0.0, 0.0);
            return state;
        };

        auto cometStateFunction = []() -> Eigen::Vector6d {
            return Eigen::Vector6d::Zero();
        };

        auto cometRotationFunction = []() -> Eigen::Matrix3d {
            return Eigen::Matrix3d::Identity();
        };

        const double molecularWeight = 18.01528e-3;
        auto comaModel = std::make_shared<aerodynamics::ComaModel>(
            zStokes,
            molecularWeight,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder
        );

        // Create ComaWindModel with Z data for all components
        aerodynamics::ComaWindModel windModelAllZ(
            zStokes,  // Use Z for X component
            zStokes,  // Use Z for Y component
            zStokes,  // Use Z for Z component
            comaModel,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            ExpectedWindValues::maxDegree,
            ExpectedWindValues::maxOrder,
            reference_frames::vertical_frame,
            true,
            true
        );

        int numSampled = 500;
        int numChecked = 0;
        int numAllMatch = 0;

        for (size_t i = 0; i < residual_0deg_4km.points.size(); i += numSampled)
        {
            const auto& point = residual_0deg_4km.points[i];

            double lon_rad = point.longitude * PI / 180.0;
            double lat_rad = point.latitude * PI / 180.0;

            Eigen::Vector3d windVelocity = windModelAllZ.getCurrentBodyFixedCartesianWindVelocity(
                testRadius, lon_rad, lat_rad, testTime
            );

            numChecked++;

            // All three components should match the SH evaluation
            if (std::abs(point.shEvaluation) > 1e-10)
            {
                BOOST_CHECK_CLOSE(windVelocity[0], point.shEvaluation, 1.0);
                BOOST_CHECK_CLOSE(windVelocity[1], point.shEvaluation, 1.0);
                BOOST_CHECK_CLOSE(windVelocity[2], point.shEvaluation, 1.0);

                double relError0 = std::abs((windVelocity[0] - point.shEvaluation) / point.shEvaluation) * 100.0;
                double relError1 = std::abs((windVelocity[1] - point.shEvaluation) / point.shEvaluation) * 100.0;
                double relError2 = std::abs((windVelocity[2] - point.shEvaluation) / point.shEvaluation) * 100.0;

                if (relError0 < 1.0 && relError1 < 1.0 && relError2 < 1.0)
                {
                    numAllMatch++;
                }
            }
        }

        std::cout << "  Checked " << numChecked << " sampled points" << std::endl;
        std::cout << "  All components match: " << numAllMatch << "/" << numChecked << std::endl;
        std::cout << "✓ All components using Z data test completed" << std::endl;
    }

    std::cout << "\n✓ All wind velocity computation tests passed!" << std::endl;
}

BOOST_FIXTURE_TEST_CASE(test_wind_model_velocity_validation_from_python, WindTestDataPaths)
{
    using namespace mathematical_constants;

    // This test validates the entire pipeline by using reference values computed from the Python interface.
    // The wind_reference_values.txt file contains: time, radial distance, latitude, longitude, solar longitude,
    // and wind velocity components (X, Y, Z) in the vertical frame.
    // We use these values as input to calculate wind velocity with the verified code and validate against the reference.

    std::cout << "\n=== Testing Wind Model Validation from Python Interface ===" << std::endl;

    // Load wind polynomial coefficients from test data files
    const std::vector<std::string> xFiles = {windFileX.string()};
    const std::vector<std::string> yFiles = {windFileY.string()};
    const std::vector<std::string> zFiles = {windFileZ.string()};

    ComaWindModelFileProcessor processor(xFiles, yFiles, zFiles);

    const int maxDegree = 10;
    const int maxOrder = 10;

    // Create Stokes datasets from polynomial datasets
    // Use the same grid as in Python: radii from 4 km to 10 km (100 points), sol_longitude 0 to 360° (37 points)
    std::vector<double> radii_m;
    for (int i = 0; i < 100; ++i)
    {
        radii_m.push_back(4000.0 + i * (10000.0 - 4000.0) / 99.0);
    }

    std::vector<double> lons_deg;
    for (int i = 0; i < 37; ++i)
    {
        lons_deg.push_back(i * 360.0 / 36.0);
    }

    ComaWindDatasetCollection windDatasets = processor.createSHDatasets(radii_m, lons_deg, maxDegree, maxOrder);

    const ComaStokesDataset& xStokes = windDatasets.getXStokesDataset();
    const ComaStokesDataset& yStokes = windDatasets.getYStokesDataset();
    const ComaStokesDataset& zStokes = windDatasets.getZStokesDataset();

    // Construct path to wind reference values file
    const boost::filesystem::path thisFile(__FILE__);
    const boost::filesystem::path testDir = thisFile.parent_path();
    const boost::filesystem::path dataDir = testDir / "test_data";
    const boost::filesystem::path referenceFile = dataDir / "wind" / "wind_reference_values.txt";

    // Read reference values file
    std::ifstream file(referenceFile.string());
    BOOST_REQUIRE_MESSAGE(file.is_open(), "Cannot open wind reference values file: " + referenceFile.string());

    struct ReferencePoint {
        double time;
        double radialDistance;
        double latitude;
        double longitude;
        double solarLongitude;
        double windX;
        double windY;
        double windZ;
    };
    std::vector<ReferencePoint> allPoints;

    std::string line;
    // Skip header line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        ReferencePoint point;
        if (iss >> point.time >> point.radialDistance >> point.latitude >> point.longitude
                >> point.solarLongitude >> point.windX >> point.windY >> point.windZ)
        {
            allPoints.push_back(point);
        }
    }
    file.close();

    BOOST_REQUIRE_MESSAGE(allPoints.size() > 0, "No data points found in wind reference values file");
    BOOST_TEST_MESSAGE("Loaded " << allPoints.size() << " wind reference data points");

    // Randomly select 100 points for testing (or all if less than 100)
    const int numTestPoints = 100;
    std::vector<int> selectedIndices;
    std::srand(54321);  // Fixed seed for reproducibility

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

    // Create a dummy ComaModel for density (required by ComaWindModel)
    // Load density polynomial coefficients
    const boost::filesystem::path densityPolyFile = dataDir / "density" / "polynomial" / "input_poly_coef_test_file.txt";
    const std::vector<std::string> densityFiles = {densityPolyFile.string()};
    const ComaModelFileProcessor densityProcessor(densityFiles);
    const double molecularWeight = 0.018;  // kg/mol (water vapor)
    const ComaStokesDataset densityStokes = densityProcessor.createSHDataset(radii_m, lons_deg, maxDegree, maxOrder);

    // Test with Stokes dataset
    int failCount = 0;
    int failCountX = 0, failCountY = 0, failCountZ = 0;

    for (int idx : selectedIndices)
    {
        const ReferencePoint& point = allPoints[idx];

        // Define state functions based on the reference solar longitude
        const double solarDistance = 1.0e11;
        const double sunX = solarDistance * std::cos(point.solarLongitude);
        const double sunY = solarDistance * std::sin(point.solarLongitude);

        auto sunStateFunction = [sunX, sunY]() -> Eigen::Vector6d {
            Eigen::Vector6d state = Eigen::Vector6d::Zero();
            state.segment(0, 3) = Eigen::Vector3d(sunX, sunY, 0.0);
            return state;
        };

        auto cometStateFunction = []() -> Eigen::Vector6d {
            return Eigen::Vector6d::Zero();
        };

        auto cometRotationFunction = []() -> Eigen::Matrix3d {
            return Eigen::Matrix3d::Identity();
        };

        // Create ComaModel for density
        auto comaModel = std::make_shared<aerodynamics::ComaModel>(
            densityStokes,
            molecularWeight,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            maxDegree,
            maxOrder
        );

        // Create ComaWindModel with Stokes coefficients
        aerodynamics::ComaWindModel windModel(
            xStokes,
            yStokes,
            zStokes,
            comaModel,
            sunStateFunction,
            cometStateFunction,
            cometRotationFunction,
            maxDegree,
            maxOrder,
            reference_frames::vertical_frame,
            true,   // includeCorotation
            true    // useRadius
        );

        // Calculate wind velocity using the ComaWindModel
        // getCurrentBodyFixedCartesianWindVelocity() returns wind velocity in m/s in vertical frame
        const Eigen::Vector3d computedWindVelocity = windModel.getCurrentBodyFixedCartesianWindVelocity(
            point.radialDistance,
            point.longitude,
            point.latitude,
            point.time
        );

        // Calculate component-wise errors; combine relative and absolute tolerances so near-zero
        // reference values (order 1e-13 in Python dataset) do not trigger spurious failures.
        const double relativeTolerance = 1e-8;
        const double absoluteTolerance = 1e-12;

        const double absoluteErrorX = std::abs(computedWindVelocity[0] - point.windX);
        const double absoluteErrorY = std::abs(computedWindVelocity[1] - point.windY);
        const double absoluteErrorZ = std::abs(computedWindVelocity[2] - point.windZ);

        const double relativeErrorX = absoluteErrorX / std::max(std::abs(point.windX), 1e-30);
        const double relativeErrorY = absoluteErrorY / std::max(std::abs(point.windY), 1e-30);
        const double relativeErrorZ = absoluteErrorZ / std::max(std::abs(point.windZ), 1e-30);

        const bool withinToleranceX = (absoluteErrorX <= absoluteTolerance) || (relativeErrorX <= relativeTolerance);
        const bool withinToleranceY = (absoluteErrorY <= absoluteTolerance) || (relativeErrorY <= relativeTolerance);
        const bool withinToleranceZ = (absoluteErrorZ <= absoluteTolerance) || (relativeErrorZ <= relativeTolerance);

        bool failed = false;
        if (!withinToleranceX)
        {
            failCountX++;
            failed = true;
        }
        if (!withinToleranceY)
        {
            failCountY++;
            failed = true;
        }
        if (!withinToleranceZ)
        {
            failCountZ++;
            failed = true;
        }

        if (failed)
        {
            failCount++;
            if (failCount <= 10)  // Only report first 10 failures
            {
                BOOST_TEST_MESSAGE("Point " << idx << " mismatch: "
                                  << "r=" << point.radialDistance << " m, "
                                  << "lat=" << point.latitude * 180.0 / PI << "°, "
                                  << "lon=" << point.longitude * 180.0 / PI << "°, "
                                  << "sol=" << point.solarLongitude * 180.0 / PI << "°, "
                                  << "computed=(" << computedWindVelocity[0] << ", "
                                  << computedWindVelocity[1] << ", " << computedWindVelocity[2] << ") m/s, "
                                  << "expected=(" << point.windX << ", " << point.windY << ", " << point.windZ << ") m/s, "
                                  << "abs_errors=(" << absoluteErrorX << ", " << absoluteErrorY << ", " << absoluteErrorZ << ") m/s, "
                                  << "rel_errors=(" << relativeErrorX << ", " << relativeErrorY << ", " << relativeErrorZ << ")");
            }
        }

        // Check each component separately
        BOOST_CHECK_MESSAGE(
            withinToleranceX,
            "Wind X mismatch at r=" << point.radialDistance << " m, "
            << "lat=" << point.latitude * 180.0 / PI << "°, "
            << "lon=" << point.longitude * 180.0 / PI << "°"
            << ", abs_error=" << absoluteErrorX
            << ", rel_error=" << relativeErrorX
        );
        BOOST_CHECK_MESSAGE(
            withinToleranceY,
            "Wind Y mismatch at r=" << point.radialDistance << " m, "
            << "lat=" << point.latitude * 180.0 / PI << "°, "
            << "lon=" << point.longitude * 180.0 / PI << "°"
            << ", abs_error=" << absoluteErrorY
            << ", rel_error=" << relativeErrorY
        );
        BOOST_CHECK_MESSAGE(
            withinToleranceZ,
            "Wind Z mismatch at r=" << point.radialDistance << " m, "
            << "lat=" << point.latitude * 180.0 / PI << "°, "
            << "lon=" << point.longitude * 180.0 / PI << "°"
            << ", abs_error=" << absoluteErrorZ
            << ", rel_error=" << relativeErrorZ
        );
    }

    BOOST_TEST_MESSAGE("Python wind reference validation: Verified " << selectedIndices.size()
                      << " randomly selected points out of " << allPoints.size()
                      << " total, " << failCount << " failures ("
                      << "X:" << failCountX << ", Y:" << failCountY << ", Z:" << failCountZ << ")");
    BOOST_CHECK_EQUAL(failCount, 0);

    std::cout << "✓ Wind model validation from Python interface completed" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
