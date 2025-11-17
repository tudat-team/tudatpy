#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <memory>
#include <functional>
#include <boost/filesystem.hpp>

// Tudat includes for the actual ComaModel classes
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"

// NOTE: The full Tudat integration requires resolving the resource header issue.
// This template shows how to use the actual ComaModel classes.
// For now, we provide a framework that demonstrates the methodology.

// When tudat/resource/resource.h is available, uncomment these:
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/aerodynamics/windModel.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"



using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;

class TudatPerformanceTest {
private:
    std::vector<double> radiiGrid_;
    std::vector<double> longitudeGrid_;

    // Test parameters
    const double minRadius_ = 4000.0;   // 4 km
    const double maxRadius_ = 20000.0;  // 20 km
    const int numRadii_ = 600;
    const int numLongitudes_ = 400;
    const double latitude_ = 0.0;       // Fixed at 0 degrees
    const double testTime_ = 490708800;   // s since J2000
    const int max_degree = 10;
    const int max_order = 10;

    // Mock state functions for ComaModel
    std::function<Eigen::Vector6d()> sunStateFunction_;
    std::function<Eigen::Vector6d()> cometStateFunction_;
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

public:
    struct TestResults {
        double polyDatasetCreationTime;
        double stokesDatasetCreationTime;
        double polyComputationTime;
        double stokesComputationTime;
        std::vector<std::vector<Eigen::Vector3d>> polyWindGrid;
        std::vector<std::vector<Eigen::Vector3d>> stokesWindGrid;
    };

    TudatPerformanceTest() {
        // Create test grids
        radiiGrid_ = createLinearGrid(minRadius_, maxRadius_, numRadii_);
        longitudeGrid_ = createLinearGrid(0.0, 360.0, numLongitudes_);

        // Setup mock state functions for ComaModel
        setupMockStateFunctions();

        std::cout << "Initialized Tudat performance test with grid: "
                  << numRadii_ << " x " << numLongitudes_
                  << " = " << (numRadii_ * numLongitudes_) << " points" << std::endl;
    }

    TestResults runTest() {
        TestResults results;

        std::cout << "\n=== Tudat Coma Wind Performance Test (Using Real ComaWindModel) ===" << std::endl;

        // Initialize result grids
        results.polyWindGrid.resize(numRadii_, std::vector<Eigen::Vector3d>(numLongitudes_));
        results.stokesWindGrid.resize(numRadii_, std::vector<Eigen::Vector3d>(numLongitudes_));

        // Test with the actual wind data files (relative to bin directory)
        std::string xWindFilePath = "../../examples/tudat/tudat/applications/coma_analysis/input/input_poly_wind_x.txt";
        std::string yWindFilePath = "../../examples/tudat/tudat/applications/coma_analysis/input/input_poly_wind_y.txt";
        std::string zWindFilePath = "../../examples/tudat/tudat/applications/coma_analysis/input/input_poly_wind_z.txt";

        // Debug: Print current working directory
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        std::cout << "Looking for wind files..." << std::endl;

        if (!fileExists(xWindFilePath) || !fileExists(yWindFilePath) || !fileExists(zWindFilePath)) {
            std::cout << "Error: Wind data files not found" << std::endl;
            return results;
        }

        std::cout << "Using wind coefficient files: " << std::endl;
        std::cout << "  X: " << xWindFilePath << std::endl;
        std::cout << "  Y: " << yWindFilePath << std::endl;
        std::cout << "  Z: " << zWindFilePath << std::endl;

        // Use ComaWindModelFileProcessor to create datasets
        std::vector<std::string> xWindFiles = {xWindFilePath};
        std::vector<std::string> yWindFiles = {yWindFilePath};
        std::vector<std::string> zWindFiles = {zWindFilePath};

        // Test 1: Create polynomial wind datasets using actual Tudat classes
        std::cout << "Creating polynomial wind datasets using ComaWindModelFileProcessor..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        ComaWindModelFileProcessor windProcessor(xWindFiles, yWindFiles, zWindFiles);
        ComaWindDatasetCollection polyWindDatasets = windProcessor.createPolyCoefDatasets();

        auto end = std::chrono::high_resolution_clock::now();
        results.polyDatasetCreationTime = std::chrono::duration<double>(end - start).count();
        std::cout << "Polynomial wind datasets creation time: " << results.polyDatasetCreationTime << " s" << std::endl;

        // Test 2: Create Stokes wind datasets using actual Tudat transformation
        std::cout << "Creating Stokes wind datasets using ComaWindModelFileProcessor..." << std::endl;
        start = std::chrono::high_resolution_clock::now();

        // Create grid for Stokes dataset
        std::vector<double> stokesRadii = radiiGrid_;
        std::vector<double> stokesLongitudes = longitudeGrid_;

        ComaWindDatasetCollection stokesWindDatasets = windProcessor.createSHDatasets(stokesRadii, stokesLongitudes);

        end = std::chrono::high_resolution_clock::now();
        results.stokesDatasetCreationTime = std::chrono::duration<double>(end - start).count();
        std::cout << "Stokes wind datasets creation time: " << results.stokesDatasetCreationTime << " s" << std::endl;

        // Test 3: Create dummy ComaModel (required by ComaWindModel constructor)
        std::cout << "Creating dummy ComaModel for wind model..." << std::endl;
        std::string densityFilePath = "../../examples/tudat/tudat/applications/coma_analysis/input/input_poly_coef_test_file.txt";

        if (!fileExists(densityFilePath)) {
            std::cout << "Error: Density data file not found at " << densityFilePath << std::endl;
            return results;
        }

        std::vector<std::string> densityFiles = {densityFilePath};
        ComaModelFileProcessor densityProcessor(densityFiles);
        ComaStokesDataset dummyStokesDataset = densityProcessor.createSHDataset(stokesRadii, stokesLongitudes);
        const double molecularWeight = 0.018015; // kg/mol for water vapor

        auto comaModel = std::make_shared<ComaModel>(dummyStokesDataset, molecularWeight,
                                                     sunStateFunction_, cometStateFunction_,
                                                     cometRotationFunction_, max_degree, max_order);

        // Test 4: Create actual ComaWindModel instances
        std::cout << "Creating ComaWindModel instances..." << std::endl;

        // Extract datasets from collections
        const ComaPolyDataset& xPolyDataset = polyWindDatasets.getXPolyDataset();
        const ComaPolyDataset& yPolyDataset = polyWindDatasets.getYPolyDataset();
        const ComaPolyDataset& zPolyDataset = polyWindDatasets.getZPolyDataset();

        const ComaStokesDataset& xStokesDataset = stokesWindDatasets.getXStokesDataset();
        const ComaStokesDataset& yStokesDataset = stokesWindDatasets.getYStokesDataset();
        const ComaStokesDataset& zStokesDataset = stokesWindDatasets.getZStokesDataset();

        ComaWindModel polyWindModel(xPolyDataset, yPolyDataset, zPolyDataset, comaModel,
                                    sunStateFunction_, cometStateFunction_,
                                    cometRotationFunction_, max_degree, max_order);

        ComaWindModel stokesWindModel(xStokesDataset, yStokesDataset, zStokesDataset, comaModel,
                                     sunStateFunction_, cometStateFunction_,
                                     cometRotationFunction_, max_degree, max_order);

        // Test 5: Measure polynomial method computation using actual ComaWindModel
        std::cout << "Computing wind velocities using polynomial ComaWindModel..." << std::endl;
        start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < numRadii_; ++i) {
            for (int j = 0; j < numLongitudes_; ++j) {
                double radius = radiiGrid_[i];
                double longitude = longitudeGrid_[j] * mathematical_constants::PI / 180.0; // Convert to radians
                double latitude = latitude_ * mathematical_constants::PI / 180.0;

                Eigen::Vector3d windVelocity = polyWindModel.getCurrentBodyFixedCartesianWindVelocity(
                    radius, longitude, latitude, testTime_);
                results.polyWindGrid[i][j] = windVelocity;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        results.polyComputationTime = std::chrono::duration<double>(end - start).count();
        std::cout << "Polynomial computation time: " << results.polyComputationTime << " s" << std::endl;

        // Test 6: Measure Stokes method computation using actual ComaWindModel
        std::cout << "Computing wind velocities using Stokes ComaWindModel..." << std::endl;
        start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < numRadii_; ++i) {
            for (int j = 0; j < numLongitudes_; ++j) {
                double radius = radiiGrid_[i];
                double longitude = longitudeGrid_[j] * mathematical_constants::PI / 180.0; // Convert to radians
                double latitude = latitude_ * mathematical_constants::PI / 180.0;

                Eigen::Vector3d windVelocity = stokesWindModel.getCurrentBodyFixedCartesianWindVelocity(
                    radius, longitude, latitude, testTime_);
                results.stokesWindGrid[i][j] = windVelocity;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        results.stokesComputationTime = std::chrono::duration<double>(end - start).count();
        std::cout << "Stokes computation time: " << results.stokesComputationTime << " s" << std::endl;

        return results;
    }

private:
    void setupMockStateFunctions() {
        // Mock sun state (position at 1 AU for solar longitude calculation)
        sunStateFunction_ = []() -> Eigen::Vector6d {
            Eigen::Vector6d sunState = Eigen::Vector6d::Zero();
            sunState[0] = 1.5e11; // 1 AU in meters (approximate)
            return sunState;
        };

        // Mock comet state (at origin)
        cometStateFunction_ = []() -> Eigen::Vector6d {
            return Eigen::Vector6d::Zero();
        };

        // Mock rotation function (will be set for each solar longitude)
        cometRotationFunction_ = []() -> Eigen::Matrix3d {
            return Eigen::Matrix3d::Identity();
        };
    }


    std::vector<double> createLinearGrid(double min, double max, int num) {
        std::vector<double> grid;
        grid.reserve(num);

        for (int i = 0; i < num; ++i) {
            double value = min + (max - min) * static_cast<double>(i) / static_cast<double>(num - 1);
            grid.push_back(value);
        }

        return grid;
    }

    bool fileExists(const std::string& filename) {
        std::ifstream file(filename);
        return file.good();
    }

public:
    void setSolarLongitude(double solarLongitudeDeg) {
        // Convert degrees to radians
        double solarLongitudeRad = solarLongitudeDeg * mathematical_constants::PI / 180.0;

        // Create rotation matrix for the solar longitude
        // This rotates around the z-axis (rotation about comet's spin axis)
        cometRotationFunction_ = [solarLongitudeRad]() -> Eigen::Matrix3d {
            Eigen::Matrix3d rotation;
            rotation << std::cos(solarLongitudeRad), -std::sin(solarLongitudeRad), 0,
                       std::sin(solarLongitudeRad),  std::cos(solarLongitudeRad), 0,
                       0,                            0,                           1;
            return rotation;
        };
    }
    void saveResults(const TestResults& results, double solarLongitudeDeg) {
        std::cout << "\nSaving results to CSV files for solar longitude " << solarLongitudeDeg << "°..." << std::endl;

        // Create subdirectory for this solar longitude (relative to bin directory)
        std::string baseOutputDir = "../../examples/tudat/tudat/applications/coma_analysis/performance_wind_analysis/output";
        std::string outputDir = baseOutputDir + "/solar_longitude_" + std::to_string(static_cast<int>(solarLongitudeDeg));

        // Create directory using system command
        std::string mkdirCmd = "mkdir -p " + outputDir;
        system(mkdirCmd.c_str());

        // Save timing results
        std::string timingPath = outputDir + "/timing_results.csv";
        std::ofstream timingFile(timingPath);
        timingFile << "Solar Longitude (deg): " << solarLongitudeDeg << std::endl;
        timingFile << "Number radii: " << numRadii_ << std::endl;
        timingFile << "Number longitude: " << numLongitudes_ << std::endl;
        timingFile << "Evaluated Points: " << (numRadii_ * numLongitudes_) << std::endl;
        timingFile << "Method,Operation,Time_seconds" << std::endl;
        timingFile << "Polynomial,Dataset_Creation," << results.polyDatasetCreationTime << std::endl;
        timingFile << "Polynomial,Wind_Computation," << results.polyComputationTime << std::endl;
        timingFile << "Stokes,Dataset_Creation," << results.stokesDatasetCreationTime << std::endl;
        timingFile << "Stokes,Wind_Computation," << results.stokesComputationTime << std::endl;
        timingFile.close();

        // Save comparison (includes both poly and stokes wind velocities)
        std::string diffPath = outputDir + "/wind_velocity_differences.csv";
        std::ofstream diffFile(diffPath);
        diffFile << "Radius_m,Longitude_deg,"
                 << "Poly_Vx,Poly_Vy,Poly_Vz,Poly_Magnitude,"
                 << "Stokes_Vx,Stokes_Vy,Stokes_Vz,Stokes_Magnitude,"
                 << "Magnitude_Diff,Magnitude_Relative_Diff_Percent,"
                 << "Vector_Diff_Norm" << std::endl;
        for (size_t i = 0; i < radiiGrid_.size(); ++i) {
            for (size_t j = 0; j < longitudeGrid_.size(); ++j) {
                Eigen::Vector3d polyWind = results.polyWindGrid[i][j];
                Eigen::Vector3d stokesWind = results.stokesWindGrid[i][j];

                double polyMag = polyWind.norm();
                double stokesMag = stokesWind.norm();
                double magDiff = std::abs(polyMag - stokesMag);
                double magRelDiff = (polyMag != 0.0) ? (magDiff / polyMag) * 100.0 : 0.0;
                double vectorDiffNorm = (polyWind - stokesWind).norm();

                diffFile << std::fixed << std::setprecision(6)
                         << radiiGrid_[i] << ","
                         << longitudeGrid_[j] << ","
                         << polyWind(0) << "," << polyWind(1) << "," << polyWind(2) << ","
                         << polyMag << ","
                         << stokesWind(0) << "," << stokesWind(1) << "," << stokesWind(2) << ","
                         << stokesMag << ","
                         << magDiff << ","
                         << magRelDiff << ","
                         << vectorDiffNorm << std::endl;
            }
        }
        diffFile.close();

        std::cout << "Results saved to " << outputDir << std::endl;
    }

    void printSummary(const TestResults& results) {
        std::cout << "\n=== TUDAT WIND PERFORMANCE SUMMARY ===" << std::endl;

        std::cout << "Grid size: " << numRadii_ << " x " << numLongitudes_
                  << " = " << (numRadii_ * numLongitudes_) << " evaluations per method" << std::endl;

        std::cout << "\nPolynomial method:" << std::endl;
        std::cout << "  Dataset creation: " << std::fixed << std::setprecision(4)
                  << results.polyDatasetCreationTime << " s" << std::endl;
        std::cout << "  Wind velocity computation: " << results.polyComputationTime << " s" << std::endl;
        std::cout << "  Per-evaluation: " << (results.polyComputationTime / (numRadii_ * numLongitudes_)) * 1000.0
                  << " ms" << std::endl;
        std::cout << "  Total: " << (results.polyDatasetCreationTime + results.polyComputationTime) << " s" << std::endl;

        std::cout << "\nStokes method:" << std::endl;
        std::cout << "  Dataset creation: " << results.stokesDatasetCreationTime << " s" << std::endl;
        std::cout << "  Wind velocity computation: " << results.stokesComputationTime << " s" << std::endl;
        std::cout << "  Per-evaluation: " << (results.stokesComputationTime / (numRadii_ * numLongitudes_)) * 1000.0
                  << " ms" << std::endl;
        std::cout << "  Total: " << (results.stokesDatasetCreationTime + results.stokesComputationTime) << " s" << std::endl;

        double polyTotal = results.polyDatasetCreationTime + results.polyComputationTime;
        double stokesTotal = results.stokesDatasetCreationTime + results.stokesComputationTime;
        double speedup = polyTotal / stokesTotal;

        std::cout << "\nComparison:" << std::endl;
        std::cout << "  Setup time ratio (Stokes/Poly): "
                  << (results.stokesDatasetCreationTime / results.polyDatasetCreationTime) << "x" << std::endl;
        std::cout << "  Computation speedup (Poly/Stokes): "
                  << (results.polyComputationTime / results.stokesComputationTime) << "x" << std::endl;
        std::cout << "  Overall speedup: " << speedup << "x ";

        if (speedup > 1.0) {
            std::cout << "(Stokes method is " << speedup << "x faster overall)" << std::endl;
        } else {
            std::cout << "(Polynomial method is " << (1.0/speedup) << "x faster overall)" << std::endl;
        }
    }
};

int main() {
    try {
        std::cout << "Tudat Coma Wind Model Performance Test" << std::endl;
        std::cout << "=======================================" << std::endl;
        std::cout << "Running tests for solar longitudes from 0° to 350° in 10° steps\n" << std::endl;

        TudatPerformanceTest test;

        // Loop through solar longitudes from 0° to 350° in 10° steps
        for (int solarLongDeg = 0; solarLongDeg <= 350; solarLongDeg += 10) {
            std::cout << "\n=== SOLAR LONGITUDE " << solarLongDeg << "° ===" << std::endl;

            // Set the solar longitude for this iteration
            test.setSolarLongitude(static_cast<double>(solarLongDeg));

            // Run the test
            auto results = test.runTest();

            // Save results to subdirectory
            test.saveResults(results, static_cast<double>(solarLongDeg));

            // Print summary for this solar longitude
            test.printSummary(results);
        }

        std::cout << "\n=== ALL TESTS COMPLETED ===" << std::endl;
        std::cout << "Results saved in examples/tudat/tudat/applications/coma_analysis/performance_wind_analysis/output/solar_longitude_*/ directories" << std::endl;
        std::cout << "Tudat wind performance test completed successfully!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}