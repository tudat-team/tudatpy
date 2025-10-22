#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <functional>
#include <regex>

// Tudat includes for the actual ComaWindModel classes
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/windModel.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;

struct DSMCDataPoint {
    double longitude;    // degrees
    double coLatitude;   // degrees (0 = north, 180 = south)
    double density;      // number density [m-3] (not used for wind comparison)
    double windSpeed;    // wind magnitude [m/s] - assumed to be radial component

    // Computed values
    double latitude;     // converted from co-latitude
    double myWindZ;      // z-component of wind from our algorithm
};

struct DSMCFile {
    std::string filename;
    double solarLongitude;  // degrees
    double cometDistance;   // km
    std::vector<DSMCDataPoint> dataPoints;
};

class DSMCWindComparison {
private:
    std::vector<DSMCFile> dsmc_files_;

    // Mock state functions for ComaWindModel
    std::function<Eigen::Vector6d()> sunStateFunction_;
    std::function<Eigen::Vector6d()> cometStateFunction_;
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

    const double testTime_ = 490708800;   // s since J2000

public:
    DSMCWindComparison() {
        setupMockStateFunctions();
    }

    bool loadDSMCFiles(const std::vector<std::string>& filePaths) {
        dsmc_files_.clear();

        for (const auto& filePath : filePaths) {
            DSMCFile dsmc_file;
            dsmc_file.filename = filePath;

            if (!parseDSMCFile(filePath, dsmc_file)) {
                std::cerr << "Failed to parse DSMC file: " << filePath << std::endl;
                return false;
            }

            dsmc_files_.push_back(dsmc_file);
            std::cout << "Loaded DSMC file: " << filePath
                      << " (Solar longitude: " << dsmc_file.solarLongitude
                      << "°, Distance: " << dsmc_file.cometDistance << " km, "
                      << dsmc_file.dataPoints.size() << " data points)" << std::endl;
        }

        return true;
    }

    void runComparison() {
        for (auto& dsmc_file : dsmc_files_) {
            std::cout << "\nProcessing file: " << dsmc_file.filename << std::endl;
            std::cout << "Solar longitude: " << dsmc_file.solarLongitude << "°" << std::endl;
            std::cout << "Comet distance: " << dsmc_file.cometDistance << " km" << std::endl;

            // Set up the coma model for this solar longitude
            setSolarLongitude(dsmc_file.solarLongitude);

            // Create wind model using the polynomial coefficient files
            std::string inputDir = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/input/";
            std::vector<std::string> xWindFiles = {inputDir + "input_poly_wind_x.txt"};
            std::vector<std::string> yWindFiles = {inputDir + "input_poly_wind_y.txt"};
            std::vector<std::string> zWindFiles = {inputDir + "input_poly_wind_z.txt"};

            // Verify files exist
            for (const auto& file : {xWindFiles[0], yWindFiles[0], zWindFiles[0]}) {
                if (!fileExists(file)) {
                    std::cerr << "Error: Wind polynomial file not found at " << file << std::endl;
                    continue;
                }
            }

            // Create processor for wind model
            ComaWindModelFileProcessor windProcessor(xWindFiles, yWindFiles, zWindFiles);

            std::cout << "Creating wind datasets..." << std::endl;
            ComaWindDatasetCollection windDatasets = windProcessor.createPolyCoefDatasets();

            // Create a dummy ComaModel (needed for ComaWindModel constructor, but we only use wind)
            // We'll create a minimal one just for the wind model to reference
            std::string polyFilePath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/input/input_poly_coef_test_file.txt";

            if (!fileExists(polyFilePath)) {
                std::cerr << "Error: Density data file not found at " << polyFilePath << std::endl;
                continue;
            }

            std::vector<std::string> polyFiles = {polyFilePath};
            ComaModelFileProcessor processor(polyFiles);
            std::vector<double> stokesRadii = {5000, dsmc_file.cometDistance * 1000.0};
            std::vector<double> stokesLongitudes = {6.0, dsmc_file.solarLongitude};
            ComaStokesDataset stokesDataset = processor.createSHDataset(stokesRadii, stokesLongitudes);
            const double molecularWeight = 0.018015; // kg/mol for water vapor

            auto comaModel = std::make_shared<ComaModel>(
                stokesDataset, molecularWeight, sunStateFunction_,
                cometStateFunction_, cometRotationFunction_, 10, 10);

            // Create ComaWindModel with the polynomial datasets
            std::cout << "Creating ComaWindModel..." << std::endl;

            // Extract datasets from the collection using getter methods
            const ComaPolyDataset& xDataset = windDatasets.getXPolyDataset();
            const ComaPolyDataset& yDataset = windDatasets.getYPolyDataset();
            const ComaPolyDataset& zDataset = windDatasets.getZPolyDataset();

            ComaWindModel windModel(
                xDataset, yDataset, zDataset,
                comaModel,
                sunStateFunction_, cometStateFunction_, cometRotationFunction_,
                10, 10
            );

            // Compute wind velocities for each DSMC data point
            std::cout << "Computing wind velocities for " << dsmc_file.dataPoints.size() << " points..." << std::endl;

            for (auto& point : dsmc_file.dataPoints) {
                // Convert coordinates
                point.latitude = 90.0 - point.coLatitude; // co-latitude to latitude

                double radius = dsmc_file.cometDistance * 1000.0; // km to m
                double longitude = point.longitude * mathematical_constants::PI / 180.0; // deg to rad
                double latitude = point.latitude * mathematical_constants::PI / 180.0; // deg to rad

                // Get wind velocity vector from the wind model
                // Note: currentAltitude parameter is actually radius from center
                Eigen::Vector3d windVelocity = windModel.getCurrentBodyFixedCartesianWindVelocity(
                    radius, longitude, latitude, testTime_);

                // Extract z-component (radial outward) for comparison with DSMC wind magnitude
                // The ComaWindModel returns wind in modified vertical frame where z is radial outward
                point.myWindZ = windVelocity(2);
            }

            // Save results to CSV
            saveComparisonResults(dsmc_file);
        }
    }

private:
    bool parseDSMCFile(const std::string& filePath, DSMCFile& dsmc_file) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        bool headerParsed = false;
        bool dataSection = false;

        while (std::getline(file, line)) {
            if (!headerParsed && line.find("sub-solar-lon") != std::string::npos) {
                // Parse header line to get solar longitude and distance
                if (!std::getline(file, line)) return false;

                std::istringstream iss(line);
                std::string date, helioDistStr, subSolarLatStr, subSolarLonStr, cometoDistStr;

                iss >> date >> helioDistStr >> subSolarLatStr >> subSolarLonStr >> cometoDistStr;

                try {
                    dsmc_file.solarLongitude = std::stod(subSolarLonStr);
                    dsmc_file.cometDistance = std::stod(cometoDistStr);
                    headerParsed = true;
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing header: " << e.what() << std::endl;
                    return false;
                }
                continue;
            }

            // Check for start of data section
            if (line.find("-----") != std::string::npos) {
                dataSection = true;
                continue;
            }

            // Parse data lines
            if (dataSection && headerParsed) {
                std::istringstream iss(line);
                DSMCDataPoint point;

                // Read all four columns: longitude, co-latitude, density, wind speed
                if (iss >> point.longitude >> point.coLatitude >> point.density >> point.windSpeed) {
                    dsmc_file.dataPoints.push_back(point);
                }
            }
        }

        return headerParsed && !dsmc_file.dataPoints.empty();
    }

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

    void setSolarLongitude(double solarLongitudeDeg) {
        // Convert degrees to radians
        double solarLongitudeRad = solarLongitudeDeg * mathematical_constants::PI / 180.0;

        // Create rotation matrix for the solar longitude
        cometRotationFunction_ = [solarLongitudeRad]() -> Eigen::Matrix3d {
            Eigen::Matrix3d rotation;
            rotation << std::cos(solarLongitudeRad),  std::sin(solarLongitudeRad), 0,
                       -std::sin(solarLongitudeRad),  std::cos(solarLongitudeRad), 0,
                       0,                            0,                           1;
            return rotation;
        };
    }

    bool fileExists(const std::string& filename) {
        std::ifstream file(filename);
        return file.good();
    }

    void saveComparisonResults(const DSMCFile& dsmc_file) {
        // Create output filename based on input filename
        std::string baseFilename = dsmc_file.filename;
        size_t lastSlash = baseFilename.find_last_of("/\\");
        if (lastSlash != std::string::npos) {
            baseFilename = baseFilename.substr(lastSlash + 1);
        }

        // Remove extension
        size_t lastDot = baseFilename.find_last_of(".");
        if (lastDot != std::string::npos) {
            baseFilename = baseFilename.substr(0, lastDot);
        }

        std::string outputPath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_Wind_comparison/output/" + baseFilename + "_wind_comparison.csv";

        // Create output directory
        std::string mkdirCmd = "mkdir -p /Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_Wind_comparison/output";
        system(mkdirCmd.c_str());

        std::ofstream outputFile(outputPath);
        if (!outputFile.is_open()) {
            std::cerr << "Failed to create output file: " << outputPath << std::endl;
            return;
        }

        // Write header
        outputFile << "# Solar longitude: " << dsmc_file.solarLongitude << " degrees\n";
        outputFile << "# Comet distance: " << dsmc_file.cometDistance << " km\n";
        outputFile << "longitude_deg,latitude_deg,dsmc_wind_m_s,computed_wind_z_m_s,difference_m_s,relative_error_percent\n";

        // Write data
        for (const auto& point : dsmc_file.dataPoints) {
            double diff = point.myWindZ - point.windSpeed;
            double relativeError = (point.windSpeed != 0.0) ? (diff / point.windSpeed) * 100.0 : 0.0;

            outputFile << std::fixed << std::setprecision(6)
                      << point.longitude << ","
                      << point.latitude << ","
                      << std::scientific << std::setprecision(6) << point.windSpeed << ","
                      << point.myWindZ << ","
                      << diff << ","
                      << std::fixed << std::setprecision(2) << relativeError << "\n";
        }

        outputFile.close();
        std::cout << "Results saved to: " << outputPath << std::endl;

        // Print some statistics
        double totalDiff = 0.0;
        double totalAbsDiff = 0.0;
        double maxAbsDiff = 0.0;
        double totalRelError = 0.0;
        double maxRelError = 0.0;

        for (const auto& point : dsmc_file.dataPoints) {
            double diff = point.myWindZ - point.windSpeed;
            double absDiff = std::abs(diff);
            double relError = (point.windSpeed != 0.0) ? std::abs(diff / point.windSpeed) * 100.0 : 0.0;

            totalDiff += diff;
            totalAbsDiff += absDiff;
            maxAbsDiff = std::max(maxAbsDiff, absDiff);
            totalRelError += relError;
            maxRelError = std::max(maxRelError, relError);
        }

        size_t numPoints = dsmc_file.dataPoints.size();
        std::cout << "Statistics for " << numPoints << " points:" << std::endl;
        std::cout << "  Average difference: " << std::scientific << (totalDiff / numPoints) << " m/s" << std::endl;
        std::cout << "  Average absolute difference: " << std::scientific << (totalAbsDiff / numPoints) << " m/s" << std::endl;
        std::cout << "  Maximum absolute difference: " << std::scientific << maxAbsDiff << " m/s" << std::endl;
        std::cout << "  Average relative error: " << std::fixed << std::setprecision(2) << (totalRelError / numPoints) << "%" << std::endl;
        std::cout << "  Maximum relative error: " << std::fixed << std::setprecision(2) << maxRelError << "%" << std::endl;
    }
};

int main() {
    try {
        std::cout << "DSMC Wind Comparison Tool for Coma Wind Model" << std::endl;
        std::cout << "==============================================" << std::endl;

        DSMCWindComparison comparison;

        // Load DSMC input files
        std::vector<std::string> dsmc_files = {
            "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_Wind_comparison/input/r_cometFixed_ep10-000_4km.txt",
            "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_Wind_comparison/input/r_cometFixed_ep10-030_10km.txt"
        };

        if (!comparison.loadDSMCFiles(dsmc_files)) {
            std::cerr << "Failed to load DSMC files" << std::endl;
            return 1;
        }

        // Run the comparison
        comparison.runComparison();

        std::cout << "\nDSMC wind comparison completed successfully!" << std::endl;
        std::cout << "Results saved in ./DSMC_Wind_comparison/output/ directory" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
