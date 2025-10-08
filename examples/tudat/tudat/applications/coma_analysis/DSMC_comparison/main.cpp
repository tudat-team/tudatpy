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

// Tudat includes for the actual ComaModel classes
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;

struct DSMCDataPoint {
    double longitude;    // degrees
    double coLatitude;   // degrees (0 = north, 180 = south)
    double density;      // number density [m-3]
    double gasSpeed;     // [m/s]

    // Computed values
    double latitude;     // converted from co-latitude
    double myDensity;    // density from our algorithm
};

struct DSMCFile {
    std::string filename;
    double solarLongitude;  // degrees
    double cometDistance;   // km
    std::vector<DSMCDataPoint> dataPoints;
};

class DSMCComparison {
private:
    std::vector<DSMCFile> dsmc_files_;

    // Mock state functions for ComaModel
    std::function<Eigen::Vector6d()> sunStateFunction_;
    std::function<Eigen::Vector6d()> cometStateFunction_;
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

    const double testTime_ = 490708800;   // s since J2000

public:
    DSMCComparison() {
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

            // Create Stokes dataset using the same approach as performance_analysis
            std::string polyFilePath = "/Users/markusreichel/PhD/tudatpy/tests/test_tudat/src/astro/aerodynamics/test_data/input_poly_coef_test_file.txt";

            if (!fileExists(polyFilePath)) {
                std::cerr << "Error: Test data file not found at " << polyFilePath << std::endl;
                continue;
            }

            // Create processor and Stokes dataset
            std::vector<std::string> polyFiles = {polyFilePath};
            ComaModelFileProcessor processor(polyFiles);

            // Create grid for Stokes dataset - use a reasonable range around the DSMC data points
            std::vector<double> stokesRadii = {dsmc_file.cometDistance * 1000.0,5000}; // Convert km to m
            std::vector<double> stokesLongitudes = {dsmc_file.solarLongitude, 10};   // in degree

            std::cout << "Creating Stokes dataset..." << std::endl;
            ComaStokesDataset stokesDataset = processor.createSHDataset(stokesRadii, stokesLongitudes);

            // Create ComaModel with Stokes dataset
            // Using molecular weight for H2O: 18.015 g/mol = 0.018015 kg/mol
            const double molecularWeight = 0.018015; // kg/mol for water vapor
            ComaModel stokesComaModel(stokesDataset, molecularWeight, sunStateFunction_, cometStateFunction_,
                                     cometRotationFunction_, 10, 10);

            // Compute densities for each DSMC data point
            std::cout << "Computing densities for " << dsmc_file.dataPoints.size() << " points..." << std::endl;

            for (auto& point : dsmc_file.dataPoints) {
                // Convert coordinates
                point.latitude = 90.0 - point.coLatitude; // co-latitude to latitude

                double radius = dsmc_file.cometDistance * 1000.0; // km to m
                double longitude = point.longitude * mathematical_constants::PI / 180.0; // deg to rad
                double latitude = point.latitude * mathematical_constants::PI / 180.0; // deg to rad

                // Get number density from our algorithm
                point.myDensity = stokesComaModel.getNumberDensity(radius, longitude, latitude, testTime_);
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

                if (iss >> point.longitude >> point.coLatitude >> point.density >> point.gasSpeed) {
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

        std::string outputPath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_comparison/output/" + baseFilename + "_comparison.csv";

        // Create output directory
        std::string mkdirCmd = "mkdir -p /Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_comparison/output";
        system(mkdirCmd.c_str());

        std::ofstream outputFile(outputPath);
        if (!outputFile.is_open()) {
            std::cerr << "Failed to create output file: " << outputPath << std::endl;
            return;
        }

        // Write header
        outputFile << "# Solar longitude: " << dsmc_file.solarLongitude << " degrees\n";
        outputFile << "# Comet distance: " << dsmc_file.cometDistance << " km\n";
        outputFile << "longitude_deg,latitude_deg,dsmc_density_linear,my_density_linear,difference_linear,";
        outputFile << "dsmc_density_log2,my_density_log2,difference_log2\n";

        // Write data
        for (const auto& point : dsmc_file.dataPoints) {
            double dsmc_log2 = std::log2(point.density);
            double my_linear = std::exp2(point.myDensity); // Convert log2 back to linear
            double diff_log2 = point.myDensity - dsmc_log2;
            double diff_linear = my_linear - point.density;

            outputFile << std::fixed << std::setprecision(6)
                      << point.longitude << ","
                      << point.latitude << ","
                      << std::scientific << std::setprecision(6) << point.density << ","
                      << my_linear << ","
                      << diff_linear << ","
                      << std::fixed << std::setprecision(6) << dsmc_log2 << ","
                      << point.myDensity << ","
                      << diff_log2 << "\n";
        }

        outputFile.close();
        std::cout << "Results saved to: " << outputPath << std::endl;

        // Print some statistics
        double totalDiffLog2 = 0.0;
        double totalDiffLinear = 0.0;
        double totalAbsDiffLog2 = 0.0;
        double totalAbsDiffLinear = 0.0;
        double maxAbsDiffLog2 = 0.0;
        double maxAbsDiffLinear = 0.0;

        for (const auto& point : dsmc_file.dataPoints) {
            double dsmc_log2 = std::log2(point.density);
            double my_linear = std::exp2(point.myDensity);
            double diff_log2 = point.myDensity - dsmc_log2;
            double diff_linear = my_linear - point.density;
            double abs_diff_log2 = std::abs(diff_log2);
            double abs_diff_linear = std::abs(diff_linear);

            totalDiffLog2 += diff_log2;
            totalDiffLinear += diff_linear;
            totalAbsDiffLog2 += abs_diff_log2;
            totalAbsDiffLinear += abs_diff_linear;
            maxAbsDiffLog2 = std::max(maxAbsDiffLog2, abs_diff_log2);
            maxAbsDiffLinear = std::max(maxAbsDiffLinear, abs_diff_linear);
        }

        size_t numPoints = dsmc_file.dataPoints.size();
        std::cout << "Statistics for " << numPoints << " points:" << std::endl;
        std::cout << "  Average difference (log2): " << (totalDiffLog2 / numPoints) << std::endl;
        std::cout << "  Average absolute difference (log2): " << (totalAbsDiffLog2 / numPoints) << std::endl;
        std::cout << "  Maximum absolute difference (log2): " << maxAbsDiffLog2 << std::endl;
        std::cout << "  Average difference (linear): " << std::scientific << (totalDiffLinear / numPoints) << std::endl;
        std::cout << "  Average absolute difference (linear): " << std::scientific << (totalAbsDiffLinear / numPoints) << std::endl;
        std::cout << "  Maximum absolute difference (linear): " << std::scientific << maxAbsDiffLinear << std::endl;
    }
};

int main() {
    try {
        std::cout << "DSMC Comparison Tool for Coma Model" << std::endl;
        std::cout << "====================================" << std::endl;

        DSMCComparison comparison;

        // Load DSMC input files
        std::vector<std::string> dsmc_files = {
            "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_comparison/input/r_cometFixed_ep10-000_4km.txt",
            "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/DSMC_comparison/input/r_cometFixed_ep10-030_10km.txt"
        };

        if (!comparison.loadDSMCFiles(dsmc_files)) {
            std::cerr << "Failed to load DSMC files" << std::endl;
            return 1;
        }

        // Run the comparison
        comparison.runComparison();

        std::cout << "\nDSMC comparison completed successfully!" << std::endl;
        std::cout << "Results saved in ./DSMC_comparison/output/ directory" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}