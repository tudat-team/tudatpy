#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <memory>
#include <functional>

// Tudat includes for the actual ComaModel classes
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;

class GridDensityAnalysis {
private:
    // Fixed parameters
    double radius_;          // meters
    const double testTime_ = 490708800;   // s since J2000

    // Grid parameters
    int numLatitudePoints_;
    int numLongitudePoints_;

    // Mock state functions for ComaModel
    std::function<Eigen::Vector6d()> sunStateFunction_;
    std::function<Eigen::Vector6d()> cometStateFunction_;
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

public:
    GridDensityAnalysis(double radiusMeters,
                       int numLatPoints = 90, int numLonPoints = 180)
        : radius_(radiusMeters),
          numLatitudePoints_(numLatPoints),
          numLongitudePoints_(numLonPoints) {

        setupMockStateFunctions();

        std::cout << "Initialized Grid Density Analysis:" << std::endl;
        std::cout << "  Radius: " << radius_ << " m" << std::endl;
        std::cout << "  Grid: " << numLatitudePoints_ << " x " << numLongitudePoints_
                  << " = " << (numLatitudePoints_ * numLongitudePoints_) << " points" << std::endl;
    }

    void runAnalysis(const std::vector<double>& solarLongitudes) {
        std::cout << "\n=== Running Grid Density Analysis for " << solarLongitudes.size()
                  << " solar longitudes ===" << std::endl;

        // Test data file path
        std::string polyFilePath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/input/input_poly_coef_test_file.txt";

        if (!fileExists(polyFilePath)) {
            std::cerr << "Error: Test data file not found at " << polyFilePath << std::endl;
            return;
        }

        std::cout << "Using poly coefficient file: " << polyFilePath << std::endl;

        // Create processor
        std::vector<std::string> polyFiles = {polyFilePath};
        ComaModelFileProcessor processor(polyFiles);

        // Create latitude and longitude grids (fixed for all solar longitudes)
        // Avoid exact poles (±90°) where Legendre polynomial derivatives have singularities
        std::vector<double> latitudes = createLinearGrid(-89.9, 89.9, numLatitudePoints_);
        std::vector<double> longitudes = createLinearGrid(1.0, 360.0, numLongitudePoints_);

        // Loop through each solar longitude
        for (const auto& solarLongitude : solarLongitudes) {
            std::cout << "\n--- Processing Solar Longitude: " << solarLongitude << "° ---" << std::endl;

            // Set solar longitude for rotation
            setSolarLongitude(solarLongitude);

            // Create grid for Stokes dataset
            std::vector<double> stokesRadii = {radius_,5000};
            std::vector<double> stokesLongitudes = {solarLongitude,5};

            std::cout << "Creating Stokes dataset..." << std::endl;
            ComaStokesDataset stokesDataset = processor.createSHDataset(stokesRadii, stokesLongitudes);

            // Create ComaModel with Stokes dataset
            // Using molecular weight for H2O: 18.015 g/mol = 0.018015 kg/mol
            const double molecularWeight = 0.018015; // kg/mol for water vapor
            ComaModel stokesComaModel(stokesDataset, molecularWeight, sunStateFunction_, cometStateFunction_,
                                     cometRotationFunction_, 10, 10);

            std::cout << "Computing densities for " << (numLatitudePoints_ * numLongitudePoints_)
                      << " grid points..." << std::endl;

            // Compute densities and save to CSV
            std::string outputPath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/grid_density_analysis/output/grid_density_r"
                                   + std::to_string(static_cast<int>(radius_)) + "m_ls"
                                   + std::to_string(static_cast<int>(solarLongitude)) + "deg.csv";

            std::ofstream outputFile(outputPath);
            if (!outputFile.is_open()) {
                std::cerr << "Failed to create output file: " << outputPath << std::endl;
                continue;
            }

            // Write header
            outputFile << "# Radius: " << radius_ << " m\n";
            outputFile << "# Solar Longitude: " << solarLongitude << " degrees\n";
            outputFile << "latitude,longitude,density_log2,density_linear\n";

            // Compute and write data
            int pointsComputed = 0;
            for (const auto& lat : latitudes) {
                for (const auto& lon : longitudes) {
                    // Convert to radians for computation
                    double latRad = lat * mathematical_constants::PI / 180.0;
                    double lonRad = lon * mathematical_constants::PI / 180.0;

                    // Get density (returns actual number density)
                    double densityLinear = stokesComaModel.getNumberDensity(radius_, lonRad, latRad, testTime_);
                    double densityLog2 = std::log2(densityLinear);

                    // Write to file
                    outputFile << std::fixed << std::setprecision(6)
                              << lat << ","
                              << lon << ","
                              << densityLog2 << ","
                              << std::scientific << std::setprecision(6) << densityLinear << "\n";

                    pointsComputed++;
                    if (pointsComputed % 1000 == 0) {
                        std::cout << "  Computed " << pointsComputed << " / "
                                  << (numLatitudePoints_ * numLongitudePoints_) << " points" << std::endl;
                    }
                }
            }

            outputFile.close();
            std::cout << "Results saved to: " << outputPath << std::endl;

            // Print some statistics
            printStatistics(outputPath);
        }
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

    void printStatistics(const std::string& filePath) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            return;
        }

        std::string line;
        // Skip header lines
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);

        double minLog2 = std::numeric_limits<double>::max();
        double maxLog2 = std::numeric_limits<double>::lowest();
        double minLinear = std::numeric_limits<double>::max();
        double maxLinear = std::numeric_limits<double>::lowest();
        double sumLog2 = 0.0;
        double sumLinear = 0.0;
        int count = 0;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string lat, lon, log2, linear;

            if (std::getline(iss, lat, ',') &&
                std::getline(iss, lon, ',') &&
                std::getline(iss, log2, ',') &&
                std::getline(iss, linear)) {

                double densLog2 = std::stod(log2);
                double densLinear = std::stod(linear);

                minLog2 = std::min(minLog2, densLog2);
                maxLog2 = std::max(maxLog2, densLog2);
                minLinear = std::min(minLinear, densLinear);
                maxLinear = std::max(maxLinear, densLinear);
                sumLog2 += densLog2;
                sumLinear += densLinear;
                count++;
            }
        }

        file.close();

        std::cout << "\n=== Density Statistics ===" << std::endl;
        std::cout << "Points analyzed: " << count << std::endl;
        std::cout << "\nLog2 scale:" << std::endl;
        std::cout << "  Min: " << std::fixed << std::setprecision(6) << minLog2 << std::endl;
        std::cout << "  Max: " << maxLog2 << std::endl;
        std::cout << "  Mean: " << (sumLog2 / count) << std::endl;
        std::cout << "\nLinear scale:" << std::endl;
        std::cout << "  Min: " << std::scientific << std::setprecision(6) << minLinear << " m^-3" << std::endl;
        std::cout << "  Max: " << maxLinear << " m^-3" << std::endl;
        std::cout << "  Mean: " << (sumLinear / count) << " m^-3" << std::endl;
    }
};

int main() {
    try {
        std::cout << "Grid Density Analysis Tool for Coma Model" << std::endl;
        std::cout << "==========================================" << std::endl;

        // Define parameters for the analysis
        double radius = 10000.0;        // 10 km in meters
        int numLatPoints = 180;          // Number of latitude points
        int numLonPoints = 360;         // Number of longitude points

        // Define array of solar longitudes to analyze
        std::vector<double> solarLongitudes = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0};

        std::cout << "Will analyze " << solarLongitudes.size() << " solar longitudes: ";
        for (size_t i = 0; i < solarLongitudes.size(); ++i) {
            std::cout << solarLongitudes[i] << "°";
            if (i < solarLongitudes.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;

        GridDensityAnalysis analysis(radius, numLatPoints, numLonPoints);
        analysis.runAnalysis(solarLongitudes);

        std::cout << "\n=== Grid density analysis completed successfully! ===" << std::endl;
        std::cout << "Results saved in grid_density_analysis/output/ directory" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
