#include <iostream>
#include <iomanip>
#include <memory>
#include <functional>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

using namespace tudat;
using namespace tudat::aerodynamics;
using namespace tudat::simulation_setup;

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "SOLAR LONGITUDE DEBUG INVESTIGATION" << std::endl;
    std::cout << "========================================" << std::endl;

    // Test parameters
    const double solarLongitudeDeg = 100.0;  // Test at 90 degrees
    const double solarLongitudeRad = solarLongitudeDeg * mathematical_constants::PI / 180.0;

    const double testRadius = 6000.0;  // 4 km (matching performance_analysis first radius)
    const double testBodyLongitude = 0.0;  // 0 degrees in body frame
    const double testBodyLatitude = 0.0;   // 0 degrees
    const double testTime = 490708800.0;   // s since J2000

    std::cout << "\nTest Configuration:" << std::endl;
    std::cout << "  Solar longitude: " << solarLongitudeDeg << " degrees (" << solarLongitudeRad << " rad)" << std::endl;
    std::cout << "  Test point: r=" << testRadius << "m, body_lon=" << testBodyLongitude << "deg, body_lat=" << testBodyLatitude << "deg" << std::endl;
    std::cout << std::endl;

    // Setup mock state functions
    std::function<Eigen::Vector6d()> sunStateFunction = []() -> Eigen::Vector6d {
        Eigen::Vector6d sunState = Eigen::Vector6d::Zero();
        sunState[0] = 1.5e11; // 1 AU in meters
        return sunState;
    };

    std::function<Eigen::Vector6d()> cometStateFunction = []() -> Eigen::Vector6d {
        return Eigen::Vector6d::Zero();
    };

    // Create rotation function for the test solar longitude
    std::function<Eigen::Matrix3d()> cometRotationFunction = [solarLongitudeRad]() -> Eigen::Matrix3d {
        Eigen::Matrix3d rotation;
        rotation << std::cos(solarLongitudeRad),  std::sin(solarLongitudeRad), 0,
                   -std::sin(solarLongitudeRad),  std::cos(solarLongitudeRad), 0,
                    0,                            0,                           1;
        return rotation;
    };

    std::cout << "Rotation matrix for solar longitude " << solarLongitudeDeg << "°:" << std::endl;
    std::cout << cometRotationFunction() << std::endl << std::endl;

    // Load polynomial coefficient file
    std::string polyFilePath = "/Users/markusreichel/PhD/tudatpy/examples/tudat/tudat/applications/coma_analysis/input/input_poly_coef_test_file.txt";
    std::vector<std::string> polyFiles = {polyFilePath};

    ComaModelFileProcessor processor(polyFiles);
    ComaPolyDataset polyDataset = processor.createPolyCoefDataset();

    std::cout << "Creating Stokes dataset with full solar longitude grid (0° to 360° in 10° steps)..." << std::endl;
    // Create a grid matching performance_analysis: 4km to 20km
    std::vector<double> stokesRadii;
    for (int i = 0; i < 20; ++i) {
        stokesRadii.push_back(4000.0 + i * 1000.0);  // 4km, 6km, 8km, ..., 22km
    }
    std::vector<double> stokesLongitudes;
    for (int i = 0; i <= 36; ++i) {
        stokesLongitudes.push_back(i * 10.0);
    }

    ComaStokesDataset stokesDataset = processor.createSHDataset(stokesRadii, stokesLongitudes);

    std::cout << "Stokes dataset created with " << stokesDataset.nLongitudes() << " solar longitudes" << std::endl;
    std::cout << "Longitude grid (in radians): ";
    const auto& lons = stokesDataset.lons();
    for (size_t i = 0; i < std::min(size_t(10), lons.size()); ++i) {
        std::cout << lons[i] << " ";
    }
    std::cout << "..." << std::endl << std::endl;

    // Create both ComaModels
    const double molecularWeight = 0.018015; // kg/mol for water vapor
    const int maxDegree = 10;
    const int maxOrder = 10;

    std::cout << "Creating Polynomial ComaModel..." << std::endl;
    ComaModel polyComaModel(polyDataset, molecularWeight, sunStateFunction, cometStateFunction,
                           cometRotationFunction, maxDegree, maxOrder);

    std::cout << "Creating Stokes ComaModel..." << std::endl;
    ComaModel stokesComaModel(stokesDataset, molecularWeight, sunStateFunction, cometStateFunction,
                             cometRotationFunction, maxDegree, maxOrder);

    std::cout << "\n========================================" << std::endl;
    std::cout << "QUERYING DENSITY" << std::endl;
    std::cout << "========================================" << std::endl;

    // Convert test point to radians
    double bodyLonRad = testBodyLongitude * mathematical_constants::PI / 180.0;
    double bodyLatRad = testBodyLatitude * mathematical_constants::PI / 180.0;

    std::cout << "\n>>> POLYNOMIAL METHOD <<<" << std::endl;
    double polyDensity = polyComaModel.getNumberDensity(testRadius, bodyLonRad, bodyLatRad, testTime);

    std::cout << "\n>>> STOKES METHOD <<<" << std::endl;
    double stokesDensity = stokesComaModel.getNumberDensity(testRadius, bodyLonRad, bodyLatRad, testTime);

    std::cout << "\n========================================" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Polynomial density: " << polyDensity << " m^-3" << std::endl;
    std::cout << "Stokes density:     " << stokesDensity << " m^-3" << std::endl;
    std::cout << "Absolute difference: " << std::abs(polyDensity - stokesDensity) << " m^-3" << std::endl;
    std::cout << "Relative difference: " << (std::abs(polyDensity - stokesDensity) / polyDensity * 100.0) << " %" << std::endl;

    return 0;
}
