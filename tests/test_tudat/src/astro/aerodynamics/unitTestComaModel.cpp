#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/astro/aerodynamics/comaModel.h"
#include <boost/test/unit_test.hpp>


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_coma_settings )

BOOST_AUTO_TEST_CASE( testComaSettingsSingleFile )
{
    using namespace std;
    using namespace tudat;

    // Create a test file with a minimal valid coma model input
    boost::filesystem::path cwd = boost::filesystem::current_path();
    const boost::filesystem::path fullPath =
            "/Users/markusreichel/PhD/tudatpy/tests/test_tudat/src/astro/aerodynamics/test_data/test_input_coma.txt";

    BOOST_REQUIRE( boost::filesystem::exists(fullPath) );

    // Construct ComaSettings object
    const std::vector< std::string > fileList = { fullPath.string( ) };
    const simulation_setup::ComaSettings comaSettings( fileList, 10, 10 );

    const simulation_setup::ComaPolyDataset polyData = comaSettings.getPolyDataset();

    BOOST_REQUIRE_EQUAL( polyData.getNumFiles( ), 1 );

    // Check polyCoefs
    // Rows and Columns are flipped internally
    const Eigen::MatrixXd polyCoefs = polyData.getPolyCoefficients( 0 );
    BOOST_CHECK_EQUAL( polyCoefs.rows(), 48 );
    BOOST_CHECK_EQUAL( polyCoefs.cols(), 121 );

    BOOST_CHECK_CLOSE( polyCoefs(0, 0), 6.262302500423528E+02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(3, 1), -4.459813047951577E-02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(47, 120), -1.049515320967208E+02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(10, 22), 1.287417812579956E-01, 1.0e-12 );



     // Check SHDegreeAndOrderIndices
     const Eigen::ArrayXXi indices = polyData.getSHDegreeAndOrderIndices( 0 );
     BOOST_REQUIRE_EQUAL( indices.rows(), 2 );
     BOOST_REQUIRE_EQUAL( indices.cols(), 121 );
     BOOST_REQUIRE_EQUAL( indices(0, 120), 10 ); // First row, last column
     BOOST_REQUIRE_EQUAL( indices(1, 120), -10 ); // Second row, last column

     // Check reference radius
     const double refRadius = polyData.getReferenceRadius( 0 );
     BOOST_CHECK_EQUAL( refRadius, 10.0 );

     // Check power inv radius
     Eigen::VectorXd pwInv = polyData.getPowersInvRadius( 0 );
     BOOST_CHECK_EQUAL( pwInv.rows(), 4 );
     BOOST_CHECK_EQUAL( pwInv[0], 0 );
     BOOST_CHECK_EQUAL( pwInv[3], 3 );

     // Check max degree and order
     BOOST_CHECK_EQUAL( comaSettings.getRequestedDegree(  ), 10 );
     BOOST_CHECK_EQUAL( comaSettings.getRequestedOrder(  ), 10 );
//
//     Eigen::MatrixXd cosineCoefficients, sineCoefficients;
//
//     double distanceToCometCentre = 6.0;
//     double solarLongitude = 30.0 * M_PI / 180.0;
//
//     // Test computation of stokes coefficients
//     aerodynamics::ComaModel::testEvaluateStokesCoefficients2D(
//             distanceToCometCentre,
//             solarLongitude,
//             polyCoefs[ 0 ],
//             SHDegreeAndOrder[ 0 ],
//             powerInvRadius[ 0 ],
//             referenceRadius[ 0 ],
//             cosineCoefficients,
//             sineCoefficients,
//             maxDegree,
//             maxOrder );
//
//     // Check output size
//     BOOST_CHECK_EQUAL( cosineCoefficients.rows(), maxDegree + 1 );
//     BOOST_CHECK_EQUAL( cosineCoefficients.cols(), maxOrder + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.rows(), maxDegree + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.cols(), maxOrder + 1 );
//
//     // Check cosine coefficients
//     BOOST_CHECK_CLOSE( cosineCoefficients(0, 0), 5.393369500951372e+01, 1.0e-12 );
//     BOOST_CHECK_CLOSE( cosineCoefficients(3, 1), 1.054997670055739e-01, 1.0e-12 );
//     BOOST_CHECK_CLOSE( cosineCoefficients(5, 4), 1.207229799736584e-02, 1.0e-12 );
//     BOOST_CHECK_CLOSE( cosineCoefficients(9, 3), -1.845804426625799e-03, 1.0e-12 );
//
//     // Check sine coefficients
//     BOOST_CHECK_CLOSE( sineCoefficients(0, 0), 0.0, 1.0e-12 );
//     BOOST_CHECK_CLOSE( sineCoefficients(6, 2), 2.139319739882320e-02, 1.0e-12 );
//     BOOST_CHECK_CLOSE( sineCoefficients(7, 5), -5.401766866307728e-02, 1.0e-12 );
//     BOOST_CHECK_CLOSE( sineCoefficients(10, 8), 1.423848196013654e-02, 1.0e-12 );
//
//     // Test computation of stoeks coefficients
//     aerodynamics::ComaModel::testEvaluateStokesCoefficients2D(
//             distanceToCometCentre,
//             solarLongitude,
//             polyCoefs[ 0 ],
//             SHDegreeAndOrder[ 0 ],
//             powerInvRadius[ 0 ],
//             referenceRadius[ 0 ],
//             cosineCoefficients,
//             sineCoefficients );
//
//     // Check output size. Should default to max. degree and order
//     BOOST_CHECK_EQUAL( cosineCoefficients.rows(), maxDegree + 1 );
//     BOOST_CHECK_EQUAL( cosineCoefficients.cols(), maxOrder + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.rows(), maxDegree + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.cols(), maxOrder + 1 );
//
//     // Test computation of stokes coefficients
//     aerodynamics::ComaModel::testEvaluateStokesCoefficients2D(
//             distanceToCometCentre,
//             solarLongitude,
//             polyCoefs[ 0 ],
//             SHDegreeAndOrder[ 0 ],
//             powerInvRadius[ 0 ],
//             referenceRadius[ 0 ],
//             cosineCoefficients,
//             sineCoefficients,
//             8,
//             3 );
//
//     // Check output size. Should default to max. degree and order
//     BOOST_CHECK_EQUAL( cosineCoefficients.rows(), 8 + 1 );
//     BOOST_CHECK_EQUAL( cosineCoefficients.cols(), 3 + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.rows(), 8 + 1 );
//     BOOST_CHECK_EQUAL( sineCoefficients.cols(), 3 + 1 );
//
//     // Check that exceeding maxDegree or maxOrder throws
//     BOOST_CHECK_THROW(
//             aerodynamics::ComaModel::testEvaluateStokesCoefficients2D(
//                 distanceToCometCentre,
//                 solarLongitude,
//                 polyCoefs[0],
//                 SHDegreeAndOrder[0],
//                 powerInvRadius[0],
//                 referenceRadius[0],
//                 cosineCoefficients,
//                 sineCoefficients,
//                 maxDegree + 5, // Invalid: 15 > 10
//                 maxOrder + 3 // Invalid: 13 > 10
//             ),
//             std::runtime_error
//             );
//
//     // Verify density
//     distanceToCometCentre = 10.0; // [km]
//     solarLongitude = 30 * M_PI / 180.0;;
//     double latitude = 94 * M_PI / 180.0;
//     double longitude = -19.8 * M_PI / 180.0;
//     // Test computation of stokes coefficients
//     aerodynamics::ComaModel::testEvaluateStokesCoefficients2D(
//             distanceToCometCentre,
//             solarLongitude,
//             polyCoefs[ 0 ],
//             SHDegreeAndOrder[ 0 ],
//             powerInvRadius[ 0 ],
//             referenceRadius[ 0 ],
//             cosineCoefficients,
//             sineCoefficients,
//             10,
//             10 );
//
//     simulation_setup::SphericalHarmonicsDensity SH( sineCoefficients,
//                                                     cosineCoefficients );
//     double density = SH.calculateSurfaceSphericalHarmonics( sineCoefficients,
//                                                             cosineCoefficients,
//                                                             latitude,
//                                                             longitude,
//                                                             10,
//                                                             10 );
//
//     std::cout << "denstiy: " << density << std::endl;
//     std::cout << "Density: " << std::pow( 2.0, density ) << std::endl;
//
    // BOOST_CHECK_CLOSE( density, 1.423848196013654e-02, 1.0e-12 );
}

BOOST_AUTO_TEST_SUITE_END( )
} // namespace unit_tests
} // namespace tudat