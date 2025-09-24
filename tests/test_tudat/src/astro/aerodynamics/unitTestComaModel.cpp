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
    boost::filesystem::path cwd = boost::filesystem::current_path( );
    const boost::filesystem::path fullPath =
            "/Users/markusreichel/PhD/tudatpy/tests/test_tudat/src/astro/aerodynamics/test_data/test_input_coma.txt";

    BOOST_REQUIRE( boost::filesystem::exists(fullPath) );

    // Construct ComaSettings object
    const std::vector< std::string > fileList = { fullPath.string( ) };
    const simulation_setup::ComaSettings comaSettings( fileList, 10, 10 );

    const simulation_setup::ComaPolyDataset polyData = comaSettings.getPolyDataset( );

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
    BOOST_CHECK_EQUAL( comaSettings.getRequestedDegree( ), 10 );
    BOOST_CHECK_EQUAL( comaSettings.getRequestedOrder( ), 10 );
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


BOOST_AUTO_TEST_CASE( testPolyProcessorAndSHDataset )
{
    using namespace tudat;

    // Input
    const boost::filesystem::path fullPath =
            "/Users/markusreichel/PhD/tudatpy/tests/test_tudat/src/astro/aerodynamics/test_data/test_input_coma.txt";
    BOOST_REQUIRE( boost::filesystem::exists(fullPath) );
    const std::vector< std::string > files = { fullPath.string( ) };

    // Processor
    simulation_setup::PolyCoefFileProcessing proc( files );

    // Build poly dataset via processor
    simulation_setup::ComaPolyDataset poly = proc.createPolyCoefDataset( );
    BOOST_CHECK_EQUAL( poly.getNumFiles(), files.size() );
    BOOST_CHECK_EQUAL( poly.getFileMeta(0).sourcePath, files[0] );

    // Check polyCoefs
    // Rows and Columns are flipped internally
    const Eigen::MatrixXd polyCoefs = poly.getPolyCoefficients( 0 );
    BOOST_CHECK_EQUAL( polyCoefs.rows(), 48 );
    BOOST_CHECK_EQUAL( polyCoefs.cols(), 121 );

    BOOST_CHECK_CLOSE( polyCoefs(0, 0), 6.262302500423528E+02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(3, 1), -4.459813047951577E-02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(47, 120), -1.049515320967208E+02, 1.0e-12 );
    BOOST_CHECK_CLOSE( polyCoefs(10, 22), 1.287417812579956E-01, 1.0e-12 );

    // Check SHDegreeAndOrderIndices
    const Eigen::ArrayXXi indices = poly.getSHDegreeAndOrderIndices( 0 );
    BOOST_REQUIRE_EQUAL( indices.rows(), 2 );
    BOOST_REQUIRE_EQUAL( indices.cols(), 121 );
    BOOST_REQUIRE_EQUAL( indices(0, 120), 10 ); // First row, last column
    BOOST_REQUIRE_EQUAL( indices(1, 120), -10 ); // Second row, last column

    // Check reference radius
    const double refRadius = poly.getReferenceRadius( 0 );
    BOOST_CHECK_EQUAL( refRadius, 10.0 );

    // Check power inv radius
    Eigen::VectorXd pwInv = poly.getPowersInvRadius( 0 );
    BOOST_CHECK_EQUAL( pwInv.rows(), 4 );
    BOOST_CHECK_EQUAL( pwInv[0], 0 );
    BOOST_CHECK_EQUAL( pwInv[3], 3 );

    // Build SH dataset

    const int maxDegree = 10;
    const int maxOrder = 10;
    const std::vector< double > radii_m = { 6000.0, 10000.0 }; // meters
    const std::vector< double > lons_deg = { 0.0, 30.0 }; // degrees

    simulation_setup::ComaStokesDataset sh = proc.createSHDataset( radii_m,
                                                                   lons_deg,
                                                                   maxDegree,
                                                                   maxOrder );

    // Shapes/metadata
    BOOST_CHECK_EQUAL( sh.nFiles(), files.size() );
    BOOST_CHECK_EQUAL( sh.nRadii(), radii_m.size() );
    BOOST_CHECK_EQUAL( sh.nLongitudes(), lons_deg.size() );
    BOOST_CHECK_EQUAL( sh.nCoeffs(), static_cast<std::size_t>((maxDegree+1)*(maxOrder+2)/2) );
    BOOST_CHECK_EQUAL( sh.nmax(), maxOrder );
    BOOST_CHECK_EQUAL( sh.files()[0].source_tag, files[0] );

    const double rel_tol = 1.0e-12;

    // 1) Check the scalar pair you specified explicitly:
    {
        auto CS = sh.getCoeff( 0, /*ri*/0, /*li*/1, /*n*/0, /*m*/0 );
        BOOST_CHECK_CLOSE( CS.first, 5.393369500951372e+01, rel_tol );
        BOOST_CHECK_CLOSE( CS.second, 0.0, rel_tol );
    }

    // 2) “Cosine coefficients” at the same grid cell (ri=1 => 6000 m, li=1 => 30 deg):
    {
        auto C_0_0 = sh.getCoeff( 0, 0, 1, 0, 0 ).first; // already checked above
        auto C_3_1 = sh.getCoeff( 0, 0, 1, 3, 1 ).first;
        auto C_5_4 = sh.getCoeff( 0, 0, 1, 5, 4 ).first;
        auto C_9_3 = sh.getCoeff( 0, 0, 1, 9, 3 ).first;

        BOOST_CHECK_CLOSE( C_0_0, 5.393369500951372e+01, rel_tol );
        BOOST_CHECK_CLOSE( C_3_1, 1.054997670055739e-01, rel_tol );
        BOOST_CHECK_CLOSE( C_5_4, 1.207229799736584e-02, rel_tol );
        BOOST_CHECK_CLOSE( C_9_3, -1.845804426625799e-03, rel_tol );
    }

    // 3) “Sine coefficients” at the same grid cell:
    {
        auto S_0_0 = sh.getCoeff( 0, 0, 1, 0, 0 ).second;
        auto S_6_2 = sh.getCoeff( 0, 0, 1, 6, 2 ).second;
        auto S_7_5 = sh.getCoeff( 0, 0, 1, 7, 5 ).second;
        auto S_10_8 = sh.getCoeff( 0, 0, 1, 10, 8 ).second;

        BOOST_CHECK_CLOSE( S_0_0, 0.0, rel_tol );
        BOOST_CHECK_CLOSE( S_6_2, 2.139319739882320e-02, rel_tol );
        BOOST_CHECK_CLOSE( S_7_5, -5.401766866307728e-02, rel_tol );
        BOOST_CHECK_CLOSE( S_10_8, 1.423848196013654e-02, rel_tol );
    }

    auto [ cosineCoefficients, sineCoefficients ] = sh.getCoefficientMatrices( 0, 0, 1 );

    BOOST_CHECK_EQUAL( cosineCoefficients.rows(), maxDegree + 1 );
    BOOST_CHECK_EQUAL( cosineCoefficients.cols(), maxOrder + 1 );
    BOOST_CHECK_EQUAL( sineCoefficients.rows(), maxDegree + 1 );
    BOOST_CHECK_EQUAL( sineCoefficients.cols(), maxOrder + 1 );

    // Example: check one value
    BOOST_CHECK_CLOSE( cosineCoefficients(3,1), 1.054997670055739e-01, 1.0e-12 );
    BOOST_CHECK_CLOSE( sineCoefficients(6,2), 2.139319739882320e-02, 1.0e-12 );
}


BOOST_AUTO_TEST_CASE( testCreateSHDataset_Defaults_And_Validation )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;

    const boost::filesystem::path fullPath =
            "/Users/markusreichel/PhD/tudatpy/tests/test_tudat/src/astro/aerodynamics/test_data/test_input_coma.txt";
    BOOST_REQUIRE( boost::filesystem::exists(fullPath) );

    const std::vector< std::string > files = { fullPath.string( ) };
    PolyCoefFileProcessing proc( files );

    // Parse once to know available maxima from the file
    ComaPolyDataset poly = proc.createPolyCoefDataset( );
    BOOST_REQUIRE_EQUAL( poly.getNumFiles(), 1 );
    const int availDeg = poly.getMaxDegreeSH( 0 );
    BOOST_CHECK_EQUAL( availDeg, 10 );
    const int availOrd = poly.getSHDegreeAndOrderIndices( 0 ).row( 1 ).abs( ).maxCoeff( );

    // Radii/longitudes used in all tests
    const std::vector< double > radii_m = { 1000.0, 6000.0 };
    const std::vector< double > lons_deg = { 0.0, 30.0, 180.0 };

    // 1) Default args (-1,-1) should select available maxima
    {
        ComaStokesDataset sh = proc.createSHDataset( radii_m, lons_deg );
        BOOST_CHECK_EQUAL( sh.nmax(), availDeg ); // dataset uses degree as nmax internally
        BOOST_CHECK_EQUAL( sh.nFiles(), files.size() );
        BOOST_CHECK_EQUAL( sh.nRadii(), radii_m.size() );
        BOOST_CHECK_EQUAL( sh.nLongitudes(), lons_deg.size() );
        // nCoeffs = (nmax+1)(nmax+2)/2
        BOOST_CHECK_EQUAL( sh.nCoeffs(), static_cast<std::size_t>((availDeg+1)*(availDeg+2)/2) );
    }

    // 2) Explicit truncation to smaller (e.g., deg=6, ord=4) should work
    {
        const int reqDeg = std::min( 6, availDeg );
        const int reqOrd = std::min( 4, availOrd );
        ComaStokesDataset sh = proc.createSHDataset( radii_m, lons_deg, reqDeg, reqOrd );
        BOOST_CHECK_EQUAL( sh.nmax(), reqDeg );
        BOOST_CHECK_EQUAL( sh.nCoeffs(), static_cast<std::size_t>((reqDeg+1)*(reqDeg+2)/2) );

        // Optional: check matrix-shaped accessor has (reqDeg+1) x (reqDeg+1),
        // entries with m>reqOrd remain zero (we can just check shape here).
        auto [ Cmat, Smat ] = sh.getCoefficientMatrices( 0, /*ri*/0, /*li*/0 );
        BOOST_CHECK_EQUAL( Cmat.rows(), reqDeg + 1 );
        BOOST_CHECK_EQUAL( Cmat.cols(), reqDeg + 1 );
        BOOST_CHECK_EQUAL( Smat.rows(), reqDeg + 1 );
        BOOST_CHECK_EQUAL( Smat.cols(), reqDeg + 1 );
    }

    // 3) Request larger degree than available -> throws
    if(availDeg >= 0)
    {
        BOOST_CHECK_THROW(
                proc.createSHDataset(radii_m, lons_deg,availDeg + 1, std::max(0, availOrd)),
                std::invalid_argument
                );
    }

    // 4) Request larger order than available -> throws
    if(availOrd >= 0)
    {
        BOOST_CHECK_THROW(
                proc.createSHDataset( radii_m, lons_deg, std::max(0, availDeg), availOrd + 1),
                std::invalid_argument
                );
    }

    // 5) Empty radii -> throws
    {
        std::vector< double > emptyR;
        BOOST_CHECK_THROW(
                proc.createSHDataset(emptyR, lons_deg),
                std::invalid_argument
                );
    }

    // 6) Empty longitudes -> throws
    {
        std::vector< double > emptyL;
        BOOST_CHECK_THROW(
                proc.createSHDataset(radii_m, emptyL),
                std::invalid_argument
                );
    }
}


BOOST_AUTO_TEST_SUITE_END( )
} // namespace unit_tests
} // namespace tudat