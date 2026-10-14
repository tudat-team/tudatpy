#include "tudat/astro/earth_orientation/highFrequencyEarthOrientationModel.h"
#include "tudat/io/readHighFrequencyEopData.h"

#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace earth_orientation
{
namespace
{
// The input library has already normalized all amplitudes to SI. Filter in the
// conventional file's units, before handing the matrices to the common summation engine.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > oceanMatrices( const input_output::HighFrequencyEopData& data,
                                                             const bool polarMotion,
                                                             const double minimumAmplitude )
{
    const int firstComponent = polarMotion ? 0 : 2;
    const int components = polarMotion ? 2 : 1;
    const double factor = polarMotion ? unit_conversions::convertArcSecondsToRadians< double >( 1.0e-6 ) : 1.0e-6;
    Eigen::MatrixXd amplitudes( data.terms.size( ), 2 * components );
    Eigen::MatrixXd multipliers( data.terms.size( ), 6 );
    int retained = 0;
    for( const auto& term : data.terms )
    {
        for( int component = 0; component < components; ++component )
        {
            amplitudes( retained, 2 * component ) = term.sineAmplitudes[ firstComponent + component ];
            amplitudes( retained, 2 * component + 1 ) = term.cosineAmplitudes[ firstComponent + component ];
        }
        if( minimumAmplitude > 0.0 && amplitudes.row( retained ).norm( ) / factor <= minimumAmplitude ) continue;
        for( int argument = 0; argument < 6; ++argument ) multipliers( retained, argument ) = term.argumentMultipliers[ argument ];
        ++retained;
    }
    amplitudes.conservativeResize( retained, Eigen::NoChange );
    multipliers.conservativeResize( retained, Eigen::NoChange );
    return { amplitudes, multipliers };
}
}  // namespace

HighFrequencyEopCorrectionCalculators createHighFrequencyEopCorrectionCalculators( const OceanTideEopModel oceanModel,
                                                                                   const double minimumPolarMotionAmplitude,
                                                                                   const double minimumUt1Amplitude )
{
    if( !std::isfinite( minimumPolarMotionAmplitude ) || minimumPolarMotionAmplitude < 0.0 || !std::isfinite( minimumUt1Amplitude ) ||
        minimumUt1Amplitude < 0.0 )
        throw std::runtime_error( "High-frequency EOP amplitude cutoffs must be finite and nonnegative" );
    if( oceanModel == OceanTideEopModel::iers_2010 )
        return { getDefaultPolarMotionCorrectionCalculator( minimumPolarMotionAmplitude ),
                 getDefaultUT1CorrectionCalculator( minimumUt1Amplitude ) };
    if( oceanModel != OceanTideEopModel::desai_sibois_2016 ) throw std::runtime_error( "Unknown ocean-tide EOP model" );

    const auto ocean = input_output::getDesaiSibois2016EopData( );
    const auto polarOcean = oceanMatrices( ocean, true, minimumPolarMotionAmplitude );
    const auto ut1Ocean = oceanMatrices( ocean, false, minimumUt1Amplitude );
    const std::string path = paths::getEarthOrientationDataFilesPath( ) + "/";
    auto polarLibration = readAmplitudesAndFundamentalArgumentMultipliers(
            path + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
            path + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
            minimumPolarMotionAmplitude );
    auto ut1Libration = readAmplitudesAndFundamentalArgumentMultipliers(
            path + "utcLibrationAmplitudes.txt", path + "utcLibrationFundamentalArgumentMultipliers.txt", minimumUt1Amplitude );
    polarLibration.first *= unit_conversions::convertArcSecondsToRadians< double >( 1.0e-6 );
    ut1Libration.first *= 1.0e-6;
    // Exactly two groups: retained libration and the replacement ocean table.
    // Published LOD amplitudes are metadata, not an extra UT1-rate contribution.
    return { std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                     std::vector< Eigen::MatrixXd >{ polarLibration.first, polarOcean.first },
                     std::vector< Eigen::MatrixXd >{ polarLibration.second, polarOcean.second },
                     std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ) ),
             std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                     std::vector< Eigen::MatrixXd >{ ut1Libration.first, ut1Ocean.first },
                     std::vector< Eigen::MatrixXd >{ ut1Libration.second, ut1Ocean.second },
                     std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ) ) };
}

}  // namespace earth_orientation
}  // namespace tudat
