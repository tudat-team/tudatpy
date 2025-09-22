#include "tudat/astro/aerodynamics/ComaModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/utilityMacros.h"
#include <stdexcept>
#include <cmath>
#include <utility>

#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{
namespace aerodynamics
{
ComaModel::ComaModel( const std::vector< Eigen::MatrixXd >& polyCoefficients,
                      const std::vector< Eigen::ArrayXXi >& SHDegreeAndOrder,
                      const std::vector< Eigen::VectorXd >& powersInvRadius,
                      const std::vector< double >& referenceRadius,
                      const std::vector< std::vector< double > >& TimePeriods,
                      const int& maximumDegree = -1,
                      const int& maximumOrder = -1 ) :
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    polyCoefficients_( polyCoefficients ),
    SHDegreeAndOrder_( SHDegreeAndOrder ),
    referenceRadius_( referenceRadius ),
    powersInvRadius_( powersInvRadius ),
    TimePeriods_( TimePeriods ),
    constantTemperature_( 300.0 ), // TODO: Placeholder
    specificGasConstant_( 5 ), //TODO: Placeholder
    ratioOfSpecificHeats_( 1.4 ) // TODO: Placeholder
{
}


double ComaModel::getDensity( const double altitude,
                              const double longitude,
                              const double latitude,
                              const double time )
{
    const int currentTimeIntervalIndex = findTimeIntervalIndex( time );
    const Eigen::MatrixXd currentPolyCoefs = polyCoefficients_[ currentTimeIntervalIndex ];
    const Eigen::ArrayXXi currentSHDegreeAndOrder = SHDegreeAndOrder_[ currentTimeIntervalIndex ];
    const Eigen::ArrayXXd currentPowersInvRadius = powersInvRadius_[ currentTimeIntervalIndex ];
    const double currentReferenceRadius = referenceRadius_[ currentTimeIntervalIndex ];

    const double raw = std::sqrt( currentPolyCoefs.cols( ) );
    const int N = static_cast< int >(std::round( raw ));
    if(N * N != currentPolyCoefs.cols( ))
    {
        std::cerr << "[ERROR] polyCoefs.cols() is not a perfect square: " << currentPolyCoefs.cols( ) << std::endl;
        std::exit( EXIT_FAILURE );
    }

    Eigen::MatrixXd cosineCoefficients, sineCoefficients;
    constexpr double distanceToCometCentre = 0.0; // TODO: Dummy. Figure out how to retrieve distance to comet centre
    constexpr double solarLongitude = 0.0; // TODO: Dummy. Figure out how to retrieve solar longitude

    evaluateStokesCoefficients2D( distanceToCometCentre,
                                  solarLongitude,
                                  currentPolyCoefs,
                                  currentSHDegreeAndOrder,
                                  currentPowersInvRadius,
                                  currentReferenceRadius,
                                  cosineCoefficients,
                                  sineCoefficients,
                                  maximumDegree_,
                                  maximumOrder_ );

    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    TUDAT_UNUSED_PARAMETER( altitude );

    return altitude; // TODO: implement actual SH-based density model
}


double ComaModel::getPressure( const double altitude,
                               const double longitude,
                               const double latitude,
                               const double time )
{
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    return altitude; // TODO: implement actual SH-based pressure model
}

double ComaModel::getTemperature( const double altitude,
                                  const double longitude,
                                  const double latitude,
                                  const double time )
{
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    return altitude; // TODO: implement actual SH-based temperature model
}

double ComaModel::getSpeedOfSound( const double altitude,
                                   const double longitude,
                                   const double latitude,
                                   const double time )
{
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    return altitude; // TODO: implement actual SH-based speed of sound model
}


int ComaModel::findTimeIntervalIndex( const double time ) const
{
    for(std::size_t i = 0; i < TimePeriods_.size( ); ++i)
    {
        const auto& period = TimePeriods_[ i ];

        // Check if the time period has exactly two values
        if(period.size( ) != 2)
        {
            std::ostringstream err;
            err << "Invalid time period at index " << i
                    << ": expected 2 values (start and end), got "
                    << period.size( ) << ".";
            throw std::runtime_error( err.str( ) );
        }

        const double startTime = period[ 0 ];
        const double endTime = period[ 1 ];

        // Check if time falls within this interval
        if(time >= startTime && time <= endTime)
        {
            return static_cast< int >(i);
        }
    }

    // If we get here, no matching time period was found
    std::ostringstream err;
    err << "Time " << time << "does not fall into any"
            << " defined time period of coma model input files.";
    throw std::runtime_error( err.str( ) );
}


void ComaModel::evaluateStokesCoefficients2D(
        const double distanceToCometCentre,
        const double solarLongitude,
        const Eigen::ArrayXXd& polyCoefficients,
        const Eigen::ArrayXXi& atDegreeAndOrder,
        const Eigen::VectorXd& atPowersInvRadius,
        const double refRadius,
        Eigen::MatrixXd& cosineCoefficients,
        Eigen::MatrixXd& sineCoefficients,
        int maxDegree,
        int maxOrder )
{
    const int maxDegAvailable = atDegreeAndOrder.row( 0 ).maxCoeff( );
    const int maxOrdAvailable = atDegreeAndOrder.row( 1 ).abs( ).maxCoeff( );

    if(maxDegree < 0) maxDegree = maxDegAvailable;
    if(maxOrder < 0) maxOrder = maxOrdAvailable;

    if(maxDegree > maxDegAvailable || maxOrder > maxOrdAvailable)
    {
        std::ostringstream err;
        err << "[FATAL] Requested maxDegree=" << maxDegree
                << ", maxOrder=" << maxOrder
                << " exceeds available (degree=" << maxDegAvailable
                << ", order=" << maxOrdAvailable << ")";
        throw std::runtime_error( err.str( ) );
    }

    const Eigen::Index numRadialTerms = atPowersInvRadius.size( );
    const Eigen::Index numIntervals = polyCoefficients.rows( ) / numRadialTerms;

    cosineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );
    sineCoefficients = Eigen::MatrixXd::Zero( maxDegree + 1, maxOrder + 1 );

    const bool usePolyvalForm = ( distanceToCometCentre <= refRadius || refRadius < 1.0e-10 );
    const double scaling = usePolyvalForm
            ? ( ( refRadius < 1.0e-10 ) ? 1.0 / distanceToCometCentre : 1.0 / distanceToCometCentre - 1.0 / refRadius )
            : ( ( refRadius < 1.0e-10 ) ? 1.0 / distanceToCometCentre : refRadius / distanceToCometCentre );

    for(int i = 0; i < polyCoefficients.cols( ); ++i)
    {
        const int l = atDegreeAndOrder( 0, i );
        const int m = atDegreeAndOrder( 1, i );
        const int absM = std::abs( m );

        if(l > maxDegree || absM > maxOrder)
            continue;

        const double value = usePolyvalForm
                ? radialPolyvalAndTemporalIFFT( scaling,
                                                solarLongitude,
                                                polyCoefficients.col( i ).reshaped( numIntervals, numRadialTerms ).matrix( ),
                                                atPowersInvRadius.array( ) )
                : radialSquaredAndTemporalIFFT( scaling,
                                                solarLongitude,
                                                polyCoefficients.block( 0, i, numIntervals, 1 ).matrix( ) );

        ( m >= 0 ? cosineCoefficients : sineCoefficients )( l, absM ) = value;
    }
}

//multilinearinterpolator


/// Evaluate general Stokes' coefficient (Cnm or Snm) using temporal Fourier IFFT-like form
double ComaModel::radialPolyvalAndTemporalIFFT( const double nucleusToSpacecraftDistance,
                                                const double solarLongitude,
                                                const Eigen::MatrixXd& polyCoefs,
                                                const Eigen::ArrayXd& atPowersInvRadius )
{
    const Eigen::Index N = polyCoefs.rows( ); // Number of angular basis functions

    // Construct alpha-dependent cosine/sine vector (FFT-like weights)
    auto alphaPart = ( Eigen::ArrayXd( N ) <<
        cos( Eigen::ArrayXd::LinSpaced( N / 2 + 1, 0.0, static_cast< double >(N) / 2.0 ) * solarLongitude ),
        sin( Eigen::ArrayXd::LinSpaced( N % 2 == 0 ? N / 2 - 1 : N / 2, 1.0, static_cast< double >(N) / 2 - ( N % 2 == 0 ? 1 : 0 ) ) *
                solarLongitude )
    ).finished( );

    return ( alphaPart.matrix( ).transpose( ) * polyCoefs * pow( nucleusToSpacecraftDistance, atPowersInvRadius ).matrix( ) ).value( ) /
            static_cast< double >(N);
}


// Reduced version: works with 1D radial coefficients (e.g. after collapsing temporal part)
double ComaModel::radialSquaredAndTemporalIFFT( const double nucleusToSpacecraftDistance,
                                                const double solarLongitude,
                                                const Eigen::VectorXd& polyCoefs )
{
    const Eigen::Index N = polyCoefs.rows( );

    auto alphaPart = ( Eigen::ArrayXd( N ) <<
        cos( Eigen::ArrayXd::LinSpaced( N / 2 + 1, 0.0, static_cast< double >(N) / 2.0 ) * solarLongitude ),
        sin( Eigen::ArrayXd::LinSpaced( N % 2 == 0 ? N / 2 - 1 : N / 2, 1.0, static_cast< double >(N) / 2 - ( N % 2 == 0 ? 1 : 0 ) ) *
                solarLongitude )
    ).finished( );

    return ( alphaPart.matrix( ).transpose( ) * polyCoefs * ( nucleusToSpacecraftDistance * nucleusToSpacecraftDistance ) ).value( ) /
            static_cast< double >(N);
}

SphericalHarmonicsDensity::SphericalHarmonicsDensity( Eigen::MatrixXd sineCoefficients,
                                                      Eigen::MatrixXd cosineCoefficients,
                                                      std::string fixedReferenceFrame ) :
    sineCoefficients_( std::move( sineCoefficients ) ),
    cosineCoefficients_( std::move( cosineCoefficients ) ),
    fixedReferenceFrame_( std::move( fixedReferenceFrame ) ),
    maximumDegree_( static_cast< int >(cosineCoefficients_.rows( )) ),
    maximumOrder_( static_cast< int >(cosineCoefficients_.cols( )) ),
    sphericalHarmonicsCache_( basic_mathematics::SphericalHarmonicsCache( true, true) )
{
    if(cosineCoefficients_.rows( ) != sineCoefficients_.rows( ) ||
        cosineCoefficients_.cols( ) != sineCoefficients_.cols( ))
    {
        throw std::runtime_error( "Spherical harmonics coefficient sizes are incompatible." );
    }

    sphericalHarmonicsCache_.resetMaximumDegreeAndOrder( maximumDegree_ + 1, maximumOrder_ + 1 );
    // TODO: can be set to highest Degree/Order
}


double SphericalHarmonicsDensity::calculateSurfaceSphericalHarmonics(
        const Eigen::MatrixXd& sineCoefficients,
        const Eigen::MatrixXd& cosineCoefficients,
        const double latitude,
        const double longitude,
        const int highestDegree,
        const int highestOrder )
{
    const double sineOfAngle = std::sin( latitude );
    sphericalHarmonicsCache_.updateAnglesOnly( sineOfAngle,
                                     longitude );

    const basic_mathematics::LegendreCache& legendreCacheReference = sphericalHarmonicsCache_.getLegendreCache( );


    double legendrePolynomial = TUDAT_NAN;
    double cosineOfOrderLongitude = TUDAT_NAN;
    double sineOfOrderLongitude = TUDAT_NAN;
    double value = 0.0;

    // Loop through all degrees.
    for (int l = 0; l <= highestDegree; ++l) {
        const int mmax = std::min(l, highestOrder);
        for (int m = 0; m <= mmax; ++m) {
            const double P = legendreCacheReference.getLegendrePolynomial(l, m); // normalized, correct argument
            const double cos_mλ = sphericalHarmonicsCache_.getCosineOfMultipleLongitude(m);

            if (m == 0) {
                value += P * cosineCoefficients(l, 0);
            } else {
                const double sin_mλ = sphericalHarmonicsCache_.getSineOfMultipleLongitude(m);
                value += P * (cosineCoefficients(l, m) * cos_mλ
                            + sineCoefficients(l, m)    * sin_mλ);
            }
        }
    }
    return value;
}





//
// double SphericalHarmonicsDensity::calculateSurfaceSphericalHarmonics(
//         const double latitude,
//         const double longitude,
//         const int highestDegree,
//         const int highestOrder )
// {
//
//     basic_mathematics::LegendreCache& legendreCacheReference = sphericalHarmonicsCache_.getLegendreCache( );
//     legendreCacheReference.update( std::sin( latitude ) );
//
//
//     double value = 0.0;
//
//     for(int degree = 0; degree < highestDegree; degree++)
//     {
//         for(int order = 0; (order <= degree) && (order < highestOrder); order++)
//         {
//             const double legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomialFromCache(
//                     degree,
//                     order,
//                     legendreCacheReference );
//
//             value += legendrePolynomial * (
//                 cosineCoefficients_( degree, order ) * std::cos( order * longitude ) +
//                 sineCoefficients_( degree, order ) * std::sin( order * longitude ) );
//         }
//     }
//     return value;
// }


// --- Unit test access wrapper ---
void ComaModel::testEvaluateStokesCoefficients2D(
        const double distanceToCometCentre,
        const double solarLongitude,
        const Eigen::ArrayXXd& polyCoefficients,
        const Eigen::ArrayXXi& atSHDegreeAndOrder,
        const Eigen::VectorXd& atPowersInvRadius,
        const double refRadius,
        Eigen::MatrixXd& cosineCoefficients,
        Eigen::MatrixXd& sineCoefficients,
        const int maxDegree,
        const int maxOrder )
{
    evaluateStokesCoefficients2D(
            distanceToCometCentre,
            solarLongitude,
            polyCoefficients,
            atSHDegreeAndOrder,
            atPowersInvRadius,
            refRadius,
            cosineCoefficients,
            sineCoefficients,
            maxDegree,
            maxOrder );
}
}
}