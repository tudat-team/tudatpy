/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/gravitation/earthSolidBodyTideGravityFieldVariations.h"

#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

namespace tudat
{
namespace gravitation
{
namespace
{
//! One IERS frequency-dependent term before its table unit scale is applied.
struct FrequencyTerm {
    //! Multipliers of the six Doodson arguments in tau, s, h, p, N-prime, ps order.
    std::array< int, 6 > multipliers;

    //! In-phase and out-of-phase Stokes amplitudes in units of 10^-12.
    double inPhase, outOfPhase;
};
#include "iers2010SolidTideFrequencyTerms.inc"

//! Evaluate one frequency-dependent term in the normalized C-iS convention.
// Adapted from the original ILR Earth solid-tide evaluator. In C-iS convention,
// eta1=-i in IERS equation 6.8b gives the positive diurnal cosine term in S21.
std::pair< double, double > calculateSingleFrequencyDependentSphericalHarmonicCorrection( const FrequencyTerm& term,
                                                                                          const Eigen::Vector6d& arguments,
                                                                                          const int order )
{
    double angle = 0.0;

    for( int i = 0; i < 6; ++i ) angle += term.multipliers[ i ] * arguments[ i ];
    const double cosine = std::cos( angle ), sine = std::sin( angle );

    if( order == 1 ) return { term.inPhase * sine + term.outOfPhase * cosine, term.inPhase * cosine - term.outOfPhase * sine };
    return { term.inPhase * cosine - term.outOfPhase * sine, order == 0 ? 0.0 : -term.inPhase * sine - term.outOfPhase * cosine };
}

//! Reject invalid physical constants before evaluating radius or gravitational-parameter ratios.
void requirePositiveFinite( const double value )
{
    if( !std::isfinite( value ) || value <= 0.0 )
        throw std::invalid_argument( "Earth gravity tide GM and radius must be positive and finite." );
}
}  // namespace

//! Retain the explicit Earth response conventions and connect the basic tide-raiser correction.
EarthSolidBodyTideGravityFieldVariations::EarthSolidBodyTideGravityFieldVariations(
        const std::function< Eigen::Vector6d( double ) >& deformedBodyState,
        const std::function< Eigen::Quaterniond( double ) >& toBodyFixedFrame,
        const std::vector< std::function< Eigen::Vector6d( double ) > >& deformingBodyStates,
        const double sourceRadius,
        const std::function< double( ) >& sourceGravitationalParameter,
        const std::vector< std::function< double( ) > >& deformingBodyGravitationalParameters,
        const std::vector< std::string >& deformingBodies,
        const std::function< Eigen::Vector6d( double ) >& doodsonArguments,
        const double fieldRadius,
        const double fieldGravitationalParameter,
        const basic_astrodynamics::TideSystem tideSystem ):
    SolidBodyTideGravityFieldVariations( deformedBodyState,
                                         toBodyFixedFrame,
                                         deformingBodyStates,
                                         sourceRadius,
                                         sourceGravitationalParameter,
                                         deformingBodyGravitationalParameters,
                                         deformingBodies,
                                         4,
                                         3 ),
    doodsonArguments_( doodsonArguments ), fieldRadius_( fieldRadius ), fieldGravitationalParameter_( fieldGravitationalParameter ),
    tideSystem_( tideSystem )
{
    requirePositiveFinite( sourceRadius );
    requirePositiveFinite( fieldRadius );
    requirePositiveFinite( fieldGravitationalParameter );
    if( !deformedBodyState || !toBodyFixedFrame || !sourceGravitationalParameter || !doodsonArguments || deformingBodyStates.empty( ) ||
        deformingBodyStates.size( ) != deformingBodyGravitationalParameters.size( ) ||
        deformingBodyStates.size( ) != deformingBodies.size( ) )
        throw std::invalid_argument( "Earth solid tide requires matching nonempty body, GM and argument callbacks." );
    std::set< std::string > names;

    for( std::size_t i = 0; i < deformingBodies.size( ); ++i )
        if( !deformingBodyStates[ i ] || !deformingBodyGravitationalParameters[ i ] || deformingBodies[ i ].empty( ) ||
            !names.insert( deformingBodies[ i ] ).second )
            throw std::invalid_argument( "Earth solid tide requires distinct tide raisers with state and GM callbacks." );
    if( tideSystem != basic_astrodynamics::TideSystem::tideFree && tideSystem != basic_astrodynamics::TideSystem::zeroTide )
        throw std::invalid_argument( "Conventional Earth solid tide requires a known tide-free or zero-tide nominal field." );
    correctionFunctions.push_back(
            [ this ]( Eigen::MatrixXd& cosine, Eigen::MatrixXd& sine ) { addBasicSolidBodyTideCorrections( cosine, sine ); } );
}

//! Add the IERS anelastic response using the current tide-raiser geometry prepared by the base model.
void EarthSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections( Eigen::MatrixXd& cosine, Eigen::MatrixXd& sine )
{
    requirePositiveFinite( radiusRatio );
    requirePositiveFinite( massRatio );
    if( !std::isfinite( sineOfLatitude ) || !std::isfinite( iLongitude.imag( ) ) )
        throw std::runtime_error( "Non-finite Earth solid-tide geometry." );
    // Dimensionless Love numbers are indexed by order; the plus terms induce degree four from degree-two forcing.
    const std::array< std::complex< double >, 3 > degreeTwoLove = { 0.30190, { 0.29830, -0.00144 }, { 0.30102, -0.00130 } };

    const std::array< double, 4 > degreeThreeLove = { 0.093, 0.093, 0.093, 0.094 };

    const std::array< double, 3 > degreeTwoLovePlus = { -0.00089, -0.00080, -0.00057 };

    for( int n = 2; n <= 3; ++n )
        for( int m = 0; m <= n; ++m )
        {
            updateTidalAmplitudeAndArgument( n, m );
            const auto correction = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                    n == 2 ? degreeTwoLove[ m ] : std::complex< double >( degreeThreeLove[ m ], 0.0 ),
                    massRatio,
                    std::pow( radiusRatio, n + 1 ),
                    tideAmplitude,
                    tideArgument,
                    n,
                    m );
            cosine( n - 2, m ) += correction.real( );
            if( m > 0 ) sine( n - 2, m ) -= correction.imag( );
        }
    // Equation 6.7 uses degree TWO forcing for degree FOUR response.
    for( int m = 0; m <= 2; ++m )
    {
        updateTidalAmplitudeAndArgument( 2, m );
        const auto correction = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                degreeTwoLovePlus[ m ], massRatio, std::pow( radiusRatio, 3 ), tideAmplitude, tideArgument, 2, m );
        cosine( 2, m ) += correction.real( );
        if( m > 0 ) sine( 2, m ) -= correction.imag( );
    }
}

//! Combine basic and frequency-dependent terms, then remove the retained permanent contribution for a zero-tide nominal field.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > EarthSolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrections(
        const double time )
{
    if( !std::isfinite( time ) ) throw std::invalid_argument( "Earth solid-tide epoch must be finite." );
    const double sourceGM = deformedBodyMass_( );

    requirePositiveFinite( sourceGM );
    auto corrections = calculateBasicSphericalHarmonicsCorrections( time );

    const Eigen::Vector6d arguments = doodsonArguments_( time );

    if( !arguments.allFinite( ) ) throw std::runtime_error( "Non-finite Earth solid-tide arguments." );
    // Published frequency amplitudes use 10^-12 Stokes units and the selected TT/UT1 arguments.
    const auto addBand = [ & ]( const auto& terms, const int order ) {
        for( const auto& term : terms )
        {
            const auto correction = calculateSingleFrequencyDependentSphericalHarmonicCorrection( term, arguments, order );

            corrections.first( 0, order ) += 1e-12 * correction.first;
            corrections.second( 0, order ) += 1e-12 * correction.second;
        }
    };

    addBand( zonalTerms, 0 );
    addBand( diurnalTerms, 1 );
    addBand( semidiurnalTerms, 2 );
    if( tideSystem_ == basic_astrodynamics::TideSystem::zeroTide ) corrections.first( 0, 0 ) -= 4.4228e-8 * ( -0.31460 ) * 0.30190;
    for( int n = 2; n <= 4; ++n )
    {
        const double scale = sourceGM / fieldGravitationalParameter_ * std::pow( deformedBodyReferenceRadius_ / fieldRadius_, n );

        corrections.first.row( n - 2 ) *= scale;
        corrections.second.row( n - 2 ) *= scale;
    }
    if( !corrections.first.allFinite( ) || !corrections.second.allFinite( ) )
        throw std::runtime_error( "Non-finite Earth solid-tide coefficients." );
    currentCosineCorrections_ = corrections.first;
    currentSineCorrections_ = corrections.second;
    return corrections;
}

//! Validate the supplied polar wobble and prepare the degree-two source-to-field reference scaling.
SolidPoleTideGravityFieldVariations::SolidPoleTideGravityFieldVariations(
        const std::function< Eigen::Vector2d( double ) >& wobbleArcseconds,
        const double sourceGravitationalParameter,
        const double sourceRadius,
        const double fieldGravitationalParameter,
        const double fieldRadius ): GravityFieldVariations( 2, 1, 2, 1 ), wobbleArcseconds_( wobbleArcseconds )
{
    requirePositiveFinite( sourceGravitationalParameter );
    requirePositiveFinite( sourceRadius );
    requirePositiveFinite( fieldGravitationalParameter );
    requirePositiveFinite( fieldRadius );
    if( !wobbleArcseconds ) throw std::invalid_argument( "Solid pole gravity tide requires a wobble callback." );
    coefficientScale_ = sourceGravitationalParameter / fieldGravitationalParameter * std::pow( sourceRadius / fieldRadius, 2 );
}

//! Apply IERS section 6.4 to the supplied wobble and return the C21/S21 correction block.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidPoleTideGravityFieldVariations::calculateSphericalHarmonicsCorrections(
        const double time )
{
    if( !std::isfinite( time ) ) throw std::invalid_argument( "Solid pole gravity tide epoch must be finite." );
    const Eigen::Vector2d wobble = wobbleArcseconds_( time );

    if( !wobble.allFinite( ) ) throw std::runtime_error( "Non-finite solid pole gravity tide wobble." );
    Eigen::MatrixXd cosine( 1, 1 ), sine( 1, 1 );

    cosine( 0, 0 ) = -1.333e-9 * ( wobble.x( ) + 0.0115 * wobble.y( ) ) * coefficientScale_;
    sine( 0, 0 ) = -1.333e-9 * ( wobble.y( ) - 0.0115 * wobble.x( ) ) * coefficientScale_;
    return { cosine, sine };
}

}  // namespace gravitation
}  // namespace tudat
