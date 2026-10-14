

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

#include "tudat/basics/utilities.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"

namespace tudat
{

namespace basic_astrodynamics
{

// Tudat adaptation of STEP2DIU/STEP2LON argument polynomials supplied by
// the IERS Conventions Center. This differently named C++ routine combines
// their argument calculations, obtains UTC from the API's TT approximation,
// and returns radians. It is neither supplied nor endorsed by the Center.
// Original unmodified routines and their revisions are identified in the
// independent test fixtures' sources.json; the IERS license is reproduced below.
Eigen::Vector6d calculateIers2010SolidEarthTideArguments( const double secondsSinceJ2000 )
{
    const double utcSecondsSinceJ2000 = sofa_interface::convertTTtoUTC( secondsSinceJ2000 );
    const double utcDaysSinceMidnight = utcSecondsSinceJ2000 / physical_constants::JULIAN_DAY + 0.5;
    const double fractionalHours = 24.0 * ( utcDaysSinceMidnight - std::floor( utcDaysSinceMidnight ) );
    const double centuries = secondsSinceJ2000 / ( physical_constants::JULIAN_DAY * 36525.0 );

    const double unprecessedLunarLongitude =
            218.31664563 + ( 481267.88194 + ( -0.0014663889 + 0.00000185139 * centuries ) * centuries ) * centuries;
    const double tau = fractionalHours * 15.0 + 280.4606184 +
            ( 36000.7700536 + ( 0.00038793 - 0.0000000258 * centuries ) * centuries ) * centuries - unprecessedLunarLongitude;
    const double precession =
            ( 1.396971278 + ( 0.000308889 + ( 0.000000021 + 0.000000007 * centuries ) * centuries ) * centuries ) * centuries;
    const double solarLongitude = 280.46645 +
            ( 36000.7697489 + ( 0.00030322222 + ( 0.000000020 - 0.00000000654 * centuries ) * centuries ) * centuries ) * centuries;
    const double lunarPerigee = 83.35324312 +
            ( 4069.01363525 + ( -0.01032172222 + ( -0.0000124991 + 0.00000005263 * centuries ) * centuries ) * centuries ) * centuries;
    const double negativeLunarNode = 234.95544499 +
            ( 1934.13626197 + ( -0.00207561111 + ( -0.00000213944 + 0.00000001650 * centuries ) * centuries ) * centuries ) * centuries;
    const double solarPerigee = 282.93734098 +
            ( 1.71945766667 + ( 0.00045688889 + ( -0.00000001778 - 0.00000000334 * centuries ) * centuries ) * centuries ) * centuries;

    Eigen::Vector6d arguments;
    arguments << tau, unprecessedLunarLongitude + precession, solarLongitude, lunarPerigee, negativeLunarNode, solarPerigee;
    for( int i = 0; i < arguments.size( ); ++i )
    {
        arguments( i ) = std::fmod( arguments( i ), 360.0 ) * mathematical_constants::PI / 180.0;
    }
    return arguments;
}

//! Calculate trigonometric expressions for correction fo site displacements, step 1 of IERS Conventions 2010, Chapter 7
std::vector< Eigen::Vector3d > calculateDisplacementCorrectionTrigonometricPart( const Eigen::Vector3d& relativeBodyState,
                                                                                 const Eigen::Vector3d& stationSphericalCoordinates,
                                                                                 const std::vector< bool >& calculateTerms )
{
    std::vector< Eigen::Vector3d > stationUnitVectors;
    stationUnitVectors.resize( 3 );
    stationUnitVectors[ 0 ] = Eigen::Vector3d::UnitX( );
    stationUnitVectors[ 1 ] = Eigen::Vector3d::UnitY( );
    stationUnitVectors[ 2 ] = Eigen::Vector3d::UnitZ( );

    double stationLongitude = stationSphericalCoordinates.z( );
    double stationLatitude = stationSphericalCoordinates.y( );

    double bodyDistance = relativeBodyState.norm( );

    double sineOfStationLongitude = std::sin( stationLongitude );
    double sineOfTwiceStationLongitude = std::sin( 2.0 * stationLongitude );
    double cosineOfTwiceStationLongitude = std::cos( 2.0 * stationLongitude );

    double cosineOfStationLongitude = std::cos( stationLongitude );
    double sineOfStationLatitude = std::sin( stationLatitude );
    double cosineOfStationLatitude = std::cos( stationLatitude );
    double cosineTwiceStationLatitude = 1.0 - 2.0 * sineOfStationLatitude * sineOfStationLatitude;

    double bodyDistanceInvSquared = 1.0 / ( bodyDistance * bodyDistance );
    double xySquareDifference = ( relativeBodyState[ 0 ] * relativeBodyState[ 0 ] - relativeBodyState[ 1 ] * relativeBodyState[ 1 ] );
    std::vector< Eigen::Vector3d > trigonometricParts;
    trigonometricParts.resize( 6 );
    assert( calculateTerms.size( ) == 6 );
    if( calculateTerms[ 0 ] == 1 )
    {
        // Eq. 7.8
        trigonometricParts[ 0 ] = 3.0 * bodyDistanceInvSquared *
                ( -sineOfStationLatitude * sineOfStationLatitude * relativeBodyState[ 2 ] *
                          ( relativeBodyState[ 0 ] * cosineOfStationLongitude + relativeBodyState[ 1 ] * sineOfStationLongitude ) *
                          stationUnitVectors[ 1 ] +
                  sineOfStationLatitude * cosineTwiceStationLatitude * relativeBodyState[ 2 ] *
                          ( relativeBodyState[ 0 ] * sineOfStationLongitude - relativeBodyState[ 1 ] * cosineOfStationLongitude ) *
                          stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 0 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 1 ] == 1 )
    {
        // Eq. 7.9
        trigonometricParts[ 1 ] = bodyDistanceInvSquared *
                ( -1.5 * sineOfStationLatitude * cosineOfStationLatitude *
                          ( xySquareDifference * cosineOfTwiceStationLongitude +
                            2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * sineOfTwiceStationLongitude ) *
                          stationUnitVectors[ 1 ] -
                  1.5 * sineOfStationLatitude * sineOfStationLatitude * cosineOfStationLatitude *
                          ( xySquareDifference * sineOfTwiceStationLongitude -
                            2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * cosineOfTwiceStationLongitude ) *
                          stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 1 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 2 ] == 1 )
    {
        // Eq. 7.10a
        trigonometricParts[ 2 ] = -3.0 * sineOfStationLatitude * cosineOfStationLatitude * relativeBodyState[ 2 ] *
                ( relativeBodyState[ 0 ] * sineOfStationLongitude - relativeBodyState[ 1 ] * cosineOfStationLongitude ) *
                bodyDistanceInvSquared * stationUnitVectors[ 2 ];
    }
    else
    {
        trigonometricParts[ 2 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 3 ] == 1 )
    {
        // Eq. 7.10b
        trigonometricParts[ 3 ] = bodyDistanceInvSquared *
                ( -3.0 * cosineTwiceStationLatitude * relativeBodyState[ 2 ] *
                          ( relativeBodyState[ 0 ] * sineOfStationLongitude - relativeBodyState[ 1 ] * cosineOfStationLongitude ) *
                          stationUnitVectors[ 1 ] -
                  3.0 * sineOfStationLatitude * relativeBodyState[ 2 ] *
                          ( relativeBodyState[ 0 ] * cosineOfStationLongitude + relativeBodyState[ 1 ] * sineOfStationLongitude ) *
                          stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 3 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 4 ] == 1 )
    {
        // Eq. 7.11a
        trigonometricParts[ 4 ] = -0.75 * cosineOfStationLatitude * cosineOfStationLatitude * bodyDistanceInvSquared *
                ( xySquareDifference * sineOfTwiceStationLongitude -
                  2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * cosineOfTwiceStationLongitude ) *
                stationUnitVectors[ 2 ];
    }
    else
    {
        trigonometricParts[ 4 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 5 ] == 1 )
    {
        // Eq. 7.11b
        trigonometricParts[ 5 ] = bodyDistanceInvSquared *
                ( 1.5 * sineOfStationLatitude * cosineOfStationLatitude *
                          ( xySquareDifference * sineOfTwiceStationLongitude -
                            2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * cosineOfTwiceStationLongitude ) *
                          stationUnitVectors[ 1 ] -
                  1.5 * cosineOfStationLatitude *
                          ( xySquareDifference * cosineOfTwiceStationLongitude +
                            2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * sineOfTwiceStationLongitude ) *
                          stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 5 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    return trigonometricParts;
}

Eigen::Vector3d calculateFirstCorrectionStepDisplacements( const Eigen::Vector3d& relativeBodyState,
                                                           const Eigen::Vector3d& stationSphericalCoordinates,
                                                           const double gravitationalParameterRatio,
                                                           const double bodyEquatorialRadius,
                                                           const std::vector< double >& correctionLoveAndShidaNumbers,
                                                           const std::vector< bool >& areFirstStepCorrectionsCalculated )
{
    std::vector< Eigen::Vector3d > trigonometricParts = calculateDisplacementCorrectionTrigonometricPart(
            relativeBodyState, stationSphericalCoordinates, areFirstStepCorrectionsCalculated );

    double distanceRatio = bodyEquatorialRadius / relativeBodyState.norm( );
    double distanceRatioSquared = distanceRatio * distanceRatio;
    double commonTerm = gravitationalParameterRatio * distanceRatioSquared * distanceRatioSquared * relativeBodyState.norm( );

    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );
    for( int i = 0; i < 6; i++ )
    {
        if( areFirstStepCorrectionsCalculated[ i ] == 1 )
        {
            switch( i )
            {
                case 0:
                    displacement += correctionLoveAndShidaNumbers[ 0 ] * trigonometricParts[ 0 ];
                    break;
                case 1:
                    displacement += correctionLoveAndShidaNumbers[ 1 ] * trigonometricParts[ 1 ];
                    break;
                case 2:
                    displacement += correctionLoveAndShidaNumbers[ 2 ] * trigonometricParts[ 2 ];
                    break;
                case 3:
                    displacement += correctionLoveAndShidaNumbers[ 3 ] * trigonometricParts[ 3 ];
                    break;
                case 4:
                    displacement += correctionLoveAndShidaNumbers[ 4 ] * trigonometricParts[ 4 ];
                    break;
                case 5:
                    displacement += correctionLoveAndShidaNumbers[ 5 ] * trigonometricParts[ 5 ];
                    break;
            }
        }
    }
    displacement *= commonTerm;

    return displacement;
}

//! Function to calculate diurnal tide-frequency dependent site displacement correction to solid Earth tide displacements.
Eigen::Vector3d calculateDiurnalFrequencyDependentDisplacementCorrection( const double tideArgument,
                                                                          const double stationLongitude,
                                                                          const double stationLatitude,
                                                                          const Eigen::Vector4d& amplitudes )
{
    // Pre-calculate sine and cosines.
    double sineOfArgument = std::sin( tideArgument + stationLongitude );
    double cosineOfArgument = std::cos( tideArgument + stationLongitude );

    // Calculate displacements in local coordinates (east, north, up), IERS 2010 Conventions 7.12a and 7.12b
    Eigen::Vector3d displacement;
    displacement << ( amplitudes( 2 ) * cosineOfArgument - amplitudes( 3 ) * sineOfArgument ) * std::sin( stationLatitude ),
            ( amplitudes( 2 ) * sineOfArgument + amplitudes( 3 ) * cosineOfArgument ) * std::cos( 2.0 * stationLatitude ),
            ( amplitudes( 0 ) * sineOfArgument + amplitudes( 1 ) * cosineOfArgument ) * std::sin( 2.0 * stationLatitude );
    return displacement;
}

//! Function to calculate long period tide-frequency dependent site displacement correction to solid Earth tide displacements.
Eigen::Vector3d calculateLongPeriodFrequencyDependentDisplacementCorrection( const double tideArgument,
                                                                             const double stationLatitude,
                                                                             const Eigen::Vector4d& amplitudes )
{
    // Pre-calculate sine and cosines.
    double sineOfArgument = std::sin( tideArgument );
    double cosineOfArgument = std::cos( tideArgument );

    // Calculate displacements in local coordinates (east, north, up), IERS 2010 Conventions 7.13a and 7.13b
    Eigen::Vector3d displacement;
    displacement << 0.0, ( amplitudes( 2 ) * cosineOfArgument + amplitudes( 3 ) * sineOfArgument ) * std::sin( 2.0 * stationLatitude ),
            ( amplitudes( 0 ) * cosineOfArgument + amplitudes( 1 ) * sineOfArgument ) *
            ( 1.5 * std::sin( stationLatitude ) * std::sin( stationLatitude ) - 0.5 );
    return displacement;
}

//! Function to calculate tide-frequency dependent site displacement corrections to solid Earth tide displacements.
Eigen::Vector3d calculateFrequencyDependentDisplacementCorrections( const Eigen::MatrixXd& doodsonMultipliers,
                                                                    const Eigen::Matrix< double, 6, 1 >& doodsonArguments,
                                                                    const Eigen::MatrixXd& amplitudes,
                                                                    const Eigen::Vector3d& stationSphericalCoordinates )
{
    // Calculate current phase of tides.
    Eigen::VectorXd tideArguments = doodsonMultipliers * doodsonArguments;

    // Initialize displacements to zero.
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );

    // Loop over all tides for which corrections are to be applied and calculate corrections.
    for( int i = 0; i < doodsonMultipliers.rows( ); i++ )
    {
        // If tide is long-period (i.e. independent of mean Lunar time), use specific long-period correction method.
        if( doodsonMultipliers( i, 0 ) == 0 )
        {
            displacement += calculateLongPeriodFrequencyDependentDisplacementCorrection(
                    tideArguments( i ), stationSphericalCoordinates.y( ), amplitudes.block( i, 0, 1, 4 ).transpose( ) );
        }

        // If tide is not long-period (i.e. not independent of mean Lunar time), use specific diurnal correction method.
        else if( doodsonMultipliers( i, 0 ) == 1 )
        {
            displacement += calculateDiurnalFrequencyDependentDisplacementCorrection( tideArguments( i ),
                                                                                      stationSphericalCoordinates.z( ),
                                                                                      stationSphericalCoordinates.y( ),
                                                                                      amplitudes.block( i, 0, 1, 4 ).transpose( ) );
        }
        else
        {
            std::cerr << "Warning, frequency dependent love number value corrections only available for diurnal and long period tides."
                      << std::endl;
        }
    }
    return displacement;
}

Iers2010EarthDeformation::Iers2010EarthDeformation(
        const std::function< Eigen::Vector6d( const double ) > deformedBodyEphemeris,
        const std::vector< std::function< Eigen::Vector6d( const double ) > >& deformingBodyEphemerides,
        const std::function< Eigen::Quaterniond( const double ) > deformedBodyRotation,
        const std::function< double( ) > gravitionalParameterOfDeformedBody,
        const std::vector< std::function< double( ) > >& gravitionalParametersOfDeformingBodies,
        const double deformedBodyEquatorialRadius,
        const std::map< int, std::pair< double, double > >& nominalDisplacementLoveNumbers,
        const std::vector< double >& firstStepLatitudeDependenceTerms,
        const std::vector< bool >& areFirstStepCorrectionsCalculated,
        const std::vector< double >& correctionLoveAndShidaNumbers,
        const std::string& diurnalFile,
        const std::string& longPeriodFile,
        const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction,
        const std::function< Eigen::Quaterniond( const EarthDeformationEpoch& ) > contextRotation,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > contextClock ):
    BasicTidalBodyDeformation( deformedBodyEphemeris,
                               deformingBodyEphemerides,
                               deformedBodyRotation,
                               gravitionalParameterOfDeformedBody,
                               gravitionalParametersOfDeformingBodies,
                               deformedBodyEquatorialRadius,
                               nominalDisplacementLoveNumbers ),
    firstStepLatitudeDependenceTerms_( firstStepLatitudeDependenceTerms ),
    areFirstStepCorrectionsCalculated_( areFirstStepCorrectionsCalculated ),
    correctionLoveAndShidaNumbers_( correctionLoveAndShidaNumbers ), doodsonArgumentFunction_( doodsonArgumentFunction ),
    contextRotation_( contextRotation ), contextClock_( contextClock ), usesNativeUtcPhase_( !doodsonArgumentFunction )
{
    if( usesNativeUtcPhase_ )
    {
        doodsonArgumentFunction_ = calculateIers2010SolidEarthTideArguments;
    }
    assert( firstStepLatitudeDependenceTerms_.size( ) == 2 );
    assert( areFirstStepCorrectionsCalculated_.size( ) == 6 );
    assert( correctionLoveAndShidaNumbers_.size( ) == 6 );

    doodsonMultipliers_.resize( 0, 6 );
    tideCorrectionAmplitudes_.resize( 0, 4 );
    areFrequencyDependentTermsCalculated_ = 0;
    bool areDiurnalTermsCalculated = 0;

    if( !( diurnalFile == "" ) )
    {
        Eigen::MatrixXd rawDiurnalFile = input_output::readMatrixFromFile( diurnalFile );
        if( !( rawDiurnalFile.rows( ) == 0 ) )
        {
            assert( rawDiurnalFile.cols( ) == 10 );
            doodsonMultipliers_.resize( rawDiurnalFile.rows( ), 6 );
            doodsonMultipliers_ = rawDiurnalFile.block( 0, 0, rawDiurnalFile.rows( ), 6 );

            tideCorrectionAmplitudes_.resize( rawDiurnalFile.rows( ), 4 );
            tideCorrectionAmplitudes_ = rawDiurnalFile.block( 0, 6, rawDiurnalFile.rows( ), 4 );

            areFrequencyDependentTermsCalculated_ = 1;
            areDiurnalTermsCalculated = 1;
        }
    }
    if( !( longPeriodFile == "" ) )
    {
        Eigen::MatrixXd rawLongPeriodFile = input_output::readMatrixFromFile( longPeriodFile );
        if( !( rawLongPeriodFile.rows( ) == 0 ) )
        {
            {
                assert( rawLongPeriodFile.cols( ) == 10 );

                int numberOfExistingRows = 0;
                if( areDiurnalTermsCalculated )
                {
                    numberOfExistingRows = doodsonMultipliers_.rows( );
                    doodsonMultipliers_.conservativeResize( numberOfExistingRows + rawLongPeriodFile.rows( ), 6 );
                    tideCorrectionAmplitudes_.conservativeResize( numberOfExistingRows + rawLongPeriodFile.rows( ), 4 );
                }
                else
                {
                    doodsonMultipliers_.resize( rawLongPeriodFile.rows( ), 6 );
                    tideCorrectionAmplitudes_.resize( rawLongPeriodFile.rows( ), 6 );
                }
                doodsonMultipliers_.block( numberOfExistingRows, 0, rawLongPeriodFile.rows( ), 6 ) =
                        rawLongPeriodFile.block( 0, 0, rawLongPeriodFile.rows( ), 6 );

                tideCorrectionAmplitudes_.block( numberOfExistingRows, 0, rawLongPeriodFile.rows( ), 4 ) =
                        rawLongPeriodFile.block( 0, 6, rawLongPeriodFile.rows( ), 4 );

                areFrequencyDependentTermsCalculated_ = 1;
            }
        }
    }

    if( areFrequencyDependentTermsCalculated_ )
    {
        tideCorrectionAmplitudes_ = tideCorrectionAmplitudes_ / 1000.0;
    }
}

void Iers2010EarthDeformation::updateBodyProperties( const EarthDeformationEpoch& epoch )
{
    // Explicit mixed-clock evaluation: user ephemeris callbacks retain TDB,
    // while a factory-owned rotation callback can use the declared clock.
    const double tdb = epoch.getTdbSeconds( );
    const Eigen::Vector3d bodyPosition = deformedBodyStateFunction_( tdb ).head< 3 >( );
    const Eigen::Quaterniond rotation = contextRotation_ ? contextRotation_( epoch ) : deformedBodyRotationFunction_( tdb );
    const double gm = gravitionalParameterOfDeformedBody_( );
    for( int i = 0; i < numberOfBodies_; ++i )
    {
        gravitationalParameterRatios_[ i ] = gravitionalParametersOfDeformingBodies_[ i ]( ) / gm;
        relativeBodyStates_[ i ] = rotation * ( deformingBodyStateFunctions_[ i ]( tdb ).head< 3 >( ) - bodyPosition );
    }
    // These values may have used precise context rotation. Invalidate only the
    // local legacy scalar cache so a later scalar call cannot reuse them.
    currentTime_ = TUDAT_NAN;
}

Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    const auto station = std::make_shared< ground_stations::GroundStationState >(
            position,
            coordinate_conversions::cartesian_position,
            std::make_shared< SphericalBodyShapeModel >( deformedBodyEquatorialRadius_ ) );
    return calculateDisplacement( epoch, station );
}
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement( const EarthDeformationEpoch& epoch,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    if( !usesNativeUtcPhase_ ) return calculateDisplacement( epoch.getTdbSeconds( ), station );
    const auto native = contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch;
    updateBodyProperties( native );
    Eigen::Vector6d arguments = doodsonArgumentFunction_( static_cast< double >( native.getTt( ) ) );
    for( int i = 0; i < 6; ++i ) arguments[ i ] = std::fmod( arguments[ i ], 2.0 * mathematical_constants::PI );
    return calculateDisplacement( station, arguments );
}
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                                         const Eigen::Vector3d& position )
{
    if( !usesNativeUtcPhase_ ) return BodyDeformationModel::calculateDisplacementVelocity( epoch, position );
    const auto native = contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch;
    native.getTtRateWithRespectToTdb( );  // An interpolated clock without rates is explicitly unsupported.
    return calculateRoutedDisplacementVelocity( native, position );
}
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacementVelocity(
        const EarthDeformationEpoch& epoch,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    if( !usesNativeUtcPhase_ ) return BodyDeformationModel::calculateDisplacementVelocity( epoch, station );
    const auto native = contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch;
    native.getTtRateWithRespectToTdb( );
    return calculateRoutedDisplacementVelocity( native, station );
}
Eigen::Vector6d Iers2010EarthDeformation::calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return usesNativeUtcPhase_ ? calculateNativeDisplacementState( epoch, position )
                               : BodyDeformationModel::calculateDisplacementState( epoch, position );
}
Eigen::Vector6d Iers2010EarthDeformation::calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                                      const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return usesNativeUtcPhase_ ? calculateNativeDisplacementState( epoch, station )
                               : BodyDeformationModel::calculateDisplacementState( epoch, station );
}
void Iers2010EarthDeformation::validateContextDerivativeInterval( const EarthDeformationEpoch&,
                                                                  const EarthDeformationEpoch& lower,
                                                                  const EarthDeformationEpoch& upper ) const
{
    if( usesNativeUtcPhase_ )
        validateEarthDeformationUtcInterval( static_cast< double >( lower.getTt( ) ), static_cast< double >( upper.getTt( ) ) );
}

void Iers2010EarthDeformation::validateDisplacementDerivativeInterval( const double, const double lower, const double upper ) const
{
    if( usesNativeUtcPhase_ )
    {
        validateEarthDeformationUtcInterval( lower, upper );
    }
}

//! Function to calculate the site displacement at a given time and site position
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement( const double ephemerisTime,
                                                                 const Eigen::Vector3d& nominalSitePositionVector )
{
    // Create and calculate spherical position and local unit vectors from site Cartesian position and put in NominalGroundStationState
    // object.
    std::shared_ptr< ground_stations::GroundStationState > nominalGroundStationState =
            std::make_shared< ground_stations::GroundStationState >(
                    nominalSitePositionVector,
                    coordinate_conversions::cartesian_position,
                    std::make_shared< basic_astrodynamics::SphericalBodyShapeModel >( deformedBodyEquatorialRadius_ ) );

    // Call second overloaded version of function for further calculation.
    return calculateDisplacement( ephemerisTime, nominalGroundStationState );
}

//! Function to calculate the site displacement at a given time and nominal site state.
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    BasicTidalBodyDeformation::updateBodyProperties( ephemerisTime );

    // If Doodson arguments are required and not provided by user, calculate them.
    Eigen::Vector6d doodsonArguments = doodsonArgumentFunction_( ephemerisTime );
    for( int i = 0; i < 6; i++ )
    {
        doodsonArguments[ i ] = fmod( doodsonArguments[ i ], 2.0 * mathematical_constants::PI );
    }

    // Call third overloaded version of function for further calculations.
    return calculateDisplacement( nominalSiteState, doodsonArguments );
}

Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement(
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState,
        const Eigen::Vector6d& doodsonArguments )
{
    std::pair< Eigen::Vector3d, Eigen::Vector3d > planetCenteredAndLocalDisplacements = calculateFirstStepDisplacements( nominalSiteState );

    Eigen::Vector3d planetCenteredDisplacement = planetCenteredAndLocalDisplacements.first;
    Eigen::Vector3d siteDisplacement = planetCenteredAndLocalDisplacements.second;

    if( areFrequencyDependentTermsCalculated_ )
    {
        siteDisplacement += calculateFrequencyDependentDisplacementCorrections(
                doodsonMultipliers_, doodsonArguments, tideCorrectionAmplitudes_, nominalSiteState->getNominalSphericalPosition( ) );
    }

    std::vector< Eigen::Vector3d > enuUnitVectors = nominalSiteState->getEnuGeocentricUnitVectors( );

    planetCenteredDisplacement += enuUnitVectors[ 0 ] * siteDisplacement[ 0 ] + enuUnitVectors[ 1 ] * siteDisplacement[ 1 ] +
            enuUnitVectors[ 2 ] * siteDisplacement[ 2 ];

    return planetCenteredDisplacement;
}

//! Calculate site displacements due to solid Earth tide deformation according to first step of Section 7.1.1 of IERS 2010 Conventions
std::pair< Eigen::Vector3d, Eigen::Vector3d > Iers2010EarthDeformation::calculateFirstStepDisplacements(
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    // Initialize displacements calculated in local frame (Step I corrections and Step II) to zero.
    Eigen::Vector3d localFrameDisplacement = Eigen::Vector3d::Zero( );

    // Calculate degree two Love and Shida numbers for first part of Step I, including latitude-dependent correction
    double sineOfStationLatitude = std::sin( nominalSiteState->getNominalLatitude( ) );
    double latitudeCorrectionTerm = ( 3.0 * sineOfStationLatitude * sineOfStationLatitude - 1.0 ) / 2.0;
    std::map< int, std::pair< double, double > > currentLoveNumbers_ = displacementLoveNumbers_;
    currentLoveNumbers_[ 2 ].first += firstStepLatitudeDependenceTerms_[ 0 ] * latitudeCorrectionTerm;
    currentLoveNumbers_[ 2 ].second += firstStepLatitudeDependenceTerms_[ 1 ] * latitudeCorrectionTerm;

    // Initialize displacements calculated in planet-centered frame (ITRS) (First part of Step I) to zero.
    Eigen::Vector3d planetCenteredDisplacement =
            calculateBasicTicalDisplacement( nominalSiteState->getNominalCartesianPosition( ), currentLoveNumbers_ );

    // Loop over contributions of each of the bodies causing deformation of Earth (typically Moon and Sun).
    for( int i = 0; i < numberOfBodies_; i++ )
    {
        // Calculate first step corrections to site displacement.
        localFrameDisplacement += calculateFirstCorrectionStepDisplacements( relativeBodyStates_.at( i ),
                                                                             nominalSiteState->getNominalSphericalPosition( ),
                                                                             gravitationalParameterRatios_.at( i ),
                                                                             deformedBodyEquatorialRadius_,
                                                                             correctionLoveAndShidaNumbers_,
                                                                             areFirstStepCorrectionsCalculated_ );
    }

    return std::pair< Eigen::Vector3d, Eigen::Vector3d >( planetCenteredDisplacement, localFrameDisplacement );
}
// std::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector6d > >
// createDoodsonArgumentInterpolator( const double intervalStartEphemerisTime,
//                                    const double intervalEndEphemerisTime,
//                                    const double timeStepEphemerisTime,
//                                    const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction )
//{
//     std::map< double, Eigen::Vector6d > doodsonArgumentMap;
//     double currentTime = intervalStartEphemerisTime;

//    while( currentTime < intervalEndEphemerisTime )
//    {
//        doodsonArgumentMap[ currentTime ] = doodsonArgumentFunction( currentTime );
//        currentTime += timeStepEphemerisTime;
//    }

//    return std::make_shared< interpolators::CubicSplineInterpolator<
//            double, Eigen::Vector6d > >( doodsonArgumentMap );
//}

}  // namespace basic_astrodynamics

}  // namespace tudat

/*
Copyright (C) 2008
IERS Conventions Center

==================================
IERS Conventions Software License
==================================

NOTICE TO USER:

BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
WHICH APPLY TO ITS USE.

1. The Software is provided by the IERS Conventions Center ("the
   Center").

2. Permission is granted to anyone to use the Software for any
   purpose, including commercial applications, free of charge,
   subject to the conditions and restrictions listed below.

3. You (the user) may adapt the Software and its algorithms for your
   own purposes and you may distribute the resulting "derived work"
   to others, provided that the derived work complies with the
   following requirements:

   a) Your work shall be clearly identified so that it cannot be
      mistaken for IERS Conventions software and that it has been
      neither distributed by nor endorsed by the Center.

   b) Your work (including source code) must contain descriptions of
      how the derived work is based upon and/or differs from the
      original Software.

   c) The name(s) of all modified routine(s) that you distribute
      shall be changed.

   d) The origin of the IERS Conventions components of your derived
      work must not be misrepresented; you must not claim that you
      wrote the original Software.

   e) The source code must be included for all routine(s) that you
      distribute.  This notice must be reproduced intact in any
      source distribution.

4. In any published work produced by the user and which includes
   results achieved by using the Software, you shall acknowledge
   that the Software was used in obtaining those results.

5. The Software is provided to the user "as is" and the Center makes
   no warranty as to its use or performance.   The Center does not
   and cannot warrant the performance or results which the user may
   obtain by using the Software.  The Center makes no warranties,
   express or implied, as to non-infringement of third party rights,
   merchantability, or fitness for any particular purpose.  In no
   event will the Center be liable to the user for any consequential,
   incidental, or special damages, including any lost profits or lost
   savings, even if a Center representative has been advised of such
   damages, or for any claim by any third party.

Correspondence concerning IERS Conventions software should be
addressed as follows:

                   Gerard Petit
   Internet email: gpetit[at]bipm.org
   Postal address: IERS Conventions Center
                   Time, frequency and gravimetry section, BIPM
                   Pavillon de Breteuil
                   92312 Sevres  FRANCE

   or

                   Brian Luzum
   Internet email: brian.luzum[at]usno.navy.mil
   Postal address: IERS Conventions Center
                   Earth Orientation Department
                   3450 Massachusetts Ave, NW
                   Washington, DC 20392


*/
