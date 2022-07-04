#include <boost/make_shared.hpp>

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

#include "tudat/basics/utilities.h"

namespace tudat
{

namespace basic_astrodynamics
{

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
        //Eq. 7.8
        trigonometricParts[ 0 ] = 3.0 * bodyDistanceInvSquared *
                ( - sineOfStationLatitude * sineOfStationLatitude * relativeBodyState[ 2 ] *
                  ( relativeBodyState[ 0 ] * cosineOfStationLongitude + relativeBodyState[ 1 ] * sineOfStationLongitude )
                  * stationUnitVectors[ 1 ] +
                  sineOfStationLatitude * cosineTwiceStationLatitude  * relativeBodyState[ 2 ] *
                  ( relativeBodyState[ 0 ] * sineOfStationLongitude - relativeBodyState[ 1 ] * cosineOfStationLongitude ) *
                  stationUnitVectors[ 0 ] );

    }
    else
    {

        trigonometricParts[ 0 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 1 ] == 1 )
    {
        //Eq. 7.9
        trigonometricParts[ 1 ] = bodyDistanceInvSquared * (
                    -1.5 * sineOfStationLatitude * cosineOfStationLatitude *
                    ( xySquareDifference * cosineOfTwiceStationLongitude + 2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] *
                      sineOfTwiceStationLongitude ) * stationUnitVectors[ 1 ] -
                    1.5 * sineOfStationLatitude * sineOfStationLatitude * cosineOfStationLatitude * (
                        xySquareDifference * sineOfTwiceStationLongitude - 2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] *
                        cosineOfTwiceStationLongitude ) * stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 1 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 2 ] == 1 )
    {
        //Eq. 7.10a
        trigonometricParts[ 2 ] = -3.0 * sineOfStationLatitude * cosineOfStationLatitude *
                relativeBodyState[ 2 ] * ( relativeBodyState[ 0 ] * sineOfStationLongitude -
                                           relativeBodyState[ 1 ] * cosineOfStationLongitude ) *
                bodyDistanceInvSquared * stationUnitVectors[ 2 ];
    }
    else
    {
        trigonometricParts[ 2 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 3 ] == 1 )
    {
        //Eq. 7.10b
        trigonometricParts[ 3 ] = bodyDistanceInvSquared * ( -3.0 * cosineTwiceStationLatitude * relativeBodyState[ 2 ] *
                                                             ( relativeBodyState[ 0 ] * sineOfStationLongitude - relativeBodyState[ 1 ] * cosineOfStationLongitude ) * stationUnitVectors[ 1 ]
                                                             -3.0 * sineOfStationLatitude * relativeBodyState[ 2 ] *
                                                             ( relativeBodyState[ 0 ] * cosineOfStationLongitude + relativeBodyState[ 1 ] * sineOfStationLongitude ) * stationUnitVectors[ 0 ] );
    }
    else
    {
        trigonometricParts[ 3 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 4 ] == 1 )
    {
        //Eq. 7.11a
        trigonometricParts[ 4 ] = -0.75 * cosineOfStationLatitude * cosineOfStationLatitude * bodyDistanceInvSquared *
                ( xySquareDifference * sineOfTwiceStationLongitude -
                  2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * cosineOfTwiceStationLongitude ) * stationUnitVectors[ 2 ];

    }
    else
    {
        trigonometricParts[ 4 ] = TUDAT_NAN * Eigen::Vector3d::Zero( );
    }

    if( calculateTerms[ 5 ] == 1 )
    {
        //Eq. 7.11b
        trigonometricParts[ 5 ] = bodyDistanceInvSquared *
                ( 1.5 * sineOfStationLatitude * cosineOfStationLatitude * (
                      xySquareDifference * sineOfTwiceStationLongitude -
                      2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * cosineOfTwiceStationLongitude )
                  * stationUnitVectors[ 1 ] -
                  1.5 * cosineOfStationLatitude * (
                      xySquareDifference * cosineOfTwiceStationLongitude +
                      2.0 * relativeBodyState[ 0 ] * relativeBodyState[ 1 ] * sineOfTwiceStationLongitude )
                  * stationUnitVectors[ 0 ] );
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
    displacement << ( amplitudes( 2 ) * cosineOfArgument - amplitudes( 3 ) * sineOfArgument ) *
                    std::sin( stationLatitude ),
            ( amplitudes( 2 ) * sineOfArgument + amplitudes( 3 ) * cosineOfArgument ) *
            std::cos( 2.0 * stationLatitude ),
            ( amplitudes( 0 ) * sineOfArgument + amplitudes( 1 ) * cosineOfArgument ) *
            std::sin( 2.0 * stationLatitude );
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
    displacement << 0.0, ( amplitudes( 2 ) * cosineOfArgument + amplitudes( 3 ) * sineOfArgument )  *
            std::sin( 2.0 * stationLatitude ), ( amplitudes( 0 ) * cosineOfArgument + amplitudes( 1 ) * sineOfArgument ) *
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
                        tideArguments( i ), stationSphericalCoordinates.y( ),
                        amplitudes.block( i, 0, 1, 4 ).transpose( ) );
        }

        // If tide is not long-period (i.e. not independent of mean Lunar time), use specific diurnal correction method.
        else if( doodsonMultipliers( i, 0 ) == 1 )
        {
            displacement += calculateDiurnalFrequencyDependentDisplacementCorrection(
                        tideArguments( i ), stationSphericalCoordinates.z( ), stationSphericalCoordinates.y( ),
                        amplitudes.block( i, 0, 1, 4 ).transpose( ) );
        }
        else
        {
            std::cerr<<"Warning, frequency dependent love number value corrections only available for diurnal and long period tides."<<std::endl;
        }
    }
    return displacement;
}

std::shared_ptr< Iers2010EarthDeformation > createDefaultEarthIers2010DeformationModel(
        const std::shared_ptr< ephemerides::Ephemeris > earthEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > lunarEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > solarEphemeris,
        const std::shared_ptr< ephemerides::RotationalEphemeris > earthRotation,
        const std::function< double( ) > gravitionalParametersOfEarth,
        const std::function< double( ) > gravitionalParametersOfMoon,
        const std::function< double( ) > gravitionalParametersOfSun)
{
    using namespace tudat::ephemerides;

    std::vector< std::function< double( ) > > gravitationalParameters;
    gravitationalParameters.resize( 2 );
    gravitationalParameters[ 0 ] = gravitionalParametersOfMoon;
    gravitationalParameters[ 1 ] = gravitionalParametersOfSun;

    std::vector< std::function< Eigen::Vector6d( const double ) > > ephemerides;
    ephemerides.resize( 2 );
    ephemerides[ 0 ] = std::bind( &Ephemeris::getCartesianState, lunarEphemeris, std::placeholders::_1 );
    ephemerides[ 1 ] = std::bind( &Ephemeris::getCartesianState, solarEphemeris, std::placeholders::_1 );

    double equatorialRadius = 6378136.6;

    using namespace iers_2010_parameters;

    std::map< int, std::pair< double, double > > nominalDisplacementLoveNumbers;
    nominalDisplacementLoveNumbers[ 2 ] = std::make_pair( PRINCIPAL_DEGREE_TWO_LOVE_NUMBER, PRINCIPAL_DEGREE_TWO_SHIDA_NUMBER );
    nominalDisplacementLoveNumbers[ 3 ] = std::make_pair( PRINCIPAL_DEGREE_THREE_LOVE_NUMBER, PRINCIPAL_DEGREE_THREE_SHIDA_NUMBER );

    std::vector< double > latitudeTerms;
    latitudeTerms.resize( 2 );
    latitudeTerms[ 0 ] = DEGREE_TWO_LATITUDE_LOVE_NUMBER;
    latitudeTerms[ 1 ] = DEGREE_TWO_LATITUDE_SHIDA_NUMBER;

    std::vector< bool > areTermsCalculated;
    areTermsCalculated.resize( 6 );
    for( int i = 0; i < 6 ; i++ )
    {
        areTermsCalculated[ i ] = 1;
    }

    std::vector< double > correctionNumbers;
    correctionNumbers.resize( 6 );
    correctionNumbers[ 2 ] = IMAGINARY_DEGREE_TWO_DIURNAL_LOVE_NUMBER;
    correctionNumbers[ 3 ] = IMAGINARY_DEGREE_TWO_DIURNAL_SHIDA_NUMBER;
    correctionNumbers[ 4 ] = IMAGINARY_DEGREE_TWO_SEMIDIURNAL_LOVE_NUMBER;
    correctionNumbers[ 5 ] = IMAGINARY_DEGREE_TWO_SEMIDIURNAL_SHIDA_NUMBER;
    correctionNumbers[ 0 ] = DEGREE_TWO_DIURNAL_TOROIDAL_LOVE_NUMBER;
    correctionNumbers[ 1 ] = DEGREE_TWO_SEMIDIURNAL_TOROIDAL_LOVE_NUMBER;

    std::string longPeriodFile = paths::getEarthOrientationDataFilesPath( ) +
            "longPeriodDisplacementFrequencyDependence.txt";
    std::string diurnalFile = paths::getEarthOrientationDataFilesPath( ) +
            "diurnalDisplacementFrequencyDependence.txt";

    std::shared_ptr< Iers2010EarthDeformation > deformationModel = std::make_shared< Iers2010EarthDeformation >
            ( std::bind( &Ephemeris::getCartesianState, earthEphemeris, std::placeholders::_1 ),
              ephemerides,
              std::bind( &RotationalEphemeris::getRotationToTargetFrame, earthRotation, std::placeholders::_1 ),
              gravitionalParametersOfEarth, gravitationalParameters,
              equatorialRadius, nominalDisplacementLoveNumbers, latitudeTerms, areTermsCalculated,
              correctionNumbers, diurnalFile, longPeriodFile );

    return deformationModel;

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
        const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
    BasicTidalBodyDeformation( deformedBodyEphemeris, deformingBodyEphemerides, deformedBodyRotation,
                               gravitionalParameterOfDeformedBody,
                               gravitionalParametersOfDeformingBodies, deformedBodyEquatorialRadius,
                               nominalDisplacementLoveNumbers ),
    firstStepLatitudeDependenceTerms_( firstStepLatitudeDependenceTerms ),
    areFirstStepCorrectionsCalculated_( areFirstStepCorrectionsCalculated ),
    correctionLoveAndShidaNumbers_( correctionLoveAndShidaNumbers ),
    doodsonArgumentFunction_( doodsonArgumentFunction )
{
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
                doodsonMultipliers_.block( numberOfExistingRows, 0, rawLongPeriodFile.rows( ), 6 )
                        = rawLongPeriodFile.block( 0, 0, rawLongPeriodFile.rows( ), 6 );

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

//! Function to calculate the site displacement at a given time and site position
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement( const double ephemerisTime,
                                                                     const Eigen::Vector3d& nominalSitePositionVector )
{

    // Create and calculate spherical position and local unit vectors from site Cartesian position and put in NominalGroundStationState object.
    std::shared_ptr< ground_stations::GroundStationState > nominalGroundStationState =
            std::make_shared< ground_stations::GroundStationState >(
                nominalSitePositionVector, coordinate_conversions::cartesian_position,
                std::make_shared< basic_astrodynamics::SphericalBodyShapeModel >( deformedBodyEquatorialRadius_ ) );

    // Call second overloaded version of function for further calculation.
    return calculateDisplacement( ephemerisTime, nominalGroundStationState );
}

//! Function to calculate the site displacement at a given time and nominal site state.
Eigen::Vector3d Iers2010EarthDeformation::calculateDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    updateBodyProperties( ephemerisTime );

    // If Doodson arguments are required and not provided by user, calculate them.
    Eigen::Vector6d doodsonArguments = doodsonArgumentFunction_( ephemerisTime );
    for( int i = 0; i < 6 ; i++ )
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

    std::pair< Eigen::Vector3d, Eigen::Vector3d > planetCenteredAndLocalDisplacements =
            calculateFirstStepDisplacements( nominalSiteState );

    Eigen::Vector3d planetCenteredDisplacement = planetCenteredAndLocalDisplacements.first;
    Eigen::Vector3d siteDisplacement = planetCenteredAndLocalDisplacements.second;

    if( areFrequencyDependentTermsCalculated_ )
    {
        siteDisplacement += calculateFrequencyDependentDisplacementCorrections(
                    doodsonMultipliers_, doodsonArguments,
                    tideCorrectionAmplitudes_, nominalSiteState->getNominalSphericalPosition( ) );
    }

    std::vector< Eigen::Vector3d > enuUnitVectors = nominalSiteState->getEnuGeocentricUnitVectors( );

    planetCenteredDisplacement += enuUnitVectors[ 0 ] * siteDisplacement[ 0 ] +
            enuUnitVectors[ 1 ] * siteDisplacement[ 1 ] +
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
    Eigen::Vector3d planetCenteredDisplacement = calculateBasicTicalDisplacement(
                nominalSiteState->getNominalCartesianPosition( ), currentLoveNumbers_ );

    // Loop over contributions of each of the bodies causing deformation of Earth (typically Moon and Sun).
    for( int i = 0; i < numberOfBodies_; i++ )
    {
        // Calculate first step corrections to site displacement.
        localFrameDisplacement += calculateFirstCorrectionStepDisplacements(
                    relativeBodyStates_.at( i ),
                    nominalSiteState->getNominalSphericalPosition( ),
                    gravitationalParameterRatios_.at( i ),
                    deformedBodyEquatorialRadius_,
                    correctionLoveAndShidaNumbers_,
                    areFirstStepCorrectionsCalculated_ );
    }

    return std::pair< Eigen::Vector3d, Eigen::Vector3d >( planetCenteredDisplacement, localFrameDisplacement );
}
//std::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector6d > >
//createDoodsonArgumentInterpolator( const double intervalStartEphemerisTime,
//                                   const double intervalEndEphemerisTime,
//                                   const double timeStepEphemerisTime,
//                                   const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction )
//{
//    std::map< double, Eigen::Vector6d > doodsonArgumentMap;
//    double currentTime = intervalStartEphemerisTime;

//    while( currentTime < intervalEndEphemerisTime )
//    {
//        doodsonArgumentMap[ currentTime ] = doodsonArgumentFunction( currentTime );
//        currentTime += timeStepEphemerisTime;
//    }

//    return std::make_shared< interpolators::CubicSplineInterpolator<
//            double, Eigen::Vector6d > >( doodsonArgumentMap );
//}

}

}
