#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{

namespace earth_orientation
{

template< >
double ShortPeriodEarthOrientationCorrectionCalculator< double >::sumCorrectionTerms( const Eigen::Vector6d& arguments )
{

    using std::cos;
    using std::sin;

    // Initialize polar motion corrections to zero.
    double currentCorrection = 0;

    double tideAngle = 0.0;

    for( unsigned int i = 0; i < argumentAmplitudes_.size( ); i++ )
    {
        Eigen::MatrixXd currentArgumentMultipliers = argumentMultipliers_.at( i );
        Eigen::MatrixXd currentAmplitudes = argumentAmplitudes_.at( i );

        // Iterte over all libration corrections.
        for( int i = 0; i < currentAmplitudes.rows( ); i++ )
        {
            // Calculate current phase angle of tide.
            tideAngle = ( arguments.transpose( ) *
                          currentArgumentMultipliers.block( i, 0, 1, 6 ).transpose( ) )( 0, 0 );

            currentCorrection += currentAmplitudes( i, 0 ) * sin( tideAngle ) +
                    currentAmplitudes( i, 1 ) * cos( tideAngle );
        }
    }


    return currentCorrection;
}

template< >
Eigen::Vector2d ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d >::sumCorrectionTerms( const Eigen::Vector6d& arguments )
{

    using std::cos;
    using std::sin;

    // Initialize polar motion corrections to zero.
    Eigen::Vector2d currentCorrection = Eigen::Vector2d::Zero( );

    double tideAngle = 0.0;

    for( unsigned int i = 0; i < argumentAmplitudes_.size( ); i++ )
    {
        Eigen::MatrixXd currentArgumentMultipliers = argumentMultipliers_.at( i );
        Eigen::MatrixXd currentAmplitudes = argumentAmplitudes_.at( i );

        // Iterte over all libration corrections.
        for( int i = 0; i < currentAmplitudes.rows( ); i++ )
        {
            // Calculate current phase angle of tide.
            tideAngle = ( arguments.transpose( ) *
                          currentArgumentMultipliers.block( i, 0, 1, 6 ).transpose( ) )( 0, 0 );

            // Calculate and add sine and cosine terms
            currentCorrection.x( ) += currentAmplitudes( i, 0 ) * sin( tideAngle ) +
                    currentAmplitudes( i, 1 ) * cos( tideAngle );
            currentCorrection.y( ) += currentAmplitudes( i, 2 ) * sin( tideAngle ) +
                    currentAmplitudes( i, 3 ) * cos( tideAngle );
        }
    }


    return currentCorrection;
}

boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > getDefaultUT1CorrectionCalculator(
        const double minimumAmplitude )
{
    return boost::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                1.0E-6, minimumAmplitude,
                std::vector< std::string >{  tudat::input_output::getDataFilesRootPath( ) +  "EarthOrientation/utcLibrationAmplitudes.txt",
                                             tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/utcOceanTidesAmplitudes.txt" },
                std::vector< std::string >{  tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/utcLibrationDoodsonMultipliers.txt",
                                             tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/utcOceanTidesDoodsonMultipliers.txt" },
                boost::bind( &sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst, _1 ) );
}

boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getDefaultPolarMotionCorrectionCalculator(
        const double minimumAmplitude )
{
    return boost::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), minimumAmplitude,
                std::vector< std::string >{ tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionLibrationAmplitudes.txt",
                                            tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionOceanTidesAmplitudes.txt", },
                std::vector< std::string >{ tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionLibrationDoodsonMultipliers.txt",
                                            tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionOceanTidesDoodsonMultipliers.txt" },
                boost::bind( &sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst, _1 ) );

}


}

}
