#include "Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "InputOutput/writeDataToFile.h"

int main( )
{
    using namespace tudat;
    using namespace tudat::earth_orientation;

    double initialTime = 1.0E7;
    double finalTime = 1.1E7;
    std::vector< double > timeSteps;
    timeSteps.push_back( 30.0 );
    timeSteps.push_back( 60.0 );
    timeSteps.push_back( 150.0 );
    timeSteps.push_back( 500.0 );
    std::vector< boost::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6,1 > > > >
            interpolators;

    for( unsigned int i = 0; i < timeSteps.size( ); i++ )
    {
        interpolators.push_back( createInterpolatorForItrsToGcrsTransformation( initialTime,
                                                                                finalTime,
                                                                                timeSteps[ i ] ) );
    }

    boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator =
            createStandardEarthOrientationCalculator( );

    std::vector< std::map< double, Eigen::MatrixXd > > interpolationResults;
    interpolationResults.resize( timeSteps.size( ) );
    std::map< double, Eigen::MatrixXd > benchmarkResults;

    double evaluationStep = 5.0;

    double currentTime = initialTime + 6.0 * timeSteps[ timeSteps.size( ) - 1 ];

    while( currentTime < finalTime - 6.0 * timeSteps[ timeSteps.size( ) - 1 ] )
    {

        benchmarkResults[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
                    currentTime,
                    tdb_scale,
                    basic_astrodynamics::JULIAN_DAY_ON_J2000 );

        for( unsigned int i = 0; i < timeSteps.size( ); i++ )
        {
            interpolationResults[ i ][ currentTime ] = interpolators[ i ]->interpolate( currentTime );
        }

        currentTime += evaluationStep;
    }

    output::writeMatrixHistoryToFile( benchmarkResults, "benchmarkEarthOrientation.dat" );
    output::writeMatrixHistoryToFile( interpolationResults[ 0 ], "interpolatedEarthOrientation0.dat" );
    output::writeMatrixHistoryToFile( interpolationResults[ 1 ], "interpolatedEarthOrientation1.dat" );
    output::writeMatrixHistoryToFile( interpolationResults[ 2 ], "interpolatedEarthOrientation2.dat" );
    output::writeMatrixHistoryToFile( interpolationResults[ 3 ], "interpolatedEarthOrientation3.dat" );

}

