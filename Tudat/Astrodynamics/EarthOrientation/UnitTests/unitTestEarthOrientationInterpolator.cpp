#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"

int main( )
{
    using namespace tudat;
    using namespace tudat::earth_orientation;

    double initialTime = 1.0E7;
    double finalTime = initialTime + 10000.0;
    std::vector< double > timeSteps;
    timeSteps.push_back( 30.0 );
    timeSteps.push_back( 60.0 );
    timeSteps.push_back( 150.0 );
    timeSteps.push_back( 500.0 );
    std::vector< boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > > >
            interpolators;

    for( unsigned int i = 0; i < timeSteps.size( ); i++ )
    {
        interpolators.push_back( createInterpolatorForItrsToGcrsAngles( initialTime,
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
                    basic_astrodynamics::tdb_scale );

        for( unsigned int i = 0; i < timeSteps.size( ); i++ )
        {
            interpolationResults[ i ][ currentTime ] = interpolators[ i ]->interpolate( currentTime );
        }

        currentTime += evaluationStep;
    }

    for( std::map< double, Eigen::MatrixXd >::const_iterator angleIterator = benchmarkResults.begin( );
         angleIterator != benchmarkResults.end( ); angleIterator++ )
    {
        std::cout<<( angleIterator->second - interpolationResults.at( 0 ).at( angleIterator->first ) ).transpose( )<<std::endl<<std::endl;
        std::cout<<( angleIterator->second - interpolationResults.at( 1 ).at( angleIterator->first ) ).transpose( )<<std::endl<<std::endl;
        std::cout<<( angleIterator->second - interpolationResults.at( 2 ).at( angleIterator->first ) ).transpose( )<<std::endl<<std::endl;
        std::cout<<( angleIterator->second - interpolationResults.at( 3 ).at( angleIterator->first ) ).transpose( )<<std::endl<<std::endl<<std::endl<<std::endl;

    }

//    output::writeMatrixHistoryToFile( benchmarkResults, "benchmarkEarthOrientation.dat" );
//    output::writeMatrixHistoryToFile( interpolationResults[ 0 ], "interpolatedEarthOrientation0.dat" );
//    output::writeMatrixHistoryToFile( interpolationResults[ 1 ], "interpolatedEarthOrientation1.dat" );
//    output::writeMatrixHistoryToFile( interpolationResults[ 2 ], "interpolatedEarthOrientation2.dat" );
//    output::writeMatrixHistoryToFile( interpolationResults[ 3 ], "interpolatedEarthOrientation3.dat" );

}

