/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110620    F.M. Engelen      File created.
 *      110721    J. Melman         Comments, variable names, and consistency modified.
 *      110722    F.M. Engelen      Removed setRelativePath function.
 *
 *    References
 *
 */

#include <sstream>
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{

//! Initialize atmosphere table reader.
void TabulatedAtmosphere::initialize( std::string atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    Eigen::MatrixXd containerOfAtmosphereTableFileData = tudat::input_output::readMatrixFromFile(
                atmosphereTableFile_, " \t", "%" );

    // Check whether data is present in the file.
    if ( containerOfAtmosphereTableFileData.rows( ) < 1 || containerOfAtmosphereTableFileData.cols( ) < 1 )
    {
        std::cerr << "The atmosphere table file is empty." << std::endl;
        std::cerr << atmosphereTableFile_ << std::endl;
    }

    // Initialize vectors.
    altitudeData_ = Eigen::VectorXd( containerOfAtmosphereTableFileData.rows( ) );
    densityData_ = Eigen::VectorXd( containerOfAtmosphereTableFileData.rows( ) );
    pressureData_= Eigen::VectorXd( containerOfAtmosphereTableFileData.rows( ) );
    temperatureData_ = Eigen::VectorXd( containerOfAtmosphereTableFileData.rows( ) );


    // Loop through all the strings stored in the container and store the data
    // in the right Eigen::VectorXd.
    for ( int i = 0; i < containerOfAtmosphereTableFileData.rows( ); i++  )
    {
        altitudeData_( i ) = containerOfAtmosphereTableFileData( i, 0 );
        densityData_( i ) = containerOfAtmosphereTableFileData( i, 1 );
        pressureData_( i ) = containerOfAtmosphereTableFileData( i, 2 );
        temperatureData_( i ) = containerOfAtmosphereTableFileData( i, 3 );
    }


    cubicSplineInterpolationForDensity_.initializeCubicSplineInterpolation(
                altitudeData_, densityData_ );
    cubicSplineInterpolationForPressure_.initializeCubicSplineInterpolation(
                altitudeData_, pressureData_ );
    cubicSplineInterpolationForTemperature_.initializeCubicSplineInterpolation(
                altitudeData_, temperatureData_ );
}

} // namespace tudat
