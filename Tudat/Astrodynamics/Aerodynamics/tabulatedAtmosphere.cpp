/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
