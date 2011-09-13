/*! \file tabulatedAtmosphere.cpp
 *    Source file that defines the tabulated atmosphere
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 Juny, 2011
 *    Last modified     : 22 July, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110721    J. Melman         Comments, variable names,
 *                                  and consistency modified.
 *      110722    F.M. Engelen      Removed setRelativePath function.
 */

// Include statements.
#include "tabulatedAtmosphere.h"

//! Initialise the atmosphere table reader.
void TabulatedAtmosphere::initialize( string atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    // Set the file reader.
    textFileReader.setRelativePath( relativePath_ );
    textFileReader.setFileName( atmosphereTableFile_ );
    textFileReader.openFile( );

    string skipcharacter = "%";
    textFileReader.skipLinesStartingWithCharacter( skipcharacter );

    // Read and store the data file.
    textFileReader.readAndStoreData( );
    containerOfAtmosphereTableFileData = textFileReader.getContainerOfData( );

    // Check whether data is present in the file.
    if ( containerOfAtmosphereTableFileData.size( ) < 1 )
    {
        std::cerr << "The atmosphere table file is empty." << std::endl;
        std::cerr << atmosphereTableFile_ << std::endl;
    }

    // Initialize vectors.
    altitudeData_ = VectorXd( containerOfAtmosphereTableFileData.size( ) );
    densityData_ = VectorXd( containerOfAtmosphereTableFileData.size( ) );
    pressureData_= VectorXd( containerOfAtmosphereTableFileData.size( ) );
    temperatureData_ = VectorXd( containerOfAtmosphereTableFileData.size( ) );
    stringstream lineStringStream( stringstream::in | stringstream::out );

    int index = 0;

    // Loop through all the strings stored in the container and store the data
    // in the right VectorXd.
    for ( iteratorContainerOfData_  = containerOfAtmosphereTableFileData.begin( );
          iteratorContainerOfData_ != containerOfAtmosphereTableFileData.end( );
          iteratorContainerOfData_++ )
    {
        lineStringStream << iteratorContainerOfData_->second;
        lineStringStream
                >> altitudeData_( index )
                >> densityData_( index )
                >> pressureData_( index )
                >> temperatureData_ ( index );
        index++;
    }

    cubicSplineInterpolationForDensity_.initializeCubicSplineInterpolation(
            altitudeData_, densityData_ );
    cubicSplineInterpolationForPressure_.initializeCubicSplineInterpolation(
            altitudeData_, pressureData_ );
    cubicSplineInterpolationForTemperature_.initializeCubicSplineInterpolation(
            altitudeData_, temperatureData_ );
}

// End of file.
