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
 *      110721    J. Melman         Comments, file names, and consistency modified.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// The provided USSA1976 table file, generated with the pascal file, has a small error
// which can be observed at the pressure at sea level. This in 101320 in the file
// but should be 101325. If this error is not acceptable, an other table file should be
// used.
// 

#ifndef TUDAT_TABULATED_ATMOSPHERE_H
#define TUDAT_TABULATED_ATMOSPHERE_H

#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

#include <Eigen/Core>
#include <iostream>
#include <string>
#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolation.h"

namespace tudat
{

//! Tabulated atmosphere class.
/*!
 * Tabulated atmospheres class, for example US1976. The default path
 * from which the files are obtained is:
 * Astrodynamics/EnvironmentModels/AtmosphereTables
 * NOTE: for the moment it only works for tables with 4 columns:
 * altitude, density, pressure and temperature.
 */
class TabulatedAtmosphere : public AtmosphereModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    TabulatedAtmosphere( ) :
        relativeDirectoryPath_( "Astrodynamics/Aerodynamics/AtmosphereTables/" ),
        atmosphereTableFile_( "" ) { }

    //! Initialize atmosphere table reader.
    /*!
     * Initializes the atmosphere table reader.
     * \param atmosphereTableFile The name of the atmosphere table.
     */
    void initialize( std::string atmosphereTableFile );

    //! Get atmosphere table file name.
    /*!
     * Returns atmosphere table file name.
     * \return The atmosphere table file.
     */
    std::string getAtmosphereTableFile( ) { return atmosphereTableFile_; }

    //! Get relative directory path.
    /*!
     * Returns relative directory path.
     * \return Relative directory path.
     */
    std::string getRelativeDirectoryPath( ) { return relativeDirectoryPath_; }

    //! Get local density.
    /*!
     * Returns the local density parameter of the atmosphere in kg per meter^3.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric density.
     */
    double getDensity( double altitude, double longitude = 0.0,
                       double latitude = 0.0, double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return cubicSplineInterpolationForDensity_.interpolate( altitude );
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric pressure.
     */
    double getPressure( double altitude, double longitude = 0.0,
                        double latitude = 0.0, double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return cubicSplineInterpolationForPressure_.interpolate( altitude );
    }

    //! Get local temperature.
    /*!
     * Returns the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return Atmospheric temperature.
     */
    double getTemperature( double altitude, double longitude = 0.0,
                           double latitude = 0.0, double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( longitude );
        TUDAT_UNUSED_PARAMETER( latitude );
        TUDAT_UNUSED_PARAMETER( time );
        return cubicSplineInterpolationForTemperature_.interpolate( altitude );
    }

protected:

private:

    //! The relative directory path.
    /*!
     *  The relative path directory path.
     */
    std::string relativeDirectoryPath_;

    //! The file name of the atmosphere table.
    /*!
     *  The file name of the atmosphere table.
     */
    std::string atmosphereTableFile_;

    //! Vector containing the altitude.
    /*!
     *  Vector containing the altitude.
     */
    Eigen::VectorXd altitudeData_;

    //! Vector containing the density data as a function of the altitude.
    /*!
     *  Vector containing the density data as a function of the altitude.
     */
    Eigen::VectorXd densityData_;

    //! Vector containing the pressure data as a function of the altitude.
    /*!
     *  Vector containing the pressure data as a function of the altitude.
     */
    Eigen::VectorXd pressureData_;

    //! Vector containing the temperature data as a function of the altitude.
    /*!
     *  Vector containing the temperature data as a function of the altitude.
     */
    Eigen::VectorXd temperatureData_;

    //! Cubic spline interpolation for density.
    /*!
     *  Cubic spline interpolation for density.
     */
    mathematics::interpolators::CubicSplineInterpolation cubicSplineInterpolationForDensity_;

    //! Cubic spline interpolation for pressure.
    /*!
     *  Cubic spline interpolation for pressure.
     */
    mathematics::interpolators::CubicSplineInterpolation cubicSplineInterpolationForPressure_;

    //! Cubic spline interpolation for temperature.
    /*!
     *  Cubic spline interpolation for temperature.
     */
    mathematics::interpolators::CubicSplineInterpolation cubicSplineInterpolationForTemperature_;
};

} // namespace tudat

#endif // TUDAT_TABULATED_ATMOSPHERE_H
