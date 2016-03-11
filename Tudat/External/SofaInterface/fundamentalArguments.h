#ifndef FUNDAMENTALARGUMENTS_H
#define FUNDAMENTALARGUMENTS_H

#include <vector>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

extern "C"
{
    #include "sofa/src/sofa.h"
    #include "sofa/src/sofam.h"
}

namespace tudat
{

namespace sofa_interface
{

//! Function to calculate mean elongation of the Moon from the Sun
/*!
 *  Function to calculate mean elongation of the Moon from the Sun fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean elongation of the Moon from the Sun
 */
double calculateMeanElongationOfMoonFromSun( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean argument of latitude of the Moon
/*!
 *  Function to calculate mean argument of latitude of the Moon fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean argument of latitude of the Moon
 */
double calculateMeanArgumentOfLatitudeOfMoon( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean anomaly of the Moon
/*!
 *  Function to calculate mean anomaly of the Moon fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean anomaly of the Moon
 */
double calculateMeanAnomalyOfMoon( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean anomaly of the Sun
/*!
 *  Function to calculate mean anomaly of the Sun fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean anomaly of the Sun
 */
double calculateMeanAnomalyOfSun( const double julianCenturiesSinceJ2000 );

//! Function to calculate longitude of the Moon's ascending node
/*!
 *  Function to calculate longitude of the Moon's ascending node fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Longitude of the Moon's ascending node
 */
double calculateLongitudeOfMoonsAscendingNode( const double julianCenturiesSinceJ2000 );

//Eigen::Vector3d calculateOceanTideArgumets( const double ephemerisTime );

Eigen::Matrix< double, 5, 1 > calculateDelaunayFundamentalArguments(
        const double ephemerisTime );

basic_mathematics::Vector6d  calculateFundamentalArguments(
        const double ephemerisTime );

basic_mathematics::Vector6d calculateDoodsonFundamentalArguments(
        const double ephemerisTime, const double terrestrialTime, const double universalTime1, const double ttAndUt1EpochShift );

basic_mathematics::Vector6d calculateDoodsonFundamentalArguments( const double ephemerisTime );

} // namespace sofa_interfaces

} // namespace tudat

#endif // FUNDAMENTALARGUMENTS_H
