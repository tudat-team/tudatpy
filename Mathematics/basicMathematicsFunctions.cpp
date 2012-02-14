/*! \file basicMathematicsFunctions.cpp
 *    Source file that defines the basicMathematicsFunctions namespace,
 *    containing all basic functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 12
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Author            : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Author            : D. Gondelach
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.gondelach@student.tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : T. Secretin
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : T.A.LeitePintoSecretin@student.tudelft.nl
 *
 *    Date created      : 3 September, 2010
 *    Last modified     : 18 January, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
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
 *      100903    K. Kumar          File created.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToExponentPower() function.
 *      102410    D. Dirkx          Minor comment changes as code check.
 *      101213    K. Kumar          Bugfix raiseToIntegerExponent(); renamed raiseToIntegerPower().
 *                                  Added computeAbsoluteValue() functions.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation().
 *      110111    J. Melman         Added computeModulo() function.
 *      110411    K. Kumar          Added convertCartesianToSpherical() function.
 *      110606    J. Melman         Removed possible singularity from
 *                                  convertCartesianToSpherical.
 *      110707    K. Kumar          Added computeSampleMean(), computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120118    D. Gondelach      Added convertCylindricalToCartesianCoordinates,
 *                                  convertCylindricalToCartesianState,
 *                                  convertCartesianToCylindricalCoordinates and
 *                                  convertCartesianToCylindricalState functions.
 *                                  Removed convertCylindricalToCartesian function.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include "Basics/basicFunctions.h"
#include "Mathematics/basicMathematicsFunctions.h"

//! Tudat library namespace.
namespace tudat
{

//! Mathematics namespace.
namespace mathematics
{

// Using declarations.
using std::map;
using std::vector;
using std::accumulate;
using std::pow;

//! Get global random number generator.
globalRandomNumberGeneratorType& getGlobalRandomNumberGenerator( )
{
  static globalRandomNumberGeneratorType globalRandomNumberGenerator(
              static_cast< unsigned int >( std::time( 0 ) ) );
  return globalRandomNumberGenerator;
}

//! Compute linear interpolation.
double computeLinearInterpolation( Eigen::VectorXd& sortedIndependentVariables,
                                   Eigen::VectorXd& associatedDependentVariables,
                                   double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestNeighbor;
    double locationTargetIndependentVariableValueInInterval;

    // Compute nearest neighbor in sorted vector of independent variables.
    // Result is always to the left of the target independent variable value.
    nearestNeighbor = basic_functions
                      ::computeNearestLeftNeighborUsingBinarySearch(
            sortedIndependentVariables,
            targetIndependentVariableValue );

    // Compute location of target independent variable value in interval
    // between nearest neighbors.
    locationTargetIndependentVariableValueInInterval
            = ( targetIndependentVariableValue
              - sortedIndependentVariables[ nearestNeighbor ] )
             / ( sortedIndependentVariables[ nearestNeighbor + 1 ]
                 - sortedIndependentVariables[ nearestNeighbor ] );

    // Return the computed value of the dependent variable.
    return ( associatedDependentVariables[ nearestNeighbor ]
             * ( 1 - locationTargetIndependentVariableValueInInterval )
             + associatedDependentVariables[ nearestNeighbor + 1 ]
             * locationTargetIndependentVariableValueInInterval );
}

//! Compute linear interpolation.
Eigen::VectorXd computeLinearInterpolation(
        std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, Eigen::VectorXd >::iterator mapIteratorIntervalLeft;
    std::map < double, Eigen::VectorXd >::iterator mapIteratorIntervalRight;

    // Compute nearest neighbor in map of data.
    // Result is always to the left of the target independent variable value.
    nearestLeftNeighbor = basic_functions::
                          computeNearestLeftNeighborUsingBinarySearch(
                                  sortedIndepedentAndDependentVariables,
                                  targetIndependentVariableValue );

    // Compute location of target independent variable value in interval
    // between nearest neighbors.
    mapIteratorIntervalLeft = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalLeft, nearestLeftNeighbor );
    mapIteratorIntervalRight = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalRight, nearestLeftNeighbor + 1 );
    locationTargetIndependentVariableValueInInterval
            = ( targetIndependentVariableValue
              - mapIteratorIntervalLeft->first )
             / ( mapIteratorIntervalRight->first
                 - mapIteratorIntervalLeft->first );

    // Return the computed value of the dependent variable.
    return ( mapIteratorIntervalLeft->second
             * ( 1 - locationTargetIndependentVariableValueInInterval )
             + mapIteratorIntervalRight->second
             * locationTargetIndependentVariableValueInInterval );
}

//! Compute linear interpolation.
State* computeLinearInterpolation(
        std::map < double, State* >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue )
{
    // Declare local variables.
    // Declare nearest neighbor.
    int nearestLeftNeighbor;

    // Pointer to state.
    State* pointerToState_ = new State;

    // Declare location of target independent variable value in interval.
    double locationTargetIndependentVariableValueInInterval;

    // Declare map iterators
    std::map < double, State* >::iterator mapIteratorIntervalLeft;
    std::map < double, State* >::iterator mapIteratorIntervalRight;

    // Compute nearest neighbor in map of data.
    // Result is always to the left of the target independent variable value.
    nearestLeftNeighbor = basic_functions::
                          computeNearestLeftNeighborUsingBinarySearch(
                                  sortedIndepedentAndDependentVariables,
                                  targetIndependentVariableValue );

    // Compute location of target independent variable value in interval
    // between nearest neighbors.
    mapIteratorIntervalLeft = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalLeft, nearestLeftNeighbor );
    mapIteratorIntervalRight = sortedIndepedentAndDependentVariables.begin( );
    advance( mapIteratorIntervalRight, nearestLeftNeighbor + 1 );
    locationTargetIndependentVariableValueInInterval
            = ( targetIndependentVariableValue
              - mapIteratorIntervalLeft->first )
             / ( mapIteratorIntervalRight->first
                 - mapIteratorIntervalLeft->first );

    // Set vector in pointer to State object to the computed value of the
    // dependent variable.

    pointerToState_->state
            = mapIteratorIntervalLeft->second->state
              * ( 1 - locationTargetIndependentVariableValueInInterval )
              + mapIteratorIntervalRight->second->state
              * locationTargetIndependentVariableValueInInterval;

    // Return pointer to State object.
    return pointerToState_;
}

//! Convert spherical to cartesian coordinates.
void convertSphericalToCartesian( double radius, double azimuthAngle, double zenithAngle,
                                  Eigen::VectorXd& cartesianCoordinates )
{
    // Declaring sine and cosine which have multiple usages to save computation time.
    double cosineOfAzimuthAngle = std::cos( azimuthAngle );
    double sineOfZenithAngle = std::sin( zenithAngle );

    // Perform transformation.
    cartesianCoordinates( 0 ) = radius * cosineOfAzimuthAngle * sineOfZenithAngle;
    cartesianCoordinates( 1 ) = radius * std::sin( azimuthAngle ) * sineOfZenithAngle;
    cartesianCoordinates( 2 ) = radius * std::cos( zenithAngle );
}

//! Convert cartesian to spherical coordinates.
void convertCartesianToSpherical( const Eigen::VectorXd& cartesianCoordinates,
                                  Eigen::VectorXd& sphericalCoordinates )
{
    // Compute transformation of Cartesian coordinates to spherical coordinates.
    sphericalCoordinates( 0 ) = cartesianCoordinates.norm( );

    // Check if coordinates are at origin.
    if ( sphericalCoordinates( 0 ) < std::numeric_limits< double >::epsilon( ) )
    {
        sphericalCoordinates( 1 ) = 0.0;
        sphericalCoordinates( 2 ) = 0.0;
    }

    // Else compute coordinates using trigonometric relationships.
    else
    {
        sphericalCoordinates( 1 ) = std::atan2( cartesianCoordinates( 1 ),
                                                cartesianCoordinates( 0 ) );
        sphericalCoordinates( 2 ) = std::acos( cartesianCoordinates( 2 )
                                               / sphericalCoordinates( 0 ) );
    }
}

//! Convert cylindrical to Cartesian coordinates.
Eigen::Vector3d convertCylindricalToCartesianCoordinates( double radius,
                                                          double azimuthAngle, double z )
{
    // Create Cartesian coordinates vector.
    Eigen::Vector3d cartesianCoordinates;

    // If radius < 0, then give warning.
    if ( radius < 0.0 )
    {
        std::cerr << "Warning: cylindrical radial coordinate is negative!\n"
                  << "This could give incorrect results!\n";
    }

    // Compute and set Cartesian coordinates.
    cartesianCoordinates << radius * std::cos( azimuthAngle ),   // x-coordinate
                            radius * std::sin( azimuthAngle ),   // y-coordinate
                            z;                                   // z-coordinate

    return cartesianCoordinates;
}

//! Convert cylindrical to Cartesian coordinates.
Eigen::Vector3d convertCylindricalToCartesianCoordinates( Eigen::Vector3d cylindricalCoordinates )
{
    // Create Cartesian coordinates vector.
    Eigen::Vector3d cartesianCoordinates;

    // If radius < 0, then give warning.
    if ( cylindricalCoordinates( 0 ) < 0.0 )
    {
        std::cerr << "Warning: cylindrical radial coordinate is negative! "
                  << "This could give incorrect results!\n";
    }

    // Compute and set Cartesian coordinates.
    cartesianCoordinates
            << cylindricalCoordinates( 0 )
               * std::cos( cylindricalCoordinates( 1 ) ), // x-coordinate
            cylindricalCoordinates( 0 )
            * std::sin( cylindricalCoordinates( 1 ) ),    // y-coordinate
            cylindricalCoordinates( 2 );                  // z-coordinate

    return cartesianCoordinates;
}

//! Convert cylindrical to Cartesian state.
Eigen::VectorXd convertCylindricalToCartesianState( Eigen::VectorXd cylindricalState )
{
    // Create Cartesian state vector, initialized with zero entries.
    Eigen::VectorXd cartesianState = Eigen::VectorXd::Zero( 6 );

    // Get azimuth angle, theta.
    double azimuthAngle = cylindricalState( 1 );

    // Compute and set Cartesian coordinates.
    cartesianState.head( 3 ) = convertCylindricalToCartesianCoordinates(
                cylindricalState.head( 3 ) );

    // If r = 0 AND Vtheta > 0, then give warning and assume Vtheta=0.
    if ( std::fabs(cylindricalState( 0 )) <= std::numeric_limits< double >::epsilon( )
         && std::fabs(cylindricalState( 4 )) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Warning: cylindrical velocity Vtheta (r*thetadot) does not equal zero while "
                  << "the radius (r) is zero! Vtheta is taken equal to zero!\n";

        // Compute and set Cartesian velocities.
        cartesianState.tail( 3 )
                << cylindricalState( 3 ) * std::cos( azimuthAngle ),   // xdot
                   cylindricalState( 3 ) * std::sin( azimuthAngle ),   // ydot
                   cylindricalState( 5 );                              // zdot
    }

    else
    {
        // Compute and set Cartesian velocities.
        cartesianState.tail( 3 )
                << cylindricalState( 3 ) * std::cos( azimuthAngle )
                   - cylindricalState( 4 ) * std::sin( azimuthAngle ),   // xdot
                   cylindricalState( 3 ) * std::sin( azimuthAngle )
                   + cylindricalState( 4 ) * std::cos( azimuthAngle ),   // ydot
                   cylindricalState( 5 );                                // zdot
    }

    return cartesianState;
}

//! Convert Cartesian to cylindrical coordinates.
Eigen::Vector3d convertCartesianToCylindricalCoordinates( Eigen::Vector3d cartesianCoordinates )
{
    // Create cylindrical coordinates vector.
    Eigen::Vector3d cylindricalCoordinates;

    // Declare new variable, the azimuth angle.
    double azimuthAngle;

    // Compute azimuth angle, theta.
    /* If x = 0, then azimuthAngle = pi/2 (y>0) or 3*pi/2 (y<0) or 0 (y=0),
       else azimuthAngle = arctan(y/x).
    */
    if ( std::fabs(cartesianCoordinates( 0 ) ) <= std::numeric_limits< double >::epsilon( ) )
    {
        azimuthAngle = computeModulo( static_cast< double >(
                                      boost::math::sign( cartesianCoordinates( 1 ) ) )
                                      * 0.5 * M_PI, 2.0 * M_PI );
    }

    else
    {
        azimuthAngle = computeModulo( std::atan2( cartesianCoordinates( 1 ),
                                                  cartesianCoordinates( 0 ) ), 2.0 * M_PI );
    }

    // Compute and set cylindrical coordinates.
    cylindricalCoordinates <<
        std::sqrt( pow( cartesianCoordinates( 0 ), 2 )
                   + pow( cartesianCoordinates( 1 ), 2 ) ), // Radius
        azimuthAngle,                                       // Azimuth angle, theta
        cartesianCoordinates( 2 );                          // z-coordinate

    return cylindricalCoordinates;
}

//! Convert Cartesian to cylindrical state.
Eigen::VectorXd convertCartesianToCylindricalState( Eigen::VectorXd cartesianState )
{
    // Create cylindrical state vector, initialized with zero entries.
    Eigen::VectorXd cylindricalState = Eigen::VectorXd::Zero( 6 );

    // Compute and set cylindrical coordinates.
    cylindricalState.head( 3 ) = convertCartesianToCylindricalCoordinates(
                cartesianState.head( 3 ) );

    // Compute and set cylindrical velocities.
    /* If radius = 0, then Vr = sqrt(xdot^2+ydot^2) and Vtheta = 0,
       else Vr = (x*xdot+y*ydot)/radius and Vtheta = (x*ydot-y*xdot)/radius.
    */
    if ( cylindricalState( 0 ) <= std::numeric_limits< double >::epsilon( ) )
    {
        cylindricalState.tail( 3 ) <<
            std::sqrt( pow( cartesianState( 3 ), 2 ) + pow( cartesianState( 4 ), 2 ) ), // Vr
            0.0,                                                                        // Vtheta
            cartesianState( 5 );                                                        // Vz
    }
    else
    {
        cylindricalState.tail( 3 ) <<
            ( cartesianState( 0 ) * cartesianState( 3 )
              + cartesianState( 1 ) * cartesianState( 4 ) ) / cylindricalState( 0 ),    // Vr
            ( cartesianState( 0 ) * cartesianState( 4 )
              - cartesianState( 1 ) * cartesianState( 3 ) ) / cylindricalState( 0 ),    // Vtheta
                cartesianState( 5 );                                                    // Vz
    }

    return cylindricalState;
}

//! Compute modulo of double.
double computeModulo( double dividend, double divisor )
{ return dividend - divisor * floor( dividend / divisor ); }

//! Compute sample mean.
double computeSampleMean( const vector< double >& sampleData )
{
    // Return sample mean.
    return accumulate( sampleData.begin( ), sampleData.end( ), 0.0 )
            / static_cast< double >( sampleData.size( ) );
}

//! Compute sample variance.
double computeSampleVariance( const vector< double >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    double sampleMean_ = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    double sumOfResidualsSquared_ = 0.0;

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared_ += pow( sampleData.at( i ) - sampleMean_, 2.0 );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared_;
}

}

}

// End of file.
