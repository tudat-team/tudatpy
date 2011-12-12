/*! \file linearAlgebra.cpp
 *    This file contains include statements for common Eigen components,
 *    typedef for vectors and a number of useful vector operation definitions.
 *
 *    Path              : /Mathematics/LinearAlgebra/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Date created      : 7 August, 2009
 *    Last modified     : 5 September, 2011
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
 *      090807    J. Melman         First creation of code.
 *      100930    D. Dirkx          Modified to comply with Tudat standards
 *      100930    J. Melman         Implemented namespace, minor comment
 *                                  changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

//! Tudat library namespace.
namespace tudat
{

//! Linear algebra namespace.
namespace linear_algebra
{

//! Determine the cosine of the angle between two vectors.
double determineCosineOfAngleBetweenVectors( const Eigen::Vector3d& vector0,
                                             const Eigen::Vector3d& vector1 )
{
    // Determine the length of the vectors.
    double normOfVector0 = vector0.norm( );
    double normOfVector1 = vector1.norm( );

    // Normalize both vectors.
    Eigen::Vector3d vector0Normalized = vector0 / normOfVector0;
    Eigen::Vector3d vector1Normalized = vector1 / normOfVector1;

    // Get the cosine of the angle by dotting the normalized vectors.
    double dotProductOfNormalizedVectors = vector0Normalized.dot( vector1Normalized );

    // Explicitly define the extreme cases, which can give problems with the acos function.
    if ( dotProductOfNormalizedVectors >= 1.0 )
    {
        return 1.0;
    }

    else if ( dotProductOfNormalizedVectors <= -1.0 )
    {
        return -1.0;
    }

    // Determine the actual angle.
    else
    {
        return dotProductOfNormalizedVectors;
    }
}

//! Determine the angle between two vectors.
double determineAngleBetweenVectors( const Eigen::Vector3d& vector0,
                                     const Eigen::Vector3d& vector1 )
{
    // Determine the cosine of the angle by using another routine.
    double dotProductOfNormalizedVectors = determineCosineOfAngleBetweenVectors( vector0,
                                                                                 vector1 );

    // Return arccosine of the above, which is effectively the angle.
    return std::acos( dotProductOfNormalizedVectors );
}

//! Determine the average of the components of a vector.
double determineAverageOfVectorComponents( const Eigen::VectorXd& vector0 )
{
    return vector0.sum( ) / vector0.rows( );
}

//! Determine the standard deviation of the components of a vector.
double determineStandardDeviationOfVectorComponents( const Eigen::VectorXd& vector0 )
{
    double varianceOfEntries = 0.0;

    // Determine average of entries.
    double averageOfEntries = determineAverageOfVectorComponents( vector0 );

    // Determine variance of entries.
    for ( int i = 0 ; i < vector0.rows( ) ; i++ )
    {
        varianceOfEntries += std::pow( ( vector0( i ) - averageOfEntries ), 2.0 );
    }

    varianceOfEntries /= vector0.rows( ) - 1;

    // Return square root of variance ( = standard deviation ).
    return std::sqrt( varianceOfEntries );
}

}

}

// End of file.
