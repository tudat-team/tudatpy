/*! \file exponentialRandomNumberGenerator.h
 *    This header file contains a class definition for generating random
 *    numbers with exponential distribution.
 *
 *    Path              : /Mathematics/RandomNumberGenerators/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 8 July, 2011
 *    Last modified     : 8 July, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *      See notes for uniform random number generator.
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
 *      110708    K. Kumar          First creation of code.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef EXPONENTIALRANDOMNUMBERGENERATOR_H
#define EXPONENTIALRANDOMNUMBERGENERATOR_H

// Include statements.
#include "Mathematics/RandomNumberGenerators/uniformRandomNumberGenerator.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Exponential random number generator class.
/*!
 * Definition of class that generates random numbers with an exponential
 * distribution.
 */
class ExponentialRandomNumberGenerator : public UniformRandomNumberGenerator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialRandomNumberGenerator( const double& exponentialRandomNumberParameter,
                                      unsigned long long seed )
        : UniformRandomNumberGenerator( seed ),
          exponentialRandomNumberParameter_( exponentialRandomNumberParameter ) { }

    //! Get exponentially-distributed, normalized, random double.
    /*!
     * Returns an exponentially-distributed, normalized, random double within
     * the interval [ 0, infinity ).
     * \return Exponentially-distributed, normalized, random double.
     */
    double getExponentiallyDistributedNormalizedRandomDouble( );

protected:

private:

    //! Exponential random number parameter.
    /*!
     * Parameter that describes random number exponential distribution.
     * The random parameter, \f$ \beta \f$, scales the exponential
     * distribution as follows:
     * \f[
     *      p( \beta * y ) = \beta * e^{ -\beta * y }
     * \f]
     */
    double exponentialRandomNumberParameter_;
};

}

#endif // EXPONENTIALRANDOMNUMBERGENERATOR_H

// End of file.
