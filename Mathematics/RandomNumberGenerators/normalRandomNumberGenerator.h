/*! \file normalRandomNumberGenerator.h
 *    This header file contains a class definition for generating random
 *    numbers with normal distribution.
 *
 *    Path              : /Mathematics/RandomNumberGenerators/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 26 July, 2011
 *    Last modified     : 10 August, 2011
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
 *      110726    K. Kumar          First creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation (equations).
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef NORMALRANDOMNUMBERGENERATOR_H
#define NORMALRANDOMNUMBERGENERATOR_H

// Include statements.
#include <cmath>
#include "Mathematics/RandomNumberGenerators/uniformRandomNumberGenerator.h"

//! Normal random number generator class.
/*!
 * Definition of class that generates random numbers with a normal
 * distribution based on the Box-Muller transformation (Press et al., 2002).
 */
class NormalRandomNumberGenerator : public UniformRandomNumberGenerator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NormalRandomNumberGenerator( const double& mean,
                                 const double& standardDeviation,
                                 unsigned long long seed );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~NormalRandomNumberGenerator( ) { }

    //! Get normally-distributed, normalized, random double.
    /*!
     * Returns an normally-distributed, normalized, random double within
     * the interval ( -infinity, infinity ).
     * \return Normally-distributed, normalized, random double.
     */
    double getNormallyDistributedNormalizedRandomDouble( );

protected:

private:

    //! Mean.
    /*!
     * Mean of normal distribution to generate random numbers.
     * The mean, \f$ \mu \f$, scales the normal distribution as follows:
     * \f[
     *      p( y ) = \frac{ 1 }{ \sqrt{ 2 * \pi } } * \sigma \beta
     *               * e^{ -\frac{ 1 }{ 2 } ( \frac{ x - \mu }{ \sigma } )^ 2 }
     * \f]
     */
    double mean_;

    //! Standard deviation.
    /*!
     * Standard deviation of normal distribution to generate random numbers.
     * The standard deviation, \f$ \sigma \f$, scales the normal distribution
     * as follows:
     * \f[
     *      p( y ) = \frac{ 1 }{ \sqrt{ 2 * \pi } } * \sigma \beta
     *               * e^{ -\frac{ 1 }{ 2 } ( \frac{ x - \mu }{ \sigma } )^ 2 }
     * \f]
     */
    double standardDeviation_;

    //! Ordinate in unit circle.
    /*!
     * Ordinate selected in unit circle as part of Box-Muller algorithm.
     */
    double ordinateInUnitCircle_;

    //! Abscissa in unit circle.
    /*!
     * Abscissa selected in unit circle as part of Box-Muller algorithm.
     */
    double abscissaInUnitCircle_;

    //! Stored value of normal random number.
    /*!
     * Stored value of normal random number.
     */
    double storedValue_;
};

#endif // NORMALRANDOMNUMBERGENERATOR_H

// End of file.
