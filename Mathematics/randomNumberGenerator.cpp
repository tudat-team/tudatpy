/*! \file randomNumberGenerator.cpp
 *    This source file contains a class implemetation for generating random numbers.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 15 October, 2010
 *    Last modified     : 25 October, 2010
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *      The random number generator implemented here is well documented in
 *      (Press et al., 2002), where the details of the algorithm steps are
 *      explained. The distinction between 64-bit and 32-bit arithmetic is
 *      based on the computer system being used to run simulations. For more
 *      information, refer to (Press et al., 2002). The use of the 'long long'
 *      variable type is to support 64-bit arithmetic in C++.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD    author        comment
 *      101015    K. Kumar      First creation of code.
 *      101020    K. Kumar      Class implemented as described by the header
 *                              file. Code comments included.
 *      101025    K. Kumar      Updated code based on D. Dirkx's codecheck
 *                              comments.
 */

// Include statements.
#include <cmath>
#include "randomNumberGenerator.h"

//! Customized constructor to pass random number generator seed.
RandomNumberGenerator::RandomNumberGenerator( unsigned long long seed )
    : randomNumberParameter2_( 4101842887655102017LL ),
    randomNumberParameter3_( 1 )
{
    randomNumberParameter1_ = seed ^ randomNumberParameter2_;
    getUniformlyDistributedRandom64BitInteger( );

    randomNumberParameter2_ = randomNumberParameter1_;
    getUniformlyDistributedRandom64BitInteger( );

    randomNumberParameter3_ = randomNumberParameter2_;
    getUniformlyDistributedRandom64BitInteger( );
}

//! Default destructor.
RandomNumberGenerator::~RandomNumberGenerator()
{
}

//! Get uniformly distributed random integer.
unsigned long long RandomNumberGenerator::
        getUniformlyDistributedRandom64BitInteger( )
{
    randomNumberParameter1_ = randomNumberParameter1_ * 2862933555777941757LL
                             + 7046029254386353087LL;

    randomNumberParameter2_ ^= randomNumberParameter2_ >> 17;
    randomNumberParameter2_ ^= randomNumberParameter2_ << 31;
    randomNumberParameter2_ ^= randomNumberParameter2_ >> 8;

    randomNumberParameter3_ = 4294957665U
                             * ( randomNumberParameter3_ & 0xffffffff )
                             + ( randomNumberParameter3_ >> 32 );

    unsigned long long randomNumberParameter4_ = randomNumberParameter1_
                            ^ ( randomNumberParameter1_ << 21 );
    randomNumberParameter4_ ^= randomNumberParameter4_ >> 35;
    randomNumberParameter4_ ^= randomNumberParameter4_ << 4;

    // Return uniformly distributed random integer using 64-bit arithmetic.
    return ( randomNumberParameter4_ + randomNumberParameter2_ )
            ^ randomNumberParameter3_;
}

//! Get uniformly distributed, normalized, random double.
double RandomNumberGenerator::
        getUniformlyDistributedNormalizedRandomDouble( )
{
    // Return uniformly distributed, normalized, random double.
    return 5.42101086242752217E-20
            * getUniformlyDistributedRandom64BitInteger( );
}

//! Get uniformly distributed random integer using 32-bit arithmetic.
unsigned int RandomNumberGenerator::
        getUniformlyDistributedRandom32BitInteger( )
{
    // Return uniformly distributed random integer using 32-bit arithmetic.
    return ( unsigned int )getUniformlyDistributedRandom64BitInteger( );
}

//! Get random plus/minus sign.
int RandomNumberGenerator::getRandomPlusMinusSign( )
{
    // Declare local variables.
    int randomPlusMinusSign_;

    // Get random integer.
    randomPlusMinusSign_ = getUniformlyDistributedRandom64BitInteger( );

    // Normalize to plus/minus 1.
    randomPlusMinusSign_ /= std::fabs( randomPlusMinusSign_ );

    // Return random plus/minus sign.
    return randomPlusMinusSign_;
}

// End of file.
