/*! \file uniformRandomNumberGenerator.cpp
 *    This source file contains a class implemetation for generating random
 *    numbers with uniform distribution.
 *
 *    Path              : /Mathematics/RandomNumberGenerators/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 15 October, 2010
 *    Last modified     : 16 May, 2011
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
 *      101015    K. Kumar          First creation of code.
 *      101020    K. Kumar          Class implemented as described by the
 *                                  header file. Code comments included.
 *      101025    K. Kumar          Updated code based on D. Dirkx's codecheck
 *                                  comments.
 *      110107    K. Kumar          Changed normalizated function to use
 *                                  climits.
 *      110117    K. Kumar          Minor layout modification; path corrected.
 *      110516    K. Kumar          Renamed file and class.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/RandomNumberGenerators/uniformRandomNumberGenerator.h"
#include "Mathematics/basicMathematicsFunctions.h"

//! Customized constructor to pass random number generator seed.
UniformRandomNumberGenerator::UniformRandomNumberGenerator( unsigned long long seed )
    : randomNumberParameter2_( 4101842887655102017LL ),
      randomNumberParameter3_( 1 )
{
    randomNumberParameter1_ = seed ^ randomNumberParameter2_;
    getUniformlyDistributedRandomInteger( );

    randomNumberParameter2_ = randomNumberParameter1_;
    getUniformlyDistributedRandomInteger( );

    randomNumberParameter3_ = randomNumberParameter2_;
    getUniformlyDistributedRandomInteger( );
}

//! Get uniformly distributed random integer.
unsigned long long UniformRandomNumberGenerator
::getUniformlyDistributedRandomInteger( )
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
double UniformRandomNumberGenerator
::getUniformlyDistributedNormalizedRandomDouble( )
{
    // Return uniformly distributed, normalized, random double.
    return  static_cast< double >( getUniformlyDistributedRandomInteger( ) )
            / ULLONG_MAX;
}

//! Get uniformly distributed random integer using 32-bit arithmetic.
unsigned int UniformRandomNumberGenerator::
        getUniformlyDistributedRandom32BitInteger( )
{
    // Return uniformly distributed random integer using 32-bit arithmetic
    return static_cast< unsigned int >(
                getUniformlyDistributedRandomInteger( ) );
}

//! Get uniformly distributed random plus/minus sign.
int UniformRandomNumberGenerator::getUniformlyDistributedRandomPlusMinusSign( )
{
    // Declare local variables.
    int randomPlusMinusSign_;

    // Get random integer.
    randomPlusMinusSign_ = static_cast< int >(
                getUniformlyDistributedRandomInteger( ) );

    // Normalize to plus/minus 1.
    randomPlusMinusSign_ /= mathematics::
                            computeAbsoluteValue( randomPlusMinusSign_ );

    // Return random plus/minus sign.
    return randomPlusMinusSign_;
}

// End of file.
