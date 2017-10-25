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
 *      120316    K. Kumar          File created
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#include <iostream>

#include <boost/assign/std/vector.hpp>

#include "Tudat/Mathematics/NumericalIntegrators/bulirschStoerRationalFunctionSequences.h"

namespace tudat
{

namespace numerical_integrators
{

//! Get rational function sequence.
const BulirschStoerRationalFunctionSequences&
BulirschStoerRationalFunctionSequences::get( RationalFunctionSequences sequence,
                                             const unsigned int lengthOfSequence )
{
    static BulirschStoerRationalFunctionSequences bulirschStoerSequence_, deufelhardSequence_;

    switch ( sequence )
    {
    case bulirschStoer:
        using namespace boost::assign;
        bulirschStoerSequence_.rationalFunctionSequence += 2, 4, 6;
        for ( unsigned int i = 3; i < lengthOfSequence; i++ )
        {
            bulirschStoerSequence_.rationalFunctionSequence.push_back(
                        2 * bulirschStoerSequence_.rationalFunctionSequence.at( i - 2 ) );
        }
        return bulirschStoerSequence_;

    case deufelhard:
        for ( unsigned int i = 0; i < lengthOfSequence; i++ )
        {
            deufelhardSequence_.rationalFunctionSequence.push_back( 2 * ( i + 1 ) );
        }
        return deufelhardSequence_;

    default: // The default case will never occur because sequence is an enum
        throw BulirschStoerRationalFunctionSequences( );
    }
}

} // namespace integrators
} // namespace tudat
