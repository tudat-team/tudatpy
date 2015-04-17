/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120716    D. Dirkx          File created.
 *      130121    K. Kumar          Added shared-ptr typedefs.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_LOOK_UP_SCHEME_H
#define TUDAT_LOOK_UP_SCHEME_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

namespace tudat
{
namespace interpolators
{

//! Enum of available lookup schemes.
/*!
 *  Enum of available lookup schemes.
 */
enum AvailableLookupScheme
{
    huntingAlgorithm,
    binarySearch
};

//! Look-up scheme class for nearest left neighbour search.
/*!
 * Look-up scheme class for nearest left neighbour search,
 * allows for different types of look-up scheme with a single interface.
 * \tparam IndependentVariableType Type of entries of vector in which lookup is to be performed.
 */
template< typename IndependentVariableType >
class LookUpScheme
{
public:

    //! Constructor, used to set data vector.
    /*!
     * Constructor, used to set data vector.
     * \param independentVariableValues vector of independent variable values in which to perform
     * lookup procedure.
     */
    LookUpScheme( const std::vector< IndependentVariableType >& independentVariableValues )
        : independentVariableValues_( independentVariableValues )
    { }

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~LookUpScheme( ) { }

    //! Find nearest left neighbour.
    /*!
     * Function finds nearest left neighbour of given value in independentVariableValues_.
     * \param valueToLookup Value of which nearest neaighbour is to be determined.
     * \return Index of entry in independentVariableValues_ vector which is nearest lower neighbour
     * to valueToLookup.
     */
    virtual int findNearestLowerNeighbour( const IndependentVariableType valueToLookup ) = 0;

protected:

    //! Vector of independent variable values in which lookup is to be performed.
    /*!
     * Vector of independent variable values in which lookup is to be performed.
     */
    std::vector< IndependentVariableType > independentVariableValues_;
};

//! Look-up scheme class for nearest left neighbour search using hunting algorithm.
/*!
 *  Look-up scheme class for nearest left neighbour search using hunting algorithm.
 *  \tparam IndependentVariableType Type of entries of vector in which lookup is to be performed.
 */
template< typename IndependentVariableType >
class HuntingAlgorithmLookupScheme: public LookUpScheme< IndependentVariableType >
{
public:

    using LookUpScheme< IndependentVariableType >::independentVariableValues_;

    //! Constructor, used to set data vector.
    /*!
     *  Constructor, used to set data vector. Initializes guess from 'previous' request to 0.
     * \param independentVariableValues vector of independent variable values in which to perform
     * lookup procedure.
     */
    HuntingAlgorithmLookupScheme( const std::vector< IndependentVariableType >&
                                  independentVariableValues )
        : LookUpScheme< IndependentVariableType >( independentVariableValues ),
          isFirstLookupDone( 0 ),
          previousNearestLowerIndex_( 0 )
    { }

    //! Default destructor
    /*!
     *  Default destructor
     */
    ~HuntingAlgorithmLookupScheme( ){ }

    //! Find nearest left neighbour.
    /*!
     * Function finds nearest left neighbour of given value in ndependentVariableValues_. If this
     * is first call of function, a binary search is used.
     * \param valueToLookup Value of which nearest neighbour is to be determined.
     * \return Index of entry in independentVariableValues_ vector which is nearest lower neighbour
     * to valueToLookup.
     */
    int findNearestLowerNeighbour( const IndependentVariableType valueToLookup )
    {
        // Initialize return value.
        int newNearestLowerIndex = 0;

        // If this is first call of function, use binary search.
        if ( !isFirstLookupDone )
        {
            newNearestLowerIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch
                    < IndependentVariableType >( independentVariableValues_, valueToLookup );
            isFirstLookupDone = 1;
        }

        else
        {
            // If requested value is in same interval, return same value as previous time.
            if ( basic_mathematics::isIndependentVariableInInterval< IndependentVariableType >
                 ( previousNearestLowerIndex_,  valueToLookup, independentVariableValues_ ) )
            {
                newNearestLowerIndex = previousNearestLowerIndex_;
            }

            // Otherwise, perform hunting algorithm.
            else
            {
                newNearestLowerIndex =
                        basic_mathematics::findNearestLeftNeighbourUsingHuntingAlgorithm<
                        IndependentVariableType >
                        (  valueToLookup, previousNearestLowerIndex_, independentVariableValues_ );
            }
        }

        // Set calculated value for use in next call.
        previousNearestLowerIndex_ = newNearestLowerIndex;

        return newNearestLowerIndex;
    }

private:

    //! Boolean to denote whether a lookup has been done.
    /*!
     * Boolean to denote whether a lookup has been done.
     */
    bool isFirstLookupDone;

    //! Nearest left index during previous call.
    /*!
     * Nearest left index during previous call
     */
    int previousNearestLowerIndex_;
};

//! Look-up scheme class for nearest left neighbour search using binary search algorithm.
/*!
 * Look-up scheme class for nearest left neighbour search using binary search algorithm.
 * \tparam IndependentVariableType Type of entries of vector in which lookup is to be performed.
 */
template< typename IndependentVariableType >
class BinarySearchLookupScheme: public LookUpScheme< IndependentVariableType >
{
public:

    using LookUpScheme< IndependentVariableType >::independentVariableValues_;

    //! Constructor, used to set data vector.
    /*!
     * Constructor, used to set data vector.
     * \param independentVariableValues vector of independent variable values in which to perform
     * lookup procedure.
     */
    BinarySearchLookupScheme(
            const std::vector< IndependentVariableType >& independentVariableValues )
        : LookUpScheme< IndependentVariableType >( independentVariableValues )
    { }

    //! Default destructor
    /*!
     *  Default destructor
     */
    ~BinarySearchLookupScheme( ){ }

    //! Find nearest left neighbour.
    /*!
     * Function finds nearest left neighbour of given value in independentVariableValues_.
     * \param valueToLookup Value of which nearest neaighbour is to be determined.
     * \return Index of entry in independentVariableValues_ vector which is nearest lower neighbour
     * to valueToLookup.
     */
    int findNearestLowerNeighbour( const IndependentVariableType valueToLookup )
    {
        return basic_mathematics::computeNearestLeftNeighborUsingBinarySearch
                < IndependentVariableType >( independentVariableValues_, valueToLookup );
    }
};

//! Typedef for shared-pointer to LookUpScheme object with double-type entries.
typedef boost::shared_ptr< LookUpScheme< double > > LookUpSchemeDoublePointer;

//! Typedef for shared-pointer to HuntingAlgorithmLookupScheme object with double-type entries.
typedef boost::shared_ptr< HuntingAlgorithmLookupScheme< double > >
HuntingAlgorithmLookupSchemeDoublePointer;

//! Typedef for shared-pointer to BinarySearchLookupScheme object with double-type entries.
typedef boost::shared_ptr< BinarySearchLookupScheme< double > >
BinarySearchLookupSchemeDoublePointer;

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_LOOK_UP_SCHEME_H
