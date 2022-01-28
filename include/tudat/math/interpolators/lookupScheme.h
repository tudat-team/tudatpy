/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_LOOK_UP_SCHEME_H
#define TUDAT_LOOK_UP_SCHEME_H

#include <vector>

#include <memory>

#include "tudat/math/basic/nearestNeighbourSearch.h"

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
    undefinedScheme,
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

    IndependentVariableType getMinimumValue( )
    {
        return independentVariableValues_.at( 0 );
    }
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
                 ( previousNearestLowerIndex_, valueToLookup, independentVariableValues_ ) )
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
typedef std::shared_ptr< LookUpScheme< double > > LookUpSchemeDoublePointer;

//! Typedef for shared-pointer to HuntingAlgorithmLookupScheme object with double-type entries.
typedef std::shared_ptr< HuntingAlgorithmLookupScheme< double > >
HuntingAlgorithmLookupSchemeDoublePointer;

//! Typedef for shared-pointer to BinarySearchLookupScheme object with double-type entries.
typedef std::shared_ptr< BinarySearchLookupScheme< double > >
BinarySearchLookupSchemeDoublePointer;

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_LOOK_UP_SCHEME_H
