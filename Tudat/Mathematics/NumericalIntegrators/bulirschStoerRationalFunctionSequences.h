/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#ifndef TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H
#define TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H

#include <vector>

namespace tudat
{

namespace numerical_integrators
{

//! Struct that defines the rational function sequence for a Bulirsch-Stoer integrator.
/*!
 * Struct that defines the rational function sequence for a Bulirsch-Stoer integrator.
 */
struct BulirschStoerRationalFunctionSequences
{

    //! Default constructor.
    /*!
     * Default constructor that initializes without setting the rational function sequence.
     */
    BulirschStoerRationalFunctionSequences( ) { }

    //! Constructor.
    /*!
     * Constructor that sets the rational function sequence.
     * \param sequence Rational function sequence.
     */
    BulirschStoerRationalFunctionSequences( const std::vector< unsigned int > sequence ) :
        rationalFunctionSequence( sequence ) { }

    //! Rational function sequence.
    /*!
     * Rational function sequence.
     */
    std::vector< unsigned int > rationalFunctionSequence;

    //! Enum of predefined rational function sequences.
    enum RationalFunctionSequences
    {
        bulirschStoer,
        deufelhard
    };

    //! Get rational function sequence.
    /*!
     * Returns requested rational function sequence.
     * \param rationalFunctionSequence The rational function sequence.
     * \return The requested rational function sequence.
     */
    static const BulirschStoerRationalFunctionSequences& get(
            RationalFunctionSequences sequence, const unsigned int lengthOfSequence = 100 );
};

} // namespace integrators
} // namespace tudat

#endif // TUDAT_BULIRSCH_STOER_RATIONAL_FUNCTION_SEQUENCES_H
