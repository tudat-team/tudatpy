/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *     Stackoverflow. C++ GCC4.4 warning: array subscript is above array bounds, 2009,
 *         http://stackoverflow.com/questions/
 *             1168525/c-gcc4-4-warning-array-subscript-is-above-array-bounds,
 *         last accessed: 27th December, 2013.
 *     GCC Mailing List. Re: How to fix 'array subscript is above array bounds' ?, 2012,
 *         http://gcc.gnu.org/ml/gcc-help/2012-04/msg00047.html, last accessed: 27th December,
 *         2013.
 *
 *    Notes
 *     Under older GCC-based compilers (4.3 and 4.4 series), it is known that this file will
 *     generate a spurious warning stating "array subscript is above array bounds". This warning
 *     can be safely ignored. It is recommended that you working with a GCC 4.5+ compiler, since
 *     the problem has been fixed in all versions that postdate 4.5. This warning has specifically
 *     been noted when compiling using the MinGW GCC 4.4.0 compiler under MS Windows. For more
 *     information on the nature of this warning, please take a look at Stackoverflow (2009) and
 *     GCC Mailing List (2012).
 *
 */

#ifndef TUDAT_MULTI_LINEAR_INTERPOLATOR_H
#define TUDAT_MULTI_LINEAR_INTERPOLATOR_H

#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

namespace tudat
{
namespace interpolators
{

//! Class for performing multi-linear interpolation for arbitrary number of independent variables.
/*!
 * Class for performing multi-linear interpolation for arbitrary number of independent variables.
 * Interpolation is calculated recursively over all dimensions of independent variables. Note
 * that the types (i.e. double, float) of all independent variables must be the same.
 * \tparam IndependentVariableType Type for independent variables.
 * \tparam DependentVariableType Type for dependent variable.
 * \tparam NumberOfDimensions Number of independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType,
          int NumberOfDimensions >
class MultiLinearInterpolator: public Interpolator< IndependentVariableType,
        DependentVariableType >
{
public:

    //! Constructor taking independent and dependent variable data.
    /*!
     * \param independentValues Vector of vectors containing data points of independent variables,
     *  each must be sorted in ascending order.
     * \param dependentData Multi-dimensional array of dependent data at each point of
     *          hyper-rectangular grid formed by independent variable points.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *          to find the nearest lower data point in the independent variables when requesting
     *          interpolation.
     */
    MultiLinearInterpolator( const std::vector< std::vector< IndependentVariableType > >
                             independentValues,
                             const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions )>
                             dependentData,
                             const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
        : independentValues_( independentValues ),
          dependentData_( dependentData )
    {
        // Check consistency of template arguments and input variables.
        if ( independentValues.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error( "Error: dimension of independent value vector provided to constructor incompatible with template parameter " );
        }

        // Check consistency of input data of dependent and independent data.
        for ( int i = 0; i < NumberOfDimensions; i++ )
        {
            if ( independentValues[ i ].size( ) != dependentData.shape( )[ i ] )
            {
                std::string errorMessage = "Warning: number of data points in dimension" +
                        std::to_string( i ) + "of independent and dependent data incompatible";
                throw std::runtime_error( errorMessage );
            }
        }

        makeLookupSchemes( selectedLookupScheme );
    }

    //! Default destructor
    /*!
     *  Default destructor
     */
    ~MultiLinearInterpolator( ){ }

    //! Function to perform interpolation.
    /*!
     *  This function performs the multilinear interpolation.
     *  \param independentValuesToInterpolate Vector of values of independent variables at which
     *  the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable in all dimensions.
     */
    DependentVariableType interpolate(
            const std::vector< IndependentVariableType >& independentValuesToInterpolate )
    {
        // Determine the nearest lower neighbours.
        std::vector< int > nearestLowerIndices;
        nearestLowerIndices.resize( NumberOfDimensions );
        for ( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            nearestLowerIndices[ i ] = lookUpSchemes_[ i ]->findNearestLowerNeighbour(
                    independentValuesToInterpolate[ i ] );
        }

        // Initialize function evaluation indices to -1 for debugging purposes.
        boost::array< int, NumberOfDimensions > interpolationIndices;
        for ( int i = 0; i < NumberOfDimensions; i++ )
        {
            interpolationIndices[ i ] = -1;
        }

        // Call first step of interpolation, this function calls itself at subsequent independent
        // variable dimensions to evaluate and properly scale dependent variable table values at
        // all 2^n grid edges.
        return performRecursiveInterpolationStep( 0, independentValuesToInterpolate,
                                                  interpolationIndices, nearestLowerIndices );
    }

    //! Function to return the number of independent variables of the interpolation.
    /*!
     *  Function to return the number of independent variables of the interpolation, i.e. size
     *  that the vector used as input for Interpolator::interpolate should be.
     *  \return Number of independent variables of the interpolation.
     */
    int getNumberOfDimensions( )
    {
        return NumberOfDimensions;
    }


private:

    //! Make the lookup scheme that is to be used.
    /*!
     * This function creates the look up scheme that is to be used in determining the interval of
     * the independent variable grid where the interpolation is to be performed. It takes the type
     * of lookup scheme as an enum and constructs the lookup scheme from the independentValues_
     * that have been set previously.
     *  \param selectedScheme Type of look-up scheme that is to be used
     */
    void makeLookupSchemes( const AvailableLookupScheme selectedScheme )
    {
        lookUpSchemes_.resize( NumberOfDimensions );
        // Find which type of scheme is used.
        switch( selectedScheme )
        {
        case binarySearch:

            for( int i = 0; i < NumberOfDimensions; i++ )
            {
                // Create binary search look up scheme.
                lookUpSchemes_[ i ] = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                        ( new BinarySearchLookupScheme< IndependentVariableType >(
                              independentValues_[ i ] ) );
            }

            break;

        case huntingAlgorithm:

            for( int i = 0; i < NumberOfDimensions; i++ )
            {
                // Create hunting scheme, which uses an intial guess from previous look-ups.
                lookUpSchemes_[ i ] = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                        ( new HuntingAlgorithmLookupScheme< IndependentVariableType >(
                              independentValues_[ i ] ) );
            }

            break;

        default:

            throw std::runtime_error( "Warning: lookup scheme not found when making scheme for 1-D interpolator" );
        }
    }

    //! Perform the step in a single dimension of the interpolation process.
    /*!
     * Function calculates single dimension of the interpolation process. Function calls itself if
     * final dimension not yet reached. Calling this function with currentVariable = 0 will result
     * in 2^{NumberOfDimensions} number of calls to the function at currentVariable =
     * NumberOfDimensions -1. As such, the complete series of calls, starting at currentVariable =
     * 0, retrieves the dependent variable values at all edges of the grid hyper-rectangle and
     * properly scales them.
     * \param currentVariable Dimension in which this interpolation step is to be performed.
     * \param independentValuesToInterpolate Vector of values of independent variables at which
     *          interpolation is to be performed.
     * \param currentArrayIndices Array of indices modified at index = currentVariable at each
     *          call of function. Variable is passed to dependentData in highest step to return
     *          data for interpolation.
     * \param nearestLowerIndices Indices in subvectors of independentValues_ vector. That is, the
     *  n-th entry of nearestLowerIndices represent the nearest lower neighbour in the n-th
     *  interpolation dimension of the independent variable vectors.
     * \return Interpolated value in a single dimension
     */
    DependentVariableType performRecursiveInterpolationStep(
            const unsigned int currentVariable,
            const std::vector< IndependentVariableType >& independentValuesToInterpolate,
            boost::array< int, NumberOfDimensions > currentArrayIndices,
            const std::vector< int >& nearestLowerIndices )
    {
        IndependentVariableType upperFraction, lowerFraction;
        DependentVariableType upperContribution, lowerContribution;

        // Calculate fractions of data points above and below independent
        // variable value to be added to interpolated value.
        upperFraction = ( independentValuesToInterpolate[ currentVariable ] -
                          independentValues_[ currentVariable ]
                          [ nearestLowerIndices[ currentVariable ] ] ) /
                ( independentValues_[ currentVariable ]
                  [ nearestLowerIndices[ currentVariable ] + 1 ] -
                  independentValues_[ currentVariable ]
                  [ nearestLowerIndices[ currentVariable ] ] );
        lowerFraction = -( independentValuesToInterpolate[ currentVariable ] -
                           independentValues_[ currentVariable ]
                           [ nearestLowerIndices[ currentVariable ] + 1 ] ) /
                ( independentValues_[ currentVariable ]
                  [ nearestLowerIndices[ currentVariable ] + 1 ] -
                  independentValues_[ currentVariable ]
                  [ nearestLowerIndices[ currentVariable ] ] );

        // If at top dimension, call dependent variable data.
        if ( currentVariable == NumberOfDimensions - 1 )
        {
            currentArrayIndices[ NumberOfDimensions - 1 ] = nearestLowerIndices[ currentVariable ];
            lowerContribution = dependentData_( currentArrayIndices );
            currentArrayIndices[ NumberOfDimensions - 1 ] = nearestLowerIndices[ currentVariable ]
                                                            + 1;
            upperContribution = dependentData_( currentArrayIndices );
        }

        // If at lower dimension, update currentArrayIndices and call function with
        // currentVariable++.
        else
        {
            currentArrayIndices[ currentVariable ] = nearestLowerIndices[ currentVariable ];
            lowerContribution = performRecursiveInterpolationStep(
                        currentVariable + 1, independentValuesToInterpolate,
                        currentArrayIndices, nearestLowerIndices );
            currentArrayIndices[ currentVariable ] = nearestLowerIndices[ currentVariable ] + 1;
            upperContribution = performRecursiveInterpolationStep(
                        currentVariable + 1, independentValuesToInterpolate,
                        currentArrayIndices, nearestLowerIndices );
        }

        // Return interpolated value.
        DependentVariableType returnValue = upperFraction * upperContribution +
                                            lowerFraction * lowerContribution;
        return returnValue;
    }

    //! Vector with pointers to look-up scheme.
    /*!
     * Pointers to the look-up schemes that is used to determine in which interval the requested
     * independent variable value falls.
     */
    std::vector< boost::shared_ptr< LookUpScheme< IndependentVariableType > > > lookUpSchemes_;

    //! Vector of vectors containing independent variables.
    /*!
     * Vector of vectors containing independent variables. The size of the outer vector is equal
     * to the number of dimensions of the interpolator.
     */
    std::vector< std::vector< IndependentVariableType > > independentValues_;

    //! Multi-dimensional array of dependent data.
    /*!
     * Multi-dimensional array of dependent data at each point of hyper-rectangular grid formed by
     * independent variable points.
     */
    boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions )> dependentData_;
};

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_MULTI_LINEAR_INTERPOLATOR_H
