/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <array>
#include <vector>

#include <boost/multi_array.hpp>
#include <memory>
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/interpolator.h"

#include "tudat/math/basic/nearestNeighbourSearch.h"
#include "tudat/math/interpolators/multiDimensionalInterpolator.h"

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
template< typename IndependentVariableType, typename DependentVariableType, unsigned int NumberOfDimensions >
class MultiLinearInterpolator : public MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >
{
public:
    // Using statements to prevent having to put 'this' everywhere in the code.
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::dependentData_;
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::independentValues_;
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::lookUpSchemes_;

    //! Default constructor taking independent and dependent variable data.
    /*!
     *  Default constructor taking independent and dependent variable data.
     *  \param independentValues Vector of vectors containing data points of independent variables,
     *      each must be sorted in ascending order.
     *  \param dependentData Multi-dimensional array of dependent data at each point of
     *      hyper-rectangular grid formed by independent variable points.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *      interpolation.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Vector of pairs of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     *  \param periods Vector of periods for periodic/angular interpolation, one per dimension. If periods[i] > 0,
     *      dimension i uses shortest-path interpolation. If periods[i] = 0 (default), standard interpolation is used.
     */
    MultiLinearInterpolator( const std::vector< std::vector< IndependentVariableType > >& independentValues,
                             const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) >& dependentData,
                             const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                             const std::vector< BoundaryInterpolationType >& boundaryHandling =
                                     std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary ),
                             const std::vector< std::pair< DependentVariableType, DependentVariableType > >& defaultExtrapolationValue =
                                     std::vector< std::pair< DependentVariableType, DependentVariableType > >(
                                             NumberOfDimensions,
                                             std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                                                             IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ),
                             const std::vector< IndependentVariableType >& periods =
                                     std::vector< IndependentVariableType >( NumberOfDimensions, IndependentVariableType( 0 ) ) ):
        MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >( boundaryHandling,
                                                                                                            defaultExtrapolationValue ),
        periods_( periods )
    {
        if( periods_.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error(
                    "Error when making multi-linear interpolator, periods vector size is incompatible "
                    "with template parameter." );
        }

        // Save (in)dependent variables
        independentValues_ = independentValues;
        dependentData_.resize(
                reinterpret_cast< std::array< size_t,
                                              boost::multi_array< DependentVariableType,
                                                                  static_cast< size_t >( NumberOfDimensions ) >::dimensionality > const& >(
                        *dependentData.shape( ) ) );  // resize dependent data container
        dependentData_ = dependentData;

        // Check consistency of template arguments and input variables.
        if( independentValues.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error(
                    "Error: dimension of independent value vector provided to constructor "
                    "incompatible with template parameter." );
        }

        // Check consistency of input data of dependent and independent data.
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            if( independentValues[ i ].size( ) != dependentData.shape( )[ i ] )
            {
                std::string errorMessage = "Error: number of data points in dimension " + std::to_string( i ) +
                        " of independent and dependent data incompatible.";
                throw std::runtime_error( errorMessage );
            }
        }

        // Create lookup scheme from independent variable data points.
        this->makeLookupSchemes( selectedLookupScheme );
    }

    //! Constructor taking independent and dependent variable data.
    /*!
     *  Constructor taking independent and dependent variable data. This constructor only requires one boundary
     *  handling method, and assumes it for each dimension.
     *  \param independentValues Vector of vectors containing data points of independent variables,
     *      each must be sorted in ascending order.
     *  \param dependentData Multi-dimensional array of dependent data at each point of
     *      hyper-rectangular grid formed by independent variable points.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *      interpolation.
     *  \param boundaryHandling Boundary handling method, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value
     *      or use_default_value_with_warning as methods for boundaryHandling.
     *  \param periods Vector of periods for periodic/angular interpolation.
     */
    MultiLinearInterpolator(
            const std::vector< std::vector< IndependentVariableType > >& independentValues,
            const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) >& dependentData,
            const AvailableLookupScheme selectedLookupScheme,
            const BoundaryInterpolationType boundaryHandling,
            const DependentVariableType& defaultExtrapolationValue = IdentityElement::getAdditionIdentity< DependentVariableType >( ),
            const std::vector< IndependentVariableType >& periods =
                    std::vector< IndependentVariableType >( NumberOfDimensions, IndependentVariableType( 0 ) ) ):
        MultiLinearInterpolator( independentValues,
                                 dependentData,
                                 selectedLookupScheme,
                                 std::vector< BoundaryInterpolationType >( NumberOfDimensions, boundaryHandling ),
                                 std::vector< std::pair< DependentVariableType, DependentVariableType > >(
                                         NumberOfDimensions,
                                         std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ),
                                 periods )
    {}

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~MultiLinearInterpolator( ) {}

    //! Function to get the periods vector for periodic interpolation.
    /*!
     *  Function to get the periods vector for periodic interpolation.
     *  \return Vector of period values (0 for non-periodic dimensions).
     */
    std::vector< IndependentVariableType > getPeriods( ) const
    {
        return periods_;
    }

    //! Function to check if a specific dimension is configured for periodic interpolation.
    /*!
     *  Function to check if a specific dimension is configured for periodic interpolation.
     *  \param dimensionIndex Index of the dimension to check.
     *  \return True if periods_[dimensionIndex] > 0, false otherwise.
     */
    bool isPeriodic( const unsigned int dimensionIndex ) const
    {
        if( dimensionIndex >= NumberOfDimensions )
        {
            throw std::runtime_error( "MultiLinearInterpolator::isPeriodic: dimensionIndex out of range" );
        }
        return periods_.at( dimensionIndex ) > IndependentVariableType( 0 );
    }

    //! Structure to cache interpolation state for batch operations
    struct InterpolationState {
        std::vector< unsigned int > nearestLowerIndices;
        std::vector< IndependentVariableType > localIndependentValues;
        bool isValid = false;

        // 2D optimization: pre-computed fractions and indices
        IndependentVariableType tx, ty, one_minus_tx, one_minus_ty;
        unsigned int i0, i1;
        unsigned int i0_next, i1_next;  // Next indices (may wrap to 0 for periodic)
    };

    //! Prepare interpolation state for batch operations (amortizes lookup cost)
    /*!
     * This function computes the grid indices and interpolation fractions once,
     * which can then be reused for multiple interpolations with different dependent data.
     * This is particularly efficient when interpolating many coefficient values at the same point.
     * \param independentValuesToInterpolate Vector of independent variable values
     * \return Interpolation state that can be passed to interpolateWithState()
     */
    InterpolationState prepareInterpolationState( const std::vector< IndependentVariableType >& independentValuesToInterpolate )
    {
        InterpolationState state;

        // Check whether size of independent variable vector is correct
        if( independentValuesToInterpolate.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error(
                    "Error in multi-dimensional interpolator. The number of independent variables "
                    "provided is incompatible with the previous definition. Provided: " +
                    std::to_string( independentValuesToInterpolate.size( ) ) + ". Needed: " + std::to_string( NumberOfDimensions ) );
        }

        // Create local copy
        state.localIndependentValues = independentValuesToInterpolate;

        // Check boundary cases
        bool useValue = false;
        DependentVariableType currentDependentVariable;
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            this->checkBoundaryCase( i, useValue, state.localIndependentValues.at( i ), currentDependentVariable );
            if( useValue )
            {
                // Cannot use state for boundary cases
                state.isValid = false;
                return state;
            }
        }

        // Determine the nearest lower neighbours
        state.nearestLowerIndices.resize( NumberOfDimensions );
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            state.nearestLowerIndices[ i ] = lookUpSchemes_[ i ]->findNearestLowerNeighbour( state.localIndependentValues[ i ] );

            // Handle last grid point: decrement for non-periodic, keep for periodic (will wrap in interpolation)
            if( state.nearestLowerIndices[ i ] == independentValues_[ i ].size( ) - 1 )
            {
                if( !this->isPeriodic( i ) )
                {
                    state.nearestLowerIndices[ i ] -= 1;
                }
            }
        }

        // Pre-compute 2D interpolation fractions if applicable
        if( NumberOfDimensions == 2 )
        {
            state.i0 = state.nearestLowerIndices[ 0 ];
            state.i1 = state.nearestLowerIndices[ 1 ];

            // Determine next indices with wrap-around for periodic dimensions
            state.i0_next = state.i0 + 1;
            state.i1_next = state.i1 + 1;

            if( this->isPeriodic( 0 ) && state.i0 == independentValues_[ 0 ].size( ) - 1 )
            {
                state.i0_next = 0;  // Wrap to first grid point
            }

            if( this->isPeriodic( 1 ) && state.i1 == independentValues_[ 1 ].size( ) - 1 )
            {
                state.i1_next = 0;  // Wrap to first grid point
            }

            // Compute interpolation fractions using shortest distance for periodic dimensions
            IndependentVariableType dx, dy;
            if( this->isPeriodic( 0 ) )
            {
                IndependentVariableType x0 = independentValues_[ 0 ][ state.i0 ];
                IndependentVariableType x1 = independentValues_[ 0 ][ state.i0_next ];
                dx = this->getShortestDistance( 0, x0, x1 );
                IndependentVariableType t_x = this->getShortestDistance( 0, x0, state.localIndependentValues[ 0 ] );
                state.tx = static_cast< IndependentVariableType >( static_cast< double >( t_x ) / static_cast< double >( dx ) );
            }
            else
            {
                dx = independentValues_[ 0 ][ state.i0_next ] - independentValues_[ 0 ][ state.i0 ];
                state.tx = ( state.localIndependentValues[ 0 ] - independentValues_[ 0 ][ state.i0 ] ) / dx;
            }

            if( this->isPeriodic( 1 ) )
            {
                IndependentVariableType y0 = independentValues_[ 1 ][ state.i1 ];
                IndependentVariableType y1 = independentValues_[ 1 ][ state.i1_next ];
                dy = this->getShortestDistance( 1, y0, y1 );
                IndependentVariableType t_y = this->getShortestDistance( 1, y0, state.localIndependentValues[ 1 ] );
                state.ty = static_cast< IndependentVariableType >( static_cast< double >( t_y ) / static_cast< double >( dy ) );
            }
            else
            {
                dy = independentValues_[ 1 ][ state.i1_next ] - independentValues_[ 1 ][ state.i1 ];
                state.ty = ( state.localIndependentValues[ 1 ] - independentValues_[ 1 ][ state.i1 ] ) / dy;
            }

            state.one_minus_tx = IndependentVariableType( 1 ) - state.tx;
            state.one_minus_ty = IndependentVariableType( 1 ) - state.ty;
        }

        state.isValid = true;
        return state;
    }

    //! Interpolate using pre-computed interpolation state (for batch operations)
    /*!
     * This function performs interpolation using a pre-computed state from prepareInterpolationState().
     * The dependent data is accessed using the cached indices and fractions, avoiding redundant lookups.
     * \param state Pre-computed interpolation state
     * \return Interpolated value of dependent variable
     */
    DependentVariableType interpolateWithState( const InterpolationState& state ) const
    {
        if( !state.isValid )
        {
            throw std::runtime_error( "Cannot interpolate with invalid interpolation state" );
        }

        // Use optimized 2D implementation if available
        if( NumberOfDimensions == 2 )
        {
            // Access the four corner values directly using pre-computed indices (with wrap-around)
            std::array< unsigned int, 2 > indices;

            indices[ 0 ] = state.i0;
            indices[ 1 ] = state.i1;
            const DependentVariableType v00 = dependentData_( indices );

            indices[ 0 ] = state.i0_next;  // May be 0 for periodic wrap
            indices[ 1 ] = state.i1;
            const DependentVariableType v10 = dependentData_( indices );

            indices[ 0 ] = state.i0;
            indices[ 1 ] = state.i1_next;  // May be 0 for periodic wrap
            const DependentVariableType v01 = dependentData_( indices );

            indices[ 0 ] = state.i0_next;  // May be 0 for periodic wrap
            indices[ 1 ] = state.i1_next;  // May be 0 for periodic wrap
            const DependentVariableType v11 = dependentData_( indices );

            // Perform bilinear interpolation using pre-computed fractions
            return state.one_minus_tx * state.one_minus_ty * v00 + state.tx * state.one_minus_ty * v10 +
                    state.one_minus_tx * state.ty * v01 + state.tx * state.ty * v11;
        }
        else
        {
            // Fall back to recursive implementation for higher dimensions
            std::array< unsigned int, NumberOfDimensions > interpolationIndices;
            for( unsigned int i = 0; i < NumberOfDimensions; i++ )
            {
                interpolationIndices[ i ] = static_cast< unsigned int >( -1 );
            }
            return performRecursiveInterpolationStep( 0, state.localIndependentValues, interpolationIndices, state.nearestLowerIndices );
        }
    }

    //! Function to perform interpolation.
    /*!
     *  This function performs the multilinear interpolation.
     *  \param independentValuesToInterpolate Vector of values of independent variables at which
     *      the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable in all dimensions.
     */
    DependentVariableType interpolate( const std::vector< IndependentVariableType >& independentValuesToInterpolate )
    {
        // Check whether size of independent variable vector is correct
        if( independentValuesToInterpolate.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error(
                    "Error in multi-dimensional interpolator. The number of independent variables "
                    "provided is incompatible with the previous definition. Provided: " +
                    std::to_string( independentValuesToInterpolate.size( ) ) + ". Needed: " + std::to_string( NumberOfDimensions ) );
        }

        // Create local copy of independent values that can be modified for boundary handling
        std::vector< IndependentVariableType > localIndependentValues = independentValuesToInterpolate;

        // Check that independent variables are in range
        bool useValue = false;
        DependentVariableType currentDependentVariable;
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            this->checkBoundaryCase( i, useValue, localIndependentValues.at( i ), currentDependentVariable );
            if( useValue )
            {
                return currentDependentVariable;
            }
        }

        // Determine the nearest lower neighbours.
        std::vector< unsigned int > nearestLowerIndices;
        nearestLowerIndices.resize( NumberOfDimensions );
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            nearestLowerIndices[ i ] = lookUpSchemes_[ i ]->findNearestLowerNeighbour( localIndependentValues[ i ] );

            // If newNearestLowerIndex is the last element of independentValues_:
            // - For non-periodic dimensions: execute extrapolation with the last and second to last elements
            // - For periodic dimensions: keep at last index (will wrap to first in interpolation functions)
            if( nearestLowerIndices[ i ] == independentValues_[ i ].size( ) - 1 )
            {
                if( !this->isPeriodic( i ) )
                {
                    nearestLowerIndices[ i ] -= 1;
                }
                // else: keep at last index for periodic wrap-around
            }
        }

        // Use optimized 2D implementation if available (avoids recursion overhead)
        if( NumberOfDimensions == 2 )
        {
            return interpolate2DOptimized( localIndependentValues, nearestLowerIndices );
        }
        else
        {
            // Initialize function evaluation indices to -1 for debugging purposes.
            std::array< unsigned int, NumberOfDimensions > interpolationIndices;
            for( unsigned int i = 0; i < NumberOfDimensions; i++ )
            {
                interpolationIndices[ i ] = -1;
            }

            // Call first step of interpolation, this function calls itself at subsequent independent
            // variable dimensions to evaluate and properly scale dependent variable table values at
            // all 2^n grid edges.
            return performRecursiveInterpolationStep( 0, localIndependentValues, interpolationIndices, nearestLowerIndices );
        }
    }

private:
    //! Function to compute the shortest angular distance between two values for a specific dimension.
    /*!
     *  Function to compute the shortest angular distance between two values, considering periodicity.
     *  For periodic dimensions (periods_[dimensionIndex] > 0), this returns the distance via the shortest path,
     *  which may wrap around the periodic boundary. For non-periodic dimensions, returns val2 - val1.
     *  \param dimensionIndex Index of the dimension.
     *  \param val1 First value.
     *  \param val2 Second value.
     *  \return Shortest distance from val1 to val2, considering periodicity.
     */
    IndependentVariableType getShortestDistance( const unsigned int dimensionIndex,
                                                 const IndependentVariableType& val1,
                                                 const IndependentVariableType& val2 ) const
    {
        if( !isPeriodic( dimensionIndex ) )
        {
            return val2 - val1;
        }

        const IndependentVariableType period = periods_.at( dimensionIndex );
        IndependentVariableType diff = val2 - val1;
        const IndependentVariableType halfPeriod = static_cast< IndependentVariableType >( static_cast< double >( period ) / 2.0 );
        const IndependentVariableType negativeHalfPeriod = static_cast< IndependentVariableType >( static_cast< double >( period ) / -2.0 );

        while( diff > halfPeriod )
        {
            diff -= period;
        }
        while( diff < negativeHalfPeriod )
        {
            diff += period;
        }

        return diff;
    }

    //! Optimized 2D bilinear interpolation (non-recursive implementation).
    /*!
     * This function performs bilinear interpolation for 2D case without recursion,
     * reducing function call overhead and improving cache performance. Supports periodic
     * interpolation for one or both dimensions.
     * \param independentValues Vector of independent variable values for interpolation
     * \param nearestLowerIndices Indices of nearest lower neighbors in each dimension
     * \return Interpolated value of dependent variable
     */
    DependentVariableType interpolate2DOptimized( const std::vector< IndependentVariableType >& independentValues,
                                                  const std::vector< unsigned int >& nearestLowerIndices )
    {
        // Pre-compute interpolation fractions for both dimensions
        const unsigned int i0 = nearestLowerIndices[ 0 ];
        const unsigned int i1 = nearestLowerIndices[ 1 ];

        // For periodic dimensions at the last grid point, wrap around to first point
        // For non-periodic dimensions, use the next point (already handled by caller)
        unsigned int i0_next = i0 + 1;
        unsigned int i1_next = i1 + 1;

        // Handle wrap-around for periodic dimension 0
        if( this->isPeriodic( 0 ) && i0 == independentValues_[ 0 ].size( ) - 1 )
        {
            i0_next = 0;  // Wrap to first grid point
        }

        // Handle wrap-around for periodic dimension 1
        if( this->isPeriodic( 1 ) && i1 == independentValues_[ 1 ].size( ) - 1 )
        {
            i1_next = 0;  // Wrap to first grid point
        }

        // Dimension 0 (x-axis)
        IndependentVariableType x0 = independentValues_[ 0 ][ i0 ];
        IndependentVariableType x1 = independentValues_[ 0 ][ i0_next ];
        IndependentVariableType x_target = independentValues[ 0 ];

        IndependentVariableType dx, t_x;
        IndependentVariableType tx;
        if( this->isPeriodic( 0 ) )
        {
            dx = this->getShortestDistance( 0, x0, x1 );
            t_x = this->getShortestDistance( 0, x0, x_target );
            // Convert to double for division to handle custom types like Time
            double txDouble = static_cast< double >( t_x ) / static_cast< double >( dx );
            tx = static_cast< IndependentVariableType >( txDouble );
        }
        else
        {
            dx = x1 - x0;
            t_x = x_target - x0;
            tx = t_x / dx;
        }

        // Dimension 1 (y-axis)
        IndependentVariableType y0 = independentValues_[ 1 ][ i1 ];
        IndependentVariableType y1 = independentValues_[ 1 ][ i1_next ];
        IndependentVariableType y_target = independentValues[ 1 ];

        IndependentVariableType dy, t_y;
        IndependentVariableType ty;
        if( this->isPeriodic( 1 ) )
        {
            dy = this->getShortestDistance( 1, y0, y1 );
            t_y = this->getShortestDistance( 1, y0, y_target );
            // Convert to double for division to handle custom types like Time
            double tyDouble = static_cast< double >( t_y ) / static_cast< double >( dy );
            ty = static_cast< IndependentVariableType >( tyDouble );
        }
        else
        {
            dy = y1 - y0;
            t_y = y_target - y0;
            ty = t_y / dy;
        }

        const IndependentVariableType one_minus_tx = IndependentVariableType( 1 ) - tx;
        const IndependentVariableType one_minus_ty = IndependentVariableType( 1 ) - ty;

        // Access the four corner values directly (avoids recursive calls)
        // For periodic dimensions, indices wrap around automatically via i0_next and i1_next
        std::array< unsigned int, 2 > indices;

        indices[ 0 ] = i0;
        indices[ 1 ] = i1;
        const DependentVariableType v00 = dependentData_( indices );

        indices[ 0 ] = i0_next;
        indices[ 1 ] = i1;
        const DependentVariableType v10 = dependentData_( indices );

        indices[ 0 ] = i0;
        indices[ 1 ] = i1_next;
        const DependentVariableType v01 = dependentData_( indices );

        indices[ 0 ] = i0_next;
        indices[ 1 ] = i1_next;
        const DependentVariableType v11 = dependentData_( indices );

        // Perform bilinear interpolation using the formula:
        // f(x,y) = (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
        return one_minus_tx * one_minus_ty * v00 + tx * one_minus_ty * v10 + one_minus_tx * ty * v01 + tx * ty * v11;
    }

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

                for( unsigned int i = 0; i < NumberOfDimensions; i++ )
                {
                    // Create binary search look up scheme.
                    lookUpSchemes_[ i ] = std::shared_ptr< LookUpScheme< IndependentVariableType > >(
                            new BinarySearchLookupScheme< IndependentVariableType >( independentValues_[ i ] ) );
                }

                break;

            case huntingAlgorithm:

                for( unsigned int i = 0; i < NumberOfDimensions; i++ )
                {
                    // Create hunting scheme, which uses an intial guess from previous look-ups.
                    lookUpSchemes_[ i ] = std::shared_ptr< LookUpScheme< IndependentVariableType > >(
                            new HuntingAlgorithmLookupScheme< IndependentVariableType >( independentValues_[ i ] ) );
                }

                break;

            default:

                throw std::runtime_error( "Warning: lookup scheme not found when making scheme for 1-D interpolator" );
        }
    }

    //! Perform the step in a single dimension of the interpolation process.
    /*!
     * Function calculates single dimension of the interpolation process. Function calls itself if
     * final dimension not yet reached. Calling this function with currentDimension = 0 will result
     * in 2^{NumberOfDimensions} number of calls to the function at currentDimension =
     * NumberOfDimensions -1. As such, the complete series of calls, starting at currentDimension =
     * 0, retrieves the dependent variable values at all edges of the grid hyper-rectangle and
     * properly scales them. For periodic dimensions, uses shortest-path interpolation.
     * \param currentDimension Dimension in which this interpolation step is to be performed.
     * \param independentValuesToInterpolate Vector of values of independent variables at which
     *          interpolation is to be performed.
     * \param currentArrayIndices Array of indices modified at index = currentDimension at each
     *          call of function. Variable is passed to dependentData in highest step to return
     *          data for interpolation.
     * \param nearestLowerIndices Indices in subvectors of independentValues_ vector. That is, the
     *  n-th entry of nearestLowerIndices represent the nearest lower neighbour in the n-th
     *  interpolation dimension of the independent variable vectors.
     * \return Interpolated value in a single dimension
     */
    DependentVariableType performRecursiveInterpolationStep( const unsigned int currentDimension,
                                                             const std::vector< IndependentVariableType >& independentValuesToInterpolate,
                                                             std::array< unsigned int, NumberOfDimensions > currentArrayIndices,
                                                             const std::vector< unsigned int >& nearestLowerIndices ) const
    {
        IndependentVariableType upperFraction, lowerFraction;
        DependentVariableType upperContribution, lowerContribution;

        // Determine the upper index: wrap around for periodic dimensions at last grid point
        unsigned int lowerIndex = nearestLowerIndices[ currentDimension ];
        unsigned int upperIndex = lowerIndex + 1;

        // Handle wrap-around for periodic dimensions
        if( this->isPeriodic( currentDimension ) && lowerIndex == independentValues_[ currentDimension ].size( ) - 1 )
        {
            upperIndex = 0;  // Wrap to first grid point
        }

        // Calculate fractions of data points above and below independent variable value
        // For periodic dimensions, use shortest angular distance; for non-periodic, use linear distance
        IndependentVariableType x0 = independentValues_[ currentDimension ][ lowerIndex ];
        IndependentVariableType x1 = independentValues_[ currentDimension ][ upperIndex ];
        IndependentVariableType x_target = independentValuesToInterpolate[ currentDimension ];

        IndependentVariableType dx;
        IndependentVariableType t_x;

        if( this->isPeriodic( currentDimension ) )
        {
            // Periodic interpolation: use shortest angular distance
            dx = this->getShortestDistance( currentDimension, x0, x1 );
            t_x = this->getShortestDistance( currentDimension, x0, x_target );

            // Compute fraction - convert to double for division to handle custom types like Time
            double fractionDouble = static_cast< double >( t_x ) / static_cast< double >( dx );
            upperFraction = static_cast< IndependentVariableType >( fractionDouble );
        }
        else
        {
            // Non-periodic interpolation: standard linear distance
            dx = x1 - x0;
            t_x = x_target - x0;
            upperFraction = t_x / dx;
        }

        lowerFraction = IndependentVariableType( 1 ) - upperFraction;

        // If at top dimension, call dependent variable data.
        if( currentDimension == NumberOfDimensions - 1 )
        {
            currentArrayIndices[ NumberOfDimensions - 1 ] = lowerIndex;
            lowerContribution = dependentData_( currentArrayIndices );
            currentArrayIndices[ NumberOfDimensions - 1 ] = upperIndex;
            upperContribution = dependentData_( currentArrayIndices );
        }

        // If at lower dimension, update currentArrayIndices and call function with
        // currentDimension++.
        else
        {
            currentArrayIndices[ currentDimension ] = lowerIndex;
            lowerContribution = performRecursiveInterpolationStep(
                    currentDimension + 1, independentValuesToInterpolate, currentArrayIndices, nearestLowerIndices );
            currentArrayIndices[ currentDimension ] = upperIndex;
            upperContribution = performRecursiveInterpolationStep(
                    currentDimension + 1, independentValuesToInterpolate, currentArrayIndices, nearestLowerIndices );
        }

        // Return interpolated value.
        DependentVariableType returnValue = upperFraction * upperContribution + lowerFraction * lowerContribution;
        return returnValue;
    }

    //! Vector of periods for periodic/angular interpolation.
    std::vector< IndependentVariableType > periods_;
};

// extern template class MultiLinearInterpolator< double, Eigen::Vector6d, 1 >;
// extern template class MultiLinearInterpolator< double, Eigen::Vector6d, 2 >;
// extern template class MultiLinearInterpolator< double, Eigen::Vector6d, 3 >;
// extern template class MultiLinearInterpolator< double, Eigen::Vector6d, 4 >;
// extern template class MultiLinearInterpolator< double, Eigen::Vector6d, 5 >;

}  // namespace interpolators

}  // namespace tudat

#endif  // TUDAT_MULTI_LINEAR_INTERPOLATOR_H
