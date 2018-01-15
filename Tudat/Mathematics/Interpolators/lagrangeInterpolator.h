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
 *      http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
 *
 */

#ifndef TUDAT_LAGRANGEINTERPOLATOR_H
#define TUDAT_LAGRANGEINTERPOLATOR_H

#include <iostream>

#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"


namespace tudat
{

namespace interpolators
{

//! Enum defining the manner in which to handle the edges of the Lagrange interpolator
enum LagrangeInterpolatorBoundaryHandling
{
    lagrange_cubic_spline_boundary_interpolation = 0,
    lagrange_no_boundary_interpolation = 1
};

//! Class to perform Lagrange polynomial interpolation
/*!
 *  Class to perform Lagrange polynomial interpolation from a set of independent and
 *  dependent values, as well as the order of the interpolation. Note that this class is optimized
 *  for many function calls to interpolate, since the denominators for
 *  the interpolations are pre-computed for all interpolation intervals.
 *  See e.g. http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html for
 *  mathematical details
 */
template< typename IndependentVariableType, typename DependentVariableType,
          typename ScalarType = IndependentVariableType >
class LagrangeInterpolator : public OneDimensionalInterpolator< IndependentVariableType,
        DependentVariableType >
{
public:

    //! Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor from vectors of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from two vectors containing the independent
     *  variables and dependent variables. A look-up scheme can be provided to override the
     *  given default.
     *  \param independentVariables Vector of values of independent variables that are used, must be
     *  sorted in ascending order.
     *  \param dependentVariables Vector of values of dependent variables that are used.
     *  \param numberOfStages Number of data points that are used to calculate the interpolating
     *  polynomial (must be even).
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     *  \param boundaryHandling Boundary handling does something.
     */
    LagrangeInterpolator( const std::vector< IndependentVariableType >& independentVariables,
                          const std::vector< DependentVariableType >& dependentVariables,
                          const int numberOfStages,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                          const LagrangeInterpolatorBoundaryHandling boundaryHandling =
            lagrange_cubic_spline_boundary_interpolation ):
        numberOfStages_( numberOfStages ), boundaryHandling_( boundaryHandling )
    {
        if( numberOfStages_ % 2 != 0 )
        {
            throw std::runtime_error(
                        "Error: Lagrange interpolator currently only handles even orders." );
        }

        // Set data vectors.
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;
        numberOfIndependentValues_ = static_cast< int >( independentValues_.size( ) );

        // Check if data is in ascending order
        if( !std::is_sorted( independentValues_.begin( ), independentValues_.end( ) ) )
        {
            throw std::runtime_error( "Error when making lagrange interpolator, input vector with independent variables should be in ascending order" );
        }

        // Verify that the initialization variables are not empty.
        if ( numberOfIndependentValues_ == 0 || dependentValues_.size( ) == 0 )
        {
            throw std::runtime_error(
                "Error: Vectors used in the Lagrange interpolator initialization are empty." );
        }

        // Check consistency of input data.
        if( static_cast< int >( dependentValues_.size( ) ) != numberOfIndependentValues_ )
        {
            throw std::runtime_error(
                "Error: indep. and dep. variables incompatible in Lagrange interpolator." );
        }

        // Define zero entry for dependent variable.
        zeroEntry_ = dependentVariables[ 0 ] - dependentVariables[ 0 ];
        if( zeroEntry_ != zeroEntry_ )
        {
            throw std::runtime_error(
                "Error: Lagrange interpolator cannot identify zero entry." );
        }

        // Create lookup scheme from independent variable values.
        this->makeLookupScheme( selectedLookupScheme );

        // Calculate denominators for each interval, to prevent recalculations dueint each
        // interpolation call.
        initializeDenominators( );
        initializeBoundaryInterpolators( selectedLookupScheme );

        // Pre-allocate cache vector for computational efficiency.
        independentVariableDifferenceCache.resize( 2 * offsetEntries_ + 2 );
    }

    //! Constructor from map of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from a map containing independent variables
     *  as key and dependent variables as value. A look-up scheme can be provided to override the
     *  given default.
     *  \param dataMap Map containing independent variables as key and dependent variables as
     *  value.
     *  \param numberOfStages Number of data points that are used to calculate the interpolating
     *  polynomial (must be even).
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     *  \param boundaryHandling Boundary Handling does something.
     */
    LagrangeInterpolator(
            const std::map< IndependentVariableType, DependentVariableType >& dataMap,
            const int numberOfStages,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const LagrangeInterpolatorBoundaryHandling boundaryHandling =
            lagrange_cubic_spline_boundary_interpolation ):
        numberOfStages_( numberOfStages ), boundaryHandling_( boundaryHandling )
    {
        if( numberOfStages_ % 2 != 0 )
        {
            throw std::runtime_error(
                        "Error: Lagrange interpolator currently only handles even orders." );
        }

        numberOfIndependentValues_ = dataMap.size( );

        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            throw std::runtime_error(
                        "The vectors used in the lagrange interpolator initialization are empty" );
        }

        // Fill data vectors with data from map.
        for( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator
             mapIterator = dataMap.begin( ); mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_.push_back( mapIterator->first );
            dependentValues_.push_back( mapIterator->second );
        }

        // Define zero entry for dependent variable.
        zeroEntry_ = dependentValues_[ 0 ] - dependentValues_[ 0 ];
        if( zeroEntry_ != zeroEntry_ )
        {
            throw std::runtime_error(
                        "Error: Lagrange interpolator cannot identify zero entry." );
        }

        // Create lookup scheme from independent variable data points.
        this->makeLookupScheme( selectedLookupScheme );

        // Calculate denominators for each interval, to prevent recalculations dueint each
        //interpolation call.
        initializeDenominators( );
        initializeBoundaryInterpolators( selectedLookupScheme );

        independentVariableDifferenceCache.resize( 2 * offsetEntries_ + 2 );
    }

    //! Destructor.
    ~LagrangeInterpolator( ){ }

    // Using statement to prevent compiler warning.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     *  Function interpolates dependent variable value at given independent variable value.
     *  The polynomial centered on the requested interval is used for the interpolation.
     *  Note that the number of data points, not the values of the independent variables
     *  are used for determining the center (in case of a non-equispaced grid). If the required
     *  interpolating polynimial goes beyond the independent variable bondaries,
     *  a cubic spline with natural boundary conditions is used.
     *  \param targetIndependentVariableValue Value of independent variable at which interpolation
     *  is to take place.
     *  \return Interpolated value of dependent variable.
     */
    DependentVariableType interpolate(
            const IndependentVariableType targetIndependentVariableValue )
    {
        using std::pow;

        if( targetIndependentVariableValue < independentValues_.at( 0 ) ||
                targetIndependentVariableValue > independentValues_.at( independentValues_.size( ) -1 ) )
        {
            std::cout << "Warning in Lagrange interpolation, outside range " <<
                       independentValues_.at( 0 ) << " " << independentValues_.at( independentValues_.size( ) -1 ) << " " <<
                       targetIndependentVariableValue << std::endl;
        }
        // Determine the lower entry in the table corresponding to the target independent variable
        // value.
        DependentVariableType interpolatedValue = zeroEntry_;

        // Find interpolation interval
        int lowerEntry = lookUpScheme_->findNearestLowerNeighbour(
                    targetIndependentVariableValue );

        // Check if requested interval is inside region in which centered lagrange interpolation
        // can be used.
        if( lowerEntry < offsetEntries_ )
        {
            if( boundaryHandling_ == lagrange_no_boundary_interpolation )
            {
                throw std::runtime_error(
                            "Error: Lagrange interpolator below allowed bounds." );
            }
            else if( numberOfStages_ > 2 )
            {
                interpolatedValue = beginInterpolator_->interpolate(
                            targetIndependentVariableValue );
            }
        }
        else if( lowerEntry >= numberOfIndependentValues_ - offsetEntries_ - 1 )
        {
            if( boundaryHandling_ == lagrange_no_boundary_interpolation )
            {
                throw std::runtime_error(
                            "Error: Lagrange interpolator above allowed bounds." );
            }
            else if( numberOfStages_ > 2 )
            {
                interpolatedValue = endInterpolator_->interpolate(
                            targetIndependentVariableValue );
            }
        }
        else
        {
            // Initialize repeated numerator to 1
            ScalarType repeatedNumerator =
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 );

            // Check if requested independent variable is equal to data point
            if( independentValues_[ lowerEntry ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry ];
            }
            else if( independentValues_[ lowerEntry + 1 ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry + 1 ];
            }
            else if( independentValues_[ lowerEntry - 1 ] == targetIndependentVariableValue )
            {
                interpolatedValue = dependentValues_[ lowerEntry - 1 ];
            }
            else
            {
                // Set up repeated numerator and cache of independent variable values from which
                // interpolant is created.
                int j = 0;
                for( int i = 0; i <= 2 * offsetEntries_ + 1; i++ )
                {
                    j = i + lowerEntry - offsetEntries_;
                    independentVariableDifferenceCache[ i ] =
                            static_cast< ScalarType >(
                                targetIndependentVariableValue - independentValues_[ j ] );

                    repeatedNumerator *= independentVariableDifferenceCache[ i ];

                }

                // Evaluate interpolating polynomial at requested data point.
                for( int i = 0; i <=  2 *offsetEntries_ + 1; i++ )
                {
                    j = i + lowerEntry - offsetEntries_;
                    interpolatedValue += dependentValues_[ j ]  *
                            ( repeatedNumerator /
                              ( independentVariableDifferenceCache[ i ] *
                                denominators[ lowerEntry ][ j - lowerEntry + offsetEntries_ ] ) );
                }
            }
        }

        return interpolatedValue;
    }

    //! Function to retrieve the number of stages of interpolator
    /*!
     * Function to retrieve the number of stages of interpolator
     * \return Number of stages of interpolator
     */
    int getNumberOfStages( )
    {
        return numberOfStages_;
    }


protected:

private:

    //! Function called at initialization which pre-computes the denominators of the
    //! interpolants at each interval.
    /*!
     *  Function called at initialization which pre-computes the denominators of the interpolants
     *  at each interval, i.e. each interval between two subsequent independent variable values.
     */
    void initializeDenominators( )
    {
        // Check validity of requested number of stages"
        if( numberOfStages_% 2 != 0 )
        {
            throw std::runtime_error(
                        "Error, Lagrange interp. only implemented for even number of stages." );
        }
        if( numberOfStages_ < 2 )
        {
            throw std::runtime_error(
                        "Error, Lagrange interplator number of stages must be greater than 2." );
        }

        // Determine offset from boundary of interpolation interval where interpolant is valid.
        offsetEntries_ = numberOfStages_ / 2 - 1;

        // Iterate over all intervals and calculate denominators
        int currentIterationStart;
        denominators.resize( numberOfIndependentValues_ );
        for( int i = offsetEntries_; i <= numberOfIndependentValues_ - offsetEntries_; i++ )
        {
            // Determine start index in independent variables for current polynomial
            currentIterationStart = i - offsetEntries_;

            denominators[ i ].resize( 2 * offsetEntries_ + 2 );

            // Calculate all denominators for single interval.
            for( int j = 0; j <= 2 * offsetEntries_ + 1; j++ )
            {
                denominators[ i ][ j ] =
                        mathematical_constants::getFloatingInteger< ScalarType >( 1 );

                for( int k = 0; k <= 2 * offsetEntries_ + 1; k++ )
                {
                    if( k != j )
                    {
                        denominators[ i ][ j ] *= static_cast< ScalarType >(
                                    independentValues_[ j + currentIterationStart ] -
                                    independentValues_[ k + currentIterationStart ] );
                    }
                }
            }
        }
    }

    //! Function called at initialization which creates the interpolators used at the boundaries
    //! of the interpolation domain.
    /*!
     *  Function called at initialization which creates the interpolators used at the boundaries
     *  of the interpolation domain.
     *  At the edges of the domain, there is no center interval available, so Runge's phenomenon
     *  can cause excessive interpolation errors, especially for higher order polynomials. In
     *  these regions, the interpolator applies any of a number of techniques, defined by the
     *  boundaryHandling_ variable.
     *  \param selectedLookupScheme Selected lookup scheme does something.
     */
    void initializeBoundaryInterpolators(
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
    {
        // Create interpolators
        if( boundaryHandling_ == lagrange_cubic_spline_boundary_interpolation )
        {
            // Ensure sufficient data points for spline.
            int cubicSplineInputSize = offsetEntries_;
            if( cubicSplineInputSize < 3 )
            {
                cubicSplineInputSize = 3;
            }

            // Set input maps for interpolators
            std::map< IndependentVariableType, DependentVariableType > startMap;
            for( int i = 0; i <= cubicSplineInputSize; i++ )
            {
                startMap[ independentValues_[ i ] ] = dependentValues_[ i ];
            }
            std::map< IndependentVariableType, DependentVariableType > endMap;
            for( int i = numberOfIndependentValues_ - cubicSplineInputSize - 1;
                 i < numberOfIndependentValues_; i++ )
            {
                endMap[ independentValues_[ i ] ] = dependentValues_[ i ];
            }

            // Create cubic spline interpolators
            beginInterpolator_ = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType, ScalarType > >( startMap );
            endInterpolator_ = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType, ScalarType > >( endMap );
        }
    }

    //! Pre-computed denominators to be used in interpolation
    std::vector< std::vector< ScalarType > > denominators;

    //! Zero entry for dependent variables
    /*!
     *  Zero entry for dependent variables, i.e. algebraic identity element for addition of
     *  dependent variables.
     */
    DependentVariableType zeroEntry_;

    //! Number of stages of interpolator
    /*!
     *  Number of stages of interpolator, i.e. order of interpolating polynomial/number of
     *  intervals used to create interpolant.
     */
    int numberOfStages_;

    //! Number of entries at edges of domain where Lagrange interpolation is not used directly.
    /*!
     *  Number of entries at edges of domain where Lagrange interpolation is not used directly.
     *  \sa initializeBoundaryInterpolators
     */
    int offsetEntries_;

    std::vector< ScalarType > independentVariableDifferenceCache;

    //! Interpolator to be used at beginning of domain.
    boost::shared_ptr< OneDimensionalInterpolator
    < IndependentVariableType, DependentVariableType > > beginInterpolator_;

    //! Interpolator to be used at end of domain.
    boost::shared_ptr< OneDimensionalInterpolator
    < IndependentVariableType, DependentVariableType > > endInterpolator_;

    //! Size of (in)dependent variable vector
    int numberOfIndependentValues_;

    //! Method to be used for handling boundaries of the interpolation domain.
    /*!
     *  Method to be used for handling boundaries of the interpolation domain.
     *  \sa initializeBoundaryInterpolators
     */
    LagrangeInterpolatorBoundaryHandling boundaryHandling_;

};

//! Typedef for LagrangeInterpolator with double as both its dependent and independent data type.
typedef LagrangeInterpolator< double, double > LagrangeInterpolatorDouble;

} // namespace interpolators

} // namespace tudat
#endif // TUDAT_LAGRANGEINTERPOLATOR_H
