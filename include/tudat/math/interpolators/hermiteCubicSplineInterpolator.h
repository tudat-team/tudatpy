/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H
#define TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H

#include <Eigen/Core>
#include <vector>

#include "tudat/math/interpolators/interpolator.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/math/basic/nearestNeighbourSearch.h"

namespace tudat
{

namespace interpolators
{

//! Hermite Cubic Spline Interpolator
template< typename IndependentVariableType, typename DependentVariableType >
class HermiteCubicSplineInterpolator : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor.
    /*!
     *  Constructor for the hermite cubic spline interpolator with (in)dependent data defined in a vector.
     *  \param independentValues Vector with the independent variable values, must be
     *      sorted in ascending order.
     *  \param dependentValues Vector with the dependent variable values.
     *  \param derivativeValues Vector with the derivative values.
     *  \param selectedLookupScheme Look-up scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    HermiteCubicSplineInterpolator(
            const std::vector< IndependentVariableType >& independentValues,
            const std::vector< DependentVariableType >& dependentValues,
            const std::vector< DependentVariableType >& derivativeValues,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
            const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        // Check consistency of input data.
        if( dependentValues.size( ) != independentValues.size( ) )
        {
            throw std::runtime_error(
                "Error: indep. and dep. variables incompatible in Hermite interpolator." );
        }

        if( dependentValues.size( ) != derivativeValues.size( ) )
        {
            throw std::runtime_error(
                "Error: derivative values incompatible in Hermite interpolator." );
        }

        independentValues_ = std::move( independentValues );
        dependentValues_ = std::move( dependentValues );
        derivativeValues_ = std::move( derivativeValues );

        // Check if data is in ascending order
        if( !std::is_sorted( independentValues_.begin( ), independentValues_.end( ) ) )
        {
            throw std::runtime_error( "Error when making hermite spline spline interpolator, input vector with independent variables should be in ascending order" );
        }

        // compute coefficients
        computeCoefficients( );

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Constructor.
    /*!
     *  Constructor for the hermite cubic spline interpolator with (in)dependent data defined in a map.
     *  \param dataMap Map with the independent variable values as keys and corresponding
     *      dependent variable values as values.
     *  \param derivativeValues
     *  \param selectedLookupScheme Look-up scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    HermiteCubicSplineInterpolator(
            const std::map< IndependentVariableType, DependentVariableType >& dataMap,
            const std::vector< DependentVariableType >& derivativeValues,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
            const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        // Check that sizes match
        if( dataMap.size( ) != derivativeValues.size( ) )
        {
            throw std::runtime_error( "Error: derivative values incompatible in Hermite interpolator." );
        }

        // Fill data vectors with data from map.
        for( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator
             mapIterator = dataMap.begin( ); mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_.push_back( std::move( mapIterator->first ) );
            dependentValues_.push_back( std::move( mapIterator->second ) );
        }

        derivativeValues_ = derivativeValues;

        // Compute coefficients
        computeCoefficients( );

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Destructor
    ~HermiteCubicSplineInterpolator( ){ }

    // Using statement to prevent compiler warning.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Get coefficients
    std::vector< std::vector< DependentVariableType > > GetCoefficients( )
    {
        return coefficients_;
    }

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     *  Function interpolates dependent variable value at given independent variable value.
     *  \param targetIndependentVariableValue Value of independent variable at which interpolation
     *      is to take place.
     *  \return Interpolated value of interpolated dependent variable.
     */
    DependentVariableType interpolate( const IndependentVariableType targetIndependentVariableValue )
    {
        // Check whether boundary handling needs to be applied, if independent variable is beyond its defined range.
        DependentVariableType targetValue;
        bool useValue = false;
        this->checkBoundaryCase( targetValue, useValue, targetIndependentVariableValue );
        if( useValue )
        {
            return targetValue;
        }

        // Determine the lower entry in the table corresponding to the target independent variable value.
        int lowerEntry_ = lookUpScheme_->findNearestLowerNeighbour( targetIndependentVariableValue );

        // Compute Hermite spline
        IndependentVariableType factor = ( targetIndependentVariableValue - independentValues_[ lowerEntry_ ] ) /
                ( independentValues_[ lowerEntry_ + 1 ] - independentValues_[ lowerEntry_ ] );
        targetValue =
                coefficients_[ 0 ][ lowerEntry_ ] * factor * factor * factor +
                coefficients_[ 1 ][ lowerEntry_ ] * factor * factor +
                coefficients_[ 2 ][ lowerEntry_ ] * factor +
                coefficients_[ 3 ][ lowerEntry_ ] ;

        return targetValue;
    }

    InterpolatorTypes getInterpolatorType( ){ return hermite_spline_interpolator; }


    std::vector< DependentVariableType > getDerivativeValues( )
    {
        return derivativeValues_;
    }
protected:

    //! Compute coefficients of the splines
    void computeCoefficients( )
    {
        // Initialize vector
        std::vector< DependentVariableType > zeroVect( independentValues_.size( ) - 1 );
        for( int i = 0 ; i < 4 ; i++ )
        {
            coefficients_.push_back( zeroVect );
        }

        // Compute coefficients for polynomials a, b, c and d so that interpolant p is:
        // p(x) = a((x-x0)/(x1-x0))^3 + b((x-x0)/(x1-x0))^2 + c((x-x0)/(x1-x0)) + d.
        for( unsigned int i = 0 ; i < ( independentValues_.size() - 1 ) ; i++ )
        {
            // Compute coefficient a
            coefficients_[ 0 ][ i ] = 2.0 * dependentValues_[ i ] - 2.0 * dependentValues_[ i + 1 ] +
                    derivativeValues_[ i ] * ( independentValues_[ i + 1 ] - independentValues_[ i ] ) +
                    derivativeValues_[ i + 1 ] * ( independentValues_[ i + 1 ] - independentValues_[ i ] );

            // Compute coefficient b
            coefficients_[ 1 ][ i ] = -3.0 * dependentValues_[ i ] + 3.0 * dependentValues_[ i + 1 ] -
                    2.0 * derivativeValues_[ i ] * ( independentValues_[ i + 1 ] - independentValues_[ i ] ) -
                    derivativeValues_[ i + 1 ] * ( independentValues_[ i + 1 ] - independentValues_[ i ] );

            // Compute coefficient c
            coefficients_[ 2 ][ i ] = derivativeValues_[ i ] * ( independentValues_[ i + 1 ] - independentValues_[ i ] );

            // Compute coefficient d
            coefficients_[ 3 ][ i ] = dependentValues_[ i ];
        }
    }



private:

    //! Derivatives of dependent variable to independent variable
    std::vector< DependentVariableType > derivativeValues_ ;

    //! Coefficients of splines
    std::vector< std::vector< DependentVariableType > > coefficients_ ;

};

//! Typede for cubic hermite spline with double (in)dependent variables.
typedef HermiteCubicSplineInterpolator< double, double > HermiteCubicSplineInterpolatorDouble;

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_HERMITE_CUBIC_SPLINE_INTERPOLATOR_H
