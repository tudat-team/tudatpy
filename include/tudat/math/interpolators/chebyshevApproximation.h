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

#ifndef TUDAT_CHEBYSHEV_APPROXIMATION_H
#define TUDAT_CHEBYSHEV_APPROXIMATION_H

#include <stdexcept>
#include <vector>
#include <map>
#include <iostream>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace interpolators
{

//! Map interpolation epoch to normalized Chebyshev variable.
/*!
 *  \param intervalStart Lower bound of interpolation interval.
 *  \param intervalEnd Upper bound of interpolation interval.
 *  \param evaluationPoint Epoch at which polynomial is evaluated.
 *  \return Normalized Chebyshev independent variable in \f$[-1,1]\f$.
 */
double determineChebyshevPolynomialIndependentVariable(
        const double intervalStart,
        const double intervalEnd,
        const double evaluationPoint );

//! Evaluate Chebyshev polynomial basis up to requested order on a physical interval.
/*!
 *  \param intervalStart Lower bound of interpolation interval.
 *  \param intervalEnd Upper bound of interpolation interval.
 *  \param evaluationPoint Epoch at which basis is evaluated.
 *  \param maximumOrder Maximum polynomial order.
 *  \return Vector of basis values \f$T_0,\ldots,T_N\f$.
 */
std::vector< double > evaluateChebyshevPolynomials(
        const double intervalStart,
        const double intervalEnd,
        const double evaluationPoint,
        const int maximumOrder );

//! Evaluate Chebyshev polynomial basis up to requested order at normalized variable.
/*!
 *  \param independentChebyshevVariable Normalized variable in \f$[-1,1]\f$.
 *  \param maximumOrder Maximum polynomial order.
 *  \return Vector of basis values \f$T_0,\ldots,T_N\f$.
 */
std::vector< double > evaluateChebyshevPolynomials(
        const double independentChebyshevVariable,
        const int maximumOrder );

//! Chebyshev interpolator class.
template< typename IndependentVariableType, typename DependentVariableType >
class ChebyshevInterpolator :
        public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor from map of Chebyshev coefficients.
    /*!
     *  \param chebyshevCoefficients Map from segment start epoch to coefficient vector.
     *  \param selectedLookupScheme Lookup scheme used to locate current interpolation segment.
     */
    ChebyshevInterpolator(
            const std::map< IndependentVariableType, std::vector< DependentVariableType > >& chebyshevCoefficients,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
    {
        if( chebyshevCoefficients.size( ) < 2 )
        {
            throw std::runtime_error(
                        "Error: ChebyshevInterpolator requires at least 2 independent values." );
        }

        independentValues_.resize( chebyshevCoefficients.size( ) );
        chebyshevCoefficients_.resize( chebyshevCoefficients.size( ) );
        dependentValues_.clear();

        int expectedOrder = -1;
        int counter = 0;
        for( const auto& entry : chebyshevCoefficients )
        {
            if( entry.second.empty( ) )
            {
                throw std::runtime_error( "Error: Chebyshev coefficient vector may not be empty." );
            }
            if( expectedOrder < 0 )
            {
                expectedOrder = static_cast< int >( entry.second.size( ) );
            }
            else if( expectedOrder != static_cast< int >( entry.second.size( ) ) )
            {
                throw std::runtime_error(
                            "Error: Inconsistent Chebyshev coefficient vector sizes across segments." );
            }
            independentValues_[ counter ] = entry.first;
            chebyshevCoefficients_[ counter ] = entry.second;
            ++counter;
        }

        chebyshevOrder_ = expectedOrder;
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Constructor from vectors of nodes and coefficients.
    /*!
     *  \param independentValues Ordered list of segment start epochs.
     *  \param chebyshevCoefficients Chebyshev coefficient vector per segment.
     *  \param selectedLookupScheme Lookup scheme used to locate current interpolation segment.
     */
    ChebyshevInterpolator(
            const std::vector< IndependentVariableType >& independentValues,
            const std::vector< std::vector< DependentVariableType > >& chebyshevCoefficients,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
    {
        if( independentValues.size( ) < 2 )
        {
            throw std::runtime_error(
                        "Error: ChebyshevInterpolator requires at least 2 independent values." );
        }
        if( chebyshevCoefficients.size( ) < independentValues.size( ) - 1 )
        {
            throw std::runtime_error(
                        "Error: Insufficient Chebyshev coefficient sets for provided independent values." );
        }
        if( chebyshevCoefficients.empty( ) )
        {
            throw std::runtime_error( "Error: empty input vectors to ChebyshevInterpolator constructor." );
        }
        for( const auto& coefficientSet : chebyshevCoefficients )
        {
            if( coefficientSet.empty( ) )
            {
                throw std::runtime_error( "Error: Chebyshev coefficient vector may not be empty." );
            }
        }
        const int expectedOrder = static_cast< int >( chebyshevCoefficients.front( ).size( ) );
        for( const auto& coefficientSet : chebyshevCoefficients )
        {
            if( static_cast< int >( coefficientSet.size( ) ) != expectedOrder )
            {
                throw std::runtime_error(
                            "Error: Inconsistent Chebyshev coefficient vector sizes across segments." );
            }
        }

        independentValues_ = independentValues;
        chebyshevCoefficients_ = chebyshevCoefficients;
        chebyshevOrder_ = expectedOrder;

        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Interpolate value.
    /*!
     *  \param independentVariableValue Epoch at which interpolated value is requested.
     *  \return Interpolated dependent-variable value.
     */
    DependentVariableType interpolate( const IndependentVariableType independentVariableValue )
    {
        DependentVariableType interpolatedValue{};

        if ( independentVariableValue <= independentValues_.front() - 0.1 * ( independentValues_[ 1 ] - independentValues_[ 0 ] ) )
        {
            throw std::runtime_error( "[ChebyshevInterpolator] extrapolation below lower limit" );
        }
        else if ( independentVariableValue >= independentValues_.back() + 
                  0.1 * ( independentValues_.back() - independentValues_[ independentValues_.size() - 2 ] ) )
        {
            throw std::runtime_error( "[ChebyshevInterpolator] extrapolation above upper limit" );
        }
        else
        {
            int nearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour( independentVariableValue );
            if( nearestLowerIndex < 0 ||
                static_cast< unsigned int >( nearestLowerIndex + 1 ) >= independentValues_.size( ) ||
                static_cast< unsigned int >( nearestLowerIndex ) >= chebyshevCoefficients_.size( ) )
            {
                throw std::runtime_error(
                            "Error: Invalid segment index encountered in ChebyshevInterpolator::interpolate." );
            }
            const auto& coeffs = chebyshevCoefficients_[ nearestLowerIndex ];
            auto chebyshevPolynomials = evaluateChebyshevPolynomials(
                        independentValues_[ nearestLowerIndex ],
                        independentValues_[ nearestLowerIndex + 1 ],
                        independentVariableValue,
                        chebyshevOrder_ );

            interpolatedValue = coeffs[ 0 ] * chebyshevPolynomials[ 0 ];
            for ( std::size_t i = 1; i < coeffs.size(); ++i )
            {
                interpolatedValue += coeffs[ i ] * chebyshevPolynomials[ i ];
            }
        }

        return interpolatedValue;
    }

    InterpolatorTypes getInterpolatorType( )
    {
        return chebyshev_interpolator;
    }

private:
    std::vector< std::vector< DependentVariableType > > chebyshevCoefficients_;
    int chebyshevOrder_;
};

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CHEBYSHEV_APPROXIMATION_H
