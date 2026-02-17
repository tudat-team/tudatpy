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
        if ( chebyshevCoefficients.empty() )
        {
            throw std::runtime_error( "Error: empty Chebyshev coefficient input map." );
        }

        independentValues_.resize( chebyshevCoefficients.size() );
        chebyshevCoefficients_.resize( chebyshevCoefficients.size() );
        dependentValues_.clear();

        int counter = 0;
        for ( const auto& entry : chebyshevCoefficients )
        {
            independentValues_[ counter ] = entry.first;
            chebyshevCoefficients_[ counter ] = entry.second;
            ++counter;
        }

        chebyshevOrder_ = static_cast< int >( chebyshevCoefficients_[ 0 ].size() );
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
        if ( independentValues.empty() || chebyshevCoefficients.empty() )
        {
            throw std::runtime_error( "Error: empty input vectors to ChebyshevInterpolator constructor." );
        }

        independentValues_ = independentValues;
        chebyshevCoefficients_ = chebyshevCoefficients;
        chebyshevOrder_ = static_cast< int >( chebyshevCoefficients_[ 0 ].size() );

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
            std::cerr << "[ChebyshevInterpolator] Warning: extrapolation below lower limit by "
                      << ( independentVariableValue - independentValues_.front() ) << "\n";
        }
        else if ( independentVariableValue >= independentValues_.back() + 
                  0.1 * ( independentValues_.back() - independentValues_[ independentValues_.size() - 2 ] ) )
        {
            std::cerr << "[ChebyshevInterpolator] Warning: extrapolation above upper limit by "
                      << ( independentVariableValue - independentValues_.back() ) << "\n";
        }
        else
        {
            int nearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour( independentVariableValue );
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
