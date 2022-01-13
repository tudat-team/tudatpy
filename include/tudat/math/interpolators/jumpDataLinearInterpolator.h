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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. statistics, Fourth Edition, Schaum's
 *          Outlines, McGraw-Hill, 2008.
 *
 */

#ifndef TUDAT_JUMPLINEARINTERPOLATOR_H
#define TUDAT_JUMPLINEARINTERPOLATOR_H

#include <map>
#include <vector>

#include <Eigen/Core>

#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/math/basic/nearestNeighbourSearch.h"

namespace tudat
{

//! interpolators namespace.
/*!
 * The interpolators namespace.
 */
namespace interpolators
{

//! Linear interpolation class with discrete jumps in data.
/*!
 *  This class is used to perform linear interpolation in a single dimension from a set of
 *  data in independent and dependent variables. The dependent variables contain discrete jumps at specified locations.
 *  The jumps are identified by a maximum allowed deviation between two dependent variable data points. The jump in data is
 *  fixed to a user-supplied value when a jump is identified.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class JumpDataLinearInterpolator : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    //! Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor from map of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from a map containing independent variables
     *  as key and dependent variables as value. A look-up scheme can be provided to override the
     *  given default.
     *  \param dataMap Map containing independent variables as key and dependent variables as value.
     *  \param maximumAllowableVariation Maximum allowable deviation between two dependent variable values, above which a
     *      jump is identified.
     *  \param jumpSize Magnitude of jumps in dependent variable data.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used to find the nearest
     *      lower data point in the independent variables when requesting interpolation.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    JumpDataLinearInterpolator( const std::map< IndependentVariableType, DependentVariableType >& dataMap,
                                const DependentVariableType maximumAllowableVariation,
                                const DependentVariableType jumpSize,
                                const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                                const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
                                const DependentVariableType& defaultExtrapolationValue = IdentityElement::getAdditionIdentity< DependentVariableType >( ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        maximumAllowableVariation_ = maximumAllowableVariation;
        jumpSize_ = jumpSize;

        // Resize data vectors of independent/dependent values.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );

        // Fill data vectors with data from map.
        int counter = 0;
        for( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator mapIterator = dataMap.begin( );
             mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = mapIterator->first;
            dependentValues_[ counter ] = mapIterator->second;
            counter++;
        }

        // Create lookup scheme from independent variable data points.
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Constructor from vectors of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from two vectors containing the independent
     *  variables and dependent variables. A look-up scheme can be provided to
     *  override the given default.
     *  \param independentValues Vector of values of independent variables that are used, must be
     *      sorted in ascending order.
     *  \param dependentValues Vector of values of dependent variables that are used.
     *  \param maximumAllowableVariation Maximum allowable deviation between two dependent variable values, above which a jump is
     *      identified.
     *  \param jumpSize Magnitude of jumps in dependent variable data.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    JumpDataLinearInterpolator( const std::vector< IndependentVariableType > independentValues,
                                const std::vector< DependentVariableType > dependentValues,
                                const DependentVariableType maximumAllowableVariation,
                                const DependentVariableType jumpSize,
                                const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                                const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
                                const DependentVariableType& defaultExtrapolationValue = IdentityElement::getAdditionIdentity< DependentVariableType >( ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        maximumAllowableVariation_ = maximumAllowableVariation;
        jumpSize_ = jumpSize;

        // Set data vectors.
        independentValues_ = independentValues;
        dependentValues_= dependentValues;

        // Create lookup scheme from independent variable values
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Destructor.
    ~JumpDataLinearInterpolator( ) { }

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     *  Function interpolates dependent variable value at given independent variable value.
     *  \param independentVariableValue Value of independent variable at which interpolation
     *  is to take place.
     *  \return Interpolated value of dependent variable.
     */
    DependentVariableType interpolate( const IndependentVariableType independentVariableValue )
    {
        // Check whether boundary handling needs to be applied, if independent variable is beyond its defined range.
        DependentVariableType interpolatedValue;
        bool useValue = false;
        this->checkBoundaryCase( interpolatedValue, useValue, independentVariableValue );
        if( useValue )
        {
            return interpolatedValue;
        }

        // Lookup nearest lower index.
        int newNearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour( independentVariableValue );


        // Check if jump occurs
        if( std::abs( dependentValues_[ newNearestLowerIndex ] - dependentValues_[ newNearestLowerIndex + 1 ] ) >
                maximumAllowableVariation_ )
        {
            double jumpSign = ( dependentValues_[ newNearestLowerIndex ] - dependentValues_[ newNearestLowerIndex + 1 ] ) /
                    std::abs( dependentValues_[ newNearestLowerIndex ] - dependentValues_[ newNearestLowerIndex + 1 ] );
            interpolatedValue = dependentValues_[ newNearestLowerIndex ] +
                    ( independentVariableValue - independentValues_[ newNearestLowerIndex ] ) /
                    ( independentValues_[ newNearestLowerIndex + 1 ] - independentValues_[ newNearestLowerIndex ] ) *
                    ( dependentValues_[ newNearestLowerIndex + 1 ] + jumpSign * jumpSize_ - dependentValues_[ newNearestLowerIndex ] );
        }
        else
        {

            interpolatedValue = dependentValues_[ newNearestLowerIndex ] +
                    ( independentVariableValue - independentValues_[ newNearestLowerIndex ] ) /
                    ( independentValues_[ newNearestLowerIndex + 1 ] - independentValues_[ newNearestLowerIndex ] ) *
                    ( dependentValues_[ newNearestLowerIndex + 1 ] - dependentValues_[ newNearestLowerIndex ] );
        }

        return interpolatedValue;
    }

    InterpolatorTypes getInterpolatorType( )
    {
        return discrete_jump_linear_interpolator;
    }
private:

    //! Maximum allowable deviation between two dependent variable values, above which a jump is identified.
    DependentVariableType maximumAllowableVariation_;

    //! Magnitude of jumps in dependent variable data.
    DependentVariableType jumpSize_;
};

} // Namespace interpolators.

} // Namespace tudat.

#endif // TUDAT_JUMPLINEARINTERPOLATOR_H


