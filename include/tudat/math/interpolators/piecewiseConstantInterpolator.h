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

#ifndef TUDAT_PIECEWISECONSTANTINTERPOLATOR_H
#define TUDAT_PIECEWISECONSTANTINTERPOLATOR_H

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace interpolators
{

//! Piecewise constant interpolator class.
/*!
 *  This class is used to perform piecewise constant interpolation in a single dimension from a set of
 *  data in independent and dependent variables. The interpolated values is equal to the dependent variable at the
 *  nearest lower neighbour at a given independent variable.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class PiecewiseConstantInterpolator : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    // Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor from vectors of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from two vectors containing the independent
     *  variables and dependent variables. A look-up scheme can be provided to
     *  override the given default.
     *  \param independentVariables Vector of values of independent variables that are used, must be
     *      sorted in ascending order.
     *  \param dependentVariables Vector of values of dependent variables that are used.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *      interpolation.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    PiecewiseConstantInterpolator( const std::vector< IndependentVariableType > independentVariables,
                                   const std::vector< DependentVariableType > dependentVariables,
                                   const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                                   const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
                                   const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;

        if( dependentValues_.size( ) != independentValues_.size( ) )
        {
            throw std::runtime_error(
                        "Warning: independent and dependent variables not of same size in piecewise constant interpolator constrcutor" );
        }

        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Constructor from map of independent/dependent data.
    /*!
     *  This constructor initializes the interpolator from a map containing independent variables
     *  as key and dependent variables as value. A look-up scheme can be provided to override the
     *  given default.
     *  \param dataMap Map containing independent variables as key and dependent variables as
     *      value.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *      interpolation.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    PiecewiseConstantInterpolator( const std::map< IndependentVariableType, DependentVariableType > dataMap,
                                   const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                                   const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
                                   const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ):
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling,
                                                                                      defaultExtrapolationValue )
    {
        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            throw std::runtime_error(
                        "The vectors used in the piecewise constant interpolator initialization are empty." );
        }

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
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Destructor
    ~PiecewiseConstantInterpolator( ){ }

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     *  Function interpolates dependent variable value at given independent variable value using piecewise constant algorithm.
     *  \param targetIndependentVariableValue Value of independent variable at which interpolation is to take place.
     *  \return Interpolated value of dependent variable.
     */
    DependentVariableType interpolate( const IndependentVariableType targetIndependentVariableValue )
    {
        // Check whether boundary handling needs to be applied, if independent variable is beyond its defined range.
        DependentVariableType interpolatedValue;
        bool useValue = false;
        this->checkBoundaryCase( interpolatedValue, useValue, targetIndependentVariableValue );
        if( useValue )
        {
            return interpolatedValue;
        }

        // Determine the lower entry in the table corresponding to the target independent variable value.
        int lowerEntry;
        if( targetIndependentVariableValue <= independentValues_.at( 0 ) )
        {
            lowerEntry = 0;
        }
        else if( targetIndependentVariableValue >= independentValues_.at( independentValues_.size( ) - 1 ) )
        {
            lowerEntry = independentValues_.size( ) - 1;
        }
        else
        {
            lowerEntry = lookUpScheme_->findNearestLowerNeighbour( targetIndependentVariableValue );
        }

        // Return interpolated value
        return dependentValues_.at( lowerEntry );
    }

    //! Function to reset the values of dependent variables used by interpolator
    /*!
     *  Function to reset the values of dependent variables used by interpolator
     *  \param dependentValues New list of values of dependent variables. List must be of same size as original list
     */
    void resetDependentValues( const std::vector< DependentVariableType >&  dependentValues )
    {
        if( dependentValues.size( ) != dependentValues_.size( ) )
        {
            throw std::runtime_error(
                        "Error when resetting dependent values for piecewise constant interpolator, sizes are inconsistent " );
        }

        dependentValues_ = dependentValues;
    }

    InterpolatorTypes getInterpolatorType( ){ return piecewise_constant_interpolator; }


protected:

private:

};

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_PIECEWISECONSTANTINTERPOLATOR_H
