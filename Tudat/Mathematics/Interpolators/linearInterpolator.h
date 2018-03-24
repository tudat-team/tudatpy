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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing, Cambridge
 *          University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's Outlines, McGraw-Hill,
 *          2008.
 *
 */

#ifndef TUDAT_LINEAR_INTERPOLATOR_H
#define TUDAT_LINEAR_INTERPOLATOR_H

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <map>
#include <vector>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"

namespace tudat
{
namespace interpolators
{

//! Linear interpolator class.
/*!
 * This class is used to perform linear interpolation in a single dimension from a set of
 * data in independent and dependent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class LinearInterpolator
        : public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    // Using statements to prevent having to put 'this' everywhere in the code.
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor from map of independent/dependent data.
    /*!
     * This constructor initializes the interpolator from a map containing independent variables
     * as key and dependent variables as value. A look-up scheme can be provided to override the
     * given default.
     * \param dataMap Map containing independent variables as key and dependent variables as
     *          value.
     * \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *          to find the nearest lower data point in the independent variables when requesting
     *          interpolation.
     */
    LinearInterpolator( const std::map< IndependentVariableType, DependentVariableType >& dataMap,
                        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
    {
        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            throw std::runtime_error(
                    "The vectors used in the linear interpolator initialization are empty." );
        }

        // Resize data vectors of independent/dependent values.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );

        // Fill data vectors with data from map.
        int counter = 0;
        for ( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator
              mapIterator = dataMap.begin( ); mapIterator != dataMap.end( ); mapIterator++ )
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
     *  sorted in ascending order.
     *  \param dependentValues Vector of values of dependent variables that are used.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     */
    LinearInterpolator( const std::vector< IndependentVariableType >& independentValues,
                        const std::vector< DependentVariableType >& dependentValues,
                        const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm )
    {
        // Verify that the initialization variables are not empty.
        if ( independentValues.size( ) == 0 || dependentValues.size( ) == 0 )
        {
            throw std::runtime_error(
                    "The vectors used in the linear interpolator initialization are empty." );
        }

        // Set data vectors.
        independentValues_ = independentValues;
        dependentValues_= dependentValues;

        // Check if data is in ascending order
        if( !std::is_sorted( independentValues_.begin( ), independentValues_.end( ) ) )
        {
            throw std::runtime_error( "Error when making linear interpolator, input vector with independent variables should be in ascending order" );
        }

        // Create lookup scheme from independent variable values
        this->makeLookupScheme( selectedLookupScheme );
    }

    //! Default destructor
    /*!
     *  Default destructor
     */
    ~LinearInterpolator( ){ }

    //! Function interpolates dependent variable value at given independent variable value.
    /*!
     * Function interpolates dependent variable value at given independent variable value.
     * \param independentVariableValue Value of independent variable at which interpolation
     *          is to take place.
     * \return Interpolated value of dependent variable.
     */
    DependentVariableType interpolate( const IndependentVariableType independentVariableValue )
    {
        // Lookup nearest lower index.
        int newNearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour(
                    independentVariableValue );

        // Perform linear interpolation.
        DependentVariableType interpolatedValue = dependentValues_[ newNearestLowerIndex ] +
                ( independentVariableValue - independentValues_[ newNearestLowerIndex ] ) /
                ( independentValues_[ newNearestLowerIndex + 1 ] -
                  independentValues_[ newNearestLowerIndex ] ) *
                ( dependentValues_[ newNearestLowerIndex + 1 ] -
                  dependentValues_[ newNearestLowerIndex ] );

        return interpolatedValue;
    }
};

//! Typedef for linear interpolator with (in)dependent variable = double.
typedef LinearInterpolator< double, double > LinearInterpolatorDouble;

//! Typedef for shared-pointer to linear interpolator with (in)dependent variable = double.
typedef boost::shared_ptr< LinearInterpolatorDouble > LinearInterpolatorDoublePointer;

//! Compute linear interpolation free function.
/*!
 * Computes linear interpolation of data provided in the form of a vector of
 * sorted indepedent variables and an associated vector of dependent variables.
 * The linear interpolation equation used is:
 * \f[
 *      y_{target} = x_{1} * ( 1 - mu ) + x_{2} * mu
 * \f]
 * where \f$ mu = \frac{ x_{target} - x_{1} } { x_{2} + x_{1} } \f$
 * and \f$ x_{2} > x_{1} \f$.
 * \param sortedIndependentVariables Vector of independent variables sorted.
 *          in ascending/descending order.
 * \param associatedDependentVariables Vector of dependent variables
 *          associated with sorted vector of independent variables.
 * \param targetIndependentVariableValue Target independent variable value
 *          in vector of sorted independent variables.
 * \return Value of dependent variable associated with target independent
 *          value in vector of sorted independent variables.
 */
double computeLinearInterpolation( const Eigen::VectorXd& sortedIndependentVariables,
                                   const Eigen::VectorXd& associatedDependentVariables,
                                   const double targetIndependentVariableValue );

//! Compute linear interpolation free function.
/*!
 * Computes linear interpolation of data provided in the form of a map of
 * independent variables and associated vectors of dependent variables.
 * The linear interpolation equation used is:
 * \f[
 *      y_{target} = x_{1} * ( 1 - mu ) + x_{2} * mu
 * \f]
 * where \f$ \mu = \frac{ x_{target} - x_{1} } { x_{2} + x_{1} } \f$
 * and \f$ x_{2} > x_{1} \f$.
 * \param sortedIndepedentAndDependentVariables Map of sorted independent
 *              variables, in ascending/descending order, and associated
 *              dependent variables.
 * \param targetIndependentVariableValue Target independent variable value
 *              in vector of sorted independent variables.
 * \return Vector of dependent variable associated with target independent
 *              value in vector of sorted independent variables.
 */
Eigen::VectorXd computeLinearInterpolation(
        const std::map< double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        const double targetIndependentVariableValue );

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_LINEAR_INTERPOLATOR_H
