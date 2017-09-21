/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100903    K. Kumar          File header and footer added.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToIntegerExponent( ) function.
 *      102410    D. Dirkx          Minor comment changes during code check.
 *      101213    K. Kumar          Modified raiseToIntegerExponent( ) function;
 *                                  renamed raiseToIntegerPower( ).
 *                                  Added computeAbsoluteValue( ) functions.
 *      110111    J. Melman         Added computeModulo( ) function.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation( ).
 *      110411    K. Kumar          Added convertCartesianToSpherical( ) function.
 *      110707    K. Kumar          Added computeSampleMean( ), computeSampleVariance( ) functions.
 *      110810    J. Leloux         Corrected doxygen documentation (equations).
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved linear interpolation functions into new Interpolators
 *                                  sub-directory.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's
 *          Outlines, McGraw-Hill, 2008.
 *
 */

#ifndef TUDAT_JUMPLINEARINTERPOLATOR_H
#define TUDAT_JUMPLINEARINTERPOLATOR_H

#include <map>
#include <vector>
#include <iomanip>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

namespace tudat
{

//! Interpolators namespace.
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
     * This constructor initializes the interpolator from a map containing independent variables
     * as key and dependent variables as value. A look-up scheme can be provided to override the
     * given default.
     * \param dataMap Map containing independent variables as key and dependent variables as
     *          value.
     * \param maximumAllowableVariation Maximum allowable deviation between two dependent variable values, above which a jump is
     *          identified.
     * \param jumpSize Magnitude of jumps in dependent variable data.
     * \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *          to find the nearest lower data point in the independent variables when requesting
     *          interpolation.
     */
    JumpDataLinearInterpolator( const std::map< IndependentVariableType, DependentVariableType >& dataMap,
                                const DependentVariableType maximumAllowableVariation,
                                const DependentVariableType jumpSize,
                                const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
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
     *  sorted in ascending order.
     *  \param dependentValues Vector of values of dependent variables that are used.
     *  \param maximumAllowableVariation Maximum allowable deviation between two dependent variable values, above which a jump is
     *          identified.
     *  \param jumpSize Magnitude of jumps in dependent variable data.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *  to find the nearest lower data point in the independent variables when requesting
     *  interpolation.
     */
    JumpDataLinearInterpolator( const std::vector< IndependentVariableType > independentValues,
                                const std::vector< DependentVariableType > dependentValues,
                                const DependentVariableType maximumAllowableVariation,
                                const DependentVariableType jumpSize,
                                const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm)
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
        // Lookup nearest lower index.
        int newNearestLowerIndex = lookUpScheme_->findNearestLowerNeighbour( independentVariableValue );

        DependentVariableType interpolatedValue;

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

private:

    //! Maximum allowable deviation between two dependent variable values, above which a jump is identified.
    DependentVariableType maximumAllowableVariation_;

    //! Magnitude of jumps in dependent variable data.
    DependentVariableType jumpSize_;
};

} // Namespace interpolators.

} // Namespace tudat.

#endif // TUDAT_JUMPLINEARINTERPOLATOR_H


