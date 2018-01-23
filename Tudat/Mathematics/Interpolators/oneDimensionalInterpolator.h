/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H
#define TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"

namespace tudat
{
namespace interpolators
{

//! Base class for interpolator with one independent independent variable.
/*!
 * Base class for the interpolators in one independent variable included in Tudat
 * \tparam IndependentVariableType Type of independent variable(s)
 * \tparam IndependentVariableType Type of dependent variable
 */
template< typename IndependentVariableType, typename DependentVariableType >
class OneDimensionalInterpolator :
        public Interpolator< IndependentVariableType, DependentVariableType >
{

public:

    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~OneDimensionalInterpolator( ) { }

    //! Function to perform interpolation.
    /*!
     * This function performs the interpolation. It calls the function that takes a single
     * independent variable value, which is to be implemented in derived classes.
     * \param independentVariableValues Vector of values of independent variables at which
     *          the value of the dependent variable is to be determined.
     * \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType
    interpolate( const std::vector< IndependentVariableType >& independentVariableValues )
    {
        // Check whether input is really 1-dimensional
        if ( independentVariableValues.size( ) != 1  )
        {
            throw std::runtime_error(
                                "Error in 1-dimensional interpolator, provided input is not 1-dimensional." );
        }

        // Call 1-dimensional interpolate function.
        return interpolate( independentVariableValues[ 0 ] );
    }

    //! Function to perform interpolation.
    /*!
     * This function performs the interpolation
     * \param independentVariableValue Independent variable value at which the value of the
     *          dependent variable is to be determined.
     * \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType
            interpolate( const IndependentVariableType independentVariableValue ) = 0;

    //! Function to perform interpolation, with non-const input argument.
    /*!
     * This function performs the interpolation, with non-const input argument. Function calls the interpolate function and is
     * included for compatibility with some function pointer binding interfaces.
     * \param independentVariableValue Independent variable value at which the value of the
     *          dependent variable is to be determined.
     * \return Interpolated value of dependent variable.
     */
    DependentVariableType
            interpolateNonConst( IndependentVariableType independentVariableValue )
    {
        return interpolate( independentVariableValue );
    }

    //! Function to return the number of independent variables of the interpolation.
    /*!
     *  Function to return the number of independent variables of the interpolation, which is always
     *  equal to 1 for this class and its derived class.
     *  \return Number of independent variables of the interpolation (1).
     */
    int getNumberOfDimensions( )
    {
        return 1;
    }

    //! Function to return the lookup scheme used by the interpolator.
    /*!
     *  Function to return the lookup scheme used by the interpolator.
     *  \return The lookup scheme used by the interpolator.
     */
    boost::shared_ptr< LookUpScheme< IndependentVariableType > > getLookUpScheme( )
    {
        return lookUpScheme_;
    }

    //! Function to return the ector with independent variables used by the interpolator.
    /*!
     *  Function to return the ector with independent variables used by the interpolator.
     *  \return Independent variables used by the interpolator.
     */
    std::vector< IndependentVariableType > getIndependentValues( )
    {
        return independentValues_;
    }

    //! Function to return the ector with dependent variables used by the interpolator.
    /*!
     *  Function to return the ector with dependent variables used by the interpolator.
     *  \return Dependent variables used by the interpolator.
     */
    std::vector< DependentVariableType > getDependentValues( )
    {
        return dependentValues_;
    }

protected:

    //! Make look-up scheme that is to be used.
    /*!
     * This function creates the look-up scheme that is to be used in determining the interval of
     * the independent variable grid where the interpolation is to be performed. It takes the type
     * of lookup scheme as an enum and constructs the look-up scheme from the independentValues_
     * that have been set previously.
     *  \param selectedScheme Type of look-up scheme that is to be used
     */
    void makeLookupScheme( const AvailableLookupScheme selectedScheme )
    {
        // Find which type of scheme is used.
        switch( selectedScheme )
        {
        case binarySearch:

            // Create binary search look up scheme.
            lookUpScheme_ = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                    ( new BinarySearchLookupScheme< IndependentVariableType >
                      ( independentValues_ ) );
            break;

        case huntingAlgorithm:

            // Create hunting scheme, which uses an intial guess from previous look-ups.
            lookUpScheme_ = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                    ( new HuntingAlgorithmLookupScheme< IndependentVariableType >
                      ( independentValues_ ) );
            break;

        default:
            throw std::runtime_error( "Warning: lookup scheme not found when making scheme for 1-D interpolator" );
        }
    }

    //! Pointer to look up scheme.
    /*!
     * Pointer to the lookup scheme that is used to determine in which interval the requested
     * independent variable value falls.
     */
    boost::shared_ptr< LookUpScheme< IndependentVariableType > > lookUpScheme_;

    //! Vector with dependent variables.
    /*!
     * Vector with dependent variables.
     */
    std::vector< DependentVariableType > dependentValues_;

    //! Vector with independent variables.
    /*!
     * Vector with independent variables.
     */
    std::vector< IndependentVariableType > independentValues_;
};

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H
