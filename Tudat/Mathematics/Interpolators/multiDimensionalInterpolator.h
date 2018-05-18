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

#ifndef TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H
#define TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H

#include <vector>
#include <iostream>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Basics/identityElements.h"

namespace tudat
{

namespace interpolators
{

//! Base class for interpolator with multiple independent variables.
/*!
 * Base class for the interpolators in multiple independent variables included in Tudat.
 * \tparam IndependentVariableType Type of independent variable(s)
 * \tparam IndependentVariableType Type of dependent variable
 * \tparam NumberOfDimensions Number of independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType, int NumberOfDimensions >
class MultiDimensionalInterpolator :
        public Interpolator< IndependentVariableType, DependentVariableType >
{

public:

    // Using statements to prevent having to put 'this' everywhere in the code.
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor.
    /*!
     * Constructor.
     */
    MultiDimensionalInterpolator(
            const std::vector< BoundaryInterpolationType > boundaryHandling =
            std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary ),
            const DependentVariableType defaultExtrapolationValue =
            IdentityElement< DependentVariableType >::getAdditionIdentity( ) ) :
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    { }

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~MultiDimensionalInterpolator( ) { }

    //! Function to perform interpolation.
    /*!
     * This function performs the interpolation. It calls the function that takes a vector of
     * independent variable values, which is to be implemented in derived classes.
     * \param independentVariableValues Vector of values of independent variables at which
     *          the value of the dependent variable is to be determined.
     * \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType
    interpolate( const std::vector< IndependentVariableType >& independentVariableValues ) = 0;

    //! Function to return the number of independent variables of the interpolation.
    /*!
     *  Function to return the number of independent variables of the interpolation, i.e. size
     *  that the vector used as input for Interpolator::interpolate should be.
     *  \return Number of independent variables of the interpolation.
     */
    int getNumberOfDimensions( )
    {
        return NumberOfDimensions;
    }

    //! Function to return the lookup scheme used by the interpolator.
    /*!
     *  Function to return the lookup scheme used by the interpolator.
     *  \return The lookup scheme used by the interpolator.
     */
    std::vector< boost::shared_ptr< LookUpScheme< IndependentVariableType > > > getLookUpScheme( )
    {
        return lookUpSchemes_;
    }

    //! Function to return the ector with independent variables used by the interpolator.
    /*!
     *  Function to return the ector with independent variables used by the interpolator.
     *  \return Independent variables used by the interpolator.
     */
    std::vector< std::vector< IndependentVariableType > > getIndependentValues( )
    {
        return independentValues_;
    }

    //! Function to return the ector with dependent variables used by the interpolator.
    /*!
     *  Function to return the ector with dependent variables used by the interpolator.
     *  \return Dependent variables used by the interpolator.
     */
    boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) > getDependentValues( )
    {
        return dependentData_;
    }

protected:

    //! Function to return the condition of the current independent variable.
    /*!
     *  Function to return the condition of the current independent variable, i.e. whether the
     *  variable is within, above or below its defined range range.
     *  \param independentVariable Value of current independent variable.
     *  \param currentVariable Value of current dimension.
     *  \return Condition with respect to boundary.
     */
    int checkInterpolationBoundary( const IndependentVariableType& independentVariable,
                                    const unsigned int& currentVariable )
    {
        int isAtBoundary = 0;
        if ( independentVariable < independentValues_.at( currentVariable ).front( ) )
        {
            isAtBoundary = -1;
        }
        else if ( independentVariable > independentValues_.at( currentVariable ).back( ) )
        {
            isAtBoundary = 1;
        }
        return isAtBoundary;
    }

    //! Function to check whether boundary handling needs to be applied, depending on method chosen.
    /*!
     *  Function to check whether boundary handling needs to be applied, depending on method chosen.
     *  If independent variable is beyond its range definition, boundary handling will be applied, depending
     *  on the method specified in boundaryHandling_.
     *  \param independentVariable Value of current independent variable.
     *  \param currentVariable Value of current dimension.
     */
    void checkBoundaryCase(
            IndependentVariableType& independentVariable,
            const unsigned int& currentVariable,
            bool& useDefault )
    {
        if ( boundaryHandling_.at( currentVariable ) != extrapolate_at_boundary )
        {
            int isAtBoundary = checkInterpolationBoundary( independentVariable, currentVariable );

            if ( isAtBoundary != 0 )
            {
                if ( boundaryHandling_.at( currentVariable ) == throw_exception_at_boundary )
                {
                    std::string errorMessage = "Error in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                            boost::lexical_cast< std::string >( currentVariable ) + " at: " +
                            boost::lexical_cast< std::string >( independentVariable ) + " but limit values are " +
                            boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).front( ) ) + " and " +
                            boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).back( ) );
                    throw std::runtime_error( errorMessage );
                }
                else if ( boundaryHandling_.at( currentVariable ) == extrapolate_at_boundary_with_warning )
                {
                    std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                            boost::lexical_cast< std::string >( currentVariable ) + " at: " +
                            boost::lexical_cast< std::string >( independentVariable ) + " but limit values are " +
                            boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).front( ) ) + " and " +
                            boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).back( ) ) + ", applying extrapolation instead.";
                    std::cerr << errorMessage << std::endl;
                }
                else if ( ( boundaryHandling_.at( currentVariable ) == use_boundary_value ) ||
                          ( boundaryHandling_.at( currentVariable ) == use_boundary_value_with_warning ) )
                {
                    if ( boundaryHandling_.at( currentVariable ) == use_boundary_value_with_warning )
                    {
                        std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                boost::lexical_cast< std::string >( currentVariable ) + " at: " +
                                boost::lexical_cast< std::string >( independentVariable ) + " but limit values are " +
                                boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).front( ) ) + " and " +
                                boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).back( ) ) + ", taking boundary value instead.";
                        std::cerr << errorMessage << std::endl;
                    }

                    if ( isAtBoundary == -1 )
                    {
                        independentVariable = independentValues_.at( currentVariable ).front( );
                    }
                    else if ( isAtBoundary == 1 )
                    {
                        independentVariable = independentValues_.at( currentVariable ).back( );
                    }
                }
                else if ( ( boundaryHandling_.at( currentVariable ) == use_default_value ) ||
                          ( boundaryHandling_.at( currentVariable ) == use_default_value_with_warning ) )
                {
                    if ( boundaryHandling_.at( currentVariable ) == use_default_value_with_warning )
                    {
                        std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                boost::lexical_cast< std::string >( currentVariable ) + " at: " +
                                boost::lexical_cast< std::string >( independentVariable ) + " but limit values are " +
                                boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).front( ) ) + " and " +
                                boost::lexical_cast< std::string >( independentValues_.at( currentVariable ).back( ) ) + ", taking default value instead.";
                        std::cerr << errorMessage << std::endl;
                    }

                    useDefault = true;
                }
                else
                {
                    throw std::runtime_error( "Error when checking interpolation boundary, boundary handling method not recognized." );
                }
            }
        }
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
        // Find which type of scheme is used.
        lookUpSchemes_.resize( NumberOfDimensions );
        switch ( selectedScheme )
        {
        case binarySearch:
        {
            for ( int i = 0; i < NumberOfDimensions; i++ )
            {
                // Create binary search look up scheme.
                lookUpSchemes_[ i ] = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                        ( new BinarySearchLookupScheme< IndependentVariableType >(
                              independentValues_[ i ] ) );
            }
            break;
        }
        case huntingAlgorithm:
        {
            for ( int i = 0; i < NumberOfDimensions; i++ )
            {
                // Create hunting scheme, which uses an intial guess from previous look-ups.
                lookUpSchemes_[ i ] = boost::shared_ptr< LookUpScheme< IndependentVariableType > >
                        ( new HuntingAlgorithmLookupScheme< IndependentVariableType >(
                              independentValues_[ i ] ) );
            }
            break;
        }
        default:
            throw std::runtime_error( "Error: lookup scheme not found when making scheme for N-D interpolator." );
        }
    }

    //! Vector with pointers to look-up scheme.
    /*!
     * Pointers to the look-up schemes that is used to determine in which interval the requested
     * independent variable value falls.
     */
    std::vector< boost::shared_ptr< LookUpScheme< IndependentVariableType > > > lookUpSchemes_;

    //! Vector of vectors containing independent variables.
    /*!
     * Vector of vectors containing independent variables. The size of the outer vector is equal
     * to the number of dimensions of the interpolator.
     */
    std::vector< std::vector< IndependentVariableType > > independentValues_;

    //! Multi-dimensional array of dependent data.
    /*!
     * Multi-dimensional array of dependent data at each point of hyper-rectangular grid formed by
     * independent variable points.
     */
    boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) > dependentData_;

    //! Behavior of interpolator when independent variable is outside range.
    /*!
     * Behavior of interpolator when independent variable is outside range.
     */
    std::vector< BoundaryInterpolationType > boundaryHandling_;

    //! Default value to be used for extrapolation.
    /*!
     * Default value to be used for extrapolation.
     */
    DependentVariableType defaultExtrapolationValue_;

};

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H
