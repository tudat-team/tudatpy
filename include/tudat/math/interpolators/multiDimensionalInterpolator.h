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

#ifndef TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H
#define TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H

#include <vector>
#include <iostream>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/interpolator.h"

#include "tudat/basics/identityElements.h"
#include "tudat/basics/tudatExceptions.h"

namespace tudat
{

namespace interpolators
{

//! Base class for interpolator with multiple independent variables.
/*!
 *  Base class for the interpolators in multiple independent variables included in Tudat.
 *  \tparam IndependentVariableType Type of independent variable(s)
 *  \tparam IndependentVariableType Type of dependent variable
 *  \tparam NumberOfDimensions Number of independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType, unsigned int NumberOfDimensions >
class MultiDimensionalInterpolator : public Interpolator< IndependentVariableType, DependentVariableType >
{
public:
    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Vector of pairs of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    MultiDimensionalInterpolator(
            const std::vector< BoundaryInterpolationType >& boundaryHandling =
                    std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary ),
            const std::vector< std::pair< DependentVariableType, DependentVariableType > >& defaultExtrapolationValue =
                    std::vector< std::pair< DependentVariableType, DependentVariableType > >(
                            NumberOfDimensions,
                            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ) ):
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    {
        // Check that the user-defined default value does not correspond to a boundary method that does not use the
        // default value
        std::pair< DependentVariableType, DependentVariableType > defaultValueForBoundaryHandling =
                std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                                IdentityElement::getAdditionIdentity< DependentVariableType >( ) );
        for( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            if( defaultExtrapolationValue_.at( i ) != defaultValueForBoundaryHandling )
            {
                if( !( ( boundaryHandling_.at( i ) == use_default_value ) ||
                       ( boundaryHandling_.at( i ) == use_default_value_with_warning ) ) )
                {
                    std::cerr << "Warning in multi-dimensional interpolator. A default value has been set (for when the "
                                 "independent variable is out-of-range) but the boundary handling method is not "
                                 "use_default_value or use_default_value_with_warning. The method selected for dimension "
                              << std::to_string( i ) << " is: " << std::to_string( boundaryHandling_.at( i ) ) << std::endl;
                }
            }
        }
    }

    //! Constructor taking single pair of default values.
    /*!
     *  Constructor taking single pair of default values.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    MultiDimensionalInterpolator( const std::vector< BoundaryInterpolationType >& boundaryHandling,
                                  const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue ):
        MultiDimensionalInterpolator(
                boundaryHandling,
                std::vector< std::pair< DependentVariableType, DependentVariableType > >( NumberOfDimensions, defaultExtrapolationValue ) )
    { }

    //! Constructor taking single default value.
    /*!
     *  Constructor taking single default value.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    MultiDimensionalInterpolator( const std::vector< BoundaryInterpolationType >& boundaryHandling,
                                  const DependentVariableType& defaultExtrapolationValue ):
        MultiDimensionalInterpolator( boundaryHandling,
                                      std::vector< std::pair< DependentVariableType, DependentVariableType > >(
                                              NumberOfDimensions,
                                              std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) )
    { }

    //! Destructor.
    /*!
     *  Destructor.
     */
    virtual ~MultiDimensionalInterpolator( ) { }

    //! Function to perform interpolation.
    /*!
     *  This function performs the interpolation. It calls the function that takes a vector of
     *  independent variable values, which is to be implemented in derived classes.
     *  \param independentVariableValues Vector of values of independent variables at which
     *          the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType interpolate( const std::vector< IndependentVariableType >& independentVariableValues ) = 0;

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
    std::vector< std::shared_ptr< LookUpScheme< IndependentVariableType > > > getLookUpScheme( )
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
     *  \param currentIndependentVariable Value of current independent variable.
     *  \param currentDimension Value of current dimension.
     *  \return Condition with respect to boundary.
     */
    int checkInterpolationBoundary( const IndependentVariableType& currentIndependentVariable, const unsigned int& currentDimension )
    {
        int isAtBoundary = 0;
        if( currentIndependentVariable < independentValues_.at( currentDimension ).front( ) )
        {
            isAtBoundary = -1;
        }
        else if( currentIndependentVariable > independentValues_.at( currentDimension ).back( ) )
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
     *  \param currentDimension Value of current dimension.
     *  \param useValue Boolean denoting whether the value given by this function needs to be used.
     *  \param currentIndependentVariable Value of current independent variable.
     *  \param dependentVariable Value of current dependent variable, in case the independent variable is out-of-range
     *      and the selected method is use_default_value or use_default_value_with_warning.
     */
    void checkBoundaryCase( const unsigned int currentDimension,
                            bool& useValue,
                            IndependentVariableType& currentIndependentVariable,
                            DependentVariableType& dependentVariable )
    {
        // If extrapolation outside domain is not allowed
        if( boundaryHandling_.at( currentDimension ) != extrapolate_at_boundary )
        {
            // If independent variable is out of range
            int isAtBoundary = checkInterpolationBoundary( currentIndependentVariable, currentDimension );
            if( isAtBoundary != 0 )
            {
                // Select course of action based on boundary handling method selected
                switch( boundaryHandling_.at( currentDimension ) )
                {
                    case throw_exception_at_boundary: {
                        throw tudat::exceptions::InterpolationOutOfBoundsError< IndependentVariableType >(
                                currentIndependentVariable,
                                independentValues_.at( currentDimension ).front( ),
                                independentValues_.at( currentDimension ).back( ),
                                currentDimension );
                        break;
                    }
                    case extrapolate_at_boundary_with_warning: {
                        // Warn user
                        std::string errorMessage =
                                "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                std::to_string( currentDimension ) + " at " + std::to_string( currentIndependentVariable ) +
                                " but limit values are " + std::to_string( independentValues_.at( currentDimension ).front( ) ) + " and " +
                                std::to_string( independentValues_.at( currentDimension ).back( ) ) + ", applying extrapolation instead.";
                        std::cerr << errorMessage << std::endl;
                        break;
                    }
                    case use_boundary_value:
                    case use_boundary_value_with_warning: {
                        // Warn user, if requested
                        if( boundaryHandling_.at( currentDimension ) == use_boundary_value_with_warning )
                        {
                            std::string errorMessage =
                                    "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                    std::to_string( currentDimension ) + " at " + std::to_string( currentIndependentVariable ) +
                                    " but limit values are " + std::to_string( independentValues_.at( currentDimension ).front( ) ) +
                                    " and " + std::to_string( independentValues_.at( currentDimension ).back( ) ) +
                                    ", taking boundary value instead.";
                            std::cerr << errorMessage << std::endl;
                        }

                        // Replace current independent variable with boundary value
                        if( isAtBoundary == -1 )
                        {
                            currentIndependentVariable = independentValues_.at( currentDimension ).front( );
                        }
                        else if( isAtBoundary == 1 )
                        {
                            currentIndependentVariable = independentValues_.at( currentDimension ).back( );
                        }
                        break;
                    }
                    case use_default_value:
                    case use_default_value_with_warning: {
                        // Warn user, if requested
                        if( boundaryHandling_.at( currentDimension ) == use_default_value_with_warning )
                        {
                            std::string errorMessage =
                                    "Warning in interpolator, requesting data point outside of boundaries, requested data of dimension " +
                                    std::to_string( currentDimension ) + " at " + std::to_string( currentIndependentVariable ) +
                                    " but limit values are " + std::to_string( independentValues_.at( currentDimension ).front( ) ) +
                                    " and " + std::to_string( independentValues_.at( currentDimension ).back( ) ) +
                                    ", taking default value instead.";
                            std::cerr << errorMessage << std::endl;
                        }

                        // Get default value
                        useValue = true;
                        if( isAtBoundary == -1 )
                        {
                            dependentVariable = defaultExtrapolationValue_.at( currentDimension ).first;
                        }
                        else if( isAtBoundary == 1 )
                        {
                            dependentVariable = defaultExtrapolationValue_.at( currentDimension ).second;
                        }
                        else
                        {
                            throw std::runtime_error( "Error when checking interpolation boundary, inconsistent data encountered" );
                        }
                        break;
                    }
                    default:
                        throw std::runtime_error( "Error when checking interpolation boundary, boundary handling method not recognized." );
                }
            }
        }
    }

    //! Make the lookup scheme that is to be used.
    /*!
     *  This function creates the look up scheme that is to be used in determining the interval of
     *  the independent variable grid where the interpolation is to be performed. It takes the type
     *  of lookup scheme as an enum and constructs the lookup scheme from the independentValues_
     *  that have been set previously.
     *  \param selectedScheme Type of look-up scheme that is to be used
     */
    void makeLookupSchemes( const AvailableLookupScheme selectedScheme )
    {
        // Find which type of scheme is used.
        lookUpSchemes_.resize( NumberOfDimensions );
        switch( selectedScheme )
        {
            case binarySearch: {
                for( unsigned int i = 0; i < NumberOfDimensions; i++ )
                {
                    // Create binary search look up scheme.
                    lookUpSchemes_[ i ] = std::shared_ptr< LookUpScheme< IndependentVariableType > >(
                            new BinarySearchLookupScheme< IndependentVariableType >( independentValues_[ i ] ) );
                }
                break;
            }
            case huntingAlgorithm: {
                for( unsigned int i = 0; i < NumberOfDimensions; i++ )
                {
                    // Create hunting scheme, which uses an intial guess from previous look-ups.
                    lookUpSchemes_[ i ] = std::shared_ptr< LookUpScheme< IndependentVariableType > >(
                            new HuntingAlgorithmLookupScheme< IndependentVariableType >( independentValues_[ i ] ) );
                }
                break;
            }
            default:
                throw std::runtime_error( "Error: lookup scheme not found when making scheme for N-D interpolator." );
        }
    }

    //! Vector with pointers to look-up scheme.
    /*!
     *  Pointers to the look-up schemes that is used to determine in which interval the requested
     *  independent variable value falls.
     */
    std::vector< std::shared_ptr< LookUpScheme< IndependentVariableType > > > lookUpSchemes_;

    //! Vector of vectors containing independent variables.
    /*!
     *  Vector of vectors containing independent variables. The size of the outer vector is equal
     *  to the number of dimensions of the interpolator.
     */
    std::vector< std::vector< IndependentVariableType > > independentValues_;

    //! Multi-dimensional array of dependent data.
    /*!
     *  Multi-dimensional array of dependent data at each point of hyper-rectangular grid formed by
     *  independent variable points.
     */
    boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) > dependentData_;

    //! Behavior of interpolator when independent variable is outside range.
    /*!
     *  Behavior of interpolator when independent variable is outside range.
     */
    std::vector< BoundaryInterpolationType > boundaryHandling_;

    //! Default value to be used for extrapolation.
    /*!
     *  Default value to be used for extrapolation.
     */
    std::vector< std::pair< DependentVariableType, DependentVariableType > > defaultExtrapolationValue_;
};

}  // namespace interpolators

}  // namespace tudat

#endif  // TUDAT_MULTI_DIMENSIONAL_INTERPOLATOR_H
