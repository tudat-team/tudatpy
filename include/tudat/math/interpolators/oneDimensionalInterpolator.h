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

#ifndef TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H
#define TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H

#include <vector>
#include <iostream>

#include <boost/lexical_cast.hpp>

#include <memory>

#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/interpolator.h"

#include "tudat/basics/identityElements.h"

namespace tudat
{

namespace interpolators
{

//! Enumeration of available interpolator types.
enum InterpolatorTypes
{
    linear_interpolator = 0,
    multi_linear_interpolator = 1,
    cubic_spline_interpolator = 2,
    lagrange_interpolator = 3,
    hermite_spline_interpolator = 4,
    piecewise_constant_interpolator = 5,
    discrete_jump_linear_interpolator = 6

};


//! Base class for interpolator with one independent variable.
/*!
 * Base class for the interpolators in one independent variable included in Tudat.
 * \tparam IndependentVariableType Type of independent variable(s)
 * \tparam IndependentVariableType Type of dependent variable
 */
template< typename IndependentVariableType, typename DependentVariableType >
class OneDimensionalInterpolator : public Interpolator< IndependentVariableType, DependentVariableType >
{
public:

    // Using statements to prevent having to put 'this' everywhere in the code.
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Constructor.
    /*!
     *  Constructor.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    OneDimensionalInterpolator(
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
            const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ):
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    { }

    //! Constructor.
    /*!throw_exception_at_boundary
     *  Constructor.
     *  \param boundaryHandling Boundary handling method in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    OneDimensionalInterpolator(
            const BoundaryInterpolationType boundaryHandling,
            const DependentVariableType& defaultExtrapolationValue ):
        OneDimensionalInterpolator( boundaryHandling, std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) )
    { }

    //! Destructor.
    /*!
     * Destructor.
     */
    virtual ~OneDimensionalInterpolator( ) { }

    //! Function to perform interpolation.
    /*!
     *  This function performs the interpolation. It calls the function that takes a single
     *  independent variable value, which is to be implemented in derived classes.
     *  \param independentVariableValues Vector of values of independent variables at which
     *      the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType
    interpolate( const std::vector< IndependentVariableType >& independentVariableValues )
    {
        // Check whether input is really 1-dimensional
        if ( independentVariableValues.size( ) != 1  )
        {
            throw std::runtime_error( "Error in 1-dimensional interpolator, provided input is not 1-dimensional." );
        }

        // Call 1-dimensional interpolate function.
        return interpolate( independentVariableValues[ 0 ] );
    }

    //! Function to perform interpolation.
    /*!
     *  This function performs the interpolation
     *  \param independentVariableValue Independent variable value at which the value of the
     *      dependent variable is to be determined.
     *  \return Interpolated value of dependent variable.
     */
    virtual DependentVariableType
    interpolate( const IndependentVariableType independentVariableValue ) = 0;

    //! Function to perform interpolation, with non-const input argument.
    /*!
     *  This function performs the interpolation, with non-const input argument. Function calls the interpolate function and is
     *  included for compatibility with some function pointer binding interfaces.
     *  \param independentVariableValue Independent variable value at which the value of the
     *          dependent variable is to be determined.
     *  \return Interpolated value of dependent variable.
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
    std::shared_ptr< LookUpScheme< IndependentVariableType > > getLookUpScheme( )
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

    BoundaryInterpolationType getBoundaryHandling( )
    {
        return boundaryHandling_;
    }

    void resetBoundaryHandling( const BoundaryInterpolationType boundaryHandling )
    {
        boundaryHandling_ = boundaryHandling;
    }

    std::pair< DependentVariableType, DependentVariableType > getDefaultExtrapolationValue( )
    {
        return defaultExtrapolationValue_;
    }

    AvailableLookupScheme getSelectedLookupScheme( )
    {
        return selectedLookupScheme_;
    }

    virtual InterpolatorTypes getInterpolatorType( ) = 0;

protected:

    //! Function to return the condition of the current independent variable.
    /*!
     *  Function to return the condition of the current independent variable, i.e. whether the
     *  variable is within, above or below its defined range range.
     *  \param targetIndependentVariable Value of independent variable (i.e., the one that is to be checked for boundary handling).
     *  \return Condition with respect to boundary.
     */
    int checkInterpolationBoundary( const IndependentVariableType& targetIndependentVariable )
    {
        int isAtBoundary = 0;
        if ( targetIndependentVariable < independentValues_.front( ) )
        {
            isAtBoundary = -1;
        }
        else if ( targetIndependentVariable > independentValues_.back( ) )
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
     *  \param dependentVariable Value of dependent variable at boundary (only used in case of use_boundary_value setting).
     *  \param useValue Boolean that specifies whether the boundary value (i.e., dependentVariable) is to be used, instead of interpolating.
     *  \param targetIndependentVariable Value of independent variable (i.e., the one that is to be checked for boundary handling).
     */
    void checkBoundaryCase(
            DependentVariableType& dependentVariable, bool& useValue,
            const IndependentVariableType& targetIndependentVariable )
    {
        // If extrapolation outside domain is not allowed
        if ( boundaryHandling_ != extrapolate_at_boundary )
        {
            // If independent variable is out of range
            int isAtBoundary = checkInterpolationBoundary( targetIndependentVariable );
            if ( isAtBoundary != 0 )
            {
                // Select course of action based on boundary handling method selected
                switch ( boundaryHandling_ )
                {
                case throw_exception_at_boundary:
                {
                    // Throw exception
                    std::string errorMessage = "Error in interpolator, requesting data point outside of boundaries, requested data at " +
                            boost::lexical_cast< std::string >( targetIndependentVariable ) + " but limit values are " +
                            boost::lexical_cast< std::string >( independentValues_.front( ) ) + " and " +
                            boost::lexical_cast< std::string >( independentValues_.back( ) );
                    throw std::runtime_error( errorMessage );
                    break;
                }
                case extrapolate_at_boundary_with_warning:
                {
                    // Warn user
                    std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data at " +
                            boost::lexical_cast< std::string >( targetIndependentVariable ) + " but limit values are " +
                            boost::lexical_cast< std::string >( independentValues_.front( ) ) + " and " +
                            boost::lexical_cast< std::string >( independentValues_.back( ) ) + ", applying extrapolation instead.";
                    std::cerr << errorMessage << std::endl;
                    break;
                }
                case use_boundary_value:
                case use_boundary_value_with_warning:
                {
                    // Warn user, if requested
                    if ( boundaryHandling_ == use_boundary_value_with_warning )
                    {
                        std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data at " +
                                boost::lexical_cast< std::string >( targetIndependentVariable ) + " but limit values are " +
                                boost::lexical_cast< std::string >( independentValues_.front( ) ) + " and " +
                                boost::lexical_cast< std::string >( independentValues_.back( ) ) + ", taking boundary value instead.";
                        std::cerr << errorMessage << std::endl;
                    }

                    // Get boundary value
                    useValue = true;
                    if ( isAtBoundary == -1 )
                    {
                        dependentVariable = dependentValues_.front( );
                    }
                    else if ( isAtBoundary == 1 )
                    {
                        dependentVariable = dependentValues_.back( );
                    }
                    else
                    {
                        throw std::runtime_error( "Error when checking interpolation boundary, inconsistent data encountered" );
                    }
                    break;
                }
                case use_default_value:
                case use_default_value_with_warning:
                {
                    // Warn user, if requested
                    if ( boundaryHandling_ == use_default_value_with_warning )
                    {
                        std::string errorMessage = "Warning in interpolator, requesting data point outside of boundaries, requested data at " +
                                boost::lexical_cast< std::string >( targetIndependentVariable ) + " but limit values are " +
                                boost::lexical_cast< std::string >( independentValues_.front( ) ) + " and " +
                                boost::lexical_cast< std::string >( independentValues_.back( ) ) + ", taking default value instead.";
                        std::cerr << errorMessage << std::endl;
                    }

                    // Get default value
                    useValue = true;
                    if ( isAtBoundary == -1 )
                    {
                        dependentVariable = defaultExtrapolationValue_.first;
                    }
                    else if ( isAtBoundary == 1 )
                    {
                        dependentVariable = defaultExtrapolationValue_.second;
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

    //! Make look-up scheme that is to be used.
    /*!
     * This function creates the look-up scheme that is to be used in determining the interval of
     * the independent variable grid where the interpolation is to be performed. It takes the type
     * of lookup scheme as an enum and constructs the look-up scheme from the independentValues_
     * that have been set previously.
     * \param selectedScheme Type of look-up scheme that is to be used
     */
    void makeLookupScheme( const AvailableLookupScheme selectedScheme )
    {
        selectedLookupScheme_ = selectedScheme;

        // Find which type of scheme is used.
        switch ( selectedLookupScheme_ )
        {
        case binarySearch:
        {
            // Create binary search look up scheme.
            lookUpScheme_ = std::make_shared< BinarySearchLookupScheme< IndependentVariableType > >
                      ( independentValues_ );
            break;
        }
        case huntingAlgorithm:
        {
            // Create hunting scheme, which uses an intial guess from previous look-ups.
            lookUpScheme_ = std::make_shared< HuntingAlgorithmLookupScheme< IndependentVariableType > >
                      ( independentValues_ );
            break;
        }
        default:
            throw std::runtime_error( "Warning: lookup scheme not found when making scheme for 1-D interpolator" );
        }
    }

    //! Pointer to look up scheme.
    /*!
     * Pointer to the lookup scheme that is used to determine in which interval the requested
     * independent variable value falls.
     */
    std::shared_ptr< LookUpScheme< IndependentVariableType > > lookUpScheme_;

    //! Vector with independent variables.
    /*!
     * Vector with independent variables.
     */
    std::vector< IndependentVariableType > independentValues_;

    //! Vector with dependent variables.
    /*!
     * Vector with dependent variables.
     */
    std::vector< DependentVariableType > dependentValues_;

    //! Behavior of interpolator when independent variable is outside range.
    /*!
     * Behavior of interpolator when independent variable is outside range.
     */
    BoundaryInterpolationType boundaryHandling_;

    AvailableLookupScheme selectedLookupScheme_ = undefinedScheme;

    //! Default value to be used for extrapolation.
    /*!
     * Default value to be used for extrapolation.
     */
    std::pair< DependentVariableType, DependentVariableType > defaultExtrapolationValue_;

};

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_ONE_DIMENSIONAL_INTERPOLATOR_H
