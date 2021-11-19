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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#ifndef TUDAT_CUBIC_SPLINE_INTERPOLATOR_H
#define TUDAT_CUBIC_SPLINE_INTERPOLATOR_H

#include <cmath>
#include <Eigen/Core>

#include <memory>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/math/basic/nearestNeighbourSearch.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace interpolators
{

//! Solve efficiently tri-diagonal matrix equation.
/*!
 *  This functions efficiently solves the matrix equation Ax = b for b, where A is tri-diagonal.
 *  The required input are the diagonal and sub/super diagonals of A, as well as the vector b.
 *  Input is given as STL vectors. The input diagonals must be at least the same size as
 *  right-hand-side. Any additional entries (at end) will be ignored in the algorithm.
 *  \tparam IndependentVariableType Type of independent variables.
 *  \tparam DependentVariableType Type of dependent variables.
 *  \param subDiagonal Sub-diagonal of matrix A.
 *  \param superDiagonal Super-diagonal of matrix A.
 *  \param diagonal Diagonal of matrix A.
 *  \param rightHandSide Right-hand-side of matrix equation.
 *  \return Solution to matrix equation.
 */
template< typename IndependentVariableType, typename DependentVariableType >
std::vector< DependentVariableType > solveTridiagonalMatrixEquation(
        const std::vector< IndependentVariableType >& subDiagonal,
        const std::vector< IndependentVariableType >& diagonal,
        const std::vector< IndependentVariableType >& superDiagonal,
        const std::vector< DependentVariableType >& rightHandSide )
{
    // Check whether input diagonals are correct size.
    unsigned int matrixSize = rightHandSide.size( );
    if ( ( diagonal.size( ) < matrixSize ) || ( superDiagonal.size( ) < matrixSize - 1 ) ||
         ( rightHandSide.size( ) < matrixSize - 1 ) )
    {
        throw std::runtime_error( "Error, input provided for diagonal and sub/super diagonals incorrect." );
    }

    // Check whether solution will not be singular.
    if ( diagonal[ 0 ] == 0.0 )
    {
        throw std::runtime_error( "Error when inverting tridiagonal system, first entry of diagonal is zero" );
    }

    std::vector< IndependentVariableType > intermediateVector( matrixSize );
    std::vector< DependentVariableType > solution( matrixSize );

    // Perform solution algorithm, from (Press W.H., et al., 2002).
    double scalingFactor = diagonal[ 0 ];
    solution[ 0 ] = rightHandSide[ 0 ] / scalingFactor;

    for ( unsigned int j = 1; j < matrixSize; j++ )
    {
        intermediateVector[ j ] = superDiagonal[ j - 1 ] / scalingFactor;
        scalingFactor = diagonal[ j ] - subDiagonal[ j - 1 ] * intermediateVector[ j ];

        // Check whether solution will not be singular.
        if ( scalingFactor == 0.0 )
        {
            throw std::runtime_error( "Error when inverting tridiagonal system, scaling factor equals zero!" );
        }
        solution[ j ] = ( rightHandSide[ j ] - subDiagonal[ j - 1 ] * solution[ j - 1 ] ) /
                scalingFactor;
    }

    for ( int j = ( matrixSize - 2 ); j >= 0 ; j-- )
    {
        solution[ j ] -= intermediateVector[ j + 1 ] * solution[ j + 1 ];
    }

    return solution;
}

//! Cubic spline interpolator, implementation from (Press W.H., et al., 2002).
/*!
 *  Cubic spline interpolator, implementation from (Press W.H., et al., 2002).
 *  Natural boundary conditions are imposed, meaning zero second derivatives (curvature) at end
 *  points. Continuity of first derivatives is imposed.
 *  \tparam IndependentVariableType Type of independent variables.
 *  \tparam DependentVariableType Type of dependent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType, typename ScalarType = IndependentVariableType >
class CubicSplineInterpolator :
        public OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >
{
public:

    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    dependentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    independentValues_;
    using OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >::
    lookUpScheme_;
    using Interpolator< IndependentVariableType, DependentVariableType >::interpolate;

    //! Cubic spline interpolator constructor.
    /*!
     *  Cubic spline interpolator constructor taking separate vectors of dependent and independent
     *  variable values.
     *  \param independentVariables Vector with the independent variable values, must be
     *      sorted in ascending order.
     *  \param dependentVariables Vector with the dependent variable values.
     *  \param selectedLookupScheme Look-up scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    CubicSplineInterpolator( const std::vector< IndependentVariableType >& independentVariables,
                             const std::vector< DependentVariableType >& dependentVariables,
                             const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                             const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary,
                             const std::pair< DependentVariableType, DependentVariableType >& defaultExtrapolationValue =
            std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                            IdentityElement::getAdditionIdentity< DependentVariableType >( ) ) ) :
        OneDimensionalInterpolator< IndependentVariableType, DependentVariableType >( boundaryHandling, defaultExtrapolationValue )
    {
        // Verify that the initialization variables are not empty.
        if ( independentVariables.size( ) == 0 || dependentVariables.size( ) == 0 )
        {
            throw std::runtime_error(
                        "The vectors used in the cubic spline interpolator initialization are empty." );
        }

        // Set dependent and independent variable values.
        independentValues_ = std::move( independentVariables );
        dependentValues_ = std::move( dependentVariables );

        // Check if data is in ascending order
        if( !std::is_sorted( independentVariables.begin( ), independentVariables.end( ) ) )
        {
            throw std::runtime_error( "Error when making cubic spline interpolator, input vector with independent variables should be in ascending order" );
        }

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );

        // Create zero value for initializing output.
        zeroValue_ = dependentVariables[ 0 ] - dependentVariables[ 0 ];

        if ( dependentValues_.size( ) != independentValues_.size( ) )
        {
            throw std::runtime_error( "Warning: independent and dependent variables not of same size in cubic spline constrcutor" );
        }

        // Calculate second derivatives of curve.
        calculateSecondDerivatives( );
    }

    //! Cubic spline interpolator constructor.
    /*!
     *  Cubic spline interpolator constructor taking separate vectors of dependent and independent
     *  variable values. This constructor takes a single default value, instead of a pair.
     *  \param independentVariables Vector with the independent variable values, must be
     *      sorted in ascending order.
     *  \param dependentVariables Vector with the dependent variable values.
     *  \param selectedLookupScheme Look-up scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    CubicSplineInterpolator( const std::vector< IndependentVariableType >& independentVariables,
                             const std::vector< DependentVariableType >& dependentVariables,
                             const AvailableLookupScheme selectedLookupScheme,
                             const BoundaryInterpolationType boundaryHandling,
                             const DependentVariableType& defaultExtrapolationValue ):
        CubicSplineInterpolator( independentVariables, dependentVariables, selectedLookupScheme, boundaryHandling,
                                 std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) )
    { }

    //! Cubic spline interpolator constructor.
    /*!
     *  Cubic spline interpolator constructor taking single map of independent and dependent
     *  variable values.
     *  \param dataMap Map with the independent variable values as keys and corresponding
     *      dependent variable values as values.
     *  \param selectedLookupScheme Lookup scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Pair of default values to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    CubicSplineInterpolator(
            const std::map< IndependentVariableType, DependentVariableType > dataMap,
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
            throw std::runtime_error( "The map used in the cubic spline interpolator initialization are empty." );
        }

        // Set data vector member variables from map.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );
        int counter = 0;
        for ( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator
              mapIterator = dataMap.begin( ); mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = std::move( mapIterator->first );
            dependentValues_[ counter ] = std::move( mapIterator->second );
            counter++;
        }

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );

        // Create zero value for initializing output.
        zeroValue_ = dependentValues_[ 0 ] - dependentValues_[ 0 ];

        // Calculate second derivatives of curve.
        calculateSecondDerivatives( );
    }

    //! Cubic spline interpolator constructor.
    /*!
     *  Cubic spline interpolator constructor taking single map of independent and dependent
     *  variable values. This constructor takes a single default value, instead of a pair.
     *  \param dataMap Map with the independent variable values as keys and corresponding
     *      dependent variable values as values.
     *  \param selectedLookupScheme Lookup scheme that is to be used when finding interval
     *      of requested independent variable value.
     *  \param boundaryHandling Boundary handling method, in case the independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case
     *      of use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    CubicSplineInterpolator(
            const std::map< IndependentVariableType, DependentVariableType > dataMap,
            const AvailableLookupScheme selectedLookupScheme,
            const BoundaryInterpolationType boundaryHandling,
            const DependentVariableType& defaultExtrapolationValue ):
        CubicSplineInterpolator( dataMap, selectedLookupScheme, boundaryHandling,
                                 std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) )
    { }

    //! Destructor.
    ~CubicSplineInterpolator( ){ }

    //! Interpolate.
    /*!
     *  Executes interpolation of data at a given target value of the independent variable, to
     *  yield an interpolated value of the dependent variable.
     *  \param targetIndependentVariableValue Target independent variable value at which point
     *      the interpolation is performed.
     *  \return Interpolated dependent variable value.
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

        // Determine the lower entry in the table corresponding to the target independent variable
        // value.
        int lowerEntry_ = lookUpScheme_->findNearestLowerNeighbour(
                    targetIndependentVariableValue );

        // Get independent variable values bounding interval in which requested value lies.
        IndependentVariableType lowerValue, upperValue;
        ScalarType squareDifference;
        lowerValue = independentValues_[ lowerEntry_ ];
        upperValue = independentValues_[ lowerEntry_ + 1 ];

        // Calculate coefficients A,B,C,D (see Numerical (Press W.H., et al., 2002))
        squareDifference = static_cast< ScalarType >( upperValue - lowerValue ) *
                static_cast< ScalarType >( upperValue - lowerValue );
        ScalarType coefficientA_ = ( upperValue - targetIndependentVariableValue )
                / static_cast< ScalarType >( upperValue - lowerValue );
        ScalarType coefficientB_ = mathematical_constants::getFloatingInteger< ScalarType >( 1.0 ) - coefficientA_;
        ScalarType coefficientC_ = ( coefficientA_ * coefficientA_ * coefficientA_ - coefficientA_ ) /
                mathematical_constants::getFloatingInteger< ScalarType >( 6.0 ) * squareDifference;
        ScalarType coefficientD_ = ( coefficientB_ * coefficientB_ * coefficientB_ - coefficientB_ ) /
                mathematical_constants::getFloatingInteger< ScalarType >( 6.0 ) * squareDifference;

        // The interpolated dependent variable value.
        return coefficientA_ * dependentValues_[ lowerEntry_ ] +
                coefficientB_ * dependentValues_[ lowerEntry_ + 1 ] +
                coefficientC_ * secondDerivativeOfCurve_[ lowerEntry_ ] +
                coefficientD_ * secondDerivativeOfCurve_[ lowerEntry_ + 1 ];
    }

    InterpolatorTypes getInterpolatorType( ){ return cubic_spline_interpolator; }

protected:

private:

    //! Calculates the second derivatives of the curve.
    /*!
     *  This function calculates the second derivatives of the curve at the nodes, assuming
     *  the first derivatives to be continuous at the nodes and imposing natural spline conditions
     *  (zero curvature at endpoints). The methodology is described in (Press W.H., et al., 2002).
     */
    void calculateSecondDerivatives( )
    {
        // Get length of vector.
        numberOfDataPoints_ = independentValues_.size( );

        // Sub-diagonal of tri-diagonal matrix.
        std::vector< ScalarType > aCoefficients_;
        aCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Diagonal of tri-diagonal matrix.
        std::vector< ScalarType > bCoefficients_;
        bCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Super-diagonal of tri-diagonal matrix.
        std::vector< ScalarType > cCoefficients_;
        cCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Right-hand side of tridiagonal matrix system
        std::vector< DependentVariableType > rCoefficients_;
        rCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Temporary value vector.
        std::vector< ScalarType > hCoefficients_;
        hCoefficients_.resize( numberOfDataPoints_ - 1 );

        // Set second derivatives of curve to zero at endpoints, i.e. impose natural spline
        // condition.
        aCoefficients_[ numberOfDataPoints_- 3 ] = mathematical_constants::getFloatingInteger< ScalarType >( 0 );
        cCoefficients_[ numberOfDataPoints_- 3 ] = mathematical_constants::getFloatingInteger< ScalarType >( 0 );

        // Compute the vectors h (temporary values),a,c,b,r.
        for ( unsigned int i = 0; i < ( numberOfDataPoints_ - 1 ); i++ )
        {
            hCoefficients_[ i ] = independentValues_[ i + 1 ] - independentValues_[ i ];
        }

        // Set tridiagonal matrix equation input.
        for ( unsigned int i = 0; i < ( numberOfDataPoints_ - 3 ); i++ )
        {
            aCoefficients_[ i ] = hCoefficients_[ i + 1 ];
            cCoefficients_[ i ] = hCoefficients_[ i + 1 ];
        }

        for ( unsigned int i = 0; i < ( numberOfDataPoints_ - 2 ); i++ )
        {
            bCoefficients_[ i ] = 2.0 * ( hCoefficients_[ i + 1 ] + hCoefficients_[ i ] );
            rCoefficients_[ i ] = 6.0 * ( ( dependentValues_[ i + 2 ]-
                                          dependentValues_[ i + 1 ] ) / hCoefficients_[ i + 1 ] -
                    ( dependentValues_[ i + 1 ] - dependentValues_[ i ] ) /
                    hCoefficients_[ i ] );
        }

        // Solve tridiagonal matrix equatuion.
        std::vector< DependentVariableType > middleSecondDerivativeOfCurvatures =
                solveTridiagonalMatrixEquation< ScalarType, DependentVariableType >
                ( aCoefficients_, bCoefficients_, cCoefficients_, rCoefficients_ );

        // Append zeros to ends of calculated second derivative values (natural spline condition).
        secondDerivativeOfCurve_.resize( numberOfDataPoints_ );
        secondDerivativeOfCurve_[ 0 ] = zeroValue_;

        for ( unsigned int i = 1; i < numberOfDataPoints_ - 1; i++ )
        {
            secondDerivativeOfCurve_[ i ] = middleSecondDerivativeOfCurvatures[ i - 1 ];
        }

        secondDerivativeOfCurve_[ numberOfDataPoints_ - 1 ] = zeroValue_;
    }

    //! Vector filled with second derivative of curvature of each point.
    std::vector< DependentVariableType > secondDerivativeOfCurve_;

    //! The number of datapoints.
    unsigned int numberOfDataPoints_;

    //! Zero value of independent variable type
    DependentVariableType zeroValue_;

};


extern template class CubicSplineInterpolator< double, Eigen::VectorXd >;
extern template class CubicSplineInterpolator< double, Eigen::Vector6d >;
extern template class CubicSplineInterpolator< double, Eigen::MatrixXd >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class CubicSplineInterpolator< Time, Eigen::VectorXd, long double >;
extern template class CubicSplineInterpolator< Time, Eigen::Vector6d, long double >;
extern template class CubicSplineInterpolator< Time, Eigen::MatrixXd, long double >;

extern template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
extern template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 6 > >;
extern template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic > >;

extern template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >;
extern template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 6 >, long double >;
extern template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic >, long double >;
#endif

//! Typedef for cubic spline interpolator with (in)dependent = double.
typedef CubicSplineInterpolator< double, double > CubicSplineInterpolatorDouble;

//! Typedef for shared-pointer to cubic spline interpolator with (in)dependent = double.
typedef std::shared_ptr< CubicSplineInterpolatorDouble > CubicSplineInterpolatorDoublePointer;

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CUBIC_SPLINE_INTERPOLATOR_H
