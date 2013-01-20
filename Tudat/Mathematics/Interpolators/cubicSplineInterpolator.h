/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110620    F.M. Engelen      File created.
 *      110707    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110714    E.A.G. Heeren     Minor spelling/lay-out corrections.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *      130114    D. Dirkx          Fixed iterator bug.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 *    Notes
 *
 */

#ifndef TUDAT_CUBIC_SPLINE_INTERPOLATOR_H
#define TUDAT_CUBIC_SPLINE_INTERPOLATOR_H

#include <cmath>
#include <iostream>

#include <Eigen/Core>

#include <boost/exception/all.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

namespace tudat
{
namespace interpolators
{

//! Solve efficiently tri-diagonal matrix equation.
/*!
 * This functions efficiently solves the matrix equation Ax = b for b, where A is tri-diagonal.
 * The required input are the diagonal and sub/super diagonals of A, as well as the vector b.
 * Input is given as STL vectors. The input diagonals must be at least the same size as
 * right-hand-side. Any additional entries (at end) will be ignored in the algorithm.
 * \tparam IndependentVariableType Type of independent variables.
 * \tparam DependentVariableType Type of dependent variables.
 * \param subDiagonal Sub-diagonal of matrix A.
 * \param superDiagonal Super-diagonal of matrix A.
 * \param diagonal Diagonal of matrix A.
 * \param diagonal Right-hand-side of matrix equation.
 * \return Solution to matrix equation.
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
        std::cerr << "Error, input provided for diagonal and sub/super "
                     " diagonals incorrect." << std::endl;
    }

    // Check whether solution will not be singular.
    if ( diagonal[ 0 ] == 0.0 )
    {
        std::cerr <<"Error when inverting tridiagonal system, "
                    "first entry of diagonal is zero" << std::endl;
    }

    std::vector< IndependentVariableType > intermediateVector( matrixSize );
    std::vector< DependentVariableType > solution( matrixSize );

    // Perform solution algorithm, from (Press W.H., et al., 2002).
    double scalingFactor = diagonal[ 0 ];
    solution[ 0 ]= rightHandSide[ 0 ] / scalingFactor;

    for ( unsigned int j = 1; j < matrixSize; j++ )
    {
        intermediateVector[ j ] = superDiagonal[ j - 1 ] / scalingFactor;
        scalingFactor = diagonal[ j ] - subDiagonal[ j - 1 ] * intermediateVector[ j ];

        // Check whether solution will not be singular.
        if ( scalingFactor == 0.0 )
        {
            std::cerr<<"Error when inverting tridiagonal system,"
                       " scaling factor equals zero!"<<std::endl;
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
 * Cubic spline interpolator, implementation from (Press W.H., et al., 2002).
 * Natural boundary conditions are imposed, meaning zero second derivatives (curvature) at end
 * points. Continuity of first derivatives is imposed.
 * \tparam IndependentVariableType Type of independent variables.
 * \tparam DependentVariableType Type of dependent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType >
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

    //! Cubic spline interpolator constructor.
    /*!
     * Cubic spline interpolator constructor taking separate vectors of dependent and independent
     * variable values.
     * \param independentVariables Vector with the independent variable values.
     * \param dependentVariables Vector with the dependent variable values.
     * \param selectedLookupScheme Look-up scheme that is to be used when finding interval
     * of requested independent variable value.
     */
    CubicSplineInterpolator( std::vector< IndependentVariableType > independentVariables,
                             std::vector< DependentVariableType > dependentVariables,
                             AvailableLookupScheme selectedLookupScheme = hunting_algorithm )

    {
        // Verify that the initialization variables are not empty.
        if ( independentVariables.size( ) == 0 || dependentVariables.size( ) == 0 )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "The vectors used in the cubic spline interpolator initialization are empty." ) ) );
        }

        // Set dependent and independent variable values.
        independentValues_ = independentVariables;
        dependentValues_ = dependentVariables;

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );

        // Create zero value for initializing output.
        zeroValue_ = dependentVariables[ 0 ] - dependentVariables[ 0 ];

        if ( dependentValues_.size( ) != independentValues_.size( ) )
        {
            std::cerr << "Warning: independent and dependent variables"
                         " not of same size in cubic spline constrcutor" << std::endl;
        }

        // Calculate second derivatives of curve.
        calculateSecondDerivatives( );
    }

    //! Cubic spline interpolator constructor.
    /*!
     * Cubic spline interpolator constructor taking single map of independent and dependent
     * variable values.
     * \param dataMap Map with the independent variable values as keys and corresponding
     * dependent variable values as values.
     * \param selectedLookupScheme Lookup scheme that is to be used when finding interval
     * of requested independent variable value.
     */
    CubicSplineInterpolator(
            const std::map< IndependentVariableType, DependentVariableType > dataMap,
            const AvailableLookupScheme selectedLookupScheme = hunting_algorithm )
    {
        // Verify that the initialization variables are not empty.
        if ( dataMap.size( ) == 0 )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
               "The vectors used in the cubic spline interpolator initialization are empty." ) ) );
        }

        // Set data vector member variables from map.
        independentValues_.resize( dataMap.size( ) );
        dependentValues_.resize( dataMap.size( ) );
        int counter = 0;
        for ( typename std::map< IndependentVariableType, DependentVariableType >::const_iterator
              mapIterator = dataMap.begin( ); mapIterator != dataMap.end( ); mapIterator++ )
        {
            independentValues_[ counter ] = mapIterator->first;
            dependentValues_[ counter ] = mapIterator->second;
            counter++;
        }

        // Create lookup scheme.
        this->makeLookupScheme( selectedLookupScheme );

        // Create zero value for initializing output.
        zeroValue_ = dependentValues_[ 0 ] - dependentValues_[ 0 ];

        // Calculate second derivatives of curve.
        calculateSecondDerivatives( );
    }

    // Statement required to prevent hiding of base class functions.
    using Interpolator< IndependentVariableType, DependentVariableType, 1 >::interpolate;

    //! Interpolate.
    /*!
     * Executes interpolation of data at a given target value of the independent variable, to
     * yield an interpolated value of the dependent variable.
     * \param targetIndependentVariableValue Target independent variable value at which point
     * the interpolation is performed.
     * \return Interpolated dependent variable value.
     */
    DependentVariableType interpolate(
            const IndependentVariableType targetIndependentVariableValue )
    {
        using std::pow;

        // Determine the lower entry in the table corresponding to the target independent variable
        // value.
        int lowerEntry_ = lookUpScheme_->findNearestLowerNeighbour(
                    targetIndependentVariableValue );

        // Get independent variable values bounding interval in which requested value lies.
        IndependentVariableType lowerValue, upperValue, squareDifference;
        lowerValue = independentValues_[ lowerEntry_ ];
        upperValue = independentValues_[ lowerEntry_ + 1 ];

        // Calculate coefficients A,B,C,D (see Numerical (Press W.H., et al., 2002))
        squareDifference = ( upperValue - lowerValue ) * ( upperValue - lowerValue );
        IndependentVariableType coefficientA_ = ( upperValue - targetIndependentVariableValue )
                / ( upperValue - lowerValue );
        IndependentVariableType coefficientB_ = 1 - coefficientA_;
        IndependentVariableType coefficientC_ = ( coefficientA_ * coefficientA_ * coefficientA_ -
                                                  coefficientA_ ) / 6.0 * squareDifference;
        IndependentVariableType coefficientD_ = ( coefficientB_ * coefficientB_ * coefficientB_ -
                                                  coefficientB_ ) / 6.0 * squareDifference;

        // The interpolated dependent variable value.
        return coefficientA_ * dependentValues_[ lowerEntry_ ] +
                coefficientB_ * dependentValues_[ lowerEntry_ + 1 ] +
                coefficientC_ * secondDerivativeOfCurve_[ lowerEntry_ ] +
                coefficientD_ * secondDerivativeOfCurve_[ lowerEntry_ + 1 ];
    }

protected:

private:

    //! Calculates the second derivatives of the curve.
    /*!
     * This function calculates the second derivatives of the curve at the nodes, assuming
     * the first derivatives to be continuous at the nodes and imposing natural spline conditions
     * (zero curvature at endpoints). The methodology is described in (Press W.H., et al., 2002).
     */
    void calculateSecondDerivatives( )
    {
        // Get length of vector.
        numberOfDataPoints_ = independentValues_.size( );

        // Sub-diagonal of tri-diagonal matrix.
        std::vector< IndependentVariableType > aCoefficients_;
        aCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Diagonal of tri-diagonal matrix.
        std::vector< IndependentVariableType > bCoefficients_;
        bCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Super-diagonal of tri-diagonal matrix.
        std::vector< IndependentVariableType > cCoefficients_;
        cCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Right-hand side of tridiagonal matrix system
        std::vector< DependentVariableType > rCoefficients_;
        rCoefficients_.resize( numberOfDataPoints_ - 2 );

        // Temporary value vector.
        std::vector< IndependentVariableType > hCoefficients_;
        hCoefficients_.resize( numberOfDataPoints_ - 1 );

        // Set second derivatives of curve to zero at endpoints, i.e. impose natural spline
        // condition.
        aCoefficients_[ numberOfDataPoints_- 3 ] = independentValues_[ 0 ] -
                                                   independentValues_[ 0 ];
        cCoefficients_[ numberOfDataPoints_- 3 ] = independentValues_[ 0 ] -
                                                   independentValues_[ 0 ];

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
                solveTridiagonalMatrixEquation< IndependentVariableType, DependentVariableType >
                ( aCoefficients_, bCoefficients_, cCoefficients_,  rCoefficients_ );

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
    /*!
     *  Vector filled with second derivative of curvature of each point.
     */
    std::vector< DependentVariableType > secondDerivativeOfCurve_;

    //! The number of datapoints.
    /*!
     * The number of datapoints.
     */
    unsigned int numberOfDataPoints_;

    //! Zero value of independent variable type
    /*!
     *  Zero value of independent variable type, computed by subtracting a value from itself.
     */
    DependentVariableType zeroValue_;
};

//! Typedef for cubic spline interpolator with (in)dependent = double.
typedef CubicSplineInterpolator< double, double > CubicSplineInterpolatorDouble;

//! Typedef for shared-pointer to cubic spline interpolator with (in)dependent = double.
typedef boost::shared_ptr< CubicSplineInterpolatorDouble > CubicSplineInterpolatorDoublePointer;

} // namespace interpolators
} // namespace tudat

#endif // TUDAT_CUBIC_SPLINE_INTERPOLATOR_H
