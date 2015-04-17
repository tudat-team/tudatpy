/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      102511    D. Dirkx          First version of file.
 *      110501    D. Dirkx          Added more comments.
 *      112701    D. Dirkx          Finalized for code check.
 *      110131    B. Romgens        Minor modifications during code check.
 *      110204    D. Dirkx          Finalized code.
 *      110615    F.M. Engelen      Made a child of Coefficient Database. Moved aerodynamic
 *                                  reference quatities to the parent class.
 *      120825    A. Ronse          Changed dataPointsOfIndependentVariables_ to array of doubles.
 *                                  Fixed bug in setMachPoint function.
 *      120912    D. Dirkx          Templatized class, adjusted to meet RAII idiom.
 *      140129    D. Dirkx          Changed Doxygen descriptions
 *      140130    T. Roegiers       Changed Doxygen descriptions
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas
 *        Aircraft Company, 1973.
 *
 *    Notes
 *
 */

#ifndef TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H
#define TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H

#include <vector>

#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{

//! Base class for aerodynamic coefficient generator.
/*!
 * Abstract base class for aerodynamic analysis method. Stores independent variable values
 * and data points of independent variables. Coefficients are stored in a multi_array of pointers.
 */
template< int NumberOfIndependentVariables, int NumberOfCoefficients = 6 >
class AerodynamicCoefficientGenerator
{
public:

    //! Default base class constructor.
    /*!
     *  Default base class constructor, sets independent variable data points and aerodynamics
     *  reference quantities.
     *  \param dataPointsOfIndependentVariables Vector of vector, where each subvector contains
     *  the data points of each of the independent variables for the coefficient generation.
     *  The number of subvectors must be equal to the number of independent variables. It is
     *  recommended that each of the subvectors is sorted in ascending order.
     *  \param referenceArea Reference area used to non-dimensionalize aerodynamic forces
     *  and moments.
     *  \param referenceLength Reference length used to non-dimensionalize aerodynamic moments.
     *  \param momentReferencePoint Reference point wrt which aerodynamic moments are calculated.
     */
    AerodynamicCoefficientGenerator(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const double referenceArea,
            const double referenceLength,
            const Eigen::Vector3d momentReferencePoint ):
        dataPointsOfIndependentVariables_( dataPointsOfIndependentVariables ),
        referenceArea_( referenceArea ),
        referenceLength_( referenceLength ),
        momentReferencePoint_( momentReferencePoint )
    {
        // Check that the size of dataPointsOfIndependentVariables matches the template parameter.
        assert( dataPointsOfIndependentVariables_.size( ) == NumberOfIndependentVariables );

        boost::array< int, NumberOfIndependentVariables > numberOfPointsPerIndependentVariables;
        for( int i = 0; i < NumberOfIndependentVariables; i++ )
        {
            numberOfPointsPerIndependentVariables[ i ] = dataPointsOfIndependentVariables_[ i ].
                                                         size( );
        }

        aerodynamicCoefficients_.resize( numberOfPointsPerIndependentVariables );
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~AerodynamicCoefficientGenerator( ) { }

    //! Get the number of points for an independent variable.
    /*!
     * Returns the number of different values of a given independent variable
     * at which coefficients are determined.
     * \param independentVariable Independent variable from which to retrieve
     * number of data points.
     * \return Number of data points for selected independent variable.
     */
    int getNumberOfValuesOfIndependentVariable( const int independentVariable ) const
    {
        return dataPointsOfIndependentVariables_ [ independentVariable ].size( );
    }

    //! Get a value of an independent variable.
    /*!
     * Returns a value of an independent variable at which aerodynamic
     * coefficients are to be determined.
     * \param independentVariable Independent variable index.
     * \param index Index from independent variable list for which to retrieve the value.
     * \return Value of independent variable at index.
     */
    double getIndependentVariablePoint( const int independentVariable, const int index ) const
    {
        return dataPointsOfIndependentVariables_[ independentVariable ][ index ];
    }

    //! Get aerodynamic coefficients.
    /*!
     * Virtual function to return the aerodynamic coefficients at a specified
     * set of independent variables.
     * \param independentVariables Array of values of independent variable indices in
     * dataPointsOfIndependentVariables_ for which to retrieve aerodynamic coefficients.
     * \return vector of coefficients at specified independent variable indices.
     */
    virtual Eigen::Matrix< double, NumberOfCoefficients, 1 > getAerodynamicCoefficients(
            const boost::array< int, NumberOfIndependentVariables > independentVariables ) = 0;

    //! Generate aerodynamic coefficients.
    /*!
     * Virtual function to generate aerodynamic coefficients for the list of independent variables
     * passed to the constructor.
     */
    virtual void generateCoefficients( ) = 0;

protected:

    //! List of pointers to VectorXds containing coefficients.
    /*!
     * List of pointers to VectorXds containing coefficients.
     */
    boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients, 1 >,
    NumberOfIndependentVariables > aerodynamicCoefficients_;

    //! Array of arrays of data points for independent variables.
    /*!
     * Array of arrays of data points for independent variables.
     */
    std::vector< std::vector< double > > dataPointsOfIndependentVariables_;

    //! Aerodynamic reference area.
    /*!
     * Reference area with which aerodynamic forces and moments are non-dimensionalized.
     */
    double referenceArea_;

    //! Aerodynamic reference length.
    /*!
     * Reference length with which aerodynamic moments are non-dimensionalized.
     */
    double referenceLength_;

    //! Aerodynamic moment reference point.
    /*!
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;
};

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H
