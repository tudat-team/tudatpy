/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas
 *        Aircraft Company, 1973.
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
 * and data points of independent variables. Currently supports Mach number, angles of attack
 * sideslip, and Reynold number as independent variables, but others can be added easily. To do so,
 * only the interface function need to be added. As an example, for adding the Mach number, only
 * the get/set numberOfMachPoints and get/set MachPoint, as well as a machIndex_, to let the code
 * know which index of the numberOfPointsPerIndependentVariables_ and
 * dataPointsOfIndependentVariables_ arrays represent the Mach number. Coefficients are stored in a
 * VectorXd of pointers which should be allocated and initialized to NULL shared_ptr by a derived
 * class. The conversion from the indices of the independent variables to the index in this
 * VectorXd is performed by the variableIndicesToListIndex function.
 */
class AerodynamicCoefficientGenerator: public AerodynamicCoefficientInterface
{
public:

    //! Default constructor.
    /*!
     * Default constructor, sets a number of independent variables to 0.
     */
    AerodynamicCoefficientGenerator( )
        : numberOfIndependentVariables_( 0 ),
          machIndex_( -0 ),
          angleOfAttackIndex_( -0 ),
          angleOfSideslipIndex_( -0 ),
          reynoldsNumberIndex_( -0 ),
          numberOfCases_( 0 )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~AerodynamicCoefficientGenerator( ) { }

    //! Set the number of independent variables
    /*!
     * Sets the number of independent variables.
     * \param numberOfVariables Number of independent variables for analysis.
     */
    void setNumberOfIndependentVariables( const int numberOfVariables );

    //! Set the number of points for Mach number.
    /*!
     * Sets the number of different Mach numbers at which coefficients are determined.
     * \param numberOfMachPoints Number of data points for Mach number.
     */
    void setNumberOfMachPoints( const int numberOfMachPoints );

    //! Sets the number of points for angle of attack.
    /*!
     * Sets the number of different angles of attack at which coefficients are determined.
     * \param numberOfAngleOfAttackPoints Number of data points for angle of attack.
     */
    void setNumberOfAngleOfAttackPoints( const int numberOfAngleOfAttackPoints );

    //! Set the number of points for angle of sideslip.
    /*!
     * Sets the number of different angles of sideslip at which coefficients are determined.
     * \param numberOfAngleOfSideslipPoints Number of data points for angle of sideslip.
     */
    void setNumberOfAngleOfSideslipPoints( const int numberOfAngleOfSideslipPoints );

    //! Set the number of points for the Reynolds Number.
    /*!
     * Sets the number of different Reynolds Number. at which coefficients are determined.
     * \param numberOfReynoldsNumberPoints Number of data points for ReynoldsNumber.
     */
    void setNumberOfReynoldsNumberPoints( const int numberOfReynoldsNumberPoints );

    //! Get the number of independent variables
    /*!
     * Returns the number of independent variables.
     * \return numberOfVariables Number of independent variables for analysis.
     */
    int getNumberOfIndependentVariables( ) { return numberOfIndependentVariables_; }

    //! Get the number of points for an independent variable.
    /*!
     * Returns the number of different values of a given independent variable
     * at which coefficients are determined.
     * \param independentVariable Independent variable from which to retrieve
     * number of data points.
     * \return Number of data points for Mach number
     */
    int getNumberOfValuesOfIndependentVariable( const int independentVariable )
    {
        return numberOfPointsPerIndependentVariables_ [ independentVariable ];
    }

    //! Get the number of points for Mach number.
    /*!
     * Returns the number of different Mach numbers at which coefficients are determined.
     * \return Number of data points for Mach number
     */
    int getNumberOfMachPoints( ) { return numberOfPointsPerIndependentVariables_ [ machIndex_ ]; }

    //! Get the number of points for angle of attack.
    /*!
     * Returns the number of different angles of attack at which coefficients are determined.
     * \return Number of data points for angle of attack.
     */
    int getNumberOfAngleOfAttackPoints ( )
    {
        return numberOfPointsPerIndependentVariables_ [ angleOfAttackIndex_ ];
    }

    //! Get the number of points for angle of sideslip.
    /*!
     * Returns the number of different angles of sideslip at which coefficients are determined.
     * \return Number of data points for angle of sideslip.
     */
    int getNumberOfAngleOfSideslipPoints ( )
    {
        return numberOfPointsPerIndependentVariables_ [ angleOfSideslipIndex_ ];
    }

    //! Get the number of points for Reynold number.
    /*!
     * Returns the number of different Reynolds numbers at which coefficients are determined.
     * \return Number of data points for Reynold number.
     */
    int getNumberOfReynoldsNumberPoints( )
    {
        return numberOfPointsPerIndependentVariables_ [ reynoldsNumberIndex_ ];
    }

    //! Set a Mach number point.
    /*!
     * Sets a Mach number at which aerodynamic coefficients are to be determined.
     * \param index Index in Mach number data point list at which to set the value.
     * \param machPoint Value of Mach number to set.
     */
    void setMachPoint( const int index, const double machPoint )
    {
        dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ] = machPoint;
    }

    //! Set an angle of attack point.
    /*!
     * Sets an angle of attack at which aerodynamic coefficients are to be determined.
     * \param index Index in angle of attack number data point list at which to set the value.
     * \param angleOfAttackPoint Value of angle of attack to set.
     */
    void setAngleOfAttackPoint( const int index, const double angleOfAttackPoint )
    {
        dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ] = angleOfAttackPoint;
    }

    //! Set an angle of sideslip point.
    /*!
     * Sets an angle of sideslip at which aerodynamic coefficients are to be determined.
     * \param index Index in angle of sideslip number data point list at which to set the value.
     * \param angleOfSideslipPoint Value of angle of sideslip to set.
     */
    void setAngleOfSideslipPoint( const int index, const double angleOfSideslipPoint )
    {
        dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ index ] = angleOfSideslipPoint;
    }

    //! Set a Reynolds Number point.
    /*!
     * Sets a ReyNolds Number at which aerodynamic coefficients are to be determined.
     * \param index Index in Reynold number data point list at which to set the value.
     * \param reynoldsNumberPoint Value of Reynold number to set.
     */
    void setReynoldsNumberPoint( const int index, const double reynoldsNumberPoint )
    {
        dataPointsOfIndependentVariables_[ reynoldsNumberIndex_ ][ index ] = reynoldsNumberPoint;
    }

    //! Get a Mach number point.
    /*!
     * Returns a Mach number at which aerodynamic coefficients are to be determined.
     * \param index Index from Mach number data point list from which to retrieve the value.
     * \return Value of Mach number at index.
     */
    double getMachPoint( const int index )
    {
        return dataPointsOfIndependentVariables_[ machIndex_ ][ index ];
    }

    //! Get an angle of attack point.
    /*!
     * Returns an angle of attack at which aerodynamic coefficients are to be determined.
     * \param index Index from angle of attack data point list from which to retrieve the value.
     * \return Value of angle of attack at index.
     */
    double getAngleOfAttackPoint( const int index )
    {
        return dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ];
    }

    //! Get an angle of sideslip point.
    /*!
     * Returns an angle of sideslip at which aerodynamic coefficients are to be determined.
     * \param index Index from angle of attack data point list from which to retrieve the value.
     * \return Value of angle of attack at index.
     */
    double getAngleOfSideslipPoint( const int index )
    {
        return dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ index ];
    }

    //! Get an Reynold Number point.
    /*!
     * Returns an Reynold Number at which aerodynamic coefficients are to be determined.
     * \param index Index from Reynold Number data point list from which to retrieve the value.
     * \return Value of Reynold Number at index.
     */
    double getReynoldsNumberPoint( const int index )
    {
        return dataPointsOfIndependentVariables_[ reynoldsNumberIndex_ ][ index ];
    }

    //! Get a value of an independent variable.
    /*!
     * Returns a value of an independent variable at which aerodynamic
     * coefficients are to be determined.
     * \param independentVariable Independent variable index.
     * \param index Index from independent variable data point list from which to retrieve the
     *                  value.
     * \return Value of angle of attack at index.
     */
    double getIndependentVariablePoint( const int independentVariable, const int index )
    {
        return dataPointsOfIndependentVariables_[ independentVariable ][ index ];
    }

    //! Get aerodynamic coefficients.
    /*!
     * Returns aerodynamic coefficients.
     * \param independentVariables Array of values of independent variable
     * indices in dataPointsOfIndependentVariables_.
     */
    virtual Eigen::VectorXd getAerodynamicCoefficients(
            const std::vector< int >& independentVariables ) = 0;

protected:

    //! List of pointers to VectorXds containing coefficients.
    /*!
     * List of pointers to VectorXds containing coefficients.
     */
    std::vector< boost::shared_ptr< Eigen::VectorXd > > vehicleCoefficients_;

    //! Convert the independent variable indices to list index in vehicleCoefficients_.
    /*!
     * Converts the independent variable indices to list index in vehicleCoefficients_.
     * \param independentVariableIndices Array of indices of independent variables.
     * \return Resulting index in vehicleCoefficients_.
     */
    int variableIndicesToListIndex( const std::vector< int >& independentVariableIndices );

    //! Number of independent variables in analysis.
    /*!
     * Number of independent variables in analysis. To be set from derived
     * class, depending on analysis type.
     */
    int numberOfIndependentVariables_;

    //! Array of number of points per independent variable in analysis.
    /*!
     * Array of number of points per independent variable. Physical meaning
     * of indices are determined by the machIndex, angleOfAttackIndex, etc.
     */
    boost::multi_array< int, 1 > numberOfPointsPerIndependentVariables_;

    //! Array of arrays of data points for independent variables.
    /*!
     * Array of arrays of data points for independent variables. Physical
     * meaning of indices are determined by the machIndex, angleOfAttackIndex, etc.
     */
    std::vector< boost::multi_array< int, 1 > > dataPointsOfIndependentVariables_;

    //! Index in independent variables arrays representing Mach number.
    /*!
     * Index in independent variables arrays representing Mach number.
     */
    int machIndex_;

    //! Index in independent variables arrays representing angle of attack.
    /*!
     * Index in independent variables arrays representing angle of attack.
     */
    int angleOfAttackIndex_;

    //! Index in independent variables arrays representing angle of sideslip.
    /*!
     * Index in independent variables arrays representing angle of sideslip.
     */
    int angleOfSideslipIndex_;

    //! Index in independent variables arrays representing Reynold Number.
    /*!
     * Index in independent variables arrays representing Reynold Number.
     */
    int reynoldsNumberIndex_;

    //! Number of data points in aerodynamic database.
    /*!
     * Number of data points in aerodynamic database, should equal product of number of data
     * points per independent variables tha are used.
     */
    int numberOfCases_;
};

} // aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H
