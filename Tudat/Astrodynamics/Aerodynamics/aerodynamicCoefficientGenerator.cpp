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
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas Aircraft
 *        Aircraft Company, 1973.
 *
 */

#include <string>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"

namespace tudat
{
namespace aerodynamics
{

//! Set the number of independent variables.
void AerodynamicCoefficientGenerator::setNumberOfIndependentVariables(
    const int numberOfVariables )
{
    // Set number of variables.
    numberOfIndependentVariables_ = numberOfVariables;

    // Allocate memory for data points.
    numberOfPointsPerIndependentVariables_.resize( boost::extents[ numberOfVariables ] );
    dataPointsOfIndependentVariables_.resize( numberOfVariables );
}

//! Set the number of points for Mach number.
void AerodynamicCoefficientGenerator::setNumberOfMachPoints( const int numberOfMachPoints )
{
    // Set value of number of mach points.
    numberOfPointsPerIndependentVariables_ [ machIndex_ ] = numberOfMachPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ machIndex_ ].resize( boost::extents[ numberOfMachPoints ] );
}

//! Set the number of points for angle of attack.
void AerodynamicCoefficientGenerator::setNumberOfAngleOfAttackPoints(
        const int numberOfAngleOfAttackPoints )
{
    // Set value of number of angle of attack points.
    numberOfPointsPerIndependentVariables_ [ angleOfAttackIndex_ ] = numberOfAngleOfAttackPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ].resize(
                boost::extents[ numberOfAngleOfAttackPoints ] );
}

//! Set the number of points for angle of sideslip.
void AerodynamicCoefficientGenerator::setNumberOfAngleOfSideslipPoints(
        const int numberOfAngleOfSideslipPoints )
{
    // Set value of number of angle of sideslip points.
    numberOfPointsPerIndependentVariables_ [ angleOfSideslipIndex_ ]
            = numberOfAngleOfSideslipPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ].resize(
                boost::extents[ numberOfAngleOfSideslipPoints ]);
}

//! Set the number of points for the Reynolds Number.
void AerodynamicCoefficientGenerator::setNumberOfReynoldsNumberPoints(
        const int numberOfReynoldsNumberPoints )
{
    // Set value of number of Reynolds Number points.
    numberOfPointsPerIndependentVariables_ [ reynoldsNumberIndex_ ] = numberOfReynoldsNumberPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ reynoldsNumberIndex_ ].resize(
                boost::extents [ numberOfReynoldsNumberPoints ] );
}

//! Convert independent variable indices to list index in vehicleCoefficients_.
int AerodynamicCoefficientGenerator::variableIndicesToListIndex(
        const std::vector< int >& independentVariableIndices )
{
    int i, j;

    // Declare requested indexs and initialize to zero.
    int coefficientsIndex_ = 0;

    // Declare variable of contribution of single index.
    int singleStepContribution_;

    // Iterate over all indices and add to coefficientsIndex.
    for ( i = 0; i < numberOfIndependentVariables_; i++ )
    {
        // Determine single step contribution.
        singleStepContribution_ = 1;

        for ( j = i + 1; j < numberOfIndependentVariables_; j++ )
        {
            singleStepContribution_ *= numberOfPointsPerIndependentVariables_[ j ];
        }

        // Add contribution to requested index.
        coefficientsIndex_ += singleStepContribution_ * independentVariableIndices[ i ];
    }

    return coefficientsIndex_;
}

} // namespace aerodynamics
} // namespace tudat
