/*!   \file aerodynamicCoefficientGenerator.cpp
 *    This file contains the definition of the aerodynamic coefficient generator
 *    base class.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : Bart Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : B.Romgens@student.tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 4 February,  2011
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas Aircraft
 *        Aircraft Company, 1973.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      102511    D. Dirkx          First version of file.
 *      110501    D. Dirkx          Added more comments.
 *      112701    D. Dirkx          Finalized for code check.
 *      110131    B. Romgens        Minor modifications during code check.
 *      110204    D. Dirkx          Finalized code.
 */

// Include statements.
#include "aerodynamicCoefficientGenerator.h"
#include "linearAlgebra.h"

//! Default constructor.
AerodynamicCoefficientGenerator::AerodynamicCoefficientGenerator( )
{
    // Set number if independent variables to 0.
    numberOfIndependentVariables_ = 0;

    // Set arrays to NULL.
    numberOfPointsPerIndependentVariables_ = NULL;
    dataPointsOfIndependentVariables_ = NULL;

    vehicleCoefficients_ = NULL;
}

//! Default destructor.
AerodynamicCoefficientGenerator::~AerodynamicCoefficientGenerator( )
{
    int i;
    // Delete data points of each independent variable and reset to NULL.
    for( i = 0; i < numberOfIndependentVariables_; i++ )
    {
        delete [ ] dataPointsOfIndependentVariables_[ i ];
        dataPointsOfIndependentVariables_[ i ] = NULL;
    }

    // Deletes arrays of data points and values of independent variables.
    delete [ ] dataPointsOfIndependentVariables_;
    delete [ ] numberOfPointsPerIndependentVariables_;
    dataPointsOfIndependentVariables_ = NULL;
    numberOfPointsPerIndependentVariables_ = NULL;
}

//! Set sreference area.
void AerodynamicCoefficientGenerator::setReferenceArea(
        const double& referenceArea )
{
    referenceArea_ = referenceArea;
}

//! Gets reference area.
double AerodynamicCoefficientGenerator::getReferenceArea( )
{
    return referenceArea_ ;
}

//! Sets reference length.
void AerodynamicCoefficientGenerator::setReferenceLength(
        double referenceLength )
{
    referenceLength_ = referenceLength;
}

//! Gets reference length.
double AerodynamicCoefficientGenerator::getReferenceLength( )
{
    return referenceLength_ ;
}

//! Sets moment reference point.
void AerodynamicCoefficientGenerator::setMomentReferencePoint(
        Vector3d momentReferencePoint)
{
    momentReferencePoint_ = momentReferencePoint;
}

//! Gets moment reference point.
VectorXd AerodynamicCoefficientGenerator::getMomentReferencePoint( )
{
    return momentReferencePoint_ ;
}

//! Sets the number of independent variables
void AerodynamicCoefficientGenerator::setNumberOfIndependentVariables(
        const int& numberOfVariables )
{
    // Set number of variables.
    numberOfIndependentVariables_ = numberOfVariables;

    // Allocate memory for data points.
    numberOfPointsPerIndependentVariables_ = new int[ numberOfVariables ];
    dataPointsOfIndependentVariables_ = new double*[ numberOfVariables ];

    // Initialize arrays of data points to NULL.
    int i;
    for( i = 0; i < numberOfVariables; i++ )
    {
        dataPointsOfIndependentVariables_[ i ] = NULL;
    }
}



//! Sets the number of points for Mach number.
void AerodynamicCoefficientGenerator::setNumberOfMachPoints(
        const int& numberOfMachPoints )
{
    // Set value of number of mach points.
    numberOfPointsPerIndependentVariables_ [ machIndex_ ] = numberOfMachPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ machIndex_ ] =
            new double[ numberOfMachPoints ];
}

//! Sets the number of points for angle of attack.
void AerodynamicCoefficientGenerator::setNumberOfAngleOfAttackPoints(
        const int& numberOfAngleOfAttackPoints )
{
    // Set value of number of angle of attack points.
    numberOfPointsPerIndependentVariables_ [ angleOfAttackIndex_ ] =
            numberOfAngleOfAttackPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ] =
            new double[ numberOfAngleOfAttackPoints ];
}

//! Sets the number of points for angle of sideslip.
void AerodynamicCoefficientGenerator::setNumberOfAngleOfSideslipPoints(
        const int& numberOfAngleOfSideslipPoints )
{
    // Set value of number of angle of sideslip points.
    numberOfPointsPerIndependentVariables_ [ angleOfSideslipIndex_ ] =
            numberOfAngleOfSideslipPoints;

    // Allocate memory for data points.
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ] =
            new double[ numberOfAngleOfSideslipPoints ];
}

//! Gets the number of independent variables
int AerodynamicCoefficientGenerator::getNumberOfIndependentVariables( )
{
    return numberOfIndependentVariables_;
}

//! Gets the number of points for an independent variable.
int AerodynamicCoefficientGenerator::getNumberOfValuesOfIndependentVariable(
        const int& independentVariable )
{
    return numberOfPointsPerIndependentVariables_ [ independentVariable ];

}

//! Gets the number of points for Mach number.
int AerodynamicCoefficientGenerator::getNumberOfMachPoints( )
{
    return numberOfPointsPerIndependentVariables_ [ machIndex_ ];
}

//! Gets the number of points for angle of attack.
int AerodynamicCoefficientGenerator::getNumberOfAngleOfAttackPoints( )
{
    return numberOfPointsPerIndependentVariables_ [ angleOfAttackIndex_ ] ;
}

//! Gets the number of points for angle of sideslip.
int AerodynamicCoefficientGenerator::getNumberOfAngleOfSideslipPoints( )
{
    return numberOfPointsPerIndependentVariables_ [ angleOfSideslipIndex_ ];
}

//! Sets a Mach number point.
void AerodynamicCoefficientGenerator::setMachPoint( const int& index,
                                                    const double& machPoint )
{
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ] = machPoint;
}

//! Sets an angle of attack point.
void AerodynamicCoefficientGenerator::setAngleOfAttackPoint(
        const int& index,
        const double& angleOfAttackPoint )
{
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ] =
            angleOfAttackPoint;
}

//! Sets an angle of sideslip point.
void AerodynamicCoefficientGenerator::setAngleOfSideslipPoint(
        const int& index,
        const double& angleOfSideslipPoint )
{
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ index ] =
            angleOfSideslipPoint;
}

//! Gets a value of an independent variable.
double AerodynamicCoefficientGenerator::getIndependentVariablePoint(
                                    const int& independentVariable,
                                    const int& index )
{
    return dataPointsOfIndependentVariables_[ independentVariable ][ index ];
}

//! Gets a Mach number point.
double AerodynamicCoefficientGenerator::getMachPoint( const int& index )
{
    return dataPointsOfIndependentVariables_[ machIndex_ ][ index ];
}

//! Gets an angle of attack point.
double AerodynamicCoefficientGenerator::getAngleOfAttackPoint( const int& index )
{
    return dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ index ];
}

//! Gets an angle of sideslip point.
double AerodynamicCoefficientGenerator::getAngleOfSideslipPoint( const int& index )
{
    return dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ index ];
}

//! Function to convert the independent variable indices to list index in
//! vehicleCoefficients_.
int AerodynamicCoefficientGenerator::variableIndicesToListIndex(
        int* independentVariableIndices)
{
    int i ,j;

    // Declare requested indexs and initialize to zero.
    int coefficientsIndex_ = 0;

    // Declare variable of contribution of single index.
    int singleStepContribution_;

    // Iterate over all indices and add to coefficientsIndex.
    for( i = 0; i < numberOfIndependentVariables_; i++ )
    {
        // Determine single step contribution.
        singleStepContribution_ = 1;
        for( j = i + 1; j < numberOfIndependentVariables_; j++ )
        {
            singleStepContribution_ *=
                    numberOfPointsPerIndependentVariables_[ j ];
        }

        // Add contribution to reqeusted index.
        coefficientsIndex_ += singleStepContribution_ *
                             independentVariableIndices[ i ];
    }

    return coefficientsIndex_;

}

// End of file.
