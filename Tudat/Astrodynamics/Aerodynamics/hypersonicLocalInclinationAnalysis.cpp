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
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *        Program, Volume II - Program Formulation, Douglas Aircraft Aircraft Company, 1973.
 *
 */

#include <iostream>
#include <string>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Geometry>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/vehicleExternalModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
{
namespace aerodynamics
{

using std::string;
using std::endl;
using tudat::mathematics::PI;

using namespace tudat::mathematics::geometric_shapes;

//! Constructor to set geometry and reference quantities.
void HypersonicLocalInclinationAnalysis::setVehicle(
        bodies::Vehicle& vehicle, std::vector< int > numberOfLines,
        std::vector< int > numberOfPoints, std::vector< bool > invertOrders )
{
    // Retrieve external surface geometry from vehicle.
    boost::shared_ptr< bodies::VehicleExternalModel > externalModel_ = vehicle.getExternalModel( );
    boost::shared_ptr< SurfaceGeometry > surface_ = externalModel_->getVehicleExternalGeometry( );

    // Set geometry if it is a single surface.
    if ( boost::dynamic_pointer_cast< SingleSurfaceGeometry > ( surface_ ) !=
         boost::shared_ptr< SingleSurfaceGeometry >( ) )
    {
        // Set number of geometries and allocate memory.
        numberOfVehicleParts_ = 1;
        vehicleParts_.resize( 1 );

        vehicleParts_[ 0 ] = boost::make_shared< LawgsPartGeometry >( );
        vehicleParts_[ 0 ]->setReversalOperator( invertOrders[ 0 ] );

        // Convert geometry to LaWGS surface mesh and set in vehicleParts_ list.
        vehicleParts_[ 0 ]->setMesh(
                boost::dynamic_pointer_cast< SingleSurfaceGeometry > ( surface_ ),
                numberOfLines[ 0 ], numberOfPoints[ 0 ] );
    }

    // Set geometry if it is a composite surface.
    else if ( boost::dynamic_pointer_cast< CompositeSurfaceGeometry >( surface_ ) !=
              boost::shared_ptr< CompositeSurfaceGeometry >( ) )
    {
        // Dynamic cast to composite surface geometry for further processing.
        boost::shared_ptr< CompositeSurfaceGeometry > compositeSurfaceGeometry_ =
                boost::dynamic_pointer_cast< CompositeSurfaceGeometry >( surface_ );

        // Set number of geometries and allocate memory.
        numberOfVehicleParts_ =
                compositeSurfaceGeometry_->getNumberOfSingleSurfaceGeometries( );
        vehicleParts_.resize( numberOfVehicleParts_ );

        // Iterate through all parts and set them in vehicleParts_ list.
        for ( int i = 0; i < numberOfVehicleParts_; i++ )
        {
            // If part is not already a LaWGS part, convert it.
            if ( boost::dynamic_pointer_cast< LawgsPartGeometry >
                 ( compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ) ) ==
                 boost::shared_ptr< LawgsPartGeometry >( ) )
            {
                vehicleParts_[ i ] = boost::make_shared< LawgsPartGeometry >( );
                vehicleParts_[ i ]->setReversalOperator( invertOrders[ i ] );

                // Convert geometry to LaWGS and set in list.
                vehicleParts_[ i ]->setMesh(
                            compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ),
                            numberOfLines[ i ], numberOfPoints[ i ] );
            }

            // Else, set geometry directly.
            else
            {
                vehicleParts_[ i ] = boost::dynamic_pointer_cast< LawgsPartGeometry >(
                            compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ) );
            }
        }
    }

    // Allocate memory for arrays of pressure coefficients, inclinations and methods.
    allocateArrays( );
}

//! Allocate pressure coefficient, inclination, and method arrays.
void HypersonicLocalInclinationAnalysis::allocateArrays( )
{
    // Allocate memory for panel inclinations and pressureCoefficient_.
    inclination_.resize( numberOfVehicleParts_ );
    pressureCoefficient_.resize( numberOfVehicleParts_ );
    for ( int i = 0 ; i < numberOfVehicleParts_; i++ )
    {
        inclination_[ i ].resize( vehicleParts_[ i ]->getNumberOfLines( ) );
        pressureCoefficient_[ i ].resize( vehicleParts_[ i ]->getNumberOfLines( ) );
        for ( int j = 0 ; j<vehicleParts_[ i ]->getNumberOfLines( ) ; j++ )
        {
            inclination_[ i ][ j ].resize( vehicleParts_[ i ]->getNumberOfPoints( ) );
            pressureCoefficient_[ i ][ j ].resize( vehicleParts_[ i ]->getNumberOfPoints( ) );
        }
    }

    // If methods have not yet been set by user, allocate memory and set
    // defaults.
    if ( selectedMethods_.num_elements( ) == 0 )
    {
        // Set memory for expansion and compression methods.
        selectedMethods_.resize( boost::extents[ 4 ][ 2 * numberOfVehicleParts_ ] );
        for ( int i = 0; i < 4 ; i++  )
        {
            for ( int j = 0; j < 2 * numberOfVehicleParts_ ; j++ )
            {
                // Sets all local inclination methods to a default of ( Modified
                // Newtonian for compression and Newtonian for expansion ).
                selectedMethods_[ i ][ j ] = 1;
            }
        }
    }
}

//! Get aerodynamic coefficients.
Eigen::VectorXd HypersonicLocalInclinationAnalysis::getAerodynamicCoefficients(
    const std::vector< int >& independentVariables )
{
    // If coefficients have not been allocated (and independent variables
    // have not been initialized), do so.
    if ( vehicleCoefficients_.size( ) == 0 )
    {
        allocateVehicleCoefficients( );
    }

    // Declare and determine index in vehicleCoefficients_ array.
    int coefficientsIndex = variableIndicesToListIndex( independentVariables );

    // If coefficients for data point have not yet been calculated, calculate them.
    if ( vehicleCoefficients_[ coefficientsIndex ] == boost::shared_ptr< Eigen::VectorXd >( ) )
    {
        determineVehicleCoefficients( independentVariables );
    }

    // Return requested coefficients.
    return *( vehicleCoefficients_[ coefficientsIndex ] );
}

//! Set local inclination methods for all parts (expansion and compression).
void HypersonicLocalInclinationAnalysis::setSelectedMethods(
        std::vector< std::vector < int > > selectedMethods )
{
    //For loops loop through input methods and set the analysis methods.
    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < numberOfVehicleParts_; j++ )
        {
            selectedMethods_[ i ][ j ] = selectedMethods[ i ][ j ];
        }
    }
}

//! Generate aerodynamic database.
void HypersonicLocalInclinationAnalysis::generateDatabase( )
{
    // If coefficients have not been allocated, do so.
    if ( vehicleCoefficients_.size( ) == 0 )
    {
        allocateVehicleCoefficients( );
    }

    int i, j, k;
    int l = 0;

    // Allocate variable to pass to coefficient determination for independent
    // variable indices.
    std::vector< int > independentVariableIndices;
    independentVariableIndices.resize( numberOfIndependentVariables_ );

    // Iterate over all combinations of independent variables.
    for ( i = 0 ; i < numberOfPointsPerIndependentVariables_[ machIndex_ ] ; i++ )
    {
        independentVariableIndices[ machIndex_ ] = i;
        for ( j = 0 ; j < numberOfPointsPerIndependentVariables_[
                angleOfAttackIndex_ ] ; j++ )
        {
            independentVariableIndices[ angleOfAttackIndex_ ] = j;
            for ( k = 0 ; k < numberOfPointsPerIndependentVariables_[
                    angleOfSideslipIndex_ ] ; k++ )
            {
                independentVariableIndices[ angleOfSideslipIndex_ ] = k;

                // If coefficient has not yet been set, calculate and set it.
                if ( vehicleCoefficients_[ l ] == boost::shared_ptr< Eigen::VectorXd >( ) )
                {
                    determineVehicleCoefficients( independentVariableIndices );
                }

                l++;
            }
        }
    }
}

//! Allocate aerodynamic coefficient array and NULL independent variables.
void HypersonicLocalInclinationAnalysis::allocateVehicleCoefficients( )
{
    // If angle of attack points have not yet been set, use defaults.
    if ( dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ].num_elements() == 0 )
    {
        setDefaultAngleOfAttackPoints( );
    }

    // If angle of sideslip points have not yet been set, use defaults.
    if ( dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ].num_elements() == 0 )
    {
        setDefaultAngleOfSideslipPoints( );
    }

    // If Mach number points have not yet been set, use defaults.
    if ( dataPointsOfIndependentVariables_[ machIndex_ ].num_elements() == 0 )
    {
        setDefaultMachPoints( );
    }

    // Determine number of combinations of independent variables.
    numberOfCases_ = numberOfPointsPerIndependentVariables_[ machIndex_ ] *
                        numberOfPointsPerIndependentVariables_[ angleOfSideslipIndex_ ] *
                        numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ];

    // Allocate memory for pointers to coefficients and initialize to NULL shared_ptrs.
    vehicleCoefficients_.resize( numberOfCases_ );
}

//! Generate aerodynamic coefficients at a single set of independent variables.
void HypersonicLocalInclinationAnalysis::determineVehicleCoefficients(
        std::vector< int > independentVariableIndices )
{
    int i;

    // Determine index in vehicleCoefficients_ at which the coefficient
    // should be set.
    int coefficientsIndex = variableIndicesToListIndex( independentVariableIndices );

    // Declare coefficients vector and initialize to zeros.
    Eigen::VectorXd coefficients = Eigen::VectorXd( 6 );
    for ( i = 0 ; i < 6 ; i++ )
    {
        coefficients[ i ] = 0;
    }

    // Declare vehicle part coefficients variable.
    Eigen::VectorXd partCoefficients;

    // Loop over all vehicle parts, calculate aerodynamic coefficients and add
    // to vehicleCoefficients_.
    for ( i = 0 ; i < numberOfVehicleParts_ ; i++ )
    {
        partCoefficients
                = determinePartCoefficients( i, independentVariableIndices );
        coefficients += partCoefficients;
    }

    // Allocate and set vehicle coefficients at given independent variables.
    vehicleCoefficients_[ coefficientsIndex ]
            = boost::make_shared< Eigen::VectorXd>( coefficients );
}

//! Determine aerodynamic coefficients of a single vehicle part.
Eigen::VectorXd HypersonicLocalInclinationAnalysis::determinePartCoefficients(
        int partNumber, std::vector< int > independentVariableIndices )
{
    // Declare and determine angles of attack and sideslip for analysis.
    double angleOfAttack = dataPointsOfIndependentVariables_[
            angleOfAttackIndex_ ][ independentVariableIndices[
                    angleOfAttackIndex_ ] ];
    double angleOfSideslip = dataPointsOfIndependentVariables_[
            angleOfSideslipIndex_ ][ independentVariableIndices[
                    angleOfSideslipIndex_ ] ];

    // Declare partCoefficient vector.
    Eigen::VectorXd partCoefficients = Eigen::VectorXd( 6 );

    // Determine panel inclinations for part.
    determineInclination( partNumber, angleOfAttack, angleOfSideslip );

    // Set pressureCoefficient_ array for given independent variables.
    determinePressureCoefficients( partNumber, independentVariableIndices );

    // Calculate force coefficients from pressure coefficients.
    partCoefficients.head( 3 ) = calculateForceCoefficients( partNumber );

    // Calculate moment coefficients from pressure coefficients.
    partCoefficients.segment( 3, 3 ) = calculateMomentCoefficients(
            partNumber );

    return partCoefficients;
}

//! Determine the pressure coefficients on a single vehicle part.
void HypersonicLocalInclinationAnalysis::determinePressureCoefficients(
        int partNumber, std::vector< int > independentVariableIndices )
{
    // Retrieve Mach number.
    double machNumber = dataPointsOfIndependentVariables_[ machIndex_ ]
                        [ independentVariableIndices[ machIndex_ ] ];

    // Determine stagnation point pressure coefficients. Value is computed once
    // here to prevent its calculation in inner loop.
    stagnationPressureCoefficient = computeStagnationPressure(
            machNumber, ratioOfSpecificHeats );
    updateCompressionPressures( machNumber, partNumber );
    updateExpansionPressures( machNumber, partNumber );
}

//! Determine force coefficients from pressure coefficients.
Eigen::VectorXd HypersonicLocalInclinationAnalysis::calculateForceCoefficients( int partNumber )
{
    int i, j;

    // Declare force coefficient vector and intialize to zeros.
    Eigen::Vector3d forceCoefficients;
    for ( i = 0 ; i < 3 ; i++)
    {
        forceCoefficients( i ) = 0.0;
    }
    
    // Loop over all panels and add pressures, scaled by panel area, to force 
    // coefficients.
    for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
    {
        for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++)
        {
            forceCoefficients = forceCoefficients -
                     pressureCoefficient_[ partNumber ][ i ][ j ] *
                     vehicleParts_[ partNumber ]->getPanelArea( i, j ) *
                     vehicleParts_[ partNumber ]->getPanelSurfaceNormal( i, j );
        }
    }

    // Normalize result by reference area.
    forceCoefficients /= referenceArea_;

    return forceCoefficients;
}

//! Determine moment coefficients from pressure coefficients.
Eigen::VectorXd HypersonicLocalInclinationAnalysis::calculateMomentCoefficients( int partNumber )
{
    int i, j;

    // Declare moment coefficient vector and intialize to zeros.
    Eigen::Vector3d momentCoefficients;
    for ( i = 0 ; i < 3 ; i++ )
    {
        momentCoefficients( i ) = 0.0;
    }

    // Declare moment arm for panel moment determination.
    Eigen::Vector3d referenceDistance ;

    // Loop over all panels and add moments due pressures.
    for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
    {
        for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
        {
            // Determine moment arm for given panel centroid.
            referenceDistance = ( vehicleParts_[ partNumber ]->
                        getPanelCentroid( i, j ) -  momentReferencePoint_ );

            momentCoefficients = momentCoefficients -
                pressureCoefficient_[ partNumber ][ i ][ j ] *
                vehicleParts_[ partNumber ]->getPanelArea( i, j ) *
                ( referenceDistance.cross( vehicleParts_[ partNumber ]->
                                           getPanelSurfaceNormal( i, j ) ) );
        }
    }

    // Scale result by reference length and area.
    momentCoefficients /= ( referenceLength_ * referenceArea_ );

    return momentCoefficients;
}

//! Determines the inclination angle of panels on a single part.
void HypersonicLocalInclinationAnalysis::determineInclination( int partNumber,
                                                               double angleOfAttack,
                                                               double angleOfSideslip )
{
    int i, j;

    // Declare free-stream velocity vector.
    Eigen::Vector3d freestreamVelocityDirection;

    // Set freestream velocity vector in body frame.
    double freestreamVelocityDirectionX = cos( angleOfAttack )* cos( angleOfSideslip );
    double freestreamVelocityDirectionY = sin( angleOfSideslip );
    double freestreamVelocityDirectionZ = sin( angleOfAttack ) * cos( angleOfSideslip );
    freestreamVelocityDirection( 0 ) = freestreamVelocityDirectionX;
    freestreamVelocityDirection( 1 ) = freestreamVelocityDirectionY;
    freestreamVelocityDirection( 2 ) = freestreamVelocityDirectionZ;

    // Declare cosine of inclination angle.
    double cosineOfInclination;

    // Loop over all panels of given vehicle part and set inclination angles.
    for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
    {
        for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
        {

            // Determine cosine of inclination angle from inner product between
            // surface normal and free-stream direction.
            cosineOfInclination = -1 * vehicleParts_[ partNumber ]->
                                  getPanelSurfaceNormal( i, j ).
                                  dot( freestreamVelocityDirection );

            // Set inclination angle.
            inclination_[ partNumber ][ i ][ j ] = PI / 2.0 - acos( cosineOfInclination );
        }
    }
}

//! Determine compression pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateCompressionPressures( double machNumber,
                                                                     int partNumber )
{
    int i, j;
    int method;

    // Determine which method to use based on Mach number regime and given
    // selectedMethods__.
    if ( machRegime_ == "Full" || machRegime_ == "High")
    {
        method = selectedMethods_[ 0 ][ partNumber ];
    }

    else
    {
        method = selectedMethods_[ 2 ][ partNumber ];
    }

    // Switch to analyze part using correct method.
    switch( method )
    {
    case 0:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // Newtonian method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                       aerodynamics::computeNewtonianPressureCoefficient(
                               inclination_[ partNumber ][ i ][ j ] );
                }
            }
        }

        break;

    case 1:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    // If panel inclination is positive, calculate pressure
                    // using Modified Newtonian method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                        aerodynamics::computeModifiedNewtonianPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ],
                            stagnationPressureCoefficient );
                }
            }
        }

        break;

    case 2:

        // Method currently disabled.
        break;

    case 3:

        // Method currently disabled.
        break;

    case 4:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // empirical Tangent Wedge method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                       aerodynamics::computeEmpiricalTangentWedgePressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber );
                }
            }
        }

        break;

    case 5:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // empirical Tangent Cone method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                       aerodynamics::computeEmpiricalTangentConePressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber );
                }
            }
        }

        break;

    case 6:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    // If panel inclination is positive, calculate pressure using
                    // Modified Dahlem Buck method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeModifiedDahlemBuckPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber );
                }
            }
        }

        break;

    case 7:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // van Dyke unified method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeVanDykeUnifiedPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber,
                            ratioOfSpecificHeats, 1 );
                }
            }
        }

        break;

    case 8:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // Smyth delta wing method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeSmythDeltaWingPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber);
                }
            }
        }

        break;

    case 9:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is positive, calculate pressure using
                // Hankey flat surface method.
                if ( inclination_[ partNumber ][ i ][ j ] > 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeHankeyFlatSurfacePressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber );
                }
            }
        }

        break;

    default:

        break;
    }
}

//! Determines expansion pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateExpansionPressures( double machNumber,
                                                                   int partNumber )
{
    int i,j;

    // Get analysis method of part to analyze.
    int method;
    if ( machRegime_ == "Full" || machRegime_ == "High" )
    {
        method = selectedMethods_[ 1 ][ partNumber ];
    }

    else
    {
        method = selectedMethods_[ 3 ][ partNumber ];
    }

    // Declare local variable.
    double freestreamPrandtlMeyerFunction;

    // Switch to analyze part using correct method.
    switch( method )
    {
    case 0:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is negative, calculate pressure using
                // vacuum method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeVacuumPressureCoefficient(
                                machNumber, ratioOfSpecificHeats );
                }
            }
        }

        break;

    case 1:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is negative, calculate pressure using
                // Newtonian expansion method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] = 0;
                }
            }
        }

        break;

    case 2:

        // Option currently disabled.
        break;

    case 3:

        // Calculate freestream Prandtl-Meyer function.
        freestreamPrandtlMeyerFunction = aerodynamics::computePrandtlMeyerFunction(
                machNumber, ratioOfSpecificHeats );

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Prandtl-Meyer expansion from freestream.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computePrandtlMeyerFreestreamPressureCoefficient(
                                inclination_[ partNumber ][ i ][ j ], machNumber,
                                ratioOfSpecificHeats, freestreamPrandtlMeyerFunction );
                }
            }
        }

        break;

    case 4:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                // If panel inclination is negative, calculate pressure using
                // high Mach base pressure method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeHighMachBasePressure( machNumber );
                }
            }
        }

        break;

    case 5:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Van Dyke unified method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeVanDykeUnifiedPressureCoefficient(
                                    inclination_[ partNumber ][ i ][ j ], machNumber,
                                    ratioOfSpecificHeats, -1 );
                }
            }
        }

        break;

    case 6:

        // Iterate over all panels on part.
        for ( i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++)
            {
                // If panel inclination is negative, calculate pressure using
                // ACM empirical method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeAcmEmpiricalPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ], machNumber );
                }
            }
        }

        break;

    default:

        break;
    }
}

//! Set the default Mach number points.
void HypersonicLocalInclinationAnalysis::setDefaultMachPoints( )
{
    // Set default points for full hypersonic analysis.
    if ( machRegime_ == "Full" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 6;
        dataPointsOfIndependentVariables_[ machIndex_ ].resize(
                    boost::extents[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ] );
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 3.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 4.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 4 ] = 10.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 5 ] = 20.0;
    }

    // Set default points for low hypersonic analysis.
    else if ( machRegime_ == "Low" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 5;
        dataPointsOfIndependentVariables_[ machIndex_ ].resize(
                    boost::extents[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ] );
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 3.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 4.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 4 ] = 10.0;
    }

    // Set default points for high hypersonic analysis.
    else if ( machRegime_ == "High" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 4;
        dataPointsOfIndependentVariables_[ machIndex_ ].resize(
                    boost::extents[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ] );
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 10.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 20.0;
    }
}

//! Set the default analysis points for angle of attack.
void HypersonicLocalInclinationAnalysis::setDefaultAngleOfAttackPoints( )
{
    // Set number of data points and allocate memory.
    numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ] = 11;
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ].resize(
                boost::extents[ numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ] ] );

    // Set default values, 0 to 40 degrees, with steps of 5 degrees.
    int i;
    for ( i = 0; i < numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ]; i++ )
    {
        dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ i ] =
                ( static_cast< double >( i ) * 5.0 * PI / 180.0 );
    }

}

//! Set the default analysis points for angle of sideslip.
void HypersonicLocalInclinationAnalysis::setDefaultAngleOfSideslipPoints( )
{
    // Set number of data points and allocate memory.
    numberOfPointsPerIndependentVariables_[ angleOfSideslipIndex_ ] = 2;
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ].resize(
                boost::extents[ numberOfPointsPerIndependentVariables_[ angleOfSideslipIndex_ ] ] );

    // Set default values, 0 and 1 degrees.
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ 0 ] = 0.0;
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ 1 ] =
            1.0 * PI / 180.0;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          HypersonicLocalInclinationAnalysis& hypersonicLocalInclinationAnalysis )
{
    stream << "This is a hypersonic local inclination analysis object."<< endl;
    stream << "The Mach regime is "
           << hypersonicLocalInclinationAnalysis.getMachRegime( ) << endl;
    stream << "The name of the vehicle being analyzed is "
           << hypersonicLocalInclinationAnalysis.getVehicleName( ) << endl;
    stream << "It contains "
           << hypersonicLocalInclinationAnalysis.getNumberOfVehicleParts( )
           << " parts in Lawgs format. " << endl;
    stream << "The names of the vehicle parts are ";

    for ( int i = 0; i < hypersonicLocalInclinationAnalysis.getNumberOfVehicleParts( ); i++ )
    {
        stream << hypersonicLocalInclinationAnalysis.getVehiclePart( i )->getName( ) << ", ";
    }

    stream << endl;

    // Return stream.
    return stream;
}

} // namespace aerodynamics
} // namespace tudat
