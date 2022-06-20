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
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *        Program, Volume II - Program Formulation, Douglas Aircraft Aircraft Company, 1973.
 *
 */

#include <string>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/pointer_cast.hpp>
#include <memory>

#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "tudat/math/geometric/compositeSurfaceGeometry.h"
#include "tudat/math/geometric/surfaceGeometry.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace aerodynamics
{

using Eigen::Vector6d;
using mathematical_constants::PI;

using namespace geometric_shapes;

//! Returns default values of mach number for use in HypersonicLocalInclinationAnalysis.
std::vector< double > getDefaultHypersonicLocalInclinationMachPoints(
        const std::string& machRegime )
{
    std::vector< double > machPoints;

    // Set default points for full hypersonic analysis.
    if ( machRegime == "Full" )
    {
        machPoints.resize( 6 );

        machPoints[ 0 ] = 3.0;
        machPoints[ 1 ] = 4.0;
        machPoints[ 2 ] = 5.0;
        machPoints[ 3 ] = 8.0;
        machPoints[ 4 ] = 10.0;
        machPoints[ 5 ] = 20.0;
    }

    // Set default points for low hypersonic analysis.
    else if ( machRegime == "Low" )
    {
        machPoints.resize( 5 );
        machPoints[ 0 ] = 3.0;
        machPoints[ 1 ] = 4.0;
        machPoints[ 2 ] = 5.0;
        machPoints[ 3 ] = 8.0;
        machPoints[ 4 ] = 10.0;
    }

    // Set default points for high hypersonic analysis.
    else if ( machRegime == "High" )
    {
        machPoints.resize( 4 );
        machPoints[ 0 ] = 5.0;
        machPoints[ 1 ] = 8.0;
        machPoints[ 2 ] = 10.0;
        machPoints[ 3 ] = 20.0;
    }
    return machPoints;
}

//! Returns default values of angle of attack for use in HypersonicLocalInclinationAnalysis.
std::vector< double > getDefaultHypersonicLocalInclinationAngleOfAttackPoints( )
{
    std::vector< double > angleOfAttackPoints;

    // Set number of data points and allocate memory.
    angleOfAttackPoints.resize( 11 );

    // Set default values, 0 to 40 degrees, with steps of 5 degrees.
    for ( int i = 0; i < 11; i++ )
    {
        angleOfAttackPoints[ i ] =
                ( static_cast< double >( i ) * 5.0 * PI / 180.0 );
    }
    return angleOfAttackPoints;
}

//! Returns default values of angle of sideslip for use in HypersonicLocalInclinationAnalysis.
std::vector< double > getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( )
{
    std::vector< double > angleOfSideslipPoints;

    // Set number of data points and allocate memory.
    angleOfSideslipPoints.resize( 2 );

    // Set default values, 0 and 1 degrees.
    angleOfSideslipPoints[ 0 ] = 0.0;
    angleOfSideslipPoints[ 1 ] = 1.0 * PI / 180.0;

    return angleOfSideslipPoints;
}

//! Function that saves the vehicle mesh data used for a HypersonicLocalInclinationAnalysis to a file
void saveVehicleMeshToFile(
        const std::shared_ptr< HypersonicLocalInclinationAnalysis > localInclinationAnalysis,
        const std::string directory,
        const std::string filePrefix )
{
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshPoints =
            localInclinationAnalysis->getMeshPoints( );
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshSurfaceNormals =
            localInclinationAnalysis->getPanelSurfaceNormals( );


//    boost::array< int, 3 > independentVariables;
//    independentVariables[ 0 ] = 0;
//    independentVariables[ 1 ] = 6;
//    independentVariables[ 2 ] = 0;

//    std::vector< std::vector< std::vector< double > > > pressureCoefficients =
//            localInclinationAnalysis->getPressureCoefficientList( independentVariables );

    int counter = 0;
    std::map< int, Eigen::Vector3d > meshPointsList;
    std::map< int, Eigen::Vector3d > surfaceNormalsList;
//    std::map< int, Eigen::Vector1d > pressureCoefficientsList;

    for( unsigned int i = 0; i < meshPoints.size( ); i++ )
    {
        for( unsigned int j = 0; j < meshPoints.at( i ).shape( )[ 0 ] - 1; j++ )
        {
            for( unsigned int k = 0; k < meshPoints.at( i ).shape( )[ 1 ] - 1; k++ )
            {
                meshPointsList[ counter ] = meshPoints[ i ][ j ][ k ];
                surfaceNormalsList[ counter ] = meshSurfaceNormals[ i ][ j ][ k ];
//                pressureCoefficientsList[ counter ] = ( Eigen::Vector1d( ) << pressureCoefficients[ i ][ j ][ k ] ).finished( );
                counter++;
            }
        }
    }

    input_output::writeDataMapToTextFile(
                meshPointsList, filePrefix + "ShapeFile.dat", directory );
    input_output::writeDataMapToTextFile(
                surfaceNormalsList, filePrefix + "SurfaceNormalFile.dat", directory );

//    input_output::writeDataMapToTextFile(
//                pressureCoefficientsList, filePrefix + "pressureCoefficientFile.dat", directory );
}


//! Default constructor.
HypersonicLocalInclinationAnalysis::HypersonicLocalInclinationAnalysis(
        const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
        const std::shared_ptr< SurfaceGeometry > inputVehicleSurface,
        const std::vector< int >& numberOfLines,
        const std::vector< int >& numberOfPoints,
        const std::vector< bool >& invertOrders,
        const std::vector< std::vector< int > >& selectedMethods,
        const double referenceArea,
        const double referenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool savePressureCoefficients )
    : AerodynamicCoefficientGenerator< 3, 6 >(
          dataPointsOfIndependentVariables, referenceLength, referenceArea, referenceLength,
          momentReferencePoint, { mach_number_dependent, angle_of_attack_dependent, angle_of_sideslip_dependent },true, false ),
      stagnationPressureCoefficient( 2.0 ),
      ratioOfSpecificHeats( 1.4 ),
      selectedMethods_( selectedMethods ),
      savePressureCoefficients_( savePressureCoefficients )
{
    // Set geometry if it is a single surface.
    if ( std::dynamic_pointer_cast< SingleSurfaceGeometry > ( inputVehicleSurface ) !=
         std::shared_ptr< SingleSurfaceGeometry >( ) )
    {
        // Set number of geometries and allocate memory.
        vehicleParts_.resize( 1 );

        vehicleParts_[ 0 ] = std::make_shared< LawgsPartGeometry >( );
        vehicleParts_[ 0 ]->setReversalOperator( invertOrders[ 0 ] );

        // If part is not already a LaWGS part, convert it.
        if ( std::dynamic_pointer_cast< LawgsPartGeometry >
             ( inputVehicleSurface ) ==
             std::shared_ptr< LawgsPartGeometry >( ) )
        {
            // Convert geometry to LaWGS surface mesh and set in vehicleParts_ list.
            vehicleParts_[ 0 ]->setMesh(
                        std::dynamic_pointer_cast< SingleSurfaceGeometry > ( inputVehicleSurface ),
                        numberOfLines[ 0 ], numberOfPoints[ 0 ] );
        }

        // Else, set geometry directly.
        else
        {
            vehicleParts_[ 0 ] = std::dynamic_pointer_cast< LawgsPartGeometry >
                    ( inputVehicleSurface );
        }
    }

    // Set geometry if it is a composite surface.
    else if ( std::dynamic_pointer_cast< CompositeSurfaceGeometry >( inputVehicleSurface ) !=
              std::shared_ptr< CompositeSurfaceGeometry >( ) )
    {
        // Dynamic cast to composite surface geometry for further processing.
        std::shared_ptr< CompositeSurfaceGeometry > compositeSurfaceGeometry_ =
                std::dynamic_pointer_cast< CompositeSurfaceGeometry >( inputVehicleSurface );

        // Set number of geometries and allocate memory.
        int numberOfVehicleParts =
                compositeSurfaceGeometry_->getNumberOfSingleSurfaceGeometries( );
        vehicleParts_.resize( numberOfVehicleParts );

        // Iterate through all parts and set them in vehicleParts_ list.
        for ( int i = 0; i < numberOfVehicleParts; i++ )
        {
            // If part is not already a LaWGS part, convert it.
            if ( std::dynamic_pointer_cast< LawgsPartGeometry >
                 ( compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ) ) ==
                 std::shared_ptr< LawgsPartGeometry >( ) )
            {
                vehicleParts_[ i ] = std::make_shared< LawgsPartGeometry >( );
                vehicleParts_[ i ]->setReversalOperator( invertOrders[ i ] );

                // Convert geometry to LaWGS and set in list.
                vehicleParts_[ i ]->setMesh(
                            compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ),
                            numberOfLines[ i ], numberOfPoints[ i ] );
            }

            // Else, set geometry directly.
            else
            {
                vehicleParts_[ i ] = std::dynamic_pointer_cast< LawgsPartGeometry >(
                            compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ) );
            }
        }
    }

    // Allocate memory for panel inclinations and pressureCoefficient_.
    inclination_.resize( vehicleParts_.size( ) );
    pressureCoefficient_.resize( vehicleParts_.size( ) );
    for ( unsigned int i = 0 ; i < vehicleParts_.size( ); i++ )
    {
        inclination_[ i ].resize( vehicleParts_[ i ]->getNumberOfLines( ) );
        pressureCoefficient_[ i ].resize( vehicleParts_[ i ]->getNumberOfLines( ) );
        for ( int j = 0 ; j < vehicleParts_[ i ]->getNumberOfLines( ) ; j++ )
        {
            inclination_[ i ][ j ].resize( vehicleParts_[ i ]->getNumberOfPoints( ) );
            pressureCoefficient_[ i ][ j ].resize( vehicleParts_[ i ]->getNumberOfPoints( ) );
        }
    }

    boost::array< int, 3 > numberOfPointsPerIndependentVariables;
    for( int i = 0; i < 3; i++ )
    {
        numberOfPointsPerIndependentVariables[ i ] =
                dataPointsOfIndependentVariables_[ i ].size( );
    }

    isCoefficientGenerated_.resize( numberOfPointsPerIndependentVariables );

    std::fill( isCoefficientGenerated_.origin( ),
               isCoefficientGenerated_.origin( ) + isCoefficientGenerated_.num_elements( ), 0 );

    generateCoefficients( );
    createInterpolator( );
}

//! Get aerodynamic coefficients.
Vector6d HypersonicLocalInclinationAnalysis::getAerodynamicCoefficientsDataPoint(
        const boost::array< int, 3 > independentVariables )
{
    if( isCoefficientGenerated_( independentVariables ) == 0 )
    {
        determineVehicleCoefficients( independentVariables );
    }

    // Return requested coefficients.
    return aerodynamicCoefficients_( independentVariables );
}

//! Generate aerodynamic database.
void HypersonicLocalInclinationAnalysis::generateCoefficients( )
{
    // Allocate variable to pass to coefficient determination for independent
    // variable indices.
    boost::array< int, 3 > independentVariableIndices;

    // Iterate over all combinations of independent variables.
    for ( unsigned int i = 0 ; i < dataPointsOfIndependentVariables_[ 0 ].size( ) ; i++ )
    {
        independentVariableIndices[ 0 ] = i;
        for ( unsigned  int j = 0 ; j < dataPointsOfIndependentVariables_[
              1 ].size( ) ; j++ )
        {
            independentVariableIndices[ 1 ] = j;
            for ( unsigned  int k = 0 ; k < dataPointsOfIndependentVariables_[
                  2 ].size( ) ; k++ )
            {
                independentVariableIndices[ 2 ] = k;

                determineVehicleCoefficients( independentVariableIndices );
            }

        }
    }
}

//! Generate aerodynamic coefficients at a single set of independent variables.
void HypersonicLocalInclinationAnalysis::determineVehicleCoefficients(
        const boost::array< int, 3 > independentVariableIndices )
{
    // Declare coefficients vector and initialize to zeros.
    Vector6d coefficients = Vector6d::Zero( );

    // Loop over all vehicle parts, calculate aerodynamic coefficients and add
    // to aerodynamicCoefficients_.
    for ( unsigned int i = 0 ; i < vehicleParts_.size( ) ; i++ )
    {
        coefficients += determinePartCoefficients( i, independentVariableIndices );
    }

    if( savePressureCoefficients_ )
    {
        pressureCoefficientList_[ independentVariableIndices ] = pressureCoefficient_;
    }

    aerodynamicCoefficients_( independentVariableIndices ) = coefficients;
    isCoefficientGenerated_( independentVariableIndices ) = 1;
}

//! Determine aerodynamic coefficients of a single vehicle part.
Vector6d HypersonicLocalInclinationAnalysis::determinePartCoefficients(
        const int partNumber, const boost::array< int, 3 > independentVariableIndices )
{
    // Declare and determine angles of attack and sideslip for analysis.
    double angleOfAttack =  dataPointsOfIndependentVariables_[ 1 ]
            [ independentVariableIndices[ 1 ] ];

    double angleOfSideslip =  dataPointsOfIndependentVariables_[ 2 ]
            [ independentVariableIndices[ 2 ] ];

    // Declare partCoefficient vector.
    Vector6d partCoefficients = Vector6d::Zero( );

    // Check whether the inclinations of the vehicle part have already been computed.
    if ( previouslyComputedInclinations_.count( std::pair< double, double >(
                                                    angleOfAttack, angleOfSideslip ) ) == 0 )
    {
        // Determine panel inclinations for part.
        determineInclinations( angleOfAttack, angleOfSideslip );

        // Add panel inclinations to container
        previouslyComputedInclinations_[ std::pair< double, double >(
                    angleOfAttack, angleOfSideslip ) ] = inclination_;
    }

    else
    {
        // Fetch inclinations from container
        inclination_ = previouslyComputedInclinations_[ std::pair< double, double >(
                    angleOfAttack, angleOfSideslip ) ];
    }

    // Set pressureCoefficient_ array for given independent variables.
    determinePressureCoefficients( partNumber, independentVariableIndices );

    // Calculate force coefficients from pressure coefficients.
    partCoefficients.segment( 0, 3 ) = calculateForceCoefficients( partNumber );

    // Calculate moment coefficients from pressure coefficients.
    partCoefficients.segment( 3, 3 ) = calculateMomentCoefficients( partNumber );

    return partCoefficients;
}

//! Determine the pressure coefficients on a single vehicle part.
void HypersonicLocalInclinationAnalysis::determinePressureCoefficients(
        const int partNumber, const boost::array< int, 3 > independentVariableIndices )
{
    // Retrieve Mach number.
    double machNumber = dataPointsOfIndependentVariables_[ 0 ]
            [ independentVariableIndices[ 0 ] ];

    // Determine stagnation point pressure coefficients. Value is computed once
    // here to prevent its calculation in inner loop.
    stagnationPressureCoefficient = computeStagnationPressure(
                machNumber, ratioOfSpecificHeats );

    updateCompressionPressures( machNumber, partNumber );
    updateExpansionPressures( machNumber, partNumber );
}

//! Determine force coefficients from pressure coefficients.
Eigen::Vector3d HypersonicLocalInclinationAnalysis::calculateForceCoefficients(
        const int partNumber )
{
    // Declare force coefficient vector and intialize to zeros.
    Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );

    // Loop over all panels and add pressures, scaled by panel area, to force
    // coefficients.
    for ( int i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
    {
        for ( int j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++)
        {
            forceCoefficients -=
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
Eigen::Vector3d HypersonicLocalInclinationAnalysis::calculateMomentCoefficients(
        const int partNumber )
{
    // Declare moment coefficient vector and intialize to zeros.
    Eigen::Vector3d momentCoefficients = Eigen::Vector3d::Zero( );

    // Declare moment arm for panel moment determination.
    Eigen::Vector3d referenceDistance;

    // Loop over all panels and add moments due pressures.
    for ( int i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
    {
        for ( int j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
        {
            // Determine moment arm for given panel centroid.
            referenceDistance = ( vehicleParts_[ partNumber ]->getPanelCentroid( i, j ) -
                                  momentReferencePoint_ );

            momentCoefficients -=
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
void HypersonicLocalInclinationAnalysis::determineInclinations( const double angleOfAttack,
                                                                const double angleOfSideslip )
{
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
    for( unsigned int k = 0; k < vehicleParts_.size( ); k++ )
    {
        for ( int i = 0 ; i < vehicleParts_[ k ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( int j = 0 ; j < vehicleParts_[ k ]->getNumberOfPoints( ) - 1 ; j++ )
            {

                // Determine cosine of inclination angle from inner product between
                // surface normal and free-stream direction.
                cosineOfInclination = vehicleParts_[ k ]->
                        getPanelSurfaceNormal( i, j ).
                        dot( freestreamVelocityDirection );

                // Set inclination angle.
                inclination_[ k ][ i ][ j ] = PI / 2.0 - acos( cosineOfInclination );
            }
        }
    }
}

//! Determine compression pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateCompressionPressures( const double machNumber,
                                                                     const int partNumber )
{
    int method = selectedMethods_[ 0 ][ partNumber ];

    std::function< double( double ) > pressureFunction;

    // Switch to analyze part using correct method.
    switch( method )
    {
    case 0:
        pressureFunction =
                std::bind( aerodynamics::computeNewtonianPressureCoefficient, std::placeholders::_1 );
        break;

    case 1:
        pressureFunction =
                std::bind( aerodynamics::computeModifiedNewtonianPressureCoefficient, std::placeholders::_1,
                           stagnationPressureCoefficient );
        break;

    case 2:
        // Method currently disabled.
        break;

    case 3:
        // Method currently disabled.
        break;

    case 4:
        pressureFunction =
                std::bind( aerodynamics::computeEmpiricalTangentWedgePressureCoefficient, std::placeholders::_1,
                           machNumber );
        break;

    case 5:
        pressureFunction =
                std::bind( aerodynamics::computeEmpiricalTangentConePressureCoefficient, std::placeholders::_1,
                           machNumber );
        break;

    case 6:
        pressureFunction =
                std::bind( aerodynamics::computeModifiedDahlemBuckPressureCoefficient, std::placeholders::_1,
                           machNumber );
        break;

    case 7:
        pressureFunction =
                std::bind( aerodynamics::computeVanDykeUnifiedPressureCoefficient, std::placeholders::_1,
                           machNumber, ratioOfSpecificHeats, 1 );
        break;

    case 8:
        pressureFunction =
                std::bind( aerodynamics::computeSmythDeltaWingPressureCoefficient, std::placeholders::_1,
                           machNumber );
        break;

    case 9:
        pressureFunction =
                std::bind( aerodynamics::computeHankeyFlatSurfacePressureCoefficient, std::placeholders::_1,
                           machNumber );
        break;

    default:
        break;
    }

    for ( int i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1; i++ )
    {
        for ( int j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
        {
            if ( inclination_[ partNumber ][ i ][ j ] > 0 )
            {
                // If panel inclination is positive, calculate pressure coefficient.
                pressureCoefficient_[ partNumber ][ i ][ j ] =
                        pressureFunction( inclination_[ partNumber ][ i ][ j ] );
            }
        }
    }
}

//! Determines expansion pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateExpansionPressures( const double machNumber,
                                                                   const int partNumber )
{
    // Get analysis method of part to analyze.
    int method = selectedMethods_[ 1 ][ partNumber ];

    if ( method == 0 || method == 1 || method == 4 )
    {
        std::function< double( ) > pressureFunction;
        switch( method )
        {
        case 0:
            pressureFunction = std::bind( &aerodynamics::computeVacuumPressureCoefficient,
                                          machNumber, ratioOfSpecificHeats );
            break;

        case 1:
            pressureFunction = [ ]( ){ return 0.0; };
            break;

        case 4:
            pressureFunction = std::bind( &aerodynamics::computeHighMachBasePressure,
                                          machNumber );
            break;

        }

        // Iterate over all panels on part.
        for ( int i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( int j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Van Dyke unified method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            pressureFunction( );
                }
            }
        }
    }

    else if( method == 3 || method == 5 || method == 6 )
    {

        std::function< double( double ) > pressureFunction;

        // Declare local variable.
        double freestreamPrandtlMeyerFunction;

        // Switch to analyze part using correct method.
        switch( method )
        {
        case 3:
            // Calculate freestream Prandtl-Meyer function.
            freestreamPrandtlMeyerFunction = aerodynamics::computePrandtlMeyerFunction(
                        machNumber, ratioOfSpecificHeats );
            pressureFunction =
                    std::bind( &aerodynamics::computePrandtlMeyerFreestreamPressureCoefficient,
                               std::placeholders::_1, machNumber, ratioOfSpecificHeats,
                               freestreamPrandtlMeyerFunction );
            break;

        case 5:
            pressureFunction =
                    std::bind( &aerodynamics::computePrandtlMeyerFreestreamPressureCoefficient,
                               std::placeholders::_1, machNumber, ratioOfSpecificHeats, -1 );
            break;

        case 6:
            pressureFunction = std::bind( &aerodynamics::computeAcmEmpiricalPressureCoefficient,
                                          std::placeholders::_1, machNumber );
            break;
        }

        // Iterate over all panels on part.
        for ( int i = 0 ; i < vehicleParts_[ partNumber ]->getNumberOfLines( ) - 1 ; i++ )
        {
            for ( int j = 0 ; j < vehicleParts_[ partNumber ]->getNumberOfPoints( ) - 1 ; j++ )
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Van Dyke unified method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            pressureFunction( inclination_[ partNumber ][ i ][ j ] );
                }
            }
        }
    }

    else
    {
        std::string errorMessage = "Error, expansion local inclination method number "
                + std::to_string( method ) + " not recognized";
        throw std::runtime_error( errorMessage );
    }
}

} // namespace aerodynamics
} // namespace tudat
