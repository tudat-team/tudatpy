/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Geometry>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Aerodynamics/rarefiedFlowAnalysis.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/InputOutput/SPARTADataReader.h"

namespace tudat
{
namespace aerodynamics
{

using Eigen::Vector6d;
using mathematical_constants::PI;

using namespace geometric_shapes;
using namespace unit_conversions;

//! Returns default values of molecular speed ratio for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowMolecularSpeedRatioPoints(
        const std::string& molecularSpeedRatioRegime )
{
    std::vector< double > molecularSpeedRatioPoints;

    // Set default points for full velocity analysis.
    if ( molecularSpeedRatioRegime == "Full" )
    {
        molecularSpeedRatioPoints.resize( 7 );
        molecularSpeedRatioPoints[ 0 ] = 1.0;
        molecularSpeedRatioPoints[ 1 ] = 2.5;
        molecularSpeedRatioPoints[ 2 ] = 5.0;
        molecularSpeedRatioPoints[ 3 ] = 10.0;
        molecularSpeedRatioPoints[ 4 ] = 25.0;
        molecularSpeedRatioPoints[ 5 ] = 50.0;
        molecularSpeedRatioPoints[ 6 ] = 100.0;
    }
    // Set default points for low velocity analysis.
    else if ( molecularSpeedRatioRegime == "Low" )
    {
        molecularSpeedRatioPoints.resize( 4 );
        molecularSpeedRatioPoints[ 0 ] = 1.0;
        molecularSpeedRatioPoints[ 1 ] = 2.5;
        molecularSpeedRatioPoints[ 2 ] = 5.0;
        molecularSpeedRatioPoints[ 3 ] = 10.0;
    }
    // Set default points for high velocity analysis.
    else if ( molecularSpeedRatioRegime == "High" )
    {
        molecularSpeedRatioPoints.resize( 4 );
        molecularSpeedRatioPoints[ 0 ] = 10.0;
        molecularSpeedRatioPoints[ 1 ] = 25.0;
        molecularSpeedRatioPoints[ 2 ] = 50.0;
        molecularSpeedRatioPoints[ 3 ] = 100.0;
    }
    return molecularSpeedRatioPoints;
}

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints(
        const std::string& angleOfAttackRegime )
{
    std::vector< double > angleOfAttackPoints;

    // Set default angles of attack
    double a = - 25;
    while ( a <= 30 )
    {
        angleOfAttackPoints.push_back( convertDegreesToRadians( a ) );
        a += 5;
    }

    // Add extra points if required
    if ( angleOfAttackRegime == "Full" )
    {
        std::vector< double > frontExtension = { convertDegreesToRadians( -75.0 ),
                                                 convertDegreesToRadians( -60.0 ),
                                                 convertDegreesToRadians( -45.0 ),
                                                 convertDegreesToRadians( -30.0 ) };
        std::vector< double > rearExtension = { convertDegreesToRadians( 30.0 ),
                                                convertDegreesToRadians( 45.0 ),
                                                convertDegreesToRadians( 60.0 ),
                                                convertDegreesToRadians( 75.0 ) };
        angleOfAttackPoints.insert( angleOfAttackPoints.begin( ), frontExtension.begin( ), frontExtension.end( ) );
        angleOfAttackPoints.insert( angleOfAttackPoints.end( ), rearExtension.begin( ), rearExtension.end( ) );
    }
    return angleOfAttackPoints;
}

//! Returns default values of angle of sideslip for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAngleOfSideslipPoints( )
{
    std::vector< double > angleOfSideslipPoints;

    // Set number of data points and allocate memory.
    angleOfSideslipPoints.resize( 1 );

    // Set default values, 0 and 1 degrees.
    angleOfSideslipPoints[ 0 ] = 0.0;

    return angleOfSideslipPoints;
}

//! Default constructor.
RarefiedFlowAnalysis::RarefiedFlowAnalysis(
        const std::string& SPARTAExecutable,
        const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
        const std::string simulationGases,
        const StandardAtmosphere& atmosphereModel,
        const std::string& geometryFileUser,
        const int referenceAxis,
        const Eigen::Vector3d& momentReferencePoint,
        const double wallTemperature = 300.0,
        const double accomodationCoefficient = 1.0 )
    : AerodynamicCoefficientGenerator< 3, 6 >(
          dataPointsOfIndependentVariables, referenceLength, referenceArea, referenceLength,
          momentReferencePoint,
          boost::assign::list_of( altitude_dependent )( mach_number_dependent )( angle_of_attack_dependent ),
          1, 0 ),
      dataPointsOfIndependentVariables_( dataPointsOfIndependentVariables ), simulationGases_( simulationGases ),
      referenceAxis_( referenceAxis ), wallTemperature_( wallTemperature ),
      accomodationCoefficient_( accomodationCoefficient )
{
    // Analyze vehicle geometry
    analyzeGeometryFile( geometryFileUser );

    // Find atmospheric conditions based on altitude
    atmosphericConditions_.resize( 6 );
    for ( unsigned int i = 0; i < simulationAltitudes.size( ); i++ )
    {
        atmosphericConditions_.at( 0 ).push_back( atmosphereModel.getDensity( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 1 ).push_back( atmosphereModel.getPressure( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 2 ).push_back( atmosphereModel.getTemperature( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 3 ).push_back( atmosphereModel.getSpecificGasConstant( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 4 ).push_back( tudat::physical_constants::MOLAR_GAS_CONSTANT /
                                                  atmosphericConditions_.at( 3 ).at( i ) );
        atmosphericConditions_.at( 5 ).push_back( tudat::physical_constants::AVOGADRO_CONSTANT *
                                                  atmosphericConditions_.at( 0 ).at( i ) /
                                                  atmosphericConditions_.at( 4 ).at( i ) );
    }

    // Get simulation conditions
    getSimulationConditions( );

    // Read SPARTA input template
    inputTemplate_ = readSPARTAInputFileTemplate( inputFileTemplate_ );

    // Copy input shape file to default name
    std::string commandString = "cp " + geometryFileUser_ + " " + geometryFileInternal_;
    std::system( commandString.c_str( ) );

    // Run SPARTA simulation
    runSPARTASimulation( SPARTAExecutable );

    // Create interpolator object
    createInterpolator( );
}

//! Get aerodynamic coefficients.
void RarefiedFlowAnalysis::analyzeGeometryFile( std::string& geometryFileUser )
{
    // Extract information on vehicle geometry
    std::pair< Eigen::Matrix< double, Eigen::Dynamic, 3 >, Eigen::Matrix< int, Eigen::Dynamic, 3 > >
            geometryData = input_output::readSPARTAGeometryFile( geometryFileUser );
    shapePoints_ = geometryData.first;
    shapeTriangles_ = geometryData.second;
    numberOfPoints_ = shapePoints_.cols( );
    numberOfTriangles_ = shapeTriangles_.cols( );

    // Get maximum and minimum values in each dimension
    maximumDimensions_ = shapePoints_.colwise( ).maxCoeff( );
    minimumDimensions_ = shapePoints_.colwise( ).minCoeff( );
    maximumDimensions_ += 0.5 * maximumDimensions_; // add extra space around shape
    minimumDimensions_ += 0.5 * minimumDimensions_; // add extra space around shape

    // Compute normal to surface elements, area of surface elements and moment arm values
    Eigen::Matrix3d currentVertices;
    Eigen::Vector3d currentNormal;
    Eigen::Vector3d currentCentroid;
    double currentNormalNorm;
    elementSurfaceNormal_.resize( 3, numberOfTriangles_ );
    elementSurfaceArea_.resize( 1, numberOfTriangles_ );
    elementMomentArm_.resize( 3, numberOfTriangles_ );
    for ( int i = 0; i < numberOfTriangles_; i++ )
    {
        // Compute properties of current surface element
        for ( unsigned int j = 0; j < 3; j++ )
        {
            currentVertices.row( j ) = shapePoints_.row( shapeTriangles_( i, j ) - 1 );
        }
        currentNormal = ( currentVertices.row( 1 ) - currentVertices.row( 0 ) ).cross(
                    currentVertices.row( 2 ) - currentVertices.row( 0 ) );
        currentNormalNorm = currentNormal.norm( );
        currentCentroid = currentVertices.colwise( ).sum( ) / 3.0;

        // Find normal, area and distance to reference point
        elementSurfaceNormal_.col( i ) = currentNormal / currentNormalNorm;
        elementSurfaceArea_( i ) = 0.5 * currentNormalNorm;
        elementMomentArm_.col( i ) = currentCentroid - momentReferencePoint;
    }

    // Compute cross-sectional area
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea_( i ) =
                0.5 * elementSurfaceNormal_.row( i ).cwiseAbs( ).dot( elementSurfaceArea_ );
    }

}

//! Generate aerodynamic database.
void RarefiedFlowAnalysis::getSimulationConditions( )
{
    // Simulation boundary and grid
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries_( 2 * i ) = minimumDimensions_( i );
        simulationBoundaries_( 2 * i + 1 ) = maximumDimensions_( i );
    }
    simulationGrid_ = ( maximumDimensions_ - minimumDimensions_ ) / gridSpacing_;

    // Convert molecular speed ratio to stream velocity and compute simulation time step and ratio of real to simulated variables
    freeStreamVelocities_.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    simulationTimeStep_.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    ratioOfRealToSimulatedParticles_.resize( simulationAltitudes.size( ), 1 );
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {
            freeStreamVelocities_( h, s ) = simulationMolecularSpeedRatios.at( s ) * std::sqrt(
                        2.0 * atmosphericConditions_.at( 3 ).at( h ) * atmosphericConditions_.at( 2 ).at( h ) );
            simulationTimeStep_( h, s ) = 0.1 * ( maximumDimensions_( std::abs( referenceAxis_ ) ) -
                                                  minimumDimensions_( std::abs( referenceAxis_ ) ) ) /
                    freeStreamVelocities_( h, s );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles_( h ) = atmosphericConditions_.at( 5 ).at( h ) *
                std::pow( gridSpacing_, 3 ) / simulatedParticlesPerCell_;
    }
}

//! Generate aerodynamic coefficients at a single set of independent variables.
void RarefiedFlowAnalysis::runSPARTASimulation( std::string SPARTAExecutable )
{
    // Generate command string for SPARTA
    std::cout << "Initiating SPARTA simulation. This may take a while." << std::endl;
    std::string runSPARTACommandString = "cd " + input_output::getSPARTADataPath( ) + "; " +
            SPARTAExecutable + " -in " + inputFile_;
    // "mpirun -np " + std::to_string( numberOfCores ) + " " +

    // Predefine variables
    std::string anglesOfAttack;
    Eigen::Vector3d velocityVector;
    std::string temporaryOutputFile;
    std::vector< std::string > outputFileExtensions = { ".400", ".600", ".800", ".1000" };
    Eigen::Matrix< double, Eigen::Dynamic, 7 > outputMatrix;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanPressureValues;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanShearValues;

    // Allocate size of array
    aerodynamicCoefficients_( boost::extents[ simulationAltitudes.size( ) ][
                simulationMolecularSpeedRatios.size( ) ][ simulationAnglesOfAttack.size( ) ] );

    // Loop over simulation parameters and run SPARTA
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {
            // Get velocity vector
            velocityVector = Eigen::Vector3d::Zero( );
            velocityVector( std::abs( referenceAxis_ ) ) = ( std::signbit( referenceAxis_ ) ? 1.0 : -1.0 ) *
                    freeStreamVelocities_( h, s );

            // Get angles of attack string
            for ( double a : simulationAnglesOfAttack )
            {
                anglesOfAttack += printToStringWithPrecision( a, 0 ) + " ";
            }

            // Print to file
            FILE * fileIdentifier = std::fopen( inputFile_.c_str( ), "w" );
            std::fprintf( fileIdentifier, inputTemplate_.c_str( ), simulationBoundaries_( 0 ), simulationBoundaries_( 1 ),
                          simulationBoundaries_( 2 ), simulationBoundaries_( 3 ), simulationBoundaries_( 4 ),
                          simulationBoundaries_( 5 ), simulationGrid_( 0 ), simulationGrid_( 1 ), simulationGrid_( 2 ),
                          atmosphericConditions_.at( 5 ).at( h ), ratioOfRealToSimulatedParticles_( h ), simulationGases_.c_str( ),
                          simulationGases_.c_str( ), velocityVector( 0 ), velocityVector( 1 ), velocityVector( 2 ),
                          simulationGases_.c_str( ), atmosphericConditions_.at( 2 ).at( h ), anglesOfAttack.c_str( ),
                          wallTemperature_, accomodationCoefficient_, simulationTimeStep_( h, s ), outputDirectory_.c_str( ) );
            std::fclose( fileIdentifier );

            // Run SPARTA
            int systemStatus = std::system( runSPARTACommandString.c_str( ) );
            if ( systemStatus != 0 )
            {
                throw std::runtime_error( "Error: SPARTA simulation failed. See the log.sparta file for more details." );
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////            PROCESS RESULTS               ////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Loop over angles of attack
            meanPressureValues.resize( 3, numberOfTriangles_ );
            meanShearValues.resize( 3, numberOfTriangles_ );
            for ( unsigned int a = 0; a < simulationAnglesOfAttack.size( ); a++ )
            {
                // Get file name
                temporaryOutputFile = outputPath_ + "/" + printToStringWithPrecision(
                            simulationAnglesOfAttack.at( a ), 0 ) + ".coeff";

                // Read output files and compute mean pressure and shear force values
                meanPressureValues.setZero( );
                meanShearValues.setZero( );
                for ( unsigned int i = 0; i < outputFileExtensions.size( ); i++ )
                {
                    outputMatrix = readMatrixFromFile( temporaryOutputFile + outputFileExtensions.at( i ), "\t ;,", "%", 9 );
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        meanPressureValues.row( j ) += outputMatrix.col( j + 1 ).transpose( );
                        meanShearValues.row( j ) += outputMatrix.col( j + 4 ).transpose( );
                    }
               }
                meanPressureValues /= outputFileExtensions.size( );
                meanShearValues /= outputFileExtensions.size( );

                // Convert pressure and shear forces to coefficients
                aerodynamicCoefficients_[ h ][ s ][ a ] = computeAerodynamicCoefficientsFromPressureShear(
                            meanPressureValues,
                            meanShearValues,
                            atmosphericConditions_.at( 0 ).at( h ), // density
                            freeStreamVelocities_( h, s ),
                            atmosphericConditions_.at( 1 ).at( h ), // pressure
                            elementSurfaceNormal,
                            elementSurfaceArea,
                            elementMomentArm,
                            shapeCrossSectionalArea( std::abs( referenceAxis_ ) ) );
                std::cout << std::endl << "Altitude: " << simulationAltitudes.at( h ) << std::endl
                          << "Speed Ratio: " << simulationMolecularSpeedRatios.at( s ) << std::endl
                          << "Angle of Attack: " << simulationAnglesOfAttack.at( a ) << std::endl
                          << "Coefficients: " << aerodynamicCoefficients_[ h ][ s ][ a ].transpose( ) << std::endl;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Here." << std::endl;
}

//! Determine the pressure coefficients on a single vehicle part.
void RarefiedFlowAnalysis::determinePressureCoefficients(
        const int partNumber, const boost::array< int, 3 > independentVariableIndices )
{

}

//! Determine force coefficients from pressure coefficients.
Eigen::Vector3d RarefiedFlowAnalysis::calculateForceCoefficients(
        const int partNumber )
{

}

//! Determine moment coefficients from pressure coefficients.
Eigen::Vector3d RarefiedFlowAnalysis::calculateMomentCoefficients(
        const int partNumber )
{

}

//! Determines the inclination angle of panels on a single part.
void RarefiedFlowAnalysis::determineInclinations( const double angleOfAttack,
                                                  const double angleOfSideslip )
{

}

//! Determine compression pressure coefficients on all parts.
void RarefiedFlowAnalysis::updateCompressionPressures( const double machNumber,
                                                       const int partNumber )
{

}

//! Determines expansion pressure coefficients on all parts.
void RarefiedFlowAnalysis::updateExpansionPressures( const double machNumber,
                                                     const int partNumber )
{

}

} // namespace aerodynamics
} // namespace tudat
