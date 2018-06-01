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
 *      Klothakis, A. and Nikolos, I., “Modeling of Rarefied Hypersonic Flows Using the Massively
 *        Parallel DSMC Kernel “SPARTA”,” in 8th GRACM International Congress on Computational Mechanics,
 *        Volos, Greece, July 2015.
 *      Dirkx, D. and Mooij, E., Conceptual Shape Optimization of Entry Vehicles. Springer, 2017.
 *
 */

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Geometry>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/spartaDataReader.h"
#include "Tudat/InputOutput/spartaInputOutput.h"

#include "Tudat/Astrodynamics/Aerodynamics/rarefiedFlowAnalysis.h"

namespace tudat
{

namespace aerodynamics
{

using namespace unit_conversions;
using namespace input_output;

//! Returns default values of altitude for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAltitudePoints(
        const std::string& targetPlanet )
{
    std::vector< double > altitudePoints;

    // Set default points for Earth.
    if ( targetPlanet == "Earth" )
    {
        altitudePoints.resize( 5 );
        altitudePoints[ 0 ] = 225.0e3;
        altitudePoints[ 1 ] = 250.0e3;
        altitudePoints[ 2 ] = 300.0e3;
        altitudePoints[ 3 ] = 400.0e3;
        altitudePoints[ 4 ] = 600.0;
    }
    // Set default points for Mars.
    else if ( targetPlanet == "Mars" )
    {
        altitudePoints.resize( 5 );
        altitudePoints[ 0 ] = 125.0e3;
        altitudePoints[ 1 ] = 150.0e3;
        altitudePoints[ 2 ] = 200.0e3;
        altitudePoints[ 3 ] = 300.0e3;
        altitudePoints[ 4 ] = 500.0e3;
    }
    // Give error otherwise.
    else
    {
        throw std::runtime_error( "Error in altitude range selection for SPARTA simulation. "
                                  "Planet not supported." );
    }
    return altitudePoints;
}

//! Returns default values of Mach number for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowMachPoints(
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

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints(
        const std::string& angleOfAttackRegime )
{
    std::vector< double > angleOfAttackPoints;

    // Set default angles of attack
    double a = - 35;
    while ( a <= 35 )
    {
        angleOfAttackPoints.push_back( convertDegreesToRadians( a ) );
        a += 5;
    }

    // Add extra points if required
    if ( angleOfAttackRegime == "Full" )
    {
        std::vector< double > frontExtension = { convertDegreesToRadians( -85.0 ),
                                                 convertDegreesToRadians( -70.0 ),
                                                 convertDegreesToRadians( -55.0 ),
                                                 convertDegreesToRadians( -40.0 ) };
        std::vector< double > rearExtension = { convertDegreesToRadians( 40.0 ),
                                                convertDegreesToRadians( 55.0 ),
                                                convertDegreesToRadians( 70.0 ),
                                                convertDegreesToRadians( 85.0 ) };
        angleOfAttackPoints.insert( angleOfAttackPoints.begin( ), frontExtension.begin( ), frontExtension.end( ) );
        angleOfAttackPoints.insert( angleOfAttackPoints.end( ), rearExtension.begin( ), rearExtension.end( ) );
    }
    return angleOfAttackPoints;
}

//! Function to sort the rows of a matrix, based on the specified column and specified order.
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > sortMatrixRows(
        const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >& matrixToBeSorted,
        const int referenceColumn, const bool descendingOrder )
{
    // Declare eventual output vector
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > sortedMatrix;
    sortedMatrix.resizeLike( matrixToBeSorted );

    // Loop over rows and assign to new matrix
    Eigen::Matrix< double, 1, Eigen::Dynamic, Eigen::RowMajor > currentRow;
    for ( unsigned int i = 0; i < matrixToBeSorted.rows( ); i++ )
    {
        // Retrieve each row in 'scrambled order' and assign it to the new matrix
        currentRow = matrixToBeSorted.row( i );

        // Based on requested order
        if ( descendingOrder )
        {
            sortedMatrix.row( matrixToBeSorted.rows( ) - currentRow[ referenceColumn ] ) = currentRow;
            // no (- 1)'s are needed since both are defined starting from 1
        }
        else
        {
            sortedMatrix.row( currentRow[ referenceColumn ] - 1 ) = currentRow;
        }
    }

    // Give output
    return sortedMatrix;
}

//! Default constructor.
RarefiedFlowAnalysis::RarefiedFlowAnalysis(
        const std::string& SPARTAExecutable,
        const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
        boost::shared_ptr< TabulatedAtmosphere > atmosphereModel,
        const std::string& simulationGases,
        const std::string& geometryFileUser,
        const double referenceArea,
        const double referenceLength,
        const int referenceAxis,
        const Eigen::Vector3d& momentReferencePoint,
        const double gridSpacing,
        const double simulatedParticlesPerCell,
        const double wallTemperature,
        const double accommodationCoefficient,
        const bool printProgressInCommandWindow,
        const std::string MPIExecutable,
        const unsigned int numberOfCores ) :
    AerodynamicCoefficientGenerator< 3, 6 >(
          dataPointsOfIndependentVariables, referenceLength, referenceArea, referenceLength,
          momentReferencePoint,
          boost::assign::list_of( altitude_dependent )( mach_number_dependent )( angle_of_attack_dependent ),
          true, true ),
      SPARTAExecutable_( SPARTAExecutable ),simulationGases_( simulationGases ), referenceAxis_( referenceAxis ),
      gridSpacing_( gridSpacing ), simulatedParticlesPerCell_( simulatedParticlesPerCell ),
      wallTemperature_( wallTemperature ), accommodationCoefficient_( accommodationCoefficient ),
      printProgressInCommandWindow_( printProgressInCommandWindow ), MPIExecutable_( MPIExecutable ),
      numberOfCores_( numberOfCores ), referenceDimension_( static_cast< unsigned int >( referenceAxis_ ) )
{
    // Check inputs
    if ( referenceDimension_ > 2 )
    {
        throw std::runtime_error( "Error in SPARTA rarefied flow analysis. Reference axis makes "
                                  "no sense for a universe with 3 spacial dimensions (i.e., our universe). "
                                  "Note that the first dimension is identified with 0." );
    }

    // Analyze vehicle geometry
    analyzeGeometryFile( geometryFileUser );

    // Find atmospheric conditions based on altitude
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        atmosphericConditions_[ density_index ].push_back(
                    atmosphereModel->getDensity( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ pressure_index ].push_back(
                    atmosphereModel->getPressure( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ temperature_index ].push_back(
                    atmosphereModel->getTemperature( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ speed_of_sound_index ].push_back(
                    atmosphereModel->getSpeedOfSound( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ number_density_index ].push_back(
                    tudat::physical_constants::AVOGADRO_CONSTANT * atmosphericConditions_[ density_index ].at( h ) /
                    atmosphereModel->getMolarMass( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
    }

    // Get simulation conditions
    getSimulationConditions( );

    // Read SPARTA input template
    inputTemplate_ = readSpartaInputFileTemplate( );

    // Copy input shape file to default name
    std::string commandString = "cp " + geometryFileUser + " " + getSpartaInternalGeometryFile( );
    std::system( commandString.c_str( ) );

    // Run SPARTA simulation
    generateCoefficients( );

    // Create interpolator object
    createInterpolator( );
}

//! Open and read geometry file for SPARTA simulation.
void RarefiedFlowAnalysis::analyzeGeometryFile( const std::string& geometryFileUser )
{
    // Extract information on vehicle geometry
    std::pair< Eigen::Matrix< double, Eigen::Dynamic, 3 >, Eigen::Matrix< int, Eigen::Dynamic, 3 > >
            geometryData = readSpartaGeometryFile( geometryFileUser );
    shapePoints_ = geometryData.first;
    shapeTriangles_ = geometryData.second;
    numberOfPoints_ = shapePoints_.rows( );
    numberOfTriangles_ = shapeTriangles_.rows( );

    // Get maximum and minimum values in each dimension
    maximumDimensions_ = shapePoints_.colwise( ).maxCoeff( );
    minimumDimensions_ = shapePoints_.colwise( ).minCoeff( );

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
        elementMomentArm_.col( i ) = currentCentroid - momentReferencePoint_;
    }

    // Compute cross-sectional area
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea_( i ) =
                0.5 * elementSurfaceNormal_.row( i ).cwiseAbs( ).dot( elementSurfaceArea_ );
    }

    // Check consistency with input dimensions
    const double tolerance = 1e-5;
    if ( std::fabs( shapeCrossSectionalArea_( referenceDimension_ ) - referenceArea_ ) > tolerance )
    {
        throw std::runtime_error( "Error in SPARTA geometry file. Input reference area does not match the "
                                  "combination of reference axis and geometry. Note that the first dimension "
                                  "is identified with 0. Tolerance set to: " + std::to_string( tolerance ) );
    }
}

//! Retrieve simulation conditions based on input and geometry.
void RarefiedFlowAnalysis::getSimulationConditions( )
{
    // Simulation boundary and grid
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries_( 2 * i ) = minimumDimensions_( i ) + 0.5 * minimumDimensions_.minCoeff( ); // add extra space around shape
        simulationBoundaries_( 2 * i + 1 ) = maximumDimensions_( i ) + 0.5 * maximumDimensions_.maxCoeff( ); // add extra space around shape
        if ( i == referenceDimension_ )
        {
            simulationBoundaries_( 2 * i ) -= 1.0; // add extra space along axis of velocity
            simulationBoundaries_( 2 * i + 1 ) += 1.0; // add extra space along axis of velocity
        }
        simulationGrid_( i ) = simulationBoundaries_( 2 * i + 1 ) - simulationBoundaries_( 2 * i );
    }
    simulationGrid_ /= gridSpacing_;

    // Convert molecular speed ratio to stream velocity and compute simulation time step and ratio of real to simulated variables
    freeStreamVelocities_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), dataPointsOfIndependentVariables_.at( 1 ).size( ) );
    simulationTimeStep_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), dataPointsOfIndependentVariables_.at( 1 ).size( ) );
    ratioOfRealToSimulatedParticles_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), 1 );
    double simulationBoxLengthAlongReferenceAxis = ( simulationBoundaries_( 2 * referenceDimension_ + 1 ) -
                                                     simulationBoundaries_( 2 * referenceDimension_ ) );
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        for ( unsigned int m = 0; m < dataPointsOfIndependentVariables_.at( 1 ).size( ); m++ )
        {
            freeStreamVelocities_( h, m ) = dataPointsOfIndependentVariables_.at( 1 ).at( m ) *
                    atmosphericConditions_[ speed_of_sound_index ].at( h );
            simulationTimeStep_( h, m ) = 0.1 * simulationBoxLengthAlongReferenceAxis /
                    freeStreamVelocities_( h, m );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles_( h ) = atmosphericConditions_[ number_density_index ].at( h ) *
                std::pow( gridSpacing_, 3 ) / simulatedParticlesPerCell_;
    }
}

//! Generate aerodynamic database.
void RarefiedFlowAnalysis::generateCoefficients( )
{
    // Inform user on progress
    std::cout << "Initiating SPARTA simulation. This may take a while." << std::endl;

    // Generate command string for SPARTA
    std::string runSPARTACommandString = "cd " + getSpartaDataPath( ) + "; ";
    if ( MPIExecutable_ != "" )
    {
        if ( numberOfCores_ < 1 )
        {
            throw std::runtime_error( "Error in SPARTA rarefied flow analysis. Number of cores needs to be "
                                      "an integer value larger or equal to one." );
        }
        runSPARTACommandString = runSPARTACommandString + MPIExecutable_ +
                " -np " + std::to_string( numberOfCores_ ) + " ";
    }
    runSPARTACommandString = runSPARTACommandString + SPARTAExecutable_ + " -echo log ";
    if ( !printProgressInCommandWindow_ )
    {
        runSPARTACommandString = runSPARTACommandString + "-screen none ";
    }
    runSPARTACommandString = runSPARTACommandString + "-in " + getSpartaInputFile( );

    // Predefine variables
    Eigen::Vector3d velocityVector;
    std::string temporaryOutputFile = getSpartaOutputPath( ) + "/coeff";
    std::vector< std::string > outputFileExtensions = { ".1000" };//{ ".400", ".600", ".800", ".1000" };
    Eigen::Matrix< double, Eigen::Dynamic, 7, Eigen::RowMajor > outputMatrix;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanPressureValues;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanShearValues;

    // Loop over simulation parameters and run SPARTA
    int systemStatus;
    // Loop over altitude
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        // Inform user on progress
        std::cout << std::endl << "Altitude: "
                  << dataPointsOfIndependentVariables_.at( 0 ).at( h ) / 1e3
                  << " km" << std::endl;

        // Loop over Mach numbers
        for ( unsigned int m = 0; m < dataPointsOfIndependentVariables_.at( 1 ).size( ); m++ )
        {
            // Inform user on progress
            std::cout << "Mach number: "
                      << dataPointsOfIndependentVariables_.at( 1 ).at( m )
                      << std::endl;

            // Loop over angles of attack
            for ( unsigned int a = 0; a < dataPointsOfIndependentVariables_.at( 2 ).size( ); a++ )
            {
                // Inform user on progress
                std::cout << "Angle of attack: "
                          << convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) )
                          << " deg" << std::endl;

                // Get velocity vector
                velocityVector = Eigen::Vector3d::Zero( );
                velocityVector( referenceDimension_ ) = ( ( referenceAxis_ >= 0 ) ? - 1.0 : 1.0 ) *
                        freeStreamVelocities_( h, m );

                // Print to file
                FILE * fileIdentifier = std::fopen( getSpartaInputFile( ).c_str( ), "w" );
                std::fprintf( fileIdentifier, inputTemplate_.c_str( ),
                              simulationBoundaries_( 0 ), simulationBoundaries_( 1 ), simulationBoundaries_( 2 ),
                              simulationBoundaries_( 3 ), simulationBoundaries_( 4 ), simulationBoundaries_( 5 ),
                              simulationGrid_( 0 ), simulationGrid_( 1 ), simulationGrid_( 2 ),
                              atmosphericConditions_[ number_density_index ].at( h ), ratioOfRealToSimulatedParticles_( h ),
                              simulationGases_.c_str( ),
                              simulationGases_.c_str( ), velocityVector( 0 ), velocityVector( 1 ), velocityVector( 2 ),
                              simulationGases_.c_str( ), atmosphericConditions_[ temperature_index ].at( h ),
                              convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) ),
                              wallTemperature_, accommodationCoefficient_,
                              simulationTimeStep_( h, m ),
                              getSpartaOutputDirectory( ).c_str( ) );
                std::fclose( fileIdentifier );

                // Run SPARTA
                systemStatus = std::system( runSPARTACommandString.c_str( ) );
                if ( systemStatus != 0 )
                {
                    throw std::runtime_error( "Error in SPARTA rarefied flow analysis. "
                                              "SPARTA Simulation failed. See the log.sparta file in "
                                              "Tudat/External/SPARTA/ for more details." );
                }

                // Loop over angles of attack
                meanPressureValues.resize( 3, numberOfTriangles_ );
                meanShearValues.resize( 3, numberOfTriangles_ );

                // Read output files and compute mean pressure and shear force values
                meanPressureValues.setZero( );
                meanShearValues.setZero( );
                for ( unsigned int i = 0; i < outputFileExtensions.size( ); i++ )
                {
                    outputMatrix = readMatrixFromFile( temporaryOutputFile + outputFileExtensions.at( i ), "\t ;,", "%", 9 );
                    outputMatrix = sortMatrixRows( outputMatrix, 0 );
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        meanPressureValues.row( j ) += outputMatrix.col( j + 1 ).transpose( );
                        meanShearValues.row( j ) += outputMatrix.col( j + 4 ).transpose( );
                    }
                }
                meanPressureValues /= outputFileExtensions.size( );
                meanShearValues /= outputFileExtensions.size( );

                // Convert pressure and shear forces to aerodynamic coefficients
                aerodynamicCoefficients_[ h ][ m ][ a ] = computeAerodynamicCoefficientsFromPressureShearForces(
                            meanPressureValues,
                            meanShearValues,
                            atmosphericConditions_[ density_index ].at( h ),
                            atmosphericConditions_[ pressure_index ].at( h ),
                            freeStreamVelocities_( h, m ),
                            elementSurfaceNormal_,
                            elementSurfaceArea_,
                            elementMomentArm_,
                            referenceArea_,
                            referenceLength_ );

                // Clean up results folder
                std::string commandString = "rm " + getSpartaOutputPath( ) + "/coeff.*";
                std::system( commandString.c_str( ) );
            }
        }
    }

    // Inform user on progress
    std::cout << std::endl << "SPARTA simulation complete." << std::endl << std::endl;
}

//! Get aerodynamic coefficients at specific conditions.
Eigen::Vector6d RarefiedFlowAnalysis::getAerodynamicCoefficientsDataPoint(
        const boost::array< int, 3 > independentVariables )
{
    // Return requested coefficients.
    return aerodynamicCoefficients_( independentVariables );
}

} // namespace aerodynamics

} // namespace tudat
