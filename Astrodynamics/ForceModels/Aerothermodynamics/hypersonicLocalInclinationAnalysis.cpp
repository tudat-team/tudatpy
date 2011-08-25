/*!   \file hypersonicLocalInclinationAnalysis.cpp
 *    This file contains the definition of the hypersonic local inclination
 *    analysis class.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : B.
 Romgens
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
#include "hypersonicLocalInclinationAnalysis.h"

// Using declarations.
using std::string;
using std::endl;

//! Default constructor.
HypersonicLocalInclinationAnalysis::HypersonicLocalInclinationAnalysis( )
{

    // Set number of independent variables.
    setNumberOfIndependentVariables( 3 );

    // Set indices of Mach number and angles of attack and sideslip in data
    // point arrays.
    machIndex_ = 0;
    angleOfAttackIndex_ = 1;
    angleOfSideslipIndex_ = 2;

    // Initialize inclinations, pressure coefficients, selected methods and
    // vehicle parts to NULL and set number of vehicle parts to 0.
    numberOfVehicleParts_ = 0;
    vehicleParts_ = NULL;
    selectedMethods_ = NULL;
    inclination_ = NULL;
    pressureCoefficient_ = NULL;

    // Set mach regime to default value.
    machRegime_ = "Full";

    //sets the ratio of specific heats, currently hard-coded, to be modified
    //in future.
    ratioOfSpecificHeats = 1.4;
}

//! Default destructor.
HypersonicLocalInclinationAnalysis::~HypersonicLocalInclinationAnalysis( )
{
    int i, j;

    // Deallocate selected methods.
    delete [ ] selectedMethods_[ 0 ];
    delete [ ] selectedMethods_[ 1 ];
    delete [ ] selectedMethods_[ 2 ];
    delete [ ] selectedMethods_[ 3 ];
    delete [ ] selectedMethods_;

    // Delete pressure coefficients and inclinations
    for( i = 0 ; i < numberOfVehicleParts_; i++)
    {
        for( j = 0 ; j < vehicleParts_[ i ].getNumberOfLines( ) - 1; j++)
        {
            delete [ ] inclination_[ i ][ j ];
            delete [ ] pressureCoefficient_[ i ][ j ];
        }
        delete [ ] inclination_[ i ];
        delete [ ] pressureCoefficient_[ i ];
    }
    delete [ ] inclination_;
    delete [ ] pressureCoefficient_;
    delete [ ] vehicleParts_;
}

//! Constructor to set geometry and reference quantities.
void HypersonicLocalInclinationAnalysis::setVehicle(
        Vehicle& vehicle,
        int* numberOfLines ,
        int* numberOfPoints,
        bool* invertOrders)
{
    int i;

    // Retrieve external surface geometry from vehicle.
    VehicleExternalModel* externalModel_ = vehicle.getPointerToExternalModel( );
    GeometricShape* surface_ = externalModel_->getVehicleExternalGeometry( );

    // Set geometry if it is a single surface.
    if ( dynamic_cast< SingleSurfaceGeometry* >( surface_ ) != NULL )
    {
        // Set number of geometries and allocate memory.
        numberOfVehicleParts_ = 1;
        vehicleParts_ = new LawgsPartGeometry[ 1 ];

        vehicleParts_->setReversalOperator( invertOrders[ 0 ] );

        // Convert geometry to LaWGS surface mesh and set in vehicleParts_ list.
        vehicleParts_[ 0 ].setMesh(
                dynamic_cast< SingleSurfaceGeometry* >( surface_ ),
                numberOfLines[ 0 ], numberOfPoints[ 0 ] );


    }
    // Set geometry if it is a composite surface.
    else if ( dynamic_cast< CompositeSurfaceGeometry* >( surface_ ) != NULL )
    {
        // Dynamic cast to composite surface geometry for further processing.
        CompositeSurfaceGeometry* compositeSurfaceGeometry_ =
                dynamic_cast< CompositeSurfaceGeometry* >( surface_ );

        // Set number of geometries and allocate memory.
        numberOfVehicleParts_ =
                compositeSurfaceGeometry_->getNumberOfSingleSurfaceGeometries( );
        vehicleParts_ = new LawgsPartGeometry[ numberOfVehicleParts_];

        // Iterate through all parts and set them in vehicleParts_ list.
        for( i = 0; i<numberOfVehicleParts_ ; i++)
        {
            // If part is not already a LaWGS part, convert it.
            if ( dynamic_cast< LawgsPartGeometry* >(
                compositeSurfaceGeometry_->getSingleSurfaceGeometry(i) ) ==
                 NULL )
            {
                vehicleParts_[ i ].setReversalOperator( invertOrders[ i ] );

                // Convert geometry to LaWGS and set in list.
                vehicleParts_[ i ].setMesh(
                        compositeSurfaceGeometry_->getSingleSurfaceGeometry(i),
                        numberOfLines[ i ],
                        numberOfPoints[ i ] );

            }
            // Else, set geometry directly.
            else
            {
                vehicleParts_[ i ] = *dynamic_cast< LawgsPartGeometry* >(
                   compositeSurfaceGeometry_->getSingleSurfaceGeometry( i ) );
            }
        }

    }
    // Allocate memory for arrays of pressure coefficients, inclinations
    // and methods.
    allocateArrays( );
}

//! Allocate pressure coefficient, inclination, and method
//! arrays.
void HypersonicLocalInclinationAnalysis::allocateArrays( )
{
    int i,j;

    // Allocate memory for panel inclinations and pressureCoefficient_.
    inclination_ = new double**[ numberOfVehicleParts_ ];
    pressureCoefficient_ = new double**[ numberOfVehicleParts_ ];
    for( i = 0 ; i < numberOfVehicleParts_ ; i++ )
    {
        inclination_[ i ] = new double*[ vehicleParts_[ i ].getNumberOfLines( ) ];
        pressureCoefficient_[ i ] =
                new double*[ vehicleParts_[ i ].getNumberOfLines( ) ];
        for ( j = 0 ; j<vehicleParts_[ i ].getNumberOfLines( ) ; j++ )
        {
            inclination_[ i ][ j ] =
                    new double[ vehicleParts_[ i ].getNumberOfPoints( ) ];
            pressureCoefficient_[ i ][ j ] =
                    new double[ vehicleParts_[ i ].getNumberOfPoints( ) ];
        }
    }

    // If methods have not yet been set by user, allocate memory and set
    // defaults.
    if( selectedMethods_ == NULL )
    {
        // Set memory for expansion and compression methods.
        selectedMethods_ = new int* [ 4 ];
        for( i = 0; i < 4 ; i++  )
        {
            // Set memory for Low and High hypersonic analysis methods for
            // each vehicle part.
            selectedMethods_[ i ] = new int[ 2 * numberOfVehicleParts_ ];
            for( j = 0; j < numberOfVehicleParts_ ; j++ )
            {
                // Sets all local inclination methods to a default of ( Modified
                // Newtonian for compression and Newtonian for expansion ).
                selectedMethods_[ i ][ j ] = 1;
            }
        }
    }
}

//! Get aerodynamic coefficients.
VectorXd HypersonicLocalInclinationAnalysis::getAerodynamicCoefficients(
        int* independentVariables )
{
    // If coefficients have not been allocated (and independent variables
    // have not been initialized), do so.
    if( vehicleCoefficients_ == NULL )
    {
        allocateVehicleCoefficients( );
    }

    // Declare and determine index in vehicleCoefficients_ array.
    int coefficientsIndex = variableIndicesToListIndex( independentVariables );

    // If coefficients for data point have not yet been calculated,
    // calculate them.
    if( vehicleCoefficients_[ coefficientsIndex ] == NULL )
    {
        determineVehicleCoefficients( independentVariables );
    }

    // Return requested coefficients.
    return *( vehicleCoefficients_[ coefficientsIndex ] );

}


//! Sets local inclination methods for all parts (expansion
//! and compression).
void HypersonicLocalInclinationAnalysis::setSelectedMethods(
        int** selectedMethods )
{
    int i,j;
    //For loops loop through input methods and set the analysis methods.
    for( i = 0; i < 4; i++ )
    {
        for( j = 0; j < numberOfVehicleParts_; j++ )
        {
            selectedMethods_[ i ][ j ] = selectedMethods[ i ][ j ];
        }
    }
}

//! Sets an analysis method on a single vehicle part.
void HypersonicLocalInclinationAnalysis::setSelectedMethod( const int& method,
                                                            const int& type,
                                                            const int& part)
{
    selectedMethods_[ type ][ part ] = method;
}

//! Generate aerodynamic database
void HypersonicLocalInclinationAnalysis::generateDatabase( )
{
    // If coefficients have not been allocated, do so.
    if( vehicleCoefficients_ == NULL )
    {
        allocateVehicleCoefficients( );
    }


    int i, j, k;
    int l = 0;

    // Allocate variable to pass to coefficient determination for independent
    // variable indices.
    int* independentVariableIndices = new int[ numberOfIndependentVariables_ ];

    // Iterate over all combinations of independent variables.
    for(i = 0 ; i < numberOfPointsPerIndependentVariables_[ machIndex_ ] ; i++ )
    {
        independentVariableIndices[ machIndex_ ] = i;
        for( j = 0 ; j < numberOfPointsPerIndependentVariables_[
                angleOfAttackIndex_ ] ; j++ )
        {
            independentVariableIndices[ angleOfAttackIndex_ ] = j;
            for( k = 0 ; k < numberOfPointsPerIndependentVariables_[
                    angleOfSideslipIndex_ ] ; k++ )
            {
                independentVariableIndices[ angleOfSideslipIndex_ ] = k;

                // If coefficient has not yet been set, calculate and set it.
                if( vehicleCoefficients_[ l ] == NULL )
                {
                    determineVehicleCoefficients( independentVariableIndices );
                }

                l++;
            }

        }
    }
    delete [ ] independentVariableIndices;
}

//! Allocates aerodynamic coefficient array and NULL independent
//! variables.
void HypersonicLocalInclinationAnalysis::allocateVehicleCoefficients( )
{
    // If angle of attack points have not yet been set, use defaults.
    if( dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ] == NULL )
    {
        setDefaultAngleOfAttackPoints( );
    }

    // If angle of sideslip points have not yet been set, use defaults.
    if( dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ] == NULL)
    {
        setDefaultAngleOfSideslipPoints( );
    }

    // If Mach number points have not yet been set, use defaults.
    if( dataPointsOfIndependentVariables_[ machIndex_ ] == NULL)
    {
        setDefaultMachPoints( );
    }

    // Determine number of combinations of independent variables.
    numberOfCases_ = numberOfPointsPerIndependentVariables_[ machIndex_ ] *
                        numberOfPointsPerIndependentVariables_[ angleOfSideslipIndex_ ] *
                        numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ];

    // Allocate memory for pointers to coefficients and initialize to NULL.
    vehicleCoefficients_ = new VectorXd*[ numberOfCases_ ];
    int i;
    for( i = 0; i < numberOfCases_ ; i++ )
    {
        vehicleCoefficients_[ i ] = NULL;
    }
}


//! Generates aerodynamic coefficients at a single set of independent
//! variables.
void HypersonicLocalInclinationAnalysis::determineVehicleCoefficients(
        int* independentVariableIndices )
{
    int i;

    // Determine index in vehicleCoefficients_ at which the coefficient
    // should be set.
    int coefficientsIndex = variableIndicesToListIndex(
            independentVariableIndices );

    // Declare coefficients vector and initialize to zeros.
    VectorXd coefficients = VectorXd( 6 );
    for ( i = 0 ; i < 6 ; i++ )
    {
        coefficients[ i ] = 0;
    }

    // Declare vehicle part coefficients variable.
    VectorXd partCoefficients;

    // Loop over all vehicle parts, calculate aerodynamic coefficients and add
    // to vehicleCoefficients_.
    for ( i = 0 ; i < numberOfVehicleParts_ ; i++ )
    {
        partCoefficients
                = determinePartCoefficients( i, independentVariableIndices );
        coefficients += partCoefficients;
    }

    // Allocate and set vehicle coefficients at given independent variables.
    vehicleCoefficients_[ coefficientsIndex ] = new VectorXd( coefficients );
}

//! Determines aerodynamic coefficients of a single vehicle part.
VectorXd HypersonicLocalInclinationAnalysis::determinePartCoefficients(
        const int& partNumber,
        int* independentVariableIndices )
{
    // Declare and determine angles of attack and sideslip for analysis.
    double angleOfAttack = dataPointsOfIndependentVariables_[
            angleOfAttackIndex_ ][ independentVariableIndices[
                    angleOfAttackIndex_ ] ];
    double angleOfSideslip = dataPointsOfIndependentVariables_[
            angleOfSideslipIndex_ ][ independentVariableIndices[
                    angleOfSideslipIndex_ ] ];

    // Declare partCoefficient vector.
    VectorXd partCoefficients = VectorXd( 6 );

    // Determine panel inclinations for part.
    determineInclination( partNumber, angleOfAttack, angleOfSideslip );

    // Set pressureCoefficient_ array for given independent variables.
    determinePressureCoefficients( partNumber, independentVariableIndices );

    // Calculate force coefficients from pressure coefficients.
    partCoefficients.start( 3 ) = calculateForceCoefficients( partNumber );

    // Calculate moment coefficients from pressure coefficients.
    partCoefficients.segment( 3, 3 ) = calculateMomentCoefficients(
            partNumber );

    return partCoefficients;
}


//! Determines the pressure coefficients on a single vehicle part
void HypersonicLocalInclinationAnalysis::determinePressureCoefficients(
        const int& partNumber,
        int* independentVariableIndices )
{
    // Retrieve Mach number.
    double machNumber = dataPointsOfIndependentVariables_[ machIndex_ ]
                        [ independentVariableIndices[ machIndex_ ] ];

    // Determine stagnation point pressure coefficients. Value is computed once
    // here to prevent its calculation in inner loop.
    stagnationPressureCoefficient = aerodynamics::computeStagnationPressure(
            machNumber, ratioOfSpecificHeats);
    updateCompressionPressures( machNumber, partNumber );
    updateExpansionPressures( machNumber, partNumber );
}

//! Determines force coefficients from pressure coefficients
VectorXd HypersonicLocalInclinationAnalysis::calculateForceCoefficients(
        const int& partNumber )
{
    int i, j;

    // Declare force coefficient vector and intialize to zeros.
    Vector3d forceCoefficients;
    for( i = 0 ; i < 3 ; i++)
    {
        forceCoefficients( i ) = 0.0;
    }
    
    // Loop over all panels and add pressures, scaled by panel area, to force 
    // coefficients.
    for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++ )
    {
        for( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
        {
            if( i == 0 && j == 0 )
            {
               //std::cout<<partNumber<<" "<<pressureCoefficient_[ partNumber ][ i ][ j ]<<"wtf"<<std::endl;
            }
            forceCoefficients = forceCoefficients -
                     pressureCoefficient_[ partNumber ][ i ][ j ] *
                     vehicleParts_[ partNumber ].getPanelArea( i, j ) *
                     vehicleParts_[ partNumber ].getPanelSurfaceNormal( i, j );
        }
    }

    // Normalize result by reference area.
    forceCoefficients /= referenceArea_;

    return forceCoefficients;
}



//! Determines moment coefficients from pressure coefficients
VectorXd HypersonicLocalInclinationAnalysis::calculateMomentCoefficients(
        const int& partNumber )
{
    int i, j;

    // Declare moment coefficient vector and intialize to zeros.
    Vector3d momentCoefficients ;
    for( i = 0 ; i < 3 ; i++ )
    {
        momentCoefficients( i ) = 0.0;
    }

    // Declare moment arm for panel moment determination.
    Vector3d referenceDistance ;

    // Loop over all panels and add moments due pressures.
    for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++ )
    {
        for( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
        {
            // Determine moment arm for given panel centroid.
            referenceDistance = ( vehicleParts_[ partNumber ].
                        getPanelCentroid( i, j ) -  momentReferencePoint_ );

            momentCoefficients = momentCoefficients -
                pressureCoefficient_[ partNumber ][ i ][ j ] *
                vehicleParts_[ partNumber ].getPanelArea( i, j ) *
                ( referenceDistance.cross( vehicleParts_[ partNumber ].
                                           getPanelSurfaceNormal( i, j ) ) );

        }
    }

    // Scale result by reference length and area.
    momentCoefficients /= ( referenceLength_ * referenceArea_ );

    return momentCoefficients;
}

//! Determines the inclination angle of panels on a single part.
void HypersonicLocalInclinationAnalysis::determineInclination(
        const int& partNumber,
        const double& angleOfAttack,
        const double& angleOfSideslip )
{
    int i, j;
    // Declare free-stream velocity vector.
    Vector3d freestreamVelocityDirection;

    // Set freestream velocity vector in body frame.
    double freestreamVelocityDirectionX = cos( angleOfAttack )*
                                          cos( angleOfSideslip );
    double freestreamVelocityDirectionY = sin( angleOfSideslip );
    double freestreamVelocityDirectionZ = sin( angleOfAttack ) *
                                          cos( angleOfSideslip );
    freestreamVelocityDirection( 0 ) = freestreamVelocityDirectionX;
    freestreamVelocityDirection( 1 ) = freestreamVelocityDirectionY;
    freestreamVelocityDirection( 2 ) = freestreamVelocityDirectionZ;

    // Declare cosine of inclination angle.
    double cosineOfInclination;

    // Loop over all panels of given vehicle part and set inclination angles.
    for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++ )
    {
        for( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
        {

            // Determine cosine of inclination angle from inner product between
            // surface normal and free-stream direction.
            cosineOfInclination = -1 * vehicleParts_[ partNumber ].
                                  getPanelSurfaceNormal( i, j ).
                                  dot( freestreamVelocityDirection );

            // Set inclination angle.
            inclination_[ partNumber ][ i ][ j ] = M_PI / 2 -
                                                   acos( cosineOfInclination );
        }
    }
}

//! Gets the number of vehicle parts.
int HypersonicLocalInclinationAnalysis::getNumberOfVehicleParts( )
{
    return numberOfVehicleParts_;
}

//! Determines compression pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateCompressionPressures(
        const double& machNumber,
        const int& partNumber)
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
    switch(method)
    {
    case 0:
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1; i++ )
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++ )
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
    }

}

//! Determines expansion pressure coefficients on all parts.
void HypersonicLocalInclinationAnalysis::updateExpansionPressures(
        const double& machNumber,
        const int& partNumber )
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
    switch(method)
    {
    case 0:
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
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
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
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
        // Option currently disabled
        break;
    case 3:
        // Calculate freestream Prandtl-Meyer function.
        freestreamPrandtlMeyerFunction = aerodynamics::computePrandtlMeyerFunction(
                machNumber, ratioOfSpecificHeats);
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Prandtl-Meyer expansion from freestream.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computePrandtlMeyerFreestreamPressureCoefficient(
                                inclination_[ partNumber ][ i ][ j ],
                                machNumber,
                                ratioOfSpecificHeats,
                                freestreamPrandtlMeyerFunction );
                }
            }
        }
        break;
    case 4:
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
            {
                // If panel inclination is negative, calculate pressure using
                // high Mach base pressure method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeHighMachBasePressure(
                            machNumber );
                }
            }
        }
        break;
    case 5:
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
            {
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    // If panel inclination is negative, calculate pressure using
                    // Van Dyke unified method.
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeVanDykeUnifiedPressureCoefficient(
                                    inclination_[ partNumber ][ i ][ j ],
                                    machNumber,
                                    ratioOfSpecificHeats,
                                    -1);

                }
            }
        }
        break;
    case 6:
        // Iterate over all panels on part.
        for( i = 0 ; i < vehicleParts_[ partNumber ].getNumberOfLines( ) - 1 ; i++)
        {
            for ( j = 0 ; j < vehicleParts_[ partNumber ].getNumberOfPoints( ) - 1 ; j++)
            {
                // If panel inclination is negative, calculate pressure using
                // ACM empirical method.
                if ( inclination_[ partNumber ][ i ][ j ] <= 0 )
                {
                    pressureCoefficient_[ partNumber ][ i ][ j ] =
                            aerodynamics::computeAcmEmpiricalPressureCoefficient(
                            inclination_[ partNumber ][ i ][ j ],
                            machNumber );
                }
            }
        }
        break;
    }

}

//! Sets the default Mach number points.
void HypersonicLocalInclinationAnalysis::setDefaultMachPoints( )
{
    // Set default points for full hypersonic analysis.
    if( machRegime_ == "Full" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 6;
        dataPointsOfIndependentVariables_[ machIndex_ ] =
            new double[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ];
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 3.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 4.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 4 ] = 10.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 5 ] = 20.0;
    }

    // Set default points for low hypersonic analysis.
    else if( machRegime_ == "Low" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 5;
        dataPointsOfIndependentVariables_[ machIndex_ ] =
            new double[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ];
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 3.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 4.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 4 ] = 10.0;
    }

    // Set default points for high hypersonic analysis.
    else if( machRegime_ == "High" )
    {
        numberOfPointsPerIndependentVariables_[ machIndex_ ] = 4;
        dataPointsOfIndependentVariables_[ machIndex_ ] =
            new double[ numberOfPointsPerIndependentVariables_[ machIndex_ ] ];
        dataPointsOfIndependentVariables_[ machIndex_ ][ 0 ] = 5.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 1 ] = 8.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 2 ] = 10.0;
        dataPointsOfIndependentVariables_[ machIndex_ ][ 3 ] = 20.0;
    }
}

//! Sets the default analysis points for angle of attack.
void HypersonicLocalInclinationAnalysis::setDefaultAngleOfAttackPoints( )
{
    //Set number of data points and allocate memory
    numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ] = 11;
    dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ] =
        new double[ numberOfPointsPerIndependentVariables_[
                angleOfAttackIndex_ ] ];

    // Set default values, 0 to 40 degrees, with steps of 5 degrees.
    int i;
    for( i = 0; i < numberOfPointsPerIndependentVariables_[ angleOfAttackIndex_ ]; i++ )
    {
        dataPointsOfIndependentVariables_[ angleOfAttackIndex_ ][ i ] =
                ( static_cast< double >( i ) * 5.0 * M_PI / 180.0 );
    }

}

//! Sets the default analysis points for angle of sideslip.
void HypersonicLocalInclinationAnalysis::setDefaultAngleOfSideslipPoints( )
{
    //Set number of data points and allocate memory
    numberOfPointsPerIndependentVariables_[ angleOfSideslipIndex_ ] = 2;
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ] =
            new double[ numberOfPointsPerIndependentVariables_[
                    angleOfSideslipIndex_ ] ];

    // Set default values, 0 and 1 degrees.
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ 0 ] = 0.0;
    dataPointsOfIndependentVariables_[ angleOfSideslipIndex_ ][ 1 ] =
            1.0 * M_PI / 180.0 ;
}

//! Gets a vehicle part.
LawgsPartGeometry HypersonicLocalInclinationAnalysis::getVehiclePart(
        const int& vehicleIndex )
{
    return vehicleParts_[ vehicleIndex ];
}

//! Function to set mach regime.
void HypersonicLocalInclinationAnalysis::setMachRegime(
        const std::string& machRegime )
{
    machRegime_ = machRegime;
}

//! Gets mach regime.
std::string HypersonicLocalInclinationAnalysis::getMachRegime( )
{
    return machRegime_;
}

//! Gets the vehicle name.
void HypersonicLocalInclinationAnalysis::setVehicleName(
        const std::string& vehicleName )
{
   vehicleName_ = vehicleName;
}

//! Gets the vehicle name.
std::string HypersonicLocalInclinationAnalysis::getVehicleName( )
{
    return vehicleName_;
}


//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          HypersonicLocalInclinationAnalysis&
                          hypersonicLocalInclinationAnalysis )
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

    for ( int i = 0;
          i < hypersonicLocalInclinationAnalysis.getNumberOfVehicleParts( );
          i++ )
    {
        stream << hypersonicLocalInclinationAnalysis
                  .getVehiclePart( i ).getName( ) << ", ";
    }

    stream << endl;

    // Return stream.
    return stream;
}

// End of file.
