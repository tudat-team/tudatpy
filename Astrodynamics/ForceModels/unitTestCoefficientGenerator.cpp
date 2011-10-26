/*! \file unitTestCoefficientGenerator.cpp
 *    This file contains the unit test of the aerodynamic coefficient generator in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : B.Romgens@student.tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas Aircraft
 *        Aircraft Company, 1973.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

// Include statements.
#include "Astrodynamics/ForceModels/hypersonicLocalInclinationAnalysis.h"
#include "Mathematics/GeometricShapes/capsule.h"
#include "Mathematics/GeometricShapes/sphereSegment.h"
#include "Output/writingOutputToFile.h"

//! Test coefficient generator.
int main( )
{
    // Declare test variable.
    bool isCoefficientGeneratorBad = 0;

    // Create test sphere.
    SphereSegment sphere = SphereSegment( );
    sphere.setRadius( 1.0 );
    sphere.setMinimumAzimuthAngle(  0.0 );
    sphere.setMaximumAzimuthAngle( 2.0 * M_PI );
    sphere.setMinimumZenithAngle( 0.0 );
    sphere.setMaximumZenithAngle( M_PI );
    VehicleExternalModel externalModel = VehicleExternalModel( );
    externalModel.setVehicleGeometry( sphere );
    Vehicle vehicle = Vehicle( );
    vehicle.setExternalModel( externalModel );

    // Create analysis object.
    HypersonicLocalInclinationAnalysis analysis = HypersonicLocalInclinationAnalysis( );

    // Set vehicle in analysis with 10,000 panels.
    int* numberOfLines = new int[ 1 ];
    int* numberOfPoints = new int[ 1 ];
    bool* invertOrder = new bool[ 1 ];
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    invertOrder[ 0 ] = 0;
    analysis.setVehicle( vehicle, numberOfLines, numberOfPoints, invertOrder );

    // Set reference quantities.
    analysis.setReferenceArea( M_PI );
    analysis.setReferenceLength( 1.0 );
    analysis.setMomentReferencePoint( Vector3d::Zero( ) );

    // Set pure Newtonian compression method for test purposes.
    analysis.setSelectedMethod( 0, 0, 0 );

    // Generate sphere database.
    analysis.generateDatabase( );

    // Allocate memory for independent variables to pass to analysis for retrieval.
    int* independentVariables = new int[ 3 ];
    independentVariables[ 0 ] = 0;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    VectorXd aerodynamicCoefficients_;
    double forceCoefficient_;

    // Iterate over all angles of attack to verify sphere coefficients.
    for ( int i = 0; i < analysis.getNumberOfMachPoints( ); i++ )
    {
        independentVariables[ 0 ] = i;

        for ( int j = 0; j < analysis.getNumberOfAngleOfAttackPoints( ); j++ )
        {
            independentVariables[ 1 ] = j;

            for ( int k = 0; k < analysis.getNumberOfAngleOfSideslipPoints( ); k++ )
            {
                independentVariables[ 2 ] = k;

                // Retrieve aerodynamic coefficients.
                aerodynamicCoefficients_ = analysis.getAerodynamicCoefficients(
                        independentVariables );
                forceCoefficient_ = sqrt( aerodynamicCoefficients_.x( )
                                        * aerodynamicCoefficients_.x( )
                                        + aerodynamicCoefficients_.y( )
                                        * aerodynamicCoefficients_.y( )
                                        + aerodynamicCoefficients_.z( )
                                        * aerodynamicCoefficients_.z( ) );

                // Check if 'total' aerodynamic coefficient is always
                // sufficiently close to zero.
                if ( fabs( forceCoefficient_ - 1.0 ) > 1.0e-2 )
                {
                    std::cerr << "Total magnitude of aerodynamic force wrong in sphere."
                              << std::endl;
                    isCoefficientGeneratorBad = true;
                }

                // Check if moment coefficients are approximately zero. Deviations
                // for pitch moment are greater due to greater range of angles of
                // attack than sideslip.
                if ( fabs( aerodynamicCoefficients_[ 3 ] ) > 1.0e-4 )
                {
                    std::cerr<<" Error, sphere roll moment coefficient not zero. "
                            <<std::endl;
                    isCoefficientGeneratorBad = true;
                }
                if ( fabs( aerodynamicCoefficients_[ 4 ] ) > 1.0e-2 )
                {
                    std::cerr << " Error, sphere pitch moment coefficient not zero. "
                              << std::endl;
                    isCoefficientGeneratorBad = true;
                }
                if ( fabs( aerodynamicCoefficients_[ 5 ] ) > 1.0E-2 )
                {
                    std::cerr << " Error, sphere yaw moment coefficient not zero. "
                              << std::endl;
                    isCoefficientGeneratorBad = true;
                }
            }
        }
    }


    // Set Apollo capsule for validation.
    Capsule capsule = Capsule( );
    capsule.setNoseRadius( 4.694 );
    capsule.setMiddleRadius( 1.956 );
    capsule.setRearAngle( -1.0 * 33.0 * M_PI / 180.0 );
    capsule.setRearLength( 2.662 );
    capsule.setSideRadius( 0.196 );
    capsule.setCapsule( );
    externalModel.setVehicleGeometry( capsule );
    vehicle.setExternalModel( externalModel );

    // Declare new analysis object
    HypersonicLocalInclinationAnalysis analysis2 = HypersonicLocalInclinationAnalysis( );

    int * numberOfLines2 = new int[ 4 ];
    int * numberOfPoints2 = new int[ 4 ];
    bool * invertOrders2 = new bool[ 4 ];

    // Set number of analysis points
    numberOfLines2[ 0 ] = 31;
    numberOfPoints2[ 0 ] = 31;
    numberOfLines2[ 1 ] = 31;
    numberOfPoints2[ 1 ] = 31;
    numberOfLines2[ 2 ] = 31;
    numberOfPoints2[ 2 ] = 2;
    numberOfLines2[ 3 ] = 11;
    numberOfPoints2[ 3 ] = 11;
    invertOrders2[ 0 ] = 1;
    invertOrders2[ 1 ] = 1;
    invertOrders2[ 2 ] = 1;
    invertOrders2[ 3 ] = 1;

    // Set capsule for analysis.
    analysis2.setVehicle( vehicle, numberOfLines2, numberOfPoints2, invertOrders2 );

    // Set reference quantities.
    analysis2.setReferenceArea( M_PI * pow( capsule.getMiddleRadius( ), 2.0 ) );
    analysis2.setReferenceLength( 3.9116 );
    VectorXd momentReference = VectorXd( 3 );
    momentReference( 0 ) = 0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;
    analysis2.setMomentReferencePoint( momentReference );

    // Set angle of attack analysis points.
    analysis2.setNumberOfAngleOfAttackPoints( 7 );
    int i;
    for ( i = 0; i < 7; i++ )
    {
        analysis2.setAngleOfAttackPoint( i, static_cast< double >( i - 6 )
                                         * 5.0 * M_PI / 180.0 );
    }

    // Generate database.
    analysis2.generateDatabase( );

    // Retrieve coefficients at zero angle of attack for comparison.
    independentVariables[ 0 ] = analysis2.getNumberOfMachPoints( ) - 1;
    independentVariables[ 1 ] = analysis2.getNumberOfMachPoints( ) - 1;
    independentVariables[ 2 ] = 0;
    aerodynamicCoefficients_ = analysis2.getAerodynamicCoefficients( independentVariables );

    // Compare values to database values.
    if ( fabs( aerodynamicCoefficients_[ 0 ] - 1.51 ) > 0.1 )
    {
        std::cerr << " Error in Apollo drag coefficient." << std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( aerodynamicCoefficients_[ 1 ] > 1.0E-15 )
    {
        std::cerr << " Error in Apollo side force coefficient." << std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( aerodynamicCoefficients_[ 2 ] > 1.0E-15 )
    {
        std::cerr << " Error in Apollo normal force coefficient." << std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( aerodynamicCoefficients_[ 3 ] > 1.0E-15 )
    {
        std::cerr<<" Error in Apollo roll moment coefficient."<<std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 4 ] +0.052 ) > 0.01 )
    {
        std::cerr << " Error in Apollo pitch moment coefficient." << std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( aerodynamicCoefficients_[ 5 ] > 1.0E-15 )
    {
        std::cerr << " Error in Apollo yaw moment coefficient." << std::endl;
        isCoefficientGeneratorBad = true;
    }

    if ( isCoefficientGeneratorBad )
    {
        std::cerr << "testCoefficientGenerator failed!" << std::endl;
    }

    return isCoefficientGeneratorBad;
}

// End of file.
