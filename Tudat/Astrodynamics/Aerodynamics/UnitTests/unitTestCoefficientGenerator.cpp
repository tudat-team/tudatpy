/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120328    D. Dirkx          Updated code to use shared_ptrs instead of raw pointers.
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas Aircraft
 *        Aircraft Company, 1973.
 *
 */

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

//! Test coefficient generator.
int main( )
{
    using namespace tudat;
    using namespace tudat::mathematics::geometric_shapes;
    using tudat::aerodynamics::HypersonicLocalInclinationAnalysis;
    using std::fabs;
    using std::vector;
    using tudat::mathematics::PI;

    // Declare test variable.
    bool isCoefficientGeneratorErroneous = false;

    // Create test sphere.
    boost::shared_ptr< SphereSegment > sphere = boost::make_shared< SphereSegment >( 1.0 );
    boost::shared_ptr< bodies::VehicleExternalModel > externalModel
            = boost::make_shared< bodies::VehicleExternalModel >( );
    externalModel->setVehicleGeometry( sphere );
    bodies::Vehicle vehicle;
    vehicle.setExternalModel( externalModel );

    // Create analysis object.
    HypersonicLocalInclinationAnalysis analysis;

    // Set vehicle in analysis with 10,000 panels.
    vector< int > numberOfLines;
    vector< int > numberOfPoints;
    vector< bool > invertOrder;
    numberOfLines.resize( 1 );
    numberOfPoints.resize( 1 );
    invertOrder.resize( 1 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    invertOrder[ 0 ] = 0;
    analysis.setVehicle( vehicle, numberOfLines, numberOfPoints, invertOrder );

    // Set reference quantities.
    analysis.setReferenceArea( PI );
    analysis.setReferenceLength( 1.0 );
    analysis.setMomentReferencePoint( Eigen::Vector3d::Zero( ) );

    // Set pure Newtonian compression method for test purposes.
    analysis.setSelectedMethod( 0, 0, 0 );

    // Generate sphere database.
    analysis.generateDatabase( );

    // Allocate memory for independent variables to pass to analysis for retrieval.
    vector< int > independentVariables;
    independentVariables.resize( 3 );
    independentVariables[ 0 ] = 0;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::VectorXd aerodynamicCoefficients_;
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
                    isCoefficientGeneratorErroneous = true;
                }

                // Check if moment coefficients are approximately zero. Deviations
                // for pitch moment are greater due to greater range of angles of
                // attack than sideslip.
                if ( fabs( aerodynamicCoefficients_[ 3 ] ) > 1.0e-4 )
                {
                    std::cerr << "Error, sphere roll moment coefficient not zero."
                              << std::endl;
                    isCoefficientGeneratorErroneous = true;
                }

                if ( fabs( aerodynamicCoefficients_[ 4 ] ) > 1.0e-2 )
                {
                    std::cerr << "Error, sphere pitch moment coefficient not zero."
                              << std::endl;
                    isCoefficientGeneratorErroneous = true;
                }

                if ( fabs( aerodynamicCoefficients_[ 5 ] ) > 1.0e-2 )
                {
                    std::cerr << "Error, sphere yaw moment coefficient not zero."
                              << std::endl;
                    isCoefficientGeneratorErroneous = true;
                }
            }
        }
    }

    // Set Apollo capsule for validation.
    boost::shared_ptr< Capsule > capsule = boost::make_shared< Capsule >(
                4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

    externalModel->setVehicleGeometry( capsule );
    vehicle.setExternalModel( externalModel );

    // Declare new analysis object
    HypersonicLocalInclinationAnalysis analysis2 = HypersonicLocalInclinationAnalysis( );

    vector< int > numberOfLines2;
    vector< int > numberOfPoints2;
    vector< bool > invertOrders2;
    numberOfLines2.resize( 4 );
    numberOfPoints2.resize( 4 );
    invertOrders2.resize( 4 );

    // Set number of analysis points
    numberOfLines2[ 0 ] = 31;
    numberOfPoints2[ 0 ] = 31;
    numberOfLines2[ 1 ] = 31;
    numberOfPoints2[ 1 ] = 31;
    numberOfLines2[ 2 ] = 31;
    numberOfPoints2[ 2 ] = 10;
    numberOfLines2[ 3 ] = 11;
    numberOfPoints2[ 3 ] = 11;
    invertOrders2[ 0 ] = 1;
    invertOrders2[ 1 ] = 1;
    invertOrders2[ 2 ] = 1;
    invertOrders2[ 3 ] = 1;

    // Set capsule for analysis.
    analysis2.setVehicle( vehicle, numberOfLines2, numberOfPoints2, invertOrders2 );

    // Set reference quantities.
    analysis2.setReferenceArea( PI * pow( capsule->getMiddleRadius( ), 2.0 ) );
    analysis2.setReferenceLength( 3.9116 );
    Eigen::VectorXd momentReference = Eigen::VectorXd( 3 );
    momentReference( 0 ) = 0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;
    analysis2.setMomentReferencePoint( momentReference );

    // Set angle of attack analysis points.
    analysis2.setNumberOfAngleOfAttackPoints( 7 );
    int i;
    for ( i = 0; i < 7; i++ )
    {
        analysis2.setAngleOfAttackPoint( i, static_cast< double >( i - 6 ) * 5.0 * PI / 180.0 );
    }

    // Generate database.
    analysis2.generateDatabase( );

    // Retrieve coefficients at zero angle of attack for comparison.
    independentVariables[ 0 ] = analysis2.getNumberOfMachPoints( ) - 1;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;
    aerodynamicCoefficients_ = analysis2.getAerodynamicCoefficients( independentVariables );

    // Compare values to database values.
    if ( fabs( aerodynamicCoefficients_[ 0 ] - 1.51 ) > 0.1 )
    {
        std::cerr << "Error in Apollo drag coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 1 ] ) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Error in Apollo side force coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 2 ] ) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Error in Apollo normal force coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 3 ] ) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Error in Apollo roll moment coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 4 ] +0.052 ) > 0.01 )
    {
        std::cerr << "Error in Apollo pitch moment coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( fabs( aerodynamicCoefficients_[ 5 ] ) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Error in Apollo yaw moment coefficient." << std::endl;
        isCoefficientGeneratorErroneous = true;
    }

    if ( isCoefficientGeneratorErroneous )
    {
        std::cerr << "testCoefficientGenerator failed!" << std::endl;
    }

    return isCoefficientGeneratorErroneous;
}
