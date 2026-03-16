#include <cstdlib>
#include <iostream>
#include <string>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/astro/gravitation/mutualForcePotential.h"

namespace tudat
{

namespace gravitation
{

namespace
{

bool isFullTwoBodyTorqueDebugEnabled( )
{
    static const bool isEnabled = [ ]( )
    {
        const char* flag = std::getenv( "TUDAT_DEBUG_FULL_TWO_BODY_TORQUE" );
        return ( flag != nullptr && std::string( flag ) != "0" );
    }( );
    return isEnabled;
}

double getExpectedC20DegreeTwoMultiplierFromDocument( const int absoluteOrder )
{
    switch( absoluteOrder )
    {
    case 0:
        return 10.0;
    case 1:
        return std::sqrt( 125.0 / 6.0 );
    case 2:
        return std::sqrt( 125.0 / 12.0 );
    default:
        return TUDAT_NAN;
    }
}

std::string debugStatusFromDifference( const double difference, const double tolerance )
{
    return ( std::fabs( difference ) <= tolerance ? "OK" : "MISMATCH" );
}

std::string debugStatusFromVectorDifference(
        const Eigen::Vector3d& difference, const double tolerance )
{
    return ( difference.norm( ) <= tolerance ? "OK" : "MISMATCH" );
}

bool isApproximatelyIdentityQuaternion( const Eigen::Quaterniond& quaternion, const double tolerance )
{
    const Eigen::Quaterniond identityQuaternion = Eigen::Quaterniond::Identity( );
    const double positiveDifferenceNorm = ( quaternion.coeffs( ) - identityQuaternion.coeffs( ) ).norm( );
    const double negativeDifferenceNorm = ( quaternion.coeffs( ) + identityQuaternion.coeffs( ) ).norm( );
    return ( std::min( positiveDifferenceNorm, negativeDifferenceNorm ) <= tolerance );
}

bool isBody1C20Only(
        const Eigen::MatrixXd& cosineCoefficientsOfBody1,
        const Eigen::MatrixXd& sineCoefficientsOfBody1,
        const double tolerance = 1.0E-30 )
{
    if( cosineCoefficientsOfBody1.rows( ) <= 2 || cosineCoefficientsOfBody1.cols( ) <= 0 ||
        sineCoefficientsOfBody1.rows( ) <= 2 || sineCoefficientsOfBody1.cols( ) <= 0 )
    {
        return false;
    }

    for( int row = 0; row < cosineCoefficientsOfBody1.rows( ); row++ )
    {
        for( int col = 0; col < cosineCoefficientsOfBody1.cols( ); col++ )
        {
            if( row == 0 && col == 0 )
            {
                continue;
            }
            if( row == 2 && col == 0 )
            {
                continue;
            }
            if( std::fabs( cosineCoefficientsOfBody1( row, col ) ) > tolerance )
            {
                return false;
            }
        }
    }

    for( int row = 0; row < sineCoefficientsOfBody1.rows( ); row++ )
    {
        for( int col = 0; col < sineCoefficientsOfBody1.cols( ); col++ )
        {
            if( std::fabs( sineCoefficientsOfBody1( row, col ) ) > tolerance )
            {
                return false;
            }
        }
    }
    return true;
}

}

//! Function to get maximum degrees of used for the spherical harmonic expansions of the two bodies
std::pair< int, int > getMaximumDegrees(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    int maximumDegree1 = 0;
    int maximumDegree2 = 0;

    unsigned int degreeOfBody1, degreeOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));

        if( degreeOfBody1 > maximumDegree1 )
        {
            maximumDegree1 = degreeOfBody1;
        }

        if( degreeOfBody2 > maximumDegree2 )
        {
            maximumDegree2 = degreeOfBody2;
        }
    }
    return std::make_pair( maximumDegree1, maximumDegree2 );
}

//! Function to compute cross-body normalization terms for mutual two-body potential
double getGammaCoefficientForMutualForcePotential(
        const int l, const int m, const int j, const int k )
{
    double gammaCoefficient = 0.0;
    
    if( ( l == 0 && m == 0 ) || ( j == 0 && k == 0 ) )
    {
        gammaCoefficient = 1.0 / std::sqrt( 4.0 * mathematical_constants::PI );
    }
    else
    {
        gammaCoefficient = std::sqrt(
                    ( 2.0 * double( l ) + 1.0 ) * ( 2.0 * double( j ) + 1 ) *
                    boost::math::factorial< double >( l + j - m - k ) * boost::math::factorial< double >( l + j + m + k ) /
                    ( boost::math::factorial< double >( l + m ) * boost::math::factorial< double >( l - m ) *
                      boost::math::factorial< double >( j + k ) * boost::math::factorial< double >( j - k ) * 4.0 * mathematical_constants::PI *
                      double( 2 * l + 2 * j + 1 ) ) );
    }
    return gammaCoefficient;
}

//! Function to compute cross-body normalization terms for mutual two-body potential, for unnormalized or fully normalized
//! coefficients
double getMutualPotentialEffectiveCoefficientMultiplier(
        const int degree1, const int order1, const int degree2, const int order2, const bool areCoefficientsNormalized )
{
    double multiplier;
    if( areCoefficientsNormalized )
    {
        multiplier =
                getGammaCoefficientForMutualForcePotential( degree1, order1, degree2, order2 ) *
                std::sqrt( 4.0 * mathematical_constants::PI * ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 )  ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 )  ) /
                           ( ( order1 == 0 && order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) ) *
                getSigmaSignFunction( order1 ) * getSigmaSignFunction( order2 ) * getSigmaSignFunction( order1 + order2 ) *
                ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) *
                ( ( degree1 % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) );
    }
    else
    {
        multiplier =
                boost::math::factorial< double >( degree1 + degree2 - std::abs( order1 + order2 ) ) /
                ( boost::math::factorial< double >( degree2 - std::abs( order2 ) ) *
                  boost::math::factorial< double >( degree1 - std::abs( order1 ) ) ) *
                getSigmaSignFunction( order1 ) * getSigmaSignFunction( order2 ) * getSigmaSignFunction( order1 + order2 ) *
                ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) *
                ( ( degree1 % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) );
    }

    if( isFullTwoBodyTorqueDebugEnabled( ) && degree1 == 2 && std::abs( order1 ) == 0 && degree2 == 2 &&
        std::abs( order2 ) <= 2 && areCoefficientsNormalized )
    {
        const double expectedMultiplier = getExpectedC20DegreeTwoMultiplierFromDocument( std::abs( order2 ) );
        const double difference = multiplier - expectedMultiplier;
        std::cout << "[FTB-DBG][STEP multiplier]"
                  << " (l1,m1,l2,m2)=(" << degree1 << "," << order1 << "," << degree2 << "," << order2 << ")"
                  << " actual=" << multiplier
                  << " expected=" << expectedMultiplier
                  << " diff=" << difference
                  << " status=" << debugStatusFromDifference( difference, 1.0E-14 )
                  << std::endl;
    }

    return multiplier;
}

double computeSingleMutualForcePotentialTerm(
        const double effectiveCosineCoefficient,
        const double effectiveSineCoefficient,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const int degreeOfBody1,
        const int orderOfBody1,
        const int degreeOfBody2,
        const int orderOfBody2 )
{
    return ( effectiveCosineCoefficient * sphericalHarmonicsCache->getCosineOfMultipleLongitude(
                 std::abs( orderOfBody1 + orderOfBody2 ) ) -
             effectiveSineCoefficient * sphericalHarmonicsCache->getSineOfMultipleLongitude(
                 std::abs( orderOfBody1 + orderOfBody2 ) ) ) *
            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 + orderOfBody2 ) );
}

double computeMutualForcePotential(
        const Eigen::Vector3d& bodyFixedPosition,
        const double effectiveGravitationalParameterOfBody1,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const int maximumDegreeOfBody1,
        const int maximumDegreeOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    
    // Determine body fixed spherical position of body udnergoing acceleration.
    Eigen::Vector3d sphericalPositon =
            coordinate_conversions::convertCartesianToSpherical( bodyFixedPosition );
    double radius = sphericalPositon.x( );
    double latitude = mathematical_constants::PI / 2.0 - sphericalPositon.y( );
    double longitude = sphericalPositon.z( );
    
    double sineOfLatitude = std::sin( latitude );
    sphericalHarmonicsCache->update( TUDAT_NAN, sineOfLatitude, longitude, TUDAT_NAN );

    double potential = 0.0;
    
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;
    
    std::vector< double > radiusRatioOfBody1List;
    double radiusRatioOfBody1 = equatorialRadiusOfBody1 / radius;
    radiusRatioOfBody1List.push_back( 1 );
    for( unsigned int i = 1; i <= static_cast< unsigned int >( maximumDegreeOfBody1 ); i++ )
    {
        radiusRatioOfBody1List.push_back( radiusRatioOfBody1List.at( i - 1 ) * radiusRatioOfBody1 );
    }
    
    std::vector< double > radiusRatioOfBody2List;
    radiusRatioOfBody2List.push_back( 1 );
    double radiusRatioOfBody2 = equatorialRadiusOfBody2 / radius;
    for( unsigned int i = 1; i <= static_cast< unsigned int >( maximumDegreeOfBody2 ); i++ )
    {
        radiusRatioOfBody2List.push_back( radiusRatioOfBody2List.at( i - 1 ) * radiusRatioOfBody2 );
    }
    
    
    double currentTerm = 0;
    for(  unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));
        
        currentTerm = 0;
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    sphericalHarmonicsCache, degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    sphericalHarmonicsCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    sphericalHarmonicsCache, degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    sphericalHarmonicsCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );
        currentTerm *= radiusRatioOfBody1List.at( degreeOfBody1 );
        currentTerm *= radiusRatioOfBody2List.at( degreeOfBody2 );
        potential += currentTerm;
    }
    
    // Multiply by central term and return
    return potential * effectiveGravitationalParameterOfBody1 / radius;
}

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
Eigen::Vector3d computeGeodesyNormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const int maximumDegree1, const int maximumDegree2, const int maximumEvaluationDegree,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
//    std::cout<<"Computing acceleration: "<<std::endl;
    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration = Eigen::Vector3d::Zero( );

    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = coordinate_conversions::
            convertCartesianToCylindrical( positionOfBodySubjectToAcceleration );

    // Compute radius coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 0 )
            = std::sqrt( cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                         + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );

    // If radius coordinate is smaller than planetary radius...
    if ( sphericalpositionOfBodySubjectToAcceleration( 0 ) < ( equatorialRadiusOfBody1 + equatorialRadiusOfBody2 ) )
    {
        // ...throw runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }

    // If radius coordinate is zero set latitude coordinate to 90 degrees.
    if ( std::fabs( cylindricalCoordinates( 0 ) ) < std::numeric_limits< double >::epsilon( ) )
    {
        sphericalpositionOfBodySubjectToAcceleration( 1 ) = mathematical_constants::PI / 2.0;
    }

    // Else compute latitude coordinate.
    else
    {
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = std::atan( cylindricalCoordinates( 2 ) / cylindricalCoordinates( 0 ) );
    }

    // Compute longitude coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 2 ) = cylindricalCoordinates( 1 );
    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    double cosineOfAngle = std::cos( sphericalpositionOfBodySubjectToAcceleration( 1 ) );

    sphericalHarmonicsCache->update( TUDAT_NAN, sineOfAngle, sphericalpositionOfBodySubjectToAcceleration( 2 ), TUDAT_NAN );

    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );

    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree, totalOrder;
    int signedOrderOfBody1, signedOrderOfBody2;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );

    bool computeTerm;

    std::vector< std::pair< double, double > > legendreTerms;

    legendreTerms.resize( ( maximumEvaluationDegree + 1 ) * ( maximumEvaluationDegree + 1 ) );
    for( unsigned int i = 0; i <= maximumEvaluationDegree; i++ )
    {
        for( unsigned int j = 0; j <= i; j++ )
        {
            // Compute geodesy-normalized Legendre polynomials.
            const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        i, j );
            const double incrementedLegendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        i, j + 1 );

            // Compute geodesy-normalized Legendre polynomial derivative.
            const double legendrePolynomialDerivative =
                    basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                        i, j, sineOfAngle,
                        legendrePolynomial, incrementedLegendrePolynomial );

            legendreTerms[ i + ( maximumEvaluationDegree + 1 ) * j ] = std::make_pair( legendrePolynomial, legendrePolynomialDerivative );
        }
    }


    // Loop through all degrees.
    std::pair< double, double > currentTerms;
    for ( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));

        totalDegree = degreeOfBody1 + degreeOfBody2;

        equatorialRadiusRatioPower = radius1Powers[ degreeOfBody1 ] * radius2Powers[ degreeOfBody2 ];

        for( unsigned j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                signedOrderOfBody1 = orderOfBody1;
                signedOrderOfBody2 = orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                computeTerm = 1;
                break;
            case 1:
                signedOrderOfBody1 = -orderOfBody1;
                signedOrderOfBody2 = orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody1 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            case 2:
                signedOrderOfBody1 = orderOfBody1;
                signedOrderOfBody2 = -orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );

                break;
            case 3:
                signedOrderOfBody1 = -orderOfBody1;
                signedOrderOfBody2 = -orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody1 == 0 || orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            }

            if( computeTerm )
            {
//                std::cout<<"Computing  "<<j<<" "<<degreeOfBody1<<" "<<orderOfBody1<<" "<<degreeOfBody2<<" "<<orderOfBody2 <<" "<<
//                           sphericalGradient.transpose( )<<" "<<
//                           sphericalpositionOfBodySubjectToAcceleration.transpose( )<<" "<<preMultiplier<<" "<<
//                           totalDegree<<" "<<totalOrder<<" "<<
//                           effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<" "<<
//                           effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<" "<<
//                           currentTerms.first<<" "<<currentTerms.second
//                        <<std::endl;

                currentTerms = legendreTerms.at( totalDegree + ( maximumEvaluationDegree + 1 ) * totalOrder );
                const double effectiveCosineCoefficient = effectiveCosineCoefficientFunction(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double effectiveSineCoefficient = effectiveSineCoefficientFunction(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                // Compute the potential gradient of a single spherical harmonic term.
                const Eigen::Vector3d termContribution = basic_mathematics::computePotentialGradient(
                            sphericalpositionOfBodySubjectToAcceleration( 0 ),
                            equatorialRadiusRatioPower,
                            sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                            sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                            cosineOfAngle,
                            preMultiplier,
                            totalDegree,
                            totalOrder,
                            effectiveCosineCoefficient,
                            effectiveSineCoefficient,
                            currentTerms.first,
                            currentTerms.second );
                sphericalGradient += termContribution;

                if( isFullTwoBodyTorqueDebugEnabled( ) &&
                    degreeOfBody1 == 2 && orderOfBody1 == 0 && degreeOfBody2 == 2 && std::abs( orderOfBody2 ) <= 2 )
                {
                    const bool expectedZeroContribution =
                            ( std::fabs( effectiveCosineCoefficient ) <= 1.0E-30 &&
                              std::fabs( effectiveSineCoefficient ) <= 1.0E-30 );
                    const std::string status = expectedZeroContribution ?
                                debugStatusFromVectorDifference( termContribution, 1.0E-30 ) :
                                "N/A";

                    std::cout << "[FTB-DBG][STEP grad_term]"
                              << " base=(l1,m1,l2,m2)=(" << degreeOfBody1 << "," << orderOfBody1
                              << "," << degreeOfBody2 << "," << orderOfBody2 << ")"
                              << " signedOrders=(" << signedOrderOfBody1 << "," << signedOrderOfBody2 << ")"
                              << " j_case=" << j
                              << " totalOrder=" << totalOrder
                              << " Ceff=" << effectiveCosineCoefficient
                              << " Seff=" << effectiveSineCoefficient
                              << " P=" << currentTerms.first
                              << " dP=" << currentTerms.second
                              << " dU=" << termContribution.transpose( )
                              << " expectedZero=" << expectedZeroContribution
                              << " status=" << status
                              << std::endl;
                }
            }
        }


    }

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.

    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
Eigen::Vector3d computeUnnormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const std::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const std::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration = Eigen::Vector3d::Zero( );
    
    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = coordinate_conversions::
            convertCartesianToCylindrical( positionOfBodySubjectToAcceleration );
    
    // Compute radius coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 0 )
            = std::sqrt( cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                         + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );
    
    // If radius coordinate is smaller than planetary radius...
    if ( sphericalpositionOfBodySubjectToAcceleration( 0 ) < ( equatorialRadiusOfBody1 + equatorialRadiusOfBody2 ) )
    {
        // ...throw runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }
    
    // If radius coordinate is zero set latitude coordinate to 90 degrees.
    if ( std::fabs( cylindricalCoordinates( 0 ) ) < std::numeric_limits< double >::epsilon( ) )
    {
        sphericalpositionOfBodySubjectToAcceleration( 1 ) = mathematical_constants::PI / 2.0;
    }
    
    // Else compute latitude coordinate.
    else
    {
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = std::atan( cylindricalCoordinates( 2 ) / cylindricalCoordinates( 0 ) );
    }
    
    // Compute longitude coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 2 ) = cylindricalCoordinates( 1 );
    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    double cosineOfAngle = std::cos( sphericalpositionOfBodySubjectToAcceleration( 1 ) );

    sphericalHarmonicsCache->update( TUDAT_NAN, sineOfAngle, sphericalpositionOfBodySubjectToAcceleration( 2 ), TUDAT_NAN );

    
    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );
    
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree, totalOrder;
    int signedOrderOfBody1, signedOrderOfBody2;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );
    
    bool computeTerm;

    // Loop through all degrees.
    for ( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));
        
        totalDegree = degreeOfBody1 + degreeOfBody2;
        
        equatorialRadiusRatioPower =
                basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody1 / sphericalpositionOfBodySubjectToAcceleration( 0 ), degreeOfBody1 ) *
                basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody2 / sphericalpositionOfBodySubjectToAcceleration( 0 ), degreeOfBody2 );
        
        for( unsigned j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                signedOrderOfBody1 = orderOfBody1;
                signedOrderOfBody2 = orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                computeTerm = 1;
                break;
            case 1:
                signedOrderOfBody1 = -orderOfBody1;
                signedOrderOfBody2 = orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody1 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            case 2:
                signedOrderOfBody1 = orderOfBody1;
                signedOrderOfBody2 = -orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );

                break;
            case 3:
                signedOrderOfBody1 = -orderOfBody1;
                signedOrderOfBody2 = -orderOfBody2;
                totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                ( orderOfBody1 == 0 || orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            }
            
            if( computeTerm )
            {
                // Compute geodesy-normalized Legendre polynomials.
                const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                            totalDegree, totalOrder );
                const double incrementedLegendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                            totalDegree, totalOrder + 1  );
                
                // Compute geodesy-normalized Legendre polynomial derivative.
                const double legendrePolynomialDerivative =
                        basic_mathematics::computeLegendrePolynomialDerivative(
                            totalOrder, sineOfAngle,
                            legendrePolynomial, incrementedLegendrePolynomial );

                // Compute the potential gradient of a single spherical harmonic term.
                sphericalGradient += basic_mathematics::computePotentialGradient(
                            sphericalpositionOfBodySubjectToAcceleration( 0 ),
                            equatorialRadiusRatioPower,
                            sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                            sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                            cosineOfAngle,
                            preMultiplier,
                            totalDegree,
                            totalOrder,
                            effectiveCosineCoefficientFunction(
                                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 ),
                            effectiveSineCoefficientFunction(
                                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 ),
                            legendrePolynomial,
                            legendrePolynomialDerivative );
            }
        }


    }

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.

    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

void computePartialDerivativesOfPotentialComponentsWrtFullCoefficients(
        std::vector< Eigen::Matrix< double, 1, 2 > >& potentialComponentsWrtFullCoefficients,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const double distance,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::function< int( const int, const int, const int, const int )> effectiveIndexFunction )
{
    int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;
    int effectiveIndex;
    double equatorialRadiusRatioPower;

    Eigen::Matrix< double, 1, 2 > currentPotentialComponentWrtFullCoefficients;
    for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));

        equatorialRadiusRatioPower = radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 ) / distance;

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = sphericalHarmonicsCache->getCosineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = sphericalHarmonicsCache->getSineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 + orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = sphericalHarmonicsCache->getCosineOfMultipleLongitude( std::abs( -orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = sphericalHarmonicsCache->getSineOfMultipleLongitude( std::abs( -orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( -orderOfBody1 + orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = sphericalHarmonicsCache->getCosineOfMultipleLongitude( std::abs( orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = sphericalHarmonicsCache->getSineOfMultipleLongitude( std::abs( orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 - orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = sphericalHarmonicsCache->getCosineOfMultipleLongitude( std::abs( -orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = sphericalHarmonicsCache->getSineOfMultipleLongitude( std::abs( -orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( -orderOfBody1 - orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

    }
}


void EffectiveMutualSphericalHarmonicsField::getCurrentEffectiveCoefficients(
        const int degree1, const int order1, const int degree2, const int order2,
        const int effectiveIndex,
        double& cosineCoefficient, double& sineCoefficient )
{
    cosineCoefficient = ( cosineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                          transformedCosineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) -
                          ( ( order1 < 0 ) ? ( -1.0 ) : ( 1.0 ) ) * ( ( order2 < 0 ) ? ( -1.0 ) : ( 1.0 ) )  *
                          sineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                          transformedSineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) );
    sineCoefficient = ( ( ( order2 < 0 ) ? ( -1.0 ) : ( 1.0 ) ) *
                        cosineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                        transformedSineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) +
                        ( ( order1 < 0 ) ? ( -1.0 ) : ( 1.0 ) )  *
                        sineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                        transformedCosineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) );

    double currentMultiplier = multipliers_.at( effectiveIndex );
    cosineCoefficient *= currentMultiplier;
    sineCoefficient *= ( ( ( order1 +  order2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) ) * currentMultiplier;

    if( isFullTwoBodyTorqueDebugEnabled( ) &&
        degree1 == 2 && std::abs( order1 ) == 0 && degree2 == 2 && std::abs( order2 ) <= 2 &&
        areCoefficientsNormalized_ && isBody1C20Only( cosineCoefficientsOfBody1_, sineCoefficientsOfBody1_ ) )
    {
        const double c20Body1 = cosineCoefficientsOfBody1_( 2, 0 );
        const double expectedCosineCoefficient =
                currentMultiplier * c20Body1 * transformedCosineCoefficientsOfBody2_( degree2, std::abs( order2 ) );
        const double expectedSineCoefficient =
                currentMultiplier * c20Body1 * transformedSineCoefficientsOfBody2_( degree2, std::abs( order2 ) );
        const double cosineDifference = cosineCoefficient - expectedCosineCoefficient;
        const double sineDifference = sineCoefficient - expectedSineCoefficient;

        std::cout << "[FTB-DBG][STEP effective_coeff]"
                  << " (l1,m1,l2,m2)=(" << degree1 << "," << order1 << "," << degree2 << "," << order2 << ")"
                  << " Ceff=" << cosineCoefficient
                  << " Ceff_expected=" << expectedCosineCoefficient
                  << " Ceff_diff=" << cosineDifference
                  << " Ceff_status=" << debugStatusFromDifference( cosineDifference, 1.0E-15 )
                  << " Seff=" << sineCoefficient
                  << " Seff_expected=" << expectedSineCoefficient
                  << " Seff_diff=" << sineDifference
                  << " Seff_status=" << debugStatusFromDifference( sineDifference, 1.0E-15 )
                  << std::endl;
    }

}

void EffectiveMutualSphericalHarmonicsField::computeCurrentEffectiveCoefficients(
        const Eigen::Quaterniond coefficientRotationQuaterion )
{
    cosineCoefficientsOfBody1_ = cosineCoefficientFunctionOfBody1_( );
    sineCoefficientsOfBody1_ = sineCoefficientFunctionOfBody1_( );
    cosineCoefficientsOfBody2_ = cosineCoefficientFunctionOfBody2_( );
    sineCoefficientsOfBody2_ = sineCoefficientFunctionOfBody2_( );

    transformationCache_->updateFromQuaternion( coefficientRotationQuaterion );
    transformationCache_->transformCoefficientsAtDegree(
                cosineCoefficientsOfBody2_,
                sineCoefficientsOfBody2_,
                transformedCosineCoefficientsOfBody2_,
                transformedSineCoefficientsOfBody2_,
                areCoefficientsNormalized_ );

    if( isFullTwoBodyTorqueDebugEnabled( ) )
    {
        const bool isIdentityRotation = isApproximatelyIdentityQuaternion( coefficientRotationQuaterion, 1.0E-15 );
        if( isIdentityRotation && cosineCoefficientsOfBody2_.rows( ) > 2 && cosineCoefficientsOfBody2_.cols( ) > 2 &&
            sineCoefficientsOfBody2_.rows( ) > 2 && sineCoefficientsOfBody2_.cols( ) > 2 )
        {
            for( int order = 0; order <= 2; order++ )
            {
                const double cosineDifference =
                        transformedCosineCoefficientsOfBody2_( 2, order ) - cosineCoefficientsOfBody2_( 2, order );
                const double sineDifference =
                        transformedSineCoefficientsOfBody2_( 2, order ) - sineCoefficientsOfBody2_( 2, order );

                std::cout << "[FTB-DBG][STEP transformed_coeff]"
                          << " (l,m)=(2," << order << ")"
                          << " C_actual=" << transformedCosineCoefficientsOfBody2_( 2, order )
                          << " C_expected=" << cosineCoefficientsOfBody2_( 2, order )
                          << " C_diff=" << cosineDifference
                          << " C_status=" << debugStatusFromDifference( cosineDifference, 1.0E-15 )
                          << " S_actual=" << transformedSineCoefficientsOfBody2_( 2, order )
                          << " S_expected=" << sineCoefficientsOfBody2_( 2, order )
                          << " S_diff=" << sineDifference
                          << " S_status=" << debugStatusFromDifference( sineDifference, 1.0E-15 )
                          << std::endl;
            }
        }
        else
        {
            std::cout << "[FTB-DBG][STEP transformed_coeff] identity_check_skipped="
                      << ( isIdentityRotation ? 0 : 1 ) << std::endl;
        }
    }

    updateEffectiveMutualPotential( );
}

void EffectiveMutualSphericalHarmonicsField::computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
        const Eigen::MatrixXd& transformedCosineCoefficients,
        const Eigen::MatrixXd& transformedSineCoefficients )
{
    transformedCosineCoefficientsOfBody2_ = transformedCosineCoefficients;
    transformedSineCoefficientsOfBody2_ = transformedSineCoefficients;

    updateEffectiveMutualPotential( );
}

void EffectiveMutualSphericalHarmonicsField::updateEffectiveMutualPotential( )
{

    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;
    int effectiveIndex;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( i ));

        effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

        getCurrentEffectiveCoefficients(
                    degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2, effectiveIndex,
                    effectiveCosineCoefficients_[ effectiveIndex ],
                    effectiveSineCoefficients_[ effectiveIndex ] );

        if( orderOfBody1 != 0 )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );

            getCurrentEffectiveCoefficients(
                        degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2, effectiveIndex,
                        effectiveCosineCoefficients_[ effectiveIndex ],
                        effectiveSineCoefficients_[ effectiveIndex ] );
        }

        if( orderOfBody2 != 0 )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );

            getCurrentEffectiveCoefficients(
                        degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2, effectiveIndex,
                        effectiveCosineCoefficients_[ effectiveIndex ],
                        effectiveSineCoefficients_[ effectiveIndex ] );
        }

        if( !( orderOfBody1 == 0 || orderOfBody2 == 0 ) )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );

            getCurrentEffectiveCoefficients(
                        degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2, effectiveIndex,
                        effectiveCosineCoefficients_[ effectiveIndex ],
                        effectiveSineCoefficients_[ effectiveIndex ] );
        }
    }
}

void EffectiveMutualSphericalHarmonicsField::computePartialsOfFullCoefficientsWrtTransformedCoefficients(
        std::vector< Eigen::Matrix2d >& fullCoefficientsWrtBody2CoefficientsList )
{
    int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;

    int effectiveIndex;

    Eigen::Matrix2d currentPartial, fullCoefficientsWrtBody2Coefficients;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( i ));

        fullCoefficientsWrtBody2Coefficients( 0, 0 ) = cosineCoefficientsOfBody1_( degreeOfBody1, orderOfBody1 );
        fullCoefficientsWrtBody2Coefficients( 0, 1 ) = -sineCoefficientsOfBody1_( degreeOfBody1, orderOfBody1 );
        fullCoefficientsWrtBody2Coefficients( 1, 0 ) = sineCoefficientsOfBody1_( degreeOfBody1, orderOfBody1 );
        fullCoefficientsWrtBody2Coefficients( 1, 1 ) = cosineCoefficientsOfBody1_( degreeOfBody1, orderOfBody1 );

        effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        {
            currentPartial = fullCoefficientsWrtBody2Coefficients * multipliers_[ effectiveIndex ];
            fullCoefficientsWrtBody2CoefficientsList[ effectiveIndex ] = currentPartial;
        }

        if( orderOfBody1 != 0 )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );

            currentPartial = fullCoefficientsWrtBody2Coefficients * multipliers_[ effectiveIndex ];
            currentPartial( 1, 0 ) *= -1.0;
            currentPartial( 0, 1 ) *= -1.0;

            if( orderOfBody1 > orderOfBody2 )
            {
                currentPartial( 1, 0 ) *= -1.0;
                currentPartial( 1, 0 ) *= -1.0;
            }
            currentPartial( 1, 0 ) *= ( ( ( -orderOfBody1 +  orderOfBody2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) );
            currentPartial( 1, 1 ) *= ( ( ( -orderOfBody1 +  orderOfBody2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) );

            fullCoefficientsWrtBody2CoefficientsList[ effectiveIndex ] = currentPartial;
        }

        if( orderOfBody2 != 0 )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPartial = fullCoefficientsWrtBody2Coefficients * multipliers_[ effectiveIndex ];
            currentPartial( 0, 1 ) *= -1.0;
            currentPartial( 1, 1 ) *= -1.0;

            currentPartial( 1, 0 ) *= ( ( ( orderOfBody1 - orderOfBody2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) );
            currentPartial( 1, 1 ) *= ( ( ( orderOfBody1 - orderOfBody2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) );

            fullCoefficientsWrtBody2CoefficientsList[ effectiveIndex ] = currentPartial;
        }

        if( orderOfBody1 != 0 && orderOfBody2 != 0 )
        {
            effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPartial = fullCoefficientsWrtBody2Coefficients * multipliers_[ effectiveIndex ];

            fullCoefficientsWrtBody2CoefficientsList[ effectiveIndex ] = currentPartial;
        }
    }
}

void EffectiveMutualSphericalHarmonicsField::initializeMultipliers( )
{
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;

    int effectiveIndex;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( i ));

        effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        multipliers_[ effectiveIndex ] =
                getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2, areCoefficientsNormalized_ );

        effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );
        multipliers_[ effectiveIndex ] =
                getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2, areCoefficientsNormalized_ );

        effectiveIndex = getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );
        multipliers_[ effectiveIndex ] =
                getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2, areCoefficientsNormalized_ );

        effectiveIndex = getEffectiveIndex( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );
        multipliers_[ effectiveIndex ]  =
                getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2, areCoefficientsNormalized_ );
    }
}


}

}
