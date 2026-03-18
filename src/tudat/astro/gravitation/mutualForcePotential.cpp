#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/astro/gravitation/mutualForcePotential.h"

namespace tudat
{

namespace gravitation
{

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
    // Implements the cross-body scaling in the effective coefficients of Dirkx et al. (2019), Eqs. (47)-(48),
    // including sigma_m sign factors from Eq. (22).
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

//! Compute full two-body mutual acceleration for normalized coefficients.
/*!
 * Evaluates the gradient of the effective potential in the body-1 frame by combining:
 * - effective coefficients from Dirkx et al. (2019), Eqs. (47)-(48);
 * - one-body-like potential summation form, Eq. (49);
 * - translational equation usage as Cartesian potential gradient, Eq. (55).
 */
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

    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );

    // Cache P_lm and dP_lm/dphi up to the requested evaluation degree for efficient Eq. (55) accumulation.
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


    // Sum contributions of each selected (l1,m1,l2,m2) pair.
    // For each pair, expand to the required signed-order combinations entering Eq. (49).
    std::pair< double, double > currentTerms;
    for ( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = std::get<0>(coefficientCombinationsToUse.at( i ));
        orderOfBody1 = std::get<1>(coefficientCombinationsToUse.at( i ));
        degreeOfBody2 = std::get<2>(coefficientCombinationsToUse.at( i ));
        orderOfBody2 = std::get<3>(coefficientCombinationsToUse.at( i ));

        totalDegree = degreeOfBody1 + degreeOfBody2;

        equatorialRadiusRatioPower = radius1Powers[ degreeOfBody1 ] * radius2Powers[ degreeOfBody2 ];

        // Radius factors correspond to the (R1/r)^l1 (R2/r)^l2 scaling in Dirkx Eq. (49), substitution Eq. (44).
        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            const bool computeTerm = getSignedOrdersForCombinationCase(
                    j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
            if( computeTerm )
            {
                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );

                currentTerms = legendreTerms.at( totalDegree + ( maximumEvaluationDegree + 1 ) * totalOrder );
                const double effectiveCosineCoefficient = effectiveCosineCoefficientFunction(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double effectiveSineCoefficient = effectiveSineCoefficientFunction(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                // Compute and accumulate one term of the Eq. (55) gradient using effective coefficients from Eqs. (47)-(48).
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
            }
        }


    }

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.

    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

//! Compute full two-body mutual acceleration for unnormalized coefficients.
/*!
 * Same flow as computeGeodesyNormalizedMutualGravitationalAccelerationSum, but using unnormalized Legendre
 * polynomials/derivatives.
 */
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
    
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );

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
        
        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            const bool computeTerm = getSignedOrdersForCombinationCase(
                    j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
            
            if( computeTerm )
            {
                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
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


//! Compute one effective coefficient pair (C_eff,S_eff) for a signed (m1,m2) combination.
/*!
 * Implements the algebraic combinations used in Dirkx et al. (2019), Eqs. (47)-(48), including signed-order
 * handling and precomputed multipliers.
 */
void EffectiveMutualSphericalHarmonicsField::getCurrentEffectiveCoefficients(
        const int degree1, const int order1, const int degree2, const int order2,
        const int effectiveIndex,
        double& cosineCoefficient, double& sineCoefficient )
{
    // Eq. (47): effective cosine coefficient combination in terms of body-1 and transformed body-2 coefficients.
    cosineCoefficient = ( cosineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                          transformedCosineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) -
                          ( ( order1 < 0 ) ? ( -1.0 ) : ( 1.0 ) ) * ( ( order2 < 0 ) ? ( -1.0 ) : ( 1.0 ) )  *
                          sineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                          transformedSineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) );
    // Eq. (48): effective sine coefficient combination with signed-order terms.
    sineCoefficient = ( ( ( order2 < 0 ) ? ( -1.0 ) : ( 1.0 ) ) *
                        cosineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                        transformedSineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) +
                        ( ( order1 < 0 ) ? ( -1.0 ) : ( 1.0 ) )  *
                        sineCoefficientsOfBody1_( degree1, std::abs( order1 ) ) *
                        transformedCosineCoefficientsOfBody2_( degree2, std::abs( order2 ) ) );

    // Eqs. (47)-(48) with Eq. (22) sign conventions captured in cached multipliers_.
    double currentMultiplier = multipliers_.at( effectiveIndex );
    cosineCoefficient *= currentMultiplier;
    sineCoefficient *= ( ( ( order1 +  order2 ) < 0 ) ? ( -1.0 ) : ( 1.0 ) ) * currentMultiplier;
}

//! Update transformed body-2 coefficients and rebuild effective coefficients for current relative orientation.
/*!
 * Applies the F2->F1 spherical-harmonic rotation and then evaluates effective coefficients used by the
 * full two-body potential/acceleration evaluation (Dirkx et al. (2019), Eqs. (47)-(49)).
 */
void EffectiveMutualSphericalHarmonicsField::computeCurrentEffectiveCoefficients(
        const Eigen::Quaterniond coefficientRotationQuaterion )
{
    // Refresh current body-fixed fields and transform body-2 coefficients from F2 to F1.
    // This is the frame transformation used before combining terms into effective coefficients (Dirkx Eq. (47)-(48)).
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

    // Populate the effective coefficients entering the Eq. (49) potential/acceleration summation.
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

//! Update all effective coefficient entries used by the full two-body potential/acceleration.
/*!
 * Iterates over selected (l1,m1,l2,m2) combinations and populates all valid signed-order cases that contribute
 * to the Eq. (49) potential summation in Dirkx et al. (2019).
 */
void EffectiveMutualSphericalHarmonicsField::updateEffectiveMutualPotential( )
{

    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get<0>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get<1>( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get<2>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get<3>( coefficientCombinationsToUse_.at( i ) );

        // Expand to signed-order combinations using centralized helper to keep acceleration/potential logic consistent.
        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            if( !getSignedOrdersForCombinationCase(
                        j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
            {
                continue;
            }
            const int effectiveIndex = getEffectiveIndex(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
            getCurrentEffectiveCoefficients(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2, effectiveIndex,
                        effectiveCosineCoefficients_[ effectiveIndex ],
                        effectiveSineCoefficients_[ effectiveIndex ] );
        }
    }
}

void EffectiveMutualSphericalHarmonicsField::computePartialsOfFullCoefficientsWrtTransformedCoefficients(
        std::vector< Eigen::Matrix2d >& fullCoefficientsWrtBody2CoefficientsList )
{
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get<0>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get<1>( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get<2>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get<3>( coefficientCombinationsToUse_.at( i ) );

        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            if( !getSignedOrdersForCombinationCase(
                        j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
            {
                continue;
            }

            const int effectiveIndex =
                    getEffectiveIndex( degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
            const double signOrderOfBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
            const double signOrderOfBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
            const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
            const double cosineCoefficientOfBody1 =
                    cosineCoefficientsOfBody1_( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
            const double sineCoefficientOfBody1 =
                    sineCoefficientsOfBody1_( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
            const double multiplier = multipliers_[ effectiveIndex ];

            Eigen::Matrix2d currentPartial = Eigen::Matrix2d::Zero( );
            currentPartial( 0, 0 ) = cosineCoefficientOfBody1 * multiplier;
            currentPartial( 0, 1 ) = -signOrderOfBody1 * signOrderOfBody2 * sineCoefficientOfBody1 * multiplier;
            currentPartial( 1, 0 ) = signOrderOfBody1 * signTotalOrder * sineCoefficientOfBody1 * multiplier;
            currentPartial( 1, 1 ) = signOrderOfBody2 * signTotalOrder * cosineCoefficientOfBody1 * multiplier;
            fullCoefficientsWrtBody2CoefficientsList[ effectiveIndex ] = currentPartial;
        }
    }
}

//! Precompute constant scaling multipliers for each effective coefficient entry.
/*!
 * Stores the degree/order dependent scalar factors of Dirkx et al. (2019), Eqs. (47)-(48), including Eq. (22)
 * sign conventions. These factors are independent of current orientation and can be cached.
 */
void EffectiveMutualSphericalHarmonicsField::initializeMultipliers( )
{
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get<0>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get<1>( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get<2>( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get<3>( coefficientCombinationsToUse_.at( i ) );

        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            if( !getSignedOrdersForCombinationCase(
                        j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
            {
                continue;
            }

            const int effectiveIndex =
                    getEffectiveIndex( degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
            // Cache Eq. (47)-(48) prefactors (including Eq. (22) sign behavior) per signed-order combination.
            multipliers_[ effectiveIndex ] = getMutualPotentialEffectiveCoefficientMultiplier(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2, areCoefficientsNormalized_ );
        }
    }
}


}

}
