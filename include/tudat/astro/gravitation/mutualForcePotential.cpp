#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/Gravitation/mutualForcePotential.h"

namespace tudat
{

namespace gravitation
{

std::pair< int, int > getMaximumDegrees(
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    int maximumDegree1 = 0;
    int maximumDegree2 = 0;

    unsigned int degreeOfBody1, degreeOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );

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


double getMutualPotentialEffectiveCoefficientMultiplier(
        const int degree1, const int order1, const int degree2, const int order2, const bool areCoefficientsNormalized )
{
    if( areCoefficientsNormalized )
    {
        return getGammaCoefficientForMutualForcePotential( degree1, order1, degree2, order2 ) *
                std::sqrt( 4.0 * mathematical_constants::PI * ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 )  ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 )  ) /
                           ( ( order1 == 0 && order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) ) *
                getSigmaSignFunction( order1 ) * getSigmaSignFunction( order2 ) * getSigmaSignFunction( order1 + order2 ) *
                ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) *
                ( ( degree1 % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) );
    }
    else
    {
        return boost::math::factorial< double >( degree1 + degree2 - std::abs( order1 + order2 ) ) /
                ( boost::math::factorial< double >( degree2 - std::abs( order2 ) ) *
                  boost::math::factorial< double >( degree1 - std::abs( order1 ) ) ) *
                getSigmaSignFunction( order1 ) * getSigmaSignFunction( order2 ) * getSigmaSignFunction( order1 + order2 ) *
                ( ( order1 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) * ( ( order2 == 0 ) ? ( 1.0 ) : ( 0.5 ) ) *
                ( ( degree1 % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) );
    }
}

double computeSingleMutualForcePotentialTerm(
        const double effectiveCosineCoefficient,
        const double effectiveSineCoefficient,
        const double polynomialParameter,
        basic_mathematics::LegendreCache* legendreCache,
        const int degreeOfBody1,
        const int orderOfBody1,
        const int degreeOfBody2,
        const int orderOfBody2 )
{
    return ( effectiveCosineCoefficient * legendreCache->getCosineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) ) -
             effectiveSineCoefficient * legendreCache->getSineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) ) ) *
            legendreCache->getOrElseUpdate( degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 + orderOfBody2 ),
                                            polynomialParameter, &basic_mathematics::computeLegendrePolynomial );
}

double computeMutualForcePotential(
        const Eigen::Vector3d& bodyFixedPosition,
        const double effectiveGravitationalParameterOfBody1,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const int maximumDegreeOfBody1,
        const int maximumDegreeOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        basic_mathematics::LegendreCache* legendreCache )
{
    
    // Determine body fixed spherical position of body udnergoing acceleration.
    Eigen::Vector3d sphericalPositon =
            coordinate_conversions::convertCartesianToSpherical( bodyFixedPosition );
    double radius = sphericalPositon.x( );
    double latitude = mathematical_constants::PI / 2.0 - sphericalPositon.y( );
    double longitude = sphericalPositon.z( );
    
    double sineOfLatitude = std::sin( latitude );
    legendreCache->updateSines( longitude );
    legendreCache->update( sineOfLatitude, &basic_mathematics::computeLegendrePolynomial );
    double potential = 0.0;
    
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;
    
    std::vector< double > radiusRatioOfBody1List;
    double radiusRatioOfBody1 = equatorialRadiusOfBody1 / radius;
    radiusRatioOfBody1List.push_back( 1 );
    for( unsigned int i = 1; i < maximumDegreeOfBody1; i++ )
    {
        radiusRatioOfBody1List.push_back( radiusRatioOfBody1List.at( i - 1 ) * radiusRatioOfBody1 );
    }
    
    std::vector< double > radiusRatioOfBody2List;
    radiusRatioOfBody2List.push_back( 1 );
    double radiusRatioOfBody2 = equatorialRadiusOfBody2 / radius;
    for( unsigned int i = 1; i < maximumDegreeOfBody2; i++ )
    {
        radiusRatioOfBody2List.push_back( radiusRatioOfBody2List.at( i - 1 ) * radiusRatioOfBody2 );
    }
    
    
    double currentTerm = 0;
    for(  unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );
        
        currentTerm = 0;
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    sineOfLatitude, legendreCache, degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    sineOfLatitude,legendreCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    sineOfLatitude, legendreCache, degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                    sineOfLatitude, legendreCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );
        currentTerm *= radiusRatioOfBody1List.at( degreeOfBody1 );
        currentTerm *= radiusRatioOfBody2List.at( degreeOfBody2 );
        
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
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const int maximumDegree1, const int maximumDegree2,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        basic_mathematics::LegendreCache* legendreCache )
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

    legendreCache->updateSines( sphericalpositionOfBodySubjectToAcceleration( 2 ) );
    legendreCache->update( sineOfAngle, &basic_mathematics::computeGeodesyLegendrePolynomial );

    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );

    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree, totalOrder;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );

    bool computeTerm;

    std::vector< std::pair< double, double > > legendreTerms;
    legendreTerms.resize( ( maximumDegree1  + maximumDegree2 + 1 ) * ( maximumDegree1  + maximumDegree2 + 1 ) );
    for( unsigned int i = 0; i <= maximumDegree1  + maximumDegree2; i++ )
    {
        for( unsigned int j = 0; j <= i; j++ )
        {
            // Compute geodesy-normalized Legendre polynomials.
            const double legendrePolynomial = legendreCache->getPolynomial(
                        i, j );
            const double incrementedLegendrePolynomial = legendreCache->getPolynomial(
                        i, j + 1 );

            // Compute geodesy-normalized Legendre polynomial derivative.
            const double legendrePolynomialDerivative =
                    basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                        i, j, sineOfAngle,
                        legendrePolynomial, incrementedLegendrePolynomial );

            legendreTerms[ i + ( maximumDegree1  + maximumDegree2 + 1 ) * j ] = std::make_pair( legendrePolynomial, legendrePolynomialDerivative );
        }
    }


    // Loop through all degrees.
    std::pair< double, double > currentTerms;
    for ( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );

        totalDegree = degreeOfBody1 + degreeOfBody2;

        equatorialRadiusRatioPower = radius1Powers[ degreeOfBody1 ] * radius2Powers[ degreeOfBody2 ];

        for( unsigned j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                totalOrder = std::abs( orderOfBody1 + orderOfBody2 );
                computeTerm = 1;
                break;
            case 1:
                totalOrder = std::abs( -orderOfBody1 + orderOfBody2 );
                ( orderOfBody1 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            case 2:
                totalOrder = std::abs( orderOfBody1 - orderOfBody2 );
                ( orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );

                break;
            case 3:
                totalOrder = std::abs( -orderOfBody1 - orderOfBody2 );
                ( orderOfBody1 == 0 || orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            }

            if( computeTerm )
            {
                currentTerms = legendreTerms.at( totalDegree + ( maximumDegree1  + maximumDegree2 + 1 ) * totalOrder );

                // Compute the potential gradient of a single spherical harmonic term.
                sphericalGradient += basic_mathematics::computePotentialGradientFromRadiusRatio(
                            sphericalpositionOfBodySubjectToAcceleration,
                            equatorialRadiusRatioPower,
                            preMultiplier,
                            cosineOfAngle,
                            totalDegree,
                            totalOrder,
                            effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                            effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                            currentTerms.first,
                            currentTerms.second,
                            legendreCache );
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
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        basic_mathematics::LegendreCache* legendreCache )
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

    legendreCache->updateSines( sphericalpositionOfBodySubjectToAcceleration( 2 ) );
    legendreCache->update( sineOfAngle, &basic_mathematics::computeLegendrePolynomial );
    
    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );
    
    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2, totalDegree, totalOrder;
    double equatorialRadiusRatioPower;
    double preMultiplier = gravitationalParameterOfBody /
            (  sphericalpositionOfBodySubjectToAcceleration( 0 ) );
    
    bool computeTerm;

    // Loop through all degrees.
    for ( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );
        
        totalDegree = degreeOfBody1 + degreeOfBody2;
        
        equatorialRadiusRatioPower =
                raiseToIntegerPower( equatorialRadiusOfBody1 / sphericalpositionOfBodySubjectToAcceleration( 0 ), degreeOfBody1 ) *
                raiseToIntegerPower( equatorialRadiusOfBody2 / sphericalpositionOfBodySubjectToAcceleration( 0 ), degreeOfBody2 );
        
        for( unsigned j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                totalOrder = std::abs( orderOfBody1 + orderOfBody2 );
                computeTerm = 1;
                break;
            case 1:
                totalOrder = std::abs( -orderOfBody1 + orderOfBody2 );
                ( orderOfBody1 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            case 2:
                totalOrder = std::abs( orderOfBody1 - orderOfBody2 );
                ( orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                
                break;
            case 3:
                totalOrder = std::abs( -orderOfBody1 - orderOfBody2 );
                ( orderOfBody1 == 0 || orderOfBody2 == 0  ) ? ( computeTerm = 0 ) : ( computeTerm = 1 );
                break;
            }
            
            if( computeTerm )
            {

                
                // Compute geodesy-normalized Legendre polynomials.
                const double legendrePolynomial = basic_mathematics::computeLegendrePolynomial(
                            totalDegree, totalOrder, sineOfAngle, legendreCache);
                const double incrementedLegendrePolynomial =
                        basic_mathematics::computeLegendrePolynomial(
                            totalDegree, totalOrder + 1, sineOfAngle, legendreCache );
                
                // Compute geodesy-normalized Legendre polynomial derivative.
                const double legendrePolynomialDerivative =
                        basic_mathematics::computeLegendrePolynomialDerivative(
                            totalOrder, sineOfAngle,
                            legendrePolynomial, incrementedLegendrePolynomial );                

                // Compute the potential gradient of a single spherical harmonic term.
                sphericalGradient += basic_mathematics::computePotentialGradientFromRadiusRatio(
                            sphericalpositionOfBodySubjectToAcceleration,
                            equatorialRadiusRatioPower,
                            preMultiplier,
                            cosineOfAngle,
                            totalDegree,
                            totalOrder,
                            effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                            effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                            legendrePolynomial,
                            legendrePolynomialDerivative,
                            legendreCache );
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
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const double distance,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        basic_mathematics::LegendreCache* legendreCache,
        const boost::function< int( const int, const int, const int, const int )> effectiveIndexFunction )
{
    int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;
    int effectiveIndex;
    double equatorialRadiusRatioPower;

    Eigen::Matrix< double, 1, 2 > currentPotentialComponentWrtFullCoefficients;
    for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );

        equatorialRadiusRatioPower = radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 ) / distance;

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = legendreCache->getCosineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = legendreCache->getSineOfMultipleLongitude( std::abs( orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= legendreCache->getPolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 + orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = legendreCache->getCosineOfMultipleLongitude( std::abs( -orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = legendreCache->getSineOfMultipleLongitude( std::abs( -orderOfBody1 + orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= legendreCache->getPolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( -orderOfBody1 + orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = legendreCache->getCosineOfMultipleLongitude( std::abs( orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = legendreCache->getSineOfMultipleLongitude( std::abs( orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= legendreCache->getPolynomial(
                        degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 - orderOfBody2 ) ) * equatorialRadiusRatioPower;
            potentialComponentsWrtFullCoefficients[ effectiveIndex ] = currentPotentialComponentWrtFullCoefficients;
        }

        {
            effectiveIndex = effectiveIndexFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );

            currentPotentialComponentWrtFullCoefficients( 0, 0 ) = legendreCache->getCosineOfMultipleLongitude( std::abs( -orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients( 0, 1 ) = legendreCache->getSineOfMultipleLongitude( std::abs( -orderOfBody1 - orderOfBody2 ) );
            currentPotentialComponentWrtFullCoefficients *= legendreCache->getPolynomial(
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

}

void EffectiveMutualSphericalHarmonicsField::computeCurrentEffectiveCoefficients(
        const double angleTheta, const double anglePhi, const double anglePsi )
{
    cosineCoefficientsOfBody1_ = cosineCoefficientFunctionOfBody1_( );
    sineCoefficientsOfBody1_ = sineCoefficientFunctionOfBody1_( );
    cosineCoefficientsOfBody2_ = cosineCoefficientFunctionOfBody2_( );
    sineCoefficientsOfBody2_ = sineCoefficientFunctionOfBody2_( );

    transformationCache_->update( angleTheta, anglePhi, anglePsi );
    transformationCache_->transformCoefficientsAtDegree(
                cosineCoefficientsOfBody2_,
                sineCoefficientsOfBody2_,
                transformedCosineCoefficientsOfBody2_,
                transformedSineCoefficientsOfBody2_,
                areCoefficientsNormalized_ );

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
        degreeOfBody1 = coefficientCombinationsToUse_.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse_.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse_.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse_.at( i ).get< 3 >( );

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
        degreeOfBody1 = coefficientCombinationsToUse_.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse_.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse_.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse_.at( i ).get< 3 >( );

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
        degreeOfBody1 = coefficientCombinationsToUse_.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse_.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse_.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse_.at( i ).get< 3 >( );

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


