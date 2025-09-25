#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/utilityMacros.h"
#include <stdexcept>
#include <cmath>
#include <utility>

#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{
namespace aerodynamics
{
ComaModel::ComaModel( const simulation_setup::ComaPolyDataset& polyDataset,
                      std::function<Eigen::Vector6d()> sunStateFunction,
                      std::function<Eigen::Vector6d()> cometStateFunction,
                      std::function<Eigen::Matrix3d()> cometRotationFunction,
                      const int& maximumDegree,
                      const int& maximumOrder ) :
    dataType_( ComaDataType::POLYNOMIAL_COEFFICIENTS ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    polyDataset_( std::make_shared<simulation_setup::ComaPolyDataset>( polyDataset ) ),
    stokesDataset_( nullptr ),
    sunStateFunction_( std::move( sunStateFunction ) ),
    cometStateFunction_( std::move( cometStateFunction ) ),
    cometRotationFunction_( std::move( cometRotationFunction ) ),
    sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
{
    if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
    {
        throw std::invalid_argument( "ComaModel: All state functions must be provided" );
    }

    if ( polyDataset_->getNumFiles() == 0 )
    {
        throw std::invalid_argument( "ComaModel: PolyDataset must contain at least one file" );
    }

    if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
    {
        throw std::invalid_argument( "ComaModel: Maximum degree and order must be >= -1" );
    }
}

ComaModel::ComaModel( const simulation_setup::ComaStokesDataset& stokesDataset,
                      std::function<Eigen::Vector6d()> sunStateFunction,
                      std::function<Eigen::Vector6d()> cometStateFunction,
                      std::function<Eigen::Matrix3d()> cometRotationFunction,
                      const int& maximumDegree,
                      const int& maximumOrder ) :
    dataType_( ComaDataType::STOKES_COEFFICIENTS ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    polyDataset_( nullptr ),
    stokesDataset_( std::make_shared<simulation_setup::ComaStokesDataset>( stokesDataset ) ),
    sunStateFunction_( std::move( sunStateFunction ) ),
    cometStateFunction_( std::move( cometStateFunction ) ),
    cometRotationFunction_( std::move( cometRotationFunction ) ),
    sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
{
    if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
    {
        throw std::invalid_argument( "ComaModel: All state functions must be provided" );
    }

    if ( stokesDataset_->nFiles() == 0 )
    {
        throw std::invalid_argument( "ComaModel: StokesDataset must contain at least one file" );
    }

    if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
    {
        throw std::invalid_argument( "ComaModel: Maximum degree and order must be >= -1" );
    }
}



double ComaModel::getDensity( const double radius,
                              const double longitude,
                              const double latitude,
                              const double time )
{
    switch ( dataType_ )
    {
        case ComaDataType::POLYNOMIAL_COEFFICIENTS:
            return computeDensityFromPolyCoefficients( radius, longitude, latitude, time );

        case ComaDataType::STOKES_COEFFICIENTS:
            return computeDensityFromStokesCoefficients( radius, longitude, latitude, time );

        default:
            throw std::runtime_error( "ComaModel: Unknown data type" );
    }
}


double ComaModel::getPressure( const double radius,
                               const double longitude,
                               const double latitude,
                               const double time )
{
    // For a coma, pressure is typically negligible compared to planetary atmospheres
    // Return a small value proportional to density
    const double density = getDensity( radius, longitude, latitude, time );

    // Assume ideal gas law: P = ρRT/M
    // Use typical values for comet coma gas (mostly water vapor)
    const double temperature = 200.0; // K, typical coma temperature
    const double molarMass = 0.018; // kg/mol, water vapor
    const double gasConstant = 8.314; // J/(mol·K)

    return density * gasConstant * temperature / molarMass;
}

double ComaModel::getTemperature( const double radius,
                                  const double longitude,
                                  const double latitude,
                                  const double time )
{
    // For comet coma, temperature varies with distance from nucleus and solar irradiation
    // Implement a simple model based on altitude and solar heating

    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );

    const double cometRadius = 1000.0; // m, typical comet radius
    const double surfaceTemperature = 200.0; // K, surface temperature
    const double spaceTemperature = 2.7; // K, cosmic background temperature

    // Simple exponential decay with altitude
    const double distanceFromSurface = radius - cometRadius;
    const double scaleHeight = 10000.0; // m, characteristic scale for temperature decay

    return spaceTemperature + ( surfaceTemperature - spaceTemperature ) *
           std::exp( -distanceFromSurface / scaleHeight );
}

double ComaModel::getSpeedOfSound( const double radius,
                                   const double longitude,
                                   const double latitude,
                                   const double time )
{
    // Speed of sound in ideal gas: c = √(γRT/M)
    // where γ is heat capacity ratio, R is gas constant, T is temperature, M is molar mass

    const double temperature = getTemperature( radius, longitude, latitude, time );

    // Typical values for water vapor (dominant component in comet coma)
    const double heatCapacityRatio = 1.33; // For water vapor
    const double molarMass = 0.018; // kg/mol
    const double gasConstant = 8.314; // J/(mol·K)

    return std::sqrt( heatCapacityRatio * gasConstant * temperature / molarMass );
}


int ComaModel::findTimeIntervalIndex( const double time ) const
{
    for ( std::size_t f = 0; f < polyDataset_->getNumFiles(); ++f )
    {
        const auto& fileMeta = polyDataset_->getFileMeta( f );
        for ( const auto& period : fileMeta.timePeriods )
        {
            if ( time >= period.first && time <= period.second )
            {
                return static_cast<int>( f );
            }
        }
    }
    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}


SphericalHarmonicsCalculator::SphericalHarmonicsCalculator( std::string fixedReferenceFrame ) :
    fixedReferenceFrame_( std::move( fixedReferenceFrame ) ),
    sphericalHarmonicsCache_( basic_mathematics::SphericalHarmonicsCache( true, true) )
{
}


double SphericalHarmonicsCalculator::calculateSurfaceSphericalHarmonics(
        const Eigen::MatrixXd& sineCoefficients,
        const Eigen::MatrixXd& cosineCoefficients,
        const double latitude,
        const double longitude,
        const int highestDegree,
        const int highestOrder )
{
    if(cosineCoefficients.rows( ) != sineCoefficients.rows( ) ||
        cosineCoefficients.cols( ) != sineCoefficients.cols( ))
    {
        throw std::runtime_error( "Spherical harmonics coefficient sizes are incompatible." );
    }

    const int maxDegree = static_cast< int >(cosineCoefficients.rows( ));
    const int maxOrder = static_cast< int >(cosineCoefficients.cols( ));

    sphericalHarmonicsCache_.resetMaximumDegreeAndOrder( maxDegree, maxOrder );

    const double sineOfAngle = std::sin( latitude );
    sphericalHarmonicsCache_.updateAnglesOnly( sineOfAngle,
                                     longitude );

    const basic_mathematics::LegendreCache& legendreCacheReference = sphericalHarmonicsCache_.getLegendreCache( );


    double legendrePolynomial = TUDAT_NAN;
    double cosineOfOrderLongitude = TUDAT_NAN;
    double sineOfOrderLongitude = TUDAT_NAN;
    double value = 0.0;

    // Loop through all degrees.
    for (int l = 0; l <= highestDegree; ++l) {
        const int mmax = std::min(l, highestOrder);
        for (int m = 0; m <= mmax; ++m) {
            const double P = legendreCacheReference.getLegendrePolynomial(l, m); // normalized, correct argument
            const double cos_mλ = sphericalHarmonicsCache_.getCosineOfMultipleLongitude(m);

            if (m == 0) {
                value += P * cosineCoefficients(l, 0);
            } else {
                const double sin_mλ = sphericalHarmonicsCache_.getSineOfMultipleLongitude(m);
                value += P * (cosineCoefficients(l, m) * cos_mλ
                            + sineCoefficients(l, m)    * sin_mλ);
            }
        }
    }
    return value;
}







double ComaModel::calculateSolarLongitude() const
{
    if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
    {
        throw std::runtime_error( "ComaModel: State functions must be initialized" );
    }

    const Eigen::Vector6d sunState = sunStateFunction_();
    const Eigen::Vector6d cometState = cometStateFunction_();
    const Eigen::Matrix3d rotationMatrix = cometRotationFunction_();

    // Calculate Sun direction in inertial frame
    const Eigen::Vector3d sunDirection = ( sunState.head<3>() - cometState.head<3>() ).normalized();

    // Transform to comet body-fixed frame
    const Eigen::Vector3d sunDirectionBodyFixed = rotationMatrix.transpose() * sunDirection;

    // Calculate solar longitude (angle from X-axis in XY plane)
    return std::atan2( sunDirectionBodyFixed.y(), sunDirectionBodyFixed.x() );
}

double ComaModel::computeDensityFromPolyCoefficients( double radius, double longitude, double latitude, double time ) const
{
    if ( !polyDataset_ )
    {
        throw std::runtime_error( "ComaModel: polyDataset_ is null" );
    }

    const int fileIndex = findTimeIntervalIndex( time );
    const auto& fileMeta = polyDataset_->getFileMeta( fileIndex );
    const auto& polyCoefficients = polyDataset_->getPolyCoefficients( fileIndex );
    const auto& shIndices = polyDataset_->getSHDegreeAndOrderIndices( fileIndex );

    const double solarLongitude = calculateSolarLongitude();

    Eigen::MatrixXd cosineCoefficients, sineCoefficients;

    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius,
        solarLongitude,
        polyCoefficients,
        shIndices,
        fileMeta.powersInvRadius,
        fileMeta.referenceRadius,
        cosineCoefficients,
        sineCoefficients,
        maximumDegree_,
        maximumOrder_ );

    return sphericalHarmonicsCalculator_->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cosineCoefficients.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cosineCoefficients.cols() - 1 );
}

double ComaModel::computeDensityFromStokesCoefficients( double radius, double longitude, double latitude, double time ) const
{
    if ( !stokesDataset_ )
    {
        throw std::runtime_error( "ComaModel: stokesDataset_ is null" );
    }

    // TODO: Implement Stokes coefficient based density calculation
    // This would involve:
    // 1. Find appropriate Stokes coefficients for the given time/position
    // 2. Interpolate if necessary
    // 3. Calculate density using spherical harmonics

    TUDAT_UNUSED_PARAMETER( radius );
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );

    throw std::runtime_error( "ComaModel: Stokes coefficient algorithm not yet implemented" );
}


}
}
