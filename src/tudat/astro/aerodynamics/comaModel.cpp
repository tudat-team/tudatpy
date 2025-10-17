#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/utilityMacros.h"
#include <stdexcept>
#include <cmath>
#include <utility>

#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"

namespace tudat
{
namespace aerodynamics
{

//=============================================================================
// ComaModel Class Implementation
//=============================================================================

//-----------------------------------------------------------------------------
// Constructors
//-----------------------------------------------------------------------------

/*!
 * \brief Constructor for ComaModel using polynomial coefficient data.
 * Initializes the coma model with polynomial coefficients for computing density
 * and other atmospheric properties using spherical harmonics expansion.
 * \param polyDataset Structured polynomial coefficient dataset
 * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
 * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
 * \param cometRotationFunction Function returning comet body-fixed rotation matrix
 * \param maximumDegree Maximum degree used to compute the coma density with SH (-1 for auto)
 * \param maximumOrder Maximum Order used to compute the coma density with SH (-1 for auto)
 */
ComaModel::ComaModel( const simulation_setup::ComaPolyDataset& polyDataset,
                      const double molecularWeight,
                      std::function<Eigen::Vector6d()> sunStateFunction,
                      std::function<Eigen::Vector6d()> cometStateFunction,
                      std::function<Eigen::Matrix3d()> cometRotationFunction,
                      const int& maximumDegree,
                      const int& maximumOrder ) :
    AtmosphereModel( false, false, true ),  // Use radius instead of altitude
    cachedSolarLongitude_( 0.0 ),
    cachedTime_( -std::numeric_limits<double>::infinity() ),
    lastFileIndex_( 0 ),
    cachedRadius_( 0.0 ),
    cachedInterpolationSolarLongitude_( 0.0 ),
    coefficientMatricesSized_( false ),
    cachedLatitude_( 0.0 ),
    cachedLongitude_( 0.0 ),
    cachedSineLatitude_( 0.0 ),
    cachedFinalDensity_( 0.0 ),
    interpolationPoint2D_( 2 ),
    interpolationPoint1D_( 1 ),
    dataType_( ComaDataType::POLYNOMIAL_COEFFICIENTS ),
    molecularWeight_( molecularWeight ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    polyDataset_( std::make_shared<simulation_setup::ComaPolyDataset>( polyDataset ) ),
    stokesDataset_( nullptr ),
    sunStateFunction_( std::move( sunStateFunction ) ),
    cometStateFunction_( std::move( cometStateFunction ) ),
    cometRotationFunction_( std::move( cometRotationFunction ) ),
    sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
{
    // Validate input parameters
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

    // Pre-allocate coefficient matrices based on maximum degree/order from dataset
    const int maxDegreeAvailable = polyDataset_->getMaxDegreeSH( 0 );
    const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : maxDegreeAvailable;
    const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : maxDegreeAvailable;
    cachedCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
    cachedSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
}

/*!
 * \brief Constructor for ComaModel using Stokes coefficient data.
 * Initializes the coma model with Stokes coefficients and pre-computes
 * interpolators for efficient density evaluation.
 * \param stokesDataset Structured Stokes coefficient dataset
 * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
 * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
 * \param cometRotationFunction Function returning comet body-fixed rotation matrix
 * \param maximumDegree Maximum degree used to compute the coma density with SH (-1 for auto)
 * \param maximumOrder Maximum Order used to compute the coma density with SH (-1 for auto)
 */
ComaModel::ComaModel( const simulation_setup::ComaStokesDataset& stokesDataset,
                      const double molecularWeight,
                      std::function<Eigen::Vector6d()> sunStateFunction,
                      std::function<Eigen::Vector6d()> cometStateFunction,
                      std::function<Eigen::Matrix3d()> cometRotationFunction,
                      const int& maximumDegree,
                      const int& maximumOrder ) :
    AtmosphereModel( false, false, true ),  // Use radius instead of altitude
    cachedSolarLongitude_( 0.0 ),
    cachedTime_( -std::numeric_limits<double>::infinity() ),
    lastFileIndex_( 0 ),
    cachedRadius_( 0.0 ),
    cachedInterpolationSolarLongitude_( 0.0 ),
    coefficientMatricesSized_( false ),
    cachedLatitude_( 0.0 ),
    cachedLongitude_( 0.0 ),
    cachedSineLatitude_( 0.0 ),
    cachedFinalDensity_( 0.0 ),
    interpolationPoint2D_( 2 ),
    interpolationPoint1D_( 1 ),
    dataType_( ComaDataType::STOKES_COEFFICIENTS ),
    molecularWeight_( molecularWeight ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    polyDataset_( nullptr ),
    stokesDataset_( std::make_shared<simulation_setup::ComaStokesDataset>( stokesDataset ) ),
    sunStateFunction_( std::move( sunStateFunction ) ),
    cometStateFunction_( std::move( cometStateFunction ) ),
    cometRotationFunction_( std::move( cometRotationFunction ) ),
    sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
{
    // Validate input parameters
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

    // Pre-allocate coefficient matrices based on nmax from dataset
    const int nmax = stokesDataset_->nmax();
    const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : nmax;
    cachedCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
    cachedSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );

    // Initialize interpolators for efficient Stokes coefficient evaluation
    initializeStokesInterpolators();
}

//-----------------------------------------------------------------------------
// Public Interface Methods (AtmosphereModel implementation)
//-----------------------------------------------------------------------------

/*!
 * \brief Compute coma density at specified location and time.
 * Uses either polynomial or Stokes coefficients depending on the data type
 * to evaluate the density field using spherical harmonics.
 * \param radius Radius from comet center at which density is to be computed [m]
 * \param longitude Longitude in comet body-fixed frame at which density is to be computed [rad]
 * \param latitude Latitude in comet body-fixed frame at which density is to be computed [rad]
 * \param time Time at which density is to be computed [s]
 * \return Coma density at specified location and time [kg/m³]
 */
double ComaModel::getDensity( const double radius,
                              const double longitude,
                              const double latitude,
                              const double time )
{
    // Validate radius against reference radius
    const int fileIndex = findTimeIntervalIndex( time );
    double referenceRadius;

    if ( dataType_ == ComaDataType::POLYNOMIAL_COEFFICIENTS )
    {
        referenceRadius = polyDataset_->getFileMeta( fileIndex ).referenceRadius;
    }
    else // ComaDataType::STOKES_COEFFICIENTS
    {
        referenceRadius = stokesDataset_->getReferenceRadius( fileIndex );
    }

    if ( radius < referenceRadius )
    {
        throw std::runtime_error( "ComaModel: Radius " + std::to_string( radius ) +
                                  " is smaller than reference radius " + std::to_string( referenceRadius ) );
    }

    // Get number density and convert to mass density
    // Convert: number_density [1/m³] × molecular_weight [kg/mol] / N_A [1/mol] = mass_density [kg/m³]
    const double numberDensity = getNumberDensity( radius, longitude, latitude, time );
    return numberDensity * molecularWeight_ / physical_constants::AVOGADRO_CONSTANT;
}

/*!
 * \brief Compute coma number density at specified location and time.
 * Returns the number density (particles per cubic meter) by evaluating the
 * spherical harmonics expansion.
 * \param radius Radius from comet center at which number density is to be computed [m]
 * \param longitude Longitude in comet body-fixed frame at which number density is to be computed [rad]
 * \param latitude Latitude in comet body-fixed frame at which number density is to be computed [rad]
 * \param time Time at which number density is to be computed [s]
 * \return Coma number density at specified location and time [m^-3]
 */
double ComaModel::getNumberDensity( const double radius,
                                    const double longitude,
                                    const double latitude,
                                    const double time )
{
    // Check if we can reuse cached final density result
    constexpr double tolerance = 1e-10;
    constexpr double toleranceSq = tolerance * tolerance;

    const double radiusDiff = radius - cachedRadius_;
    const double lonDiff = longitude - cachedLongitude_;
    const double latDiff = latitude - cachedLatitude_;
    const double timeDiff = time - cachedTime_;

    if ( cacheFlags_.densityValid &&
         radiusDiff * radiusDiff < toleranceSq &&
         lonDiff * lonDiff < toleranceSq &&
         latDiff * latDiff < toleranceSq &&
         timeDiff * timeDiff < toleranceSq )
    {
        return cachedFinalDensity_;
    }

    // Internal compute functions return log2 of number density
    // Convert to actual number density before returning
    double numberDensityLog2;

    switch ( dataType_ )
    {
        case ComaDataType::POLYNOMIAL_COEFFICIENTS:
            numberDensityLog2 = computeNumberDensityFromPolyCoefficients( radius, longitude, latitude, time );
            break;

        case ComaDataType::STOKES_COEFFICIENTS:
            numberDensityLog2 = computeNumberDensityFromStokesCoefficients( radius, longitude, latitude, time );
            break;

        default:
            throw std::runtime_error( "ComaModel: Unknown data type" );
    }

    // Convert log2(number_density) to actual number_density and cache it
    cachedFinalDensity_ = std::exp2( numberDensityLog2 );
    cacheFlags_.densityValid = true;

    return cachedFinalDensity_;
}

/*!
 * \brief Compute coma pressure at specified location and time.
 * Estimates pressure using ideal gas law with typical coma temperature
 * and composition (primarily water vapor).
 * \param radius Radius from comet center at which pressure is to be computed [m]
 * \param longitude Longitude in comet body-fixed frame at which pressure is to be computed [rad]
 * \param latitude Latitude in comet body-fixed frame at which pressure is to be computed [rad]
 * \param time Time at which pressure is to be computed [s]
 * \return Coma pressure at specified location and time [N/m²]
 */
double ComaModel::getPressure( const double radius,
                               const double longitude,
                               const double latitude,
                               const double time )
{
    // For a coma, pressure is typically negligible compared to planetary atmospheres
    // Return a small value proportional to density using ideal gas law: P = ρRT/M
    const double density = getDensity( radius, longitude, latitude, time );

    // Use typical values for comet coma gas (mostly water vapor)
    const double temperature = 200.0;      // K, typical coma temperature
    const double molarMass = 0.018;        // kg/mol, water vapor
    const double gasConstant = 8.314;      // J/(mol·K)

    return density * gasConstant * temperature / molarMass;
}

/*!
 * \brief Compute coma temperature at specified location and time.
 * Uses a simple exponential decay model with distance from the comet surface.
 * \param radius Radius from comet center at which temperature is to be computed [m]
 * \param longitude Longitude in comet body-fixed frame at which temperature is to be computed [rad]
 * \param latitude Latitude in comet body-fixed frame at which temperature is to be computed [rad]
 * \param time Time at which temperature is to be computed [s]
 * \return Coma temperature at specified location and time [K]
 */
double ComaModel::getTemperature( const double radius,
                                  const double longitude,
                                  const double latitude,
                                  const double time )
{
    // For comet coma, temperature varies with distance from nucleus and solar irradiation
    // Implement a simple exponential decay model based on altitude

    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );

    const double cometRadius = 1000.0;         // m, typical comet radius
    const double surfaceTemperature = 200.0;   // K, surface temperature
    const double spaceTemperature = 2.7;       // K, cosmic background temperature
    const double scaleHeight = 10000.0;        // m, characteristic scale for temperature decay

    // Simple exponential decay with altitude
    const double distanceFromSurface = radius - cometRadius;
    return spaceTemperature + ( surfaceTemperature - spaceTemperature ) *
           std::exp( -distanceFromSurface / scaleHeight );
}

/*!
 * \brief Compute speed of sound in coma at specified location and time.
 * Calculates speed of sound using ideal gas formula with local temperature
 * and typical water vapor properties.
 * \param radius Radius from comet center at which speed of sound is to be computed [m]
 * \param longitude Longitude in comet body-fixed frame at which speed of sound is to be computed [rad]
 * \param latitude Latitude in comet body-fixed frame at which speed of sound is to be computed [rad]
 * \param time Time at which speed of sound is to be computed [s]
 * \return Coma speed of sound at specified location and time [m/s]
 */
double ComaModel::getSpeedOfSound( const double radius,
                                   const double longitude,
                                   const double latitude,
                                   const double time )
{
    // Speed of sound in ideal gas: c = √(γRT/M)
    // where γ is heat capacity ratio, R is gas constant, T is temperature, M is molar mass
    const double temperature = getTemperature( radius, longitude, latitude, time );

    // Typical values for water vapor (dominant component in comet coma)
    const double heatCapacityRatio = 1.33;     // For water vapor
    const double molarMass = 0.018;            // kg/mol
    const double gasConstant = 8.314;          // J/(mol·K)

    return std::sqrt( heatCapacityRatio * gasConstant * temperature / molarMass );
}

//-----------------------------------------------------------------------------
// Private Helper Methods - Coordinate and Time Utilities
//-----------------------------------------------------------------------------

/*!
 * \brief Find the time interval index containing the specified time.
 * Searches through dataset files to determine which time period contains
 * the given time value. Uses binary search for Stokes datasets and caches
 * the last used file index as a hint for the next search.
 * \param time The time for which the corresponding time interval index is to be found [s]
 * \return The index of the time period in which the given time falls
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaModel::findTimeIntervalIndex( const double time ) const
{
    if ( dataType_ == ComaDataType::POLYNOMIAL_COEFFICIENTS )
    {
        // Polynomial coefficients: check last used file first (cache optimization)
        const std::size_t numFiles = polyDataset_->getNumFiles();
        if ( lastFileIndex_ >= 0 && lastFileIndex_ < static_cast<int>( numFiles ) )
        {
            const auto& fileMeta = polyDataset_->getFileMeta( lastFileIndex_ );
            for ( const auto& period : fileMeta.timePeriods )
            {
                if ( time >= period.first && time <= period.second )
                {
                    return lastFileIndex_;
                }
            }
        }

        // Cache miss - search all files
        for ( std::size_t fileIndex = 0; fileIndex < numFiles; ++fileIndex )
        {
            const auto& fileMeta = polyDataset_->getFileMeta( fileIndex );
            for ( const auto& period : fileMeta.timePeriods )
            {
                if ( time >= period.first && time <= period.second )
                {
                    lastFileIndex_ = static_cast<int>( fileIndex );
                    return lastFileIndex_;
                }
            }
        }
    }
    else if ( dataType_ == ComaDataType::STOKES_COEFFICIENTS )
    {
        const std::size_t numFiles = stokesDataset_->nFiles();

        // Check last used file first (temporal locality optimization)
        if ( lastFileIndex_ >= 0 && lastFileIndex_ < static_cast<int>( numFiles ) )
        {
            const auto& fileMeta = stokesDataset_->files()[lastFileIndex_];
            if ( time >= fileMeta.start_epoch && time <= fileMeta.end_epoch )
            {
                return lastFileIndex_;
            }
        }

        // Binary search for time interval (assumes files are ordered by time)
        int left = 0;
        int right = static_cast<int>( numFiles ) - 1;

        while ( left <= right )
        {
            const int mid = left + ( right - left ) / 2;
            const auto& fileMeta = stokesDataset_->files()[mid];

            if ( time < fileMeta.start_epoch )
            {
                right = mid - 1;
            }
            else if ( time > fileMeta.end_epoch )
            {
                left = mid + 1;
            }
            else
            {
                // Found the interval
                lastFileIndex_ = mid;
                return mid;
            }
        }
    }

    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}

/*!
 * \brief Calculate solar longitude in comet body-fixed frame with caching.
 * Computes the angle from X-axis to Sun direction projected onto the XY plane,
 * which is used to account for solar heating variations. Results are cached
 * to avoid redundant state function calls when time hasn't changed.
 * \param time Time at which to compute solar longitude [s]
 * \return Solar longitude angle from X-axis in XY plane [rad]
 * \throws std::runtime_error If state functions are not initialized
 */
double ComaModel::calculateSolarLongitude( const double time ) const
{
    constexpr double timeTolerance = 1e-10;
    constexpr double timeToleranceSq = timeTolerance * timeTolerance;
    const double timeDiff = time - cachedTime_;

    // Check if cached value is still valid (time hasn't changed significantly)
    if ( cacheFlags_.solarLongitudeValid && timeDiff * timeDiff < timeToleranceSq )
    {
        return cachedSolarLongitude_;
    }

    if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
    {
        throw std::runtime_error( "ComaModel: State functions must be initialized" );
    }

    cachedSunState_ = sunStateFunction_();
    cachedCometState_ = cometStateFunction_();
    cachedRotationMatrix_ = cometRotationFunction_();

    // Calculate Sun direction in inertial frame
    const Eigen::Vector3d sunDirection = ( cachedSunState_.head<3>() - cachedCometState_.head<3>() ).normalized();

    // Transform to comet body-fixed frame
    const Eigen::Vector3d sunDirectionBodyFixed = cachedRotationMatrix_.transpose() * sunDirection;

    // Calculate solar longitude (angle from X-axis in XY plane)
    cachedSolarLongitude_ = std::atan2( sunDirectionBodyFixed.y(), sunDirectionBodyFixed.x() );
    cachedTime_ = time;
    cacheFlags_.solarLongitudeValid = true;
    cacheFlags_.stateValid = true;

    return cachedSolarLongitude_;
}

//-----------------------------------------------------------------------------
// Private Helper Methods - Density Computation
//-----------------------------------------------------------------------------

/*!
 * \brief Compute density from polynomial coefficients.
 * Evaluates polynomial coefficients at the given radius and solar longitude,
 * then uses spherical harmonics to compute density at the surface location.
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute density [s]
 * \return Coma density [kg/m³]
 * \throws std::runtime_error If dataset is null or time is out of range
 */
double ComaModel::computeNumberDensityFromPolyCoefficients( double radius, double longitude, double latitude, double time ) const
{
    if ( !polyDataset_ )
    {
        throw std::runtime_error( "ComaModel: polyDataset_ is null" );
    }

    // Get time-dependent data
    const int fileIndex = findTimeIntervalIndex( time );
    const auto& fileMeta = polyDataset_->getFileMeta( fileIndex );
    const auto& polyCoefficients = polyDataset_->getPolyCoefficients( fileIndex );
    const auto& shIndices = polyDataset_->getSHDegreeAndOrderIndices( fileIndex );

    // Calculate current solar longitude for heliocentric dependence (with caching)
    const double solarLongitude = calculateSolarLongitude( time );

    // Only resize/zero coefficient matrices if dimensions changed
    const int maxDegreeAvailable = fileMeta.maxDegreeSH;
    const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : maxDegreeAvailable;
    const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : maxDegreeAvailable;

    // Check sizing flag first to avoid repeated size comparisons
    if ( !coefficientMatricesSized_ ||
         cachedCosineCoefficients_.rows() != effectiveMaxDegree + 1 ||
         cachedCosineCoefficients_.cols() != effectiveMaxOrder + 1 )
    {
        cachedCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        coefficientMatricesSized_ = true;
    }
    else
    {
        // Same size - just zero the contents (faster than reallocation)
        cachedCosineCoefficients_.setZero();
        cachedSineCoefficients_.setZero();
    }

    // Evaluate polynomial coefficients to get spherical harmonic coefficients
    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius,
        solarLongitude,
        polyCoefficients,
        shIndices,
        fileMeta.powersInvRadius,
        fileMeta.referenceRadius,
        cachedCosineCoefficients_,
        cachedSineCoefficients_,
        maximumDegree_,
        maximumOrder_ );

    // Update cached position values for density cache
    cachedLatitude_ = latitude;
    cachedLongitude_ = longitude;

    // Compute number density using spherical harmonics expansion (returns log2 of number density)
    return sphericalHarmonicsCalculator_->calculateSurfaceSphericalHarmonics(
        cachedSineCoefficients_, cachedCosineCoefficients_,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cachedCosineCoefficients_.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cachedCosineCoefficients_.cols() - 1 );
}

/*!
 * \brief Compute density from Stokes coefficients using interpolators.
 * Uses pre-initialized interpolators to efficiently evaluate Stokes coefficients
 * and applies distance-dependent decay for radii beyond the reference radius.
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute density [s]
 * \return Coma density [kg/m³]
 * \throws std::runtime_error If dataset is null or time is out of range
 */
double ComaModel::computeNumberDensityFromStokesCoefficients( double radius, double longitude, double latitude, double time ) const
{
    if ( !stokesDataset_ )
    {
        throw std::runtime_error( "ComaModel: stokesDataset_ is null" );
    }

    // Step 1: Get time-dependent properties
    const int fileIndex = findTimeIntervalIndex( time );
    const double solarLongitude = calculateSolarLongitude( time );

    // Step 2: Get dataset properties
    const int nmax = stokesDataset_->nmax();
    const double referenceRadius = stokesDataset_->getReferenceRadius(fileIndex);

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Check if coefficient interpolation can be skipped (same radius/solar longitude)
    constexpr double radiusTolerance = 1e-10;
    constexpr double solarLongitudeTolerance = 1e-10;
    constexpr double radiusToleranceSq = radiusTolerance * radiusTolerance;
    constexpr double solarLongitudeToleranceSq = solarLongitudeTolerance * solarLongitudeTolerance;

    const double radiusDiff = radius - cachedRadius_;
    const double solarLongitudeDiff = solarLongitude - cachedInterpolationSolarLongitude_;

    const bool sameRadius = cacheFlags_.interpolationValid &&
                           radiusDiff * radiusDiff < radiusToleranceSq;
    const bool sameSolarLongitude = cacheFlags_.interpolationValid &&
                                   solarLongitudeDiff * solarLongitudeDiff < solarLongitudeToleranceSq;

    if ( !sameRadius || !sameSolarLongitude )
    {
        // Cache miss - need to recompute coefficients
        // Step 3: Only zero coefficient matrices if dimensions match
        if ( cachedCosineCoefficients_.rows() == effectiveMaxDegree + 1 &&
             cachedCosineCoefficients_.cols() == effectiveMaxOrder + 1 )
        {
            cachedCosineCoefficients_.setZero();
            cachedSineCoefficients_.setZero();
        }
        else
        {
            cachedCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
            cachedSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        }

    // Step 4: Choose interpolation strategy based on radius
    if (radius <= referenceRadius)
    {
        interpolationPoint2D_[0] = radius;
        interpolationPoint2D_[1] = solarLongitude;

        // Prepare interpolation state once for all coefficients (batch mode)
        // Check if we have any interpolators to use for state preparation
        if ( !stokesInterpolators_[fileIndex].empty() )
        {
            // Use the first available interpolator to prepare the shared interpolation state
            auto firstInterpolator = stokesInterpolators_[fileIndex].begin()->second.first.get();
            auto interpolationState = firstInterpolator->prepareInterpolationState(interpolationPoint2D_);

            // Now use the pre-computed state for all coefficient interpolations
            for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
            {
                for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                {
                    std::pair<int,int> degreeOrderPair = {degree, order};
                    auto it = stokesInterpolators_[fileIndex].find(degreeOrderPair);

                    if ( it != stokesInterpolators_[fileIndex].end() )
                    {
                        // Use batch interpolation with pre-computed state (much faster!)
                        cachedCosineCoefficients_(degree, order) = it->second.first->interpolateWithState(interpolationState);
                        if ( order > 0 )  // Sine coefficients only exist for m > 0
                        {
                            cachedSineCoefficients_(degree, order) = it->second.second->interpolateWithState(interpolationState);
                        }
                    }
                    else
                    {
                        // Fallback: create 2D interpolator on-the-fly using helper method
                        createFallback2DInterpolator( fileIndex, degree, order,
                                                      cachedCosineCoefficients_(degree, order),
                                                      cachedSineCoefficients_(degree, order),
                                                      radius, solarLongitude );
                    }
                }
            }
        }
        else
        {
            // No interpolators available - use fallback for all coefficients
            for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
            {
                for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                {
                    createFallback2DInterpolator( fileIndex, degree, order,
                                                  cachedCosineCoefficients_(degree, order),
                                                  cachedSineCoefficients_(degree, order),
                                                  radius, solarLongitude );
                }
            }
        }
    }
    else
    {
        interpolationPoint1D_[0] = solarLongitude;

        // Prepare interpolation state once for all reduced coefficients (batch mode)
        if ( !reducedStokesInterpolators_[fileIndex].empty() )
        {
            // Use the first available interpolator to prepare the shared interpolation state
            auto firstReducedInterpolator = reducedStokesInterpolators_[fileIndex].begin()->second.first.get();
            auto reducedInterpolationState = firstReducedInterpolator->prepareInterpolationState(interpolationPoint1D_);

            // Now use the pre-computed state for all coefficient interpolations
            for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
            {
                for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                {
                    std::pair<int,int> degreeOrderPair = {degree, order};
                    auto it = reducedStokesInterpolators_[fileIndex].find(degreeOrderPair);

                    if ( it != reducedStokesInterpolators_[fileIndex].end() )
                    {
                        // Use batch interpolation with pre-computed state (much faster!)
                        cachedCosineCoefficients_(degree, order) = it->second.first->interpolateWithState(reducedInterpolationState);
                        if ( order > 0 )  // Sine coefficients only exist for m > 0
                        {
                            cachedSineCoefficients_(degree, order) = it->second.second->interpolateWithState(reducedInterpolationState);
                        }
                    }
                    else
                    {
                        // Fallback: create 1D reduced interpolator on-the-fly using helper method
                        createFallback1DInterpolator( fileIndex, degree, order,
                                                      cachedCosineCoefficients_(degree, order),
                                                      cachedSineCoefficients_(degree, order),
                                                      solarLongitude );
                    }
                }
            }
        }
        else
        {
            // No interpolators available - use fallback for all coefficients
            for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
            {
                for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                {
                    createFallback1DInterpolator( fileIndex, degree, order,
                                                  cachedCosineCoefficients_(degree, order),
                                                  cachedSineCoefficients_(degree, order),
                                                  solarLongitude );
                }
            }
        }

        // Apply decay term to the reduced coefficients
        simulation_setup::StokesCoefficientsEvaluator::applyDecayTerm(cachedCosineCoefficients_, radius, referenceRadius);
    }

        // Update interpolation cache after successful coefficient computation
        cachedRadius_ = radius;
        cachedInterpolationSolarLongitude_ = solarLongitude;
        cacheFlags_.interpolationValid = true;
    }
    // else: Cache hit - reuse previously computed coefficient matrices

    // Step 5: Compute number density using spherical harmonics expansion
    // Note: This step always runs (can't cache because lat/lon may change)
    // Update cached position values for density cache
    cachedLatitude_ = latitude;
    cachedLongitude_ = longitude;

    return sphericalHarmonicsCalculator_->calculateSurfaceSphericalHarmonics(
        cachedSineCoefficients_, cachedCosineCoefficients_,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );
}

//-----------------------------------------------------------------------------
// Private Helper Methods - Initialization
//-----------------------------------------------------------------------------

/*!
 * \brief Initialize interpolators for Stokes coefficients.
 * Pre-computes multilinear interpolators for all spherical harmonic coefficient
 * pairs to significantly improve performance during density evaluations.
 */
void ComaModel::initializeStokesInterpolators()
{
    // This function is only called from the Stokes coefficients constructor,
    // so dataType_ is guaranteed to be STOKES_COEFFICIENTS and stokesDataset_ is guaranteed to be non-null

    // Get dataset properties
    const auto& radiiGrid = stokesDataset_->radii();
    const auto& longitudeGrid = stokesDataset_->lons();
    const int nmax = stokesDataset_->nmax();
    const std::size_t nFiles = stokesDataset_->nFiles();

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Resize interpolator vectors to accommodate all files
    stokesInterpolators_.resize( nFiles );
    reducedStokesInterpolators_.resize( nFiles );
    fallbackStokesInterpolators_.resize( nFiles );
    fallbackReducedStokesInterpolators_.resize( nFiles );

    // Set up interpolation grids (shared for all interpolators)
    std::vector<std::vector<double>> independentGrids(2);
    independentGrids[0] = radiiGrid;     // Radius grid
    independentGrids[1] = longitudeGrid; // Solar longitude grid

    // Initialize interpolators for each file
    for ( std::size_t fileIndex = 0; fileIndex < nFiles; ++fileIndex )
    {
        // Initialize interpolators for each (n,m) pair for this file
        for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
        {
            for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
            {
                // Create 2D grids for this coefficient
                const std::size_t numRadii = radiiGrid.size();
                const std::size_t numLongitudes = longitudeGrid.size();

                boost::multi_array<double, 2> cosineGrid(boost::extents[numRadii][numLongitudes]);
                boost::multi_array<double, 2> sineGrid(boost::extents[numRadii][numLongitudes]);

                // Fill grids with coefficient values for this file
                for ( std::size_t radiusIndex = 0; radiusIndex < numRadii; ++radiusIndex )
                {
                    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
                    {
                        auto coeffs = stokesDataset_->getCoeff(fileIndex, radiusIndex, longitudeIndex, degree, order);
                        cosineGrid[radiusIndex][longitudeIndex] = coeffs.first;  // Cosine coefficient
                        sineGrid[radiusIndex][longitudeIndex] = coeffs.second;   // Sine coefficient
                    }
                }

                // Create and store interpolators for this coefficient and file
                auto cosineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 2>>(
                    independentGrids, cosineGrid,
                    interpolators::huntingAlgorithm,
                    interpolators::extrapolate_at_boundary
                );

                auto sineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 2>>(
                    independentGrids, sineGrid,
                    interpolators::huntingAlgorithm,
                    interpolators::extrapolate_at_boundary
                );

                std::pair<int,int> degreeOrderPair = {degree, order};
                stokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
            }
        }

        // Initialize reduced interpolators for coefficients beyond reference radius
        // These are 1D interpolators (solar longitude only)
        std::vector<std::vector<double>> reducedIndependentGrids(1);
        reducedIndependentGrids[0] = longitudeGrid; // Solar longitude grid only

        for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
        {
            for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
            {
                const std::size_t numLongitudes = longitudeGrid.size();

                boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
                boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

                // Fill grids with reduced coefficient values for this file
                for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
                {
                    auto reducedCoeffs = stokesDataset_->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
                    reducedCosineGrid[longitudeIndex] = reducedCoeffs.first;  // Cosine coefficient
                    reducedSineGrid[longitudeIndex] = reducedCoeffs.second;   // Sine coefficient
                }

                // Create and store reduced interpolators for this coefficient and file
                auto reducedCosineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 1>>(
                    reducedIndependentGrids, reducedCosineGrid,
                    interpolators::huntingAlgorithm,
                    interpolators::extrapolate_at_boundary
                );

                auto reducedSineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 1>>(
                    reducedIndependentGrids, reducedSineGrid,
                    interpolators::huntingAlgorithm,
                    interpolators::extrapolate_at_boundary
                );

                std::pair<int,int> degreeOrderPair = {degree, order};
                reducedStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Private Helper Methods - Fallback Interpolator Creation
//-----------------------------------------------------------------------------

/*!
 * \brief Helper to create 2D interpolator for Stokes coefficients on-the-fly.
 * This is a fallback method used when pre-initialized interpolators are not available.
 * Created interpolators are cached for reuse to avoid redundant construction.
 */
void ComaModel::createFallback2DInterpolator( const int fileIndex, const int degree, const int order,
                                              double& cosineCoeff, double& sineCoeff,
                                              const double radius, const double solarLongitude ) const
{
    std::pair<int,int> degreeOrderPair = {degree, order};

    // Check if we already created this interpolator
    auto it = fallbackStokesInterpolators_[fileIndex].find(degreeOrderPair);
    if ( it != fallbackStokesInterpolators_[fileIndex].end() )
    {
        // Cache hit - reuse existing interpolator
        interpolationPoint2D_[0] = radius;
        interpolationPoint2D_[1] = solarLongitude;
        cosineCoeff = it->second.first->interpolate(interpolationPoint2D_);
        if ( order > 0 )
        {
            sineCoeff = it->second.second->interpolate(interpolationPoint2D_);
        }
        else
        {
            sineCoeff = 0.0;
        }
        return;
    }

    // Cache miss - create and store the interpolator
    const auto& radiiGrid = stokesDataset_->radii();
    const auto& longitudeGrid = stokesDataset_->lons();
    const std::size_t numRadii = radiiGrid.size();
    const std::size_t numLongitudes = longitudeGrid.size();

    std::vector<std::vector<double>> independentGrids(2);
    independentGrids[0] = radiiGrid;
    independentGrids[1] = longitudeGrid;

    boost::multi_array<double, 2> cosineGrid(boost::extents[numRadii][numLongitudes]);
    boost::multi_array<double, 2> sineGrid(boost::extents[numRadii][numLongitudes]);

    for ( std::size_t radiusIndex = 0; radiusIndex < numRadii; ++radiusIndex )
    {
        for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
        {
            auto coeffs = stokesDataset_->getCoeff(fileIndex, radiusIndex, longitudeIndex, degree, order);
            cosineGrid[radiusIndex][longitudeIndex] = coeffs.first;
            sineGrid[radiusIndex][longitudeIndex] = coeffs.second;
        }
    }

    auto cosineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 2>>(
        independentGrids, cosineGrid,
        interpolators::huntingAlgorithm,
        interpolators::extrapolate_at_boundary
    );

    auto sineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 2>>(
        independentGrids, sineGrid,
        interpolators::huntingAlgorithm,
        interpolators::extrapolate_at_boundary
    );

    interpolationPoint2D_[0] = radius;
    interpolationPoint2D_[1] = solarLongitude;
    cosineCoeff = cosineInterpolator->interpolate(interpolationPoint2D_);
    if ( order > 0 )
    {
        sineCoeff = sineInterpolator->interpolate(interpolationPoint2D_);
    }
    else
    {
        sineCoeff = 0.0;
    }

    // Store interpolators in cache for future use
    fallbackStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
}

/*!
 * \brief Helper to create 1D reduced interpolator for Stokes coefficients on-the-fly.
 * This is a fallback method used when pre-initialized reduced interpolators are not available.
 * Created interpolators are cached for reuse to avoid redundant construction.
 */
void ComaModel::createFallback1DInterpolator( const int fileIndex, const int degree, const int order,
                                              double& cosineCoeff, double& sineCoeff,
                                              const double solarLongitude ) const
{
    std::pair<int,int> degreeOrderPair = {degree, order};

    // Check if we already created this interpolator
    auto it = fallbackReducedStokesInterpolators_[fileIndex].find(degreeOrderPair);
    if ( it != fallbackReducedStokesInterpolators_[fileIndex].end() )
    {
        // Cache hit - reuse existing interpolator
        interpolationPoint1D_[0] = solarLongitude;
        cosineCoeff = it->second.first->interpolate(interpolationPoint1D_);
        if ( order > 0 )
        {
            sineCoeff = it->second.second->interpolate(interpolationPoint1D_);
        }
        else
        {
            sineCoeff = 0.0;
        }
        return;
    }

    // Cache miss - create and store the interpolator
    const auto& longitudeGrid = stokesDataset_->lons();
    const std::size_t numLongitudes = longitudeGrid.size();

    std::vector<std::vector<double>> reducedIndependentGrids(1);
    reducedIndependentGrids[0] = longitudeGrid;

    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
    boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
    {
        auto reducedCoeffs = stokesDataset_->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
        reducedCosineGrid[longitudeIndex] = reducedCoeffs.first;
        reducedSineGrid[longitudeIndex] = reducedCoeffs.second;
    }

    auto reducedCosineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 1>>(
        reducedIndependentGrids, reducedCosineGrid,
        interpolators::huntingAlgorithm,
        interpolators::extrapolate_at_boundary
    );

    auto reducedSineInterpolator = std::make_unique<interpolators::MultiLinearInterpolator<double, double, 1>>(
        reducedIndependentGrids, reducedSineGrid,
        interpolators::huntingAlgorithm,
        interpolators::extrapolate_at_boundary
    );

    interpolationPoint1D_[0] = solarLongitude;
    cosineCoeff = reducedCosineInterpolator->interpolate(interpolationPoint1D_);
    if ( order > 0 )
    {
        sineCoeff = reducedSineInterpolator->interpolate(interpolationPoint1D_);
    }
    else
    {
        sineCoeff = 0.0;
    }

    // Store interpolators in cache for future use
    fallbackReducedStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
}

//=============================================================================
// SphericalHarmonicsCalculator Class Implementation
//=============================================================================

/*!
 * \brief Constructor for spherical harmonics calculator.
 * Initializes the calculator with caching enabled for both Legendre polynomials
 * and trigonometric functions to optimize repeated evaluations.
 * \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional)
 */
SphericalHarmonicsCalculator::SphericalHarmonicsCalculator( std::string fixedReferenceFrame ) :
    fixedReferenceFrame_( std::move( fixedReferenceFrame ) ),
    sphericalHarmonicsCache_( basic_mathematics::SphericalHarmonicsCache( true, true) )
{
}

/*!
 * \brief Calculate spherical harmonics field value at surface location.
 * Evaluates the spherical harmonics expansion using provided coefficients
 * and cached Legendre polynomials for efficient computation.
 * \param sineCoefficients Sine coefficients matrix
 * \param cosineCoefficients Cosine coefficients matrix
 * \param latitude Geocentric latitude [rad]
 * \param longitude Longitude [rad]
 * \param highestDegree Highest spherical harmonic degree to include in the evaluation
 * \param highestOrder Highest spherical harmonic order to include in the evaluation
 * \return Evaluated surface field value at the given latitude and longitude
 * \throws std::runtime_error If coefficient matrices have incompatible sizes
 */
double SphericalHarmonicsCalculator::calculateSurfaceSphericalHarmonics(
        const Eigen::MatrixXd& sineCoefficients,
        const Eigen::MatrixXd& cosineCoefficients,
        const double latitude,
        const double longitude,
        const int highestDegree,
        const int highestOrder )
{
    // Validate input coefficient matrices
    if(cosineCoefficients.rows( ) != sineCoefficients.rows( ) ||
        cosineCoefficients.cols( ) != sineCoefficients.cols( ))
    {
        throw std::runtime_error( "Spherical harmonics coefficient sizes are incompatible." );
    }

    // Set up cache for spherical harmonics computation
    // Only reset if degree/order changed to avoid unnecessary cache invalidation
    const int maxDegree = static_cast< int >(cosineCoefficients.rows( ));
    const int maxOrder = static_cast< int >(cosineCoefficients.cols( ));
    if ( maxDegree != lastMaxDegree_ || maxOrder != lastMaxOrder_ )
    {
        sphericalHarmonicsCache_.resetMaximumDegreeAndOrder( maxDegree, maxOrder );
        lastMaxDegree_ = maxDegree;
        lastMaxOrder_ = maxOrder;
    }
    // else: Cache is still valid - skip expensive reset operation

    // Cache trigonometric computations to avoid redundant std::sin calls
    constexpr double angleToleranceSq = 1e-20; // (1e-10)^2
    const double latDiff = latitude - lastLatitude_;
    const double lonDiff = longitude - lastLongitude_;

    double sineOfAngle;
    if ( latDiff * latDiff < angleToleranceSq && lonDiff * lonDiff < angleToleranceSq )
    {
        // Cache hit - reuse previously computed sine
        sineOfAngle = lastSineLatitude_;
        // Longitude cache is handled internally by sphericalHarmonicsCache_
    }
    else
    {
        // Cache miss - compute and cache sine of latitude
        sineOfAngle = std::sin( latitude );
        lastLatitude_ = latitude;
        lastLongitude_ = longitude;
        lastSineLatitude_ = sineOfAngle;
    }

    // Update cache with current position
    sphericalHarmonicsCache_.updateAnglesOnly( sineOfAngle, longitude );
    const basic_mathematics::LegendreCache& legendreCacheReference = sphericalHarmonicsCache_.getLegendreCache( );

    // Compute spherical harmonics expansion with reduced function call overhead
    double value = 0.0;

    // Process m=0 terms separately (only cosine coefficients)
    for (int degree = 0; degree <= highestDegree; ++degree)
    {
        const double legendrePolynomial = legendreCacheReference.getLegendrePolynomial(degree, 0);
        value += legendrePolynomial * cosineCoefficients(degree, 0);
    }

    // Process m>0 terms with pre-fetched trigonometric values
    for (int order = 1; order <= highestOrder; ++order)
    {
        const double cosineMultipleLongitude = sphericalHarmonicsCache_.getCosineOfMultipleLongitude(order);
        const double sineMultipleLongitude = sphericalHarmonicsCache_.getSineOfMultipleLongitude(order);

        // Process all degrees for this order (degree >= order)
        for (int degree = order; degree <= highestDegree; ++degree)
        {
            const double legendrePolynomial = legendreCacheReference.getLegendrePolynomial(degree, order);
            value += legendrePolynomial * (cosineCoefficients(degree, order) * cosineMultipleLongitude +
                                          sineCoefficients(degree, order) * sineMultipleLongitude);
        }
    }

    return value;
}

} // namespace aerodynamics
} // namespace tudat