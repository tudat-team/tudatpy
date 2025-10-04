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
 * the given time value.
 * \param time The time for which the corresponding time interval index is to be found [s]
 * \return The index of the time period in which the given time falls
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaModel::findTimeIntervalIndex( const double time ) const
{
    if ( dataType_ == ComaDataType::POLYNOMIAL_COEFFICIENTS )
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
    }
    else if ( dataType_ == ComaDataType::STOKES_COEFFICIENTS )
    {
        for ( std::size_t f = 0; f < stokesDataset_->nFiles(); ++f )
        {
            const auto& fileMeta = stokesDataset_->files()[f];
            if ( time >= fileMeta.start_epoch && time <= fileMeta.end_epoch )
            {
                return static_cast<int>( f );
            }
        }
    }

    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}

/*!
 * \brief Calculate solar longitude in comet body-fixed frame.
 * Computes the angle from X-axis to Sun direction projected onto the XY plane,
 * which is used to account for solar heating variations.
 * \return Solar longitude angle from X-axis in XY plane [rad]
 * \throws std::runtime_error If state functions are not initialized
 */
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
double ComaModel::computeDensityFromPolyCoefficients( double radius, double longitude, double latitude, double time ) const
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

    // Calculate current solar longitude for heliocentric dependence
    const double solarLongitude = calculateSolarLongitude();

    // Evaluate polynomial coefficients to get spherical harmonic coefficients
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

    // Compute density using spherical harmonics expansion
    return sphericalHarmonicsCalculator_->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cosineCoefficients.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cosineCoefficients.cols() - 1 );
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
double ComaModel::computeDensityFromStokesCoefficients( double radius, double longitude, double latitude, double time ) const
{
    if ( !stokesDataset_ )
    {
        throw std::runtime_error( "ComaModel: stokesDataset_ is null" );
    }

    // Step 1: Get time-dependent properties
    const int fileIndex = findTimeIntervalIndex( time );
    const double solarLongitude = calculateSolarLongitude();

    // Step 2: Get dataset properties
    const int nmax = stokesDataset_->nmax();
    const double referenceRadius = stokesDataset_->getReferenceRadius(fileIndex);

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Step 3: Initialize coefficient matrices
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero(effectiveMaxDegree + 1, effectiveMaxOrder + 1);
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero(effectiveMaxDegree + 1, effectiveMaxOrder + 1);

    // Step 4: Choose interpolation strategy based on radius
    if (radius <= referenceRadius)
    {
        // Use regular 2D interpolators (radius, solar longitude)
        const auto& radiiGrid = stokesDataset_->radii();
        std::vector<double> interpolationPoint = {radius, solarLongitude};

        for ( int n = 0; n <= effectiveMaxDegree; ++n )
        {
            for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
            {
                std::pair<int,int> nmPair = {n, m};
                auto it = stokesInterpolators_.find(nmPair);

                if ( it != stokesInterpolators_.end() )
                {
                    // Use pre-initialized 2D interpolators
                    cosineCoefficients(n, m) = it->second.first->interpolate(interpolationPoint);
                    if ( m > 0 )  // Sine coefficients only exist for m > 0
                    {
                        sineCoefficients(n, m) = it->second.second->interpolate(interpolationPoint);
                    }
                }
                else
                {
                    // Fallback: create 2D interpolator on-the-fly
                    const auto& longitudeGrid = stokesDataset_->lons();
                    const std::size_t nRadii = radiiGrid.size();
                    const std::size_t nLons = longitudeGrid.size();

                    std::vector<std::vector<double>> independentGrids(2);
                    independentGrids[0] = radiiGrid;
                    independentGrids[1] = longitudeGrid;

                    boost::multi_array<double, 2> cosineGrid(boost::extents[nRadii][nLons]);
                    boost::multi_array<double, 2> sineGrid(boost::extents[nRadii][nLons]);

                    for ( std::size_t r = 0; r < nRadii; ++r )
                    {
                        for ( std::size_t l = 0; l < nLons; ++l )
                        {
                            auto coeffs = stokesDataset_->getCoeff(fileIndex, r, l, n, m);
                            cosineGrid[r][l] = coeffs.first;
                            sineGrid[r][l] = coeffs.second;
                        }
                    }

                    interpolators::MultiLinearInterpolator<double, double, 2> cosineInterpolator(
                        independentGrids, cosineGrid,
                        interpolators::huntingAlgorithm,
                        interpolators::extrapolate_at_boundary
                    );

                    interpolators::MultiLinearInterpolator<double, double, 2> sineInterpolator(
                        independentGrids, sineGrid,
                        interpolators::huntingAlgorithm,
                        interpolators::extrapolate_at_boundary
                    );

                    cosineCoefficients(n, m) = cosineInterpolator.interpolate(interpolationPoint);
                    if ( m > 0 )
                    {
                        sineCoefficients(n, m) = sineInterpolator.interpolate(interpolationPoint);
                    }
                }
            }
        }
    }
    else
    {
        // Use reduced 1D interpolators (solar longitude only) + apply decay term
        std::vector<double> reducedInterpolationPoint = {solarLongitude};

        for ( int n = 0; n <= effectiveMaxDegree; ++n )
        {
            for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
            {
                std::pair<int,int> nmPair = {n, m};
                auto it = reducedStokesInterpolators_.find(nmPair);

                if ( it != reducedStokesInterpolators_.end() )
                {
                    // Use pre-initialized 1D reduced interpolators
                    cosineCoefficients(n, m) = it->second.first->interpolate(reducedInterpolationPoint);
                    if ( m > 0 )  // Sine coefficients only exist for m > 0
                    {
                        sineCoefficients(n, m) = it->second.second->interpolate(reducedInterpolationPoint);
                    }
                }
                else
                {
                    // Fallback: create 1D reduced interpolator on-the-fly
                    const auto& longitudeGrid = stokesDataset_->lons();
                    const std::size_t nLons = longitudeGrid.size();

                    std::vector<std::vector<double>> reducedIndependentGrids(1);
                    reducedIndependentGrids[0] = longitudeGrid;

                    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[nLons]);
                    boost::multi_array<double, 1> reducedSineGrid(boost::extents[nLons]);

                    for ( std::size_t l = 0; l < nLons; ++l )
                    {
                        auto reducedCoeffs = stokesDataset_->getReducedCoeff(fileIndex, l, n, m);
                        reducedCosineGrid[l] = reducedCoeffs.first;
                        reducedSineGrid[l] = reducedCoeffs.second;
                    }

                    interpolators::MultiLinearInterpolator<double, double, 1> reducedCosineInterpolator(
                        reducedIndependentGrids, reducedCosineGrid,
                        interpolators::huntingAlgorithm,
                        interpolators::extrapolate_at_boundary
                    );

                    interpolators::MultiLinearInterpolator<double, double, 1> reducedSineInterpolator(
                        reducedIndependentGrids, reducedSineGrid,
                        interpolators::huntingAlgorithm,
                        interpolators::extrapolate_at_boundary
                    );

                    cosineCoefficients(n, m) = reducedCosineInterpolator.interpolate(reducedInterpolationPoint);
                    if ( m > 0 )
                    {
                        sineCoefficients(n, m) = reducedSineInterpolator.interpolate(reducedInterpolationPoint);
                    }
                }
            }
        }

        // Apply decay term to the reduced coefficients
        simulation_setup::StokesCoefficientsEvaluator::applyDecayTerm(cosineCoefficients, radius, referenceRadius);
    }

    // Step 5: Compute density using spherical harmonics expansion
    return sphericalHarmonicsCalculator_->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
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

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Set up interpolation grids (shared for all interpolators)
    std::vector<std::vector<double>> independentGrids(2);
    independentGrids[0] = radiiGrid;     // Radius grid
    independentGrids[1] = longitudeGrid; // Solar longitude grid

    // Initialize interpolators for each (n,m) pair
    for ( int n = 0; n <= effectiveMaxDegree; ++n )
    {
        for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
        {
            // Create 2D grids for this coefficient (using first file for now)
            const std::size_t nRadii = radiiGrid.size();
            const std::size_t nLons = longitudeGrid.size();
            const int fileIndex = 0; // TODO: Could be extended for multi-file support

            boost::multi_array<double, 2> cosineGrid(boost::extents[nRadii][nLons]);
            boost::multi_array<double, 2> sineGrid(boost::extents[nRadii][nLons]);

            // Fill grids with coefficient values
            for ( std::size_t r = 0; r < nRadii; ++r )
            {
                for ( std::size_t l = 0; l < nLons; ++l )
                {
                    auto coeffs = stokesDataset_->getCoeff(fileIndex, r, l, n, m);
                    cosineGrid[r][l] = coeffs.first;  // Cosine coefficient
                    sineGrid[r][l] = coeffs.second;   // Sine coefficient
                }
            }

            // Create and store interpolators for this coefficient
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

            std::pair<int,int> nmPair = {n, m};
            stokesInterpolators_[nmPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
        }
    }

    // Initialize reduced interpolators for coefficients beyond reference radius
    // These are 1D interpolators (solar longitude only)
    std::vector<std::vector<double>> reducedIndependentGrids(1);
    reducedIndependentGrids[0] = longitudeGrid; // Solar longitude grid only

    for ( int n = 0; n <= effectiveMaxDegree; ++n )
    {
        for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
        {
            const std::size_t nLons = longitudeGrid.size();
            const int fileIndex = 0; // TODO: Could be extended for multi-file support

            boost::multi_array<double, 1> reducedCosineGrid(boost::extents[nLons]);
            boost::multi_array<double, 1> reducedSineGrid(boost::extents[nLons]);

            // Fill grids with reduced coefficient values
            for ( std::size_t l = 0; l < nLons; ++l )
            {
                auto reducedCoeffs = stokesDataset_->getReducedCoeff(fileIndex, l, n, m);
                reducedCosineGrid[l] = reducedCoeffs.first;  // Cosine coefficient
                reducedSineGrid[l] = reducedCoeffs.second;   // Sine coefficient
            }

            // Create and store reduced interpolators for this coefficient
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

            std::pair<int,int> nmPair = {n, m};
            reducedStokesInterpolators_[nmPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
        }
    }
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
    const int maxDegree = static_cast< int >(cosineCoefficients.rows( ));
    const int maxOrder = static_cast< int >(cosineCoefficients.cols( ));
    sphericalHarmonicsCache_.resetMaximumDegreeAndOrder( maxDegree, maxOrder );

    // Update cache with current position
    const double sineOfAngle = std::sin( latitude );
    sphericalHarmonicsCache_.updateAnglesOnly( sineOfAngle, longitude );
    const basic_mathematics::LegendreCache& legendreCacheReference = sphericalHarmonicsCache_.getLegendreCache( );

    // Compute spherical harmonics expansion
    double value = 0.0;
    for (int l = 0; l <= highestDegree; ++l)
    {
        const int mmax = std::min(l, highestOrder);
        for (int m = 0; m <= mmax; ++m)
        {
            const double P = legendreCacheReference.getLegendrePolynomial(l, m);
            const double cos_mλ = sphericalHarmonicsCache_.getCosineOfMultipleLongitude(m);

            if (m == 0)
            {
                value += P * cosineCoefficients(l, 0);
            }
            else
            {
                const double sin_mλ = sphericalHarmonicsCache_.getSineOfMultipleLongitude(m);
                value += P * (cosineCoefficients(l, m) * cos_mλ + sineCoefficients(l, m) * sin_mλ);
            }
        }
    }

    return value;
}

} // namespace aerodynamics
} // namespace tudat