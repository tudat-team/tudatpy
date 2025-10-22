#include "tudat/astro/aerodynamics/windModel.h"
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

/*!
 * \brief Constructor for ComaWindModel using polynomial coefficient datasets.
 * Initializes the wind model with separate datasets for x, y, and z wind components
 * using polynomial coefficients. The model can optionally share a spherical harmonics
 * calculator with the provided ComaModel for efficiency.
 * \param xPolyDataset Dataset for x-component wind (polynomial coefficients)
 * \param yPolyDataset Dataset for y-component wind (polynomial coefficients)
 * \param zPolyDataset Dataset for z-component wind (polynomial coefficients)
 * \param comaModel Shared pointer to the ComaModel for accessing state functions
 * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
 * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
 * \param cometRotationFunction Function returning comet body-fixed rotation matrix
 * \param maximumDegree Maximum degree used for computation (-1 for auto)
 * \param maximumOrder Maximum order used for computation (-1 for auto)
 * \param associatedFrame Reference frame for the wind model
 * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
 */
ComaWindModel::ComaWindModel( const simulation_setup::ComaPolyDataset& xPolyDataset,
                   const simulation_setup::ComaPolyDataset& yPolyDataset,
                   const simulation_setup::ComaPolyDataset& zPolyDataset,
                   const std::shared_ptr<ComaModel>& comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree,
                   const int& maximumOrder,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame,
                   const bool includeCorotation,
                   const bool useRadius ) :
        WindModel( associatedFrame, includeCorotation, useRadius ),
        dataType_( 0 ), // POLYNOMIAL_COEFFICIENTS
        maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder ),
        xPolyDataset_( std::make_shared<simulation_setup::ComaPolyDataset>( xPolyDataset ) ),
        yPolyDataset_( std::make_shared<simulation_setup::ComaPolyDataset>( yPolyDataset ) ),
        zPolyDataset_( std::make_shared<simulation_setup::ComaPolyDataset>( zPolyDataset ) ),
        xStokesDataset_( nullptr ),
        yStokesDataset_( nullptr ),
        zStokesDataset_( nullptr ),
        comaModel_( comaModel ),
        sunStateFunction_( std::move( sunStateFunction ) ),
        cometStateFunction_( std::move( cometStateFunction ) ),
        cometRotationFunction_( std::move( cometRotationFunction ) ),
        sphericalHarmonicsCalculator_( nullptr ),
        sharedSphericalHarmonicsCalculator_( nullptr ),
        cachedSolarLongitude_( 0.0 ),
        cachedTime_( -std::numeric_limits<double>::infinity() ),
        lastXFileIndex_( 0 ),
        lastYFileIndex_( 0 ),
        lastZFileIndex_( 0 ),
        cachedRadius_( 0.0 ),
        cachedEffectiveRadius_( 0.0 ),
        cachedInterpolationSolarLongitude_( 0.0 ),
        coefficientMatricesSized_( false ),
        cachedLatitude_( 0.0 ),
        cachedLongitude_( 0.0 ),
        cachedFinalWindVector_( Eigen::Vector3d::Zero() ),
        interpolationPoint2D_( 2 ),
        interpolationPoint1D_( 1 )
    {
        if ( !comaModel_ )
        {
            throw std::invalid_argument( "ComaWindModel: ComaModel must be provided" );
        }

        if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
        {
            throw std::invalid_argument( "ComaWindModel: All state functions must be provided" );
        }

        if ( xPolyDataset_->getNumFiles() == 0 || yPolyDataset_->getNumFiles() == 0 || zPolyDataset_->getNumFiles() == 0 )
        {
            throw std::invalid_argument( "ComaWindModel: All PolyDatasets must contain at least one file" );
        }

        if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
        {
            throw std::invalid_argument( "ComaWindModel: Maximum degree and order must be >= -1" );
        }

        // Try to share the spherical harmonics calculator with the ComaModel
        if ( comaModel_ && comaModel_->getSphericalHarmonicsCalculator() != nullptr )
        {
            // Share the calculator from ComaModel (non-owning pointer)
            sharedSphericalHarmonicsCalculator_ = comaModel_->getSphericalHarmonicsCalculator();
        }
        else
        {
            // Create our own calculator
            sphericalHarmonicsCalculator_ = std::make_unique<SphericalHarmonicsCalculator>();
        }

        // Pre-allocate coefficient matrices based on maximum degree/order from dataset
        const int maxDegreeAvailable = xPolyDataset_->getMaxDegreeSH( 0 );
        const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : maxDegreeAvailable;
        const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : maxDegreeAvailable;
        cachedXCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedXSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
    }

/*!
 * \brief Constructor for ComaWindModel using Stokes coefficient datasets.
 * Initializes the wind model with separate datasets for x, y, and z wind components
 * using Stokes coefficients. Pre-initializes interpolators for efficient evaluation
 * and can share a spherical harmonics calculator with the provided ComaModel.
 * \param xStokesDataset Dataset for x-component wind (Stokes coefficients)
 * \param yStokesDataset Dataset for y-component wind (Stokes coefficients)
 * \param zStokesDataset Dataset for z-component wind (Stokes coefficients)
 * \param comaModel Shared pointer to the ComaModel for accessing state functions
 * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
 * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
 * \param cometRotationFunction Function returning comet body-fixed rotation matrix
 * \param maximumDegree Maximum degree used for computation (-1 for auto)
 * \param maximumOrder Maximum order used for computation (-1 for auto)
 * \param associatedFrame Reference frame for the wind model
 * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
 */
ComaWindModel::ComaWindModel( const simulation_setup::ComaStokesDataset& xStokesDataset,
                   const simulation_setup::ComaStokesDataset& yStokesDataset,
                   const simulation_setup::ComaStokesDataset& zStokesDataset,
                   std::shared_ptr<ComaModel> comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree,
                   const int& maximumOrder,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame,
                   const bool includeCorotation,
                   const bool useRadius ) :
        WindModel( associatedFrame, includeCorotation, useRadius ),
        dataType_( 1 ), // STOKES_COEFFICIENTS
        maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder ),
        xPolyDataset_( nullptr ),
        yPolyDataset_( nullptr ),
        zPolyDataset_( nullptr ),
        xStokesDataset_( std::make_shared<simulation_setup::ComaStokesDataset>( xStokesDataset ) ),
        yStokesDataset_( std::make_shared<simulation_setup::ComaStokesDataset>( yStokesDataset ) ),
        zStokesDataset_( std::make_shared<simulation_setup::ComaStokesDataset>( zStokesDataset ) ),
        comaModel_( comaModel ),
        sunStateFunction_( std::move( sunStateFunction ) ),
        cometStateFunction_( std::move( cometStateFunction ) ),
        cometRotationFunction_( std::move( cometRotationFunction ) ),
        sphericalHarmonicsCalculator_( nullptr ),
        sharedSphericalHarmonicsCalculator_( nullptr ),
        cachedSolarLongitude_( 0.0 ),
        cachedTime_( -std::numeric_limits<double>::infinity() ),
        lastXFileIndex_( 0 ),
        lastYFileIndex_( 0 ),
        lastZFileIndex_( 0 ),
        cachedRadius_( 0.0 ),
        cachedEffectiveRadius_( 0.0 ),
        cachedInterpolationSolarLongitude_( 0.0 ),
        coefficientMatricesSized_( false ),
        cachedLatitude_( 0.0 ),
        cachedLongitude_( 0.0 ),
        cachedFinalWindVector_( Eigen::Vector3d::Zero() ),
        interpolationPoint2D_( 2 ),
        interpolationPoint1D_( 1 )
    {
        if ( !comaModel_ )
        {
            throw std::invalid_argument( "ComaWindModel: ComaModel must be provided" );
        }

        if ( !sunStateFunction_ || !cometStateFunction_ || !cometRotationFunction_ )
        {
            throw std::invalid_argument( "ComaWindModel: All state functions must be provided" );
        }

        if ( xStokesDataset_->nFiles() == 0 || yStokesDataset_->nFiles() == 0 || zStokesDataset_->nFiles() == 0 )
        {
            throw std::invalid_argument( "ComaWindModel: All StokesDatasets must contain at least one file" );
        }

        if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
        {
            throw std::invalid_argument( "ComaWindModel: Maximum degree and order must be >= -1" );
        }

        // Try to share the spherical harmonics calculator with the ComaModel
        if ( comaModel_ && comaModel_->getSphericalHarmonicsCalculator() != nullptr )
        {
            // Share the calculator from ComaModel (non-owning pointer)
            sharedSphericalHarmonicsCalculator_ = comaModel_->getSphericalHarmonicsCalculator();
        }
        else
        {
            // Create our own calculator
            sphericalHarmonicsCalculator_ = std::make_unique<SphericalHarmonicsCalculator>();
        }

        // Pre-allocate coefficient matrices based on nmax from dataset
        const int nmax = xStokesDataset_->nmax();
        const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : nmax;
        const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : nmax;
        cachedXCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedXSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );

        // Initialize interpolators for Stokes coefficients
        initializeStokesInterpolators();
    }

/*!
 * \brief Compute the current wind velocity vector in body-fixed coordinates with caching.
 * This method evaluates the wind velocity at the specified location and time
 * using either polynomial or Stokes coefficient data depending on the model type.
 * Results are cached to avoid redundant computations when inputs haven't changed.
 * \param currentAltitude Altitude at which wind vector is to be retrieved [m]
 * \param currentLongitude Longitude at which wind vector is to be retrieved [rad]
 * \param currentLatitude Latitude at which wind vector is to be retrieved [rad]
 * \param currentTime Time at which wind vector is to be retrieved [s]
 * \return Wind velocity vector in body-fixed, body-centered frame [m/s]
 */
Eigen::Vector3d ComaWindModel::getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime )
    {
        // Convert altitude to radius (assuming currentAltitude is already radius from center)
        const double radius = currentAltitude;

        // Check if we can reuse cached final wind vector result
        constexpr double tolerance = 1e-10;
        constexpr double toleranceSq = tolerance * tolerance;

        const double radiusDiff = radius - cachedRadius_;
        const double lonDiff = currentLongitude - cachedLongitude_;
        const double latDiff = currentLatitude - cachedLatitude_;
        const double timeDiff = currentTime - cachedTime_;

        if ( cacheFlags_.windVectorValid &&
             radiusDiff * radiusDiff < toleranceSq &&
             lonDiff * lonDiff < toleranceSq &&
             latDiff * latDiff < toleranceSq &&
             timeDiff * timeDiff < toleranceSq )
        {
            return cachedFinalWindVector_;
        }

        // Compute wind velocity vector using optimized vectorized methods
        switch ( dataType_ )
        {
            case 0: // POLYNOMIAL_COEFFICIENTS
                cachedFinalWindVector_ = computeWindVectorFromPolyCoefficients( radius, currentLongitude, currentLatitude, currentTime );
                break;

            case 1: // STOKES_COEFFICIENTS
                cachedFinalWindVector_ = computeWindVectorFromStokesCoefficients( radius, currentLongitude, currentLatitude, currentTime );
                break;

            default:
                throw std::runtime_error( "ComaWindModel: Unknown data type" );
        }

        // Cache the input parameters for future validation
        cachedRadius_ = radius;
        cachedLongitude_ = currentLongitude;
        cachedLatitude_ = currentLatitude;
        cachedTime_ = currentTime;
        cacheFlags_.windVectorValid = true;

        return cachedFinalWindVector_;
    }

/*!
 * \brief Compute complete wind vector from Stokes coefficients (vectorized, all 3 components).
 * This optimized method computes all wind components in a single pass, using batch interpolation
 * and sharing caches for maximum performance. Eliminates redundant computations and cache misses.
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute the wind vector [s]
 * \return Complete wind velocity vector [m/s]
 * \throws std::runtime_error If datasets are null or time is out of range
 */
Eigen::Vector3d ComaWindModel::computeWindVectorFromStokesCoefficients(
    double radius, double longitude, double latitude, double time ) const
{
    if ( !xStokesDataset_ || !yStokesDataset_ || !zStokesDataset_ )
    {
        throw std::runtime_error( "ComaWindModel: Stokes datasets are null" );
    }

    // Step 1: Get time-dependent properties (shared for all components)
    const int fileIndex = findTimeIntervalIndex( time, xStokesDataset_ );  // Assume all datasets have same time structure
    const double solarLongitude = calculateSolarLongitude( time );

    // Step 2: Get dataset properties
    const int nmax = xStokesDataset_->nmax();
    const double referenceRadius = xStokesDataset_->getReferenceRadius(fileIndex);

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Clamp radius to reference radius - wind velocity is constant beyond reference radius
    const double effectiveRadius = std::min(radius, referenceRadius);

    // Check if coefficient interpolation can be skipped (same radius/solar longitude)
    constexpr double radiusTolerance = 1e-10;
    constexpr double solarLongitudeTolerance = 1e-10;
    constexpr double radiusToleranceSq = radiusTolerance * radiusTolerance;
    constexpr double solarLongitudeToleranceSq = solarLongitudeTolerance * solarLongitudeTolerance;

    const double radiusDiff = effectiveRadius - cachedEffectiveRadius_;
    const double solarLongitudeDiff = solarLongitude - cachedInterpolationSolarLongitude_;

    const bool sameRadius = cacheFlags_.interpolationValid &&
                           radiusDiff * radiusDiff < radiusToleranceSq;
    const bool sameSolarLongitude = cacheFlags_.interpolationValid &&
                                   solarLongitudeDiff * solarLongitudeDiff < solarLongitudeToleranceSq;

    if ( !sameRadius || !sameSolarLongitude )
    {
        // Cache miss - need to recompute coefficients for all 3 components
        // Step 3: Ensure all coefficient matrices have correct dimensions
        for (auto* cosineCoeffs : {&cachedXCosineCoefficients_, &cachedYCosineCoefficients_, &cachedZCosineCoefficients_})
        {
            if ( cosineCoeffs->rows() == effectiveMaxDegree + 1 &&
                 cosineCoeffs->cols() == effectiveMaxOrder + 1 )
            {
                cosineCoeffs->setZero();
            }
            else
            {
                *cosineCoeffs = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
            }
        }
        for (auto* sineCoeffs : {&cachedXSineCoefficients_, &cachedYSineCoefficients_, &cachedZSineCoefficients_})
        {
            if ( sineCoeffs->rows() == effectiveMaxDegree + 1 &&
                 sineCoeffs->cols() == effectiveMaxOrder + 1 )
            {
                sineCoeffs->setZero();
            }
            else
            {
                *sineCoeffs = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
            }
        }

        // Step 4: Interpolate coefficients using 2D interpolation (radius is clamped to referenceRadius)
        {
            interpolationPoint2D_[0] = effectiveRadius;
            interpolationPoint2D_[1] = solarLongitude;

            // Prepare interpolation state once for all coefficients (batch mode)
            if ( !xStokesInterpolators_.empty() )
            {
                // Use the first available interpolator to prepare the shared interpolation state
                auto firstInterpolator = xStokesInterpolators_.begin()->second.first.get();
                auto interpolationState = firstInterpolator->prepareInterpolationState(interpolationPoint2D_);

                // Now interpolate ALL THREE components using the same pre-computed state
                for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
                {
                    for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                    {
                        std::pair<int,int> degreeOrderPair = {degree, order};

                        // X-component
                        auto itX = xStokesInterpolators_.find(degreeOrderPair);
                        if ( itX != xStokesInterpolators_.end() )
                        {
                            cachedXCosineCoefficients_(degree, order) = itX->second.first->interpolateWithState(interpolationState);
                            if ( order > 0 )
                            {
                                cachedXSineCoefficients_(degree, order) = itX->second.second->interpolateWithState(interpolationState);
                            }
                        }
                        else
                        {
                            createFallback2DInterpolator( xStokesDataset_, xFallbackStokesInterpolators_, fileIndex, degree, order,
                                                          cachedXCosineCoefficients_(degree, order),
                                                          cachedXSineCoefficients_(degree, order),
                                                          effectiveRadius, solarLongitude );
                        }

                        // Y-component
                        auto itY = yStokesInterpolators_.find(degreeOrderPair);
                        if ( itY != yStokesInterpolators_.end() )
                        {
                            cachedYCosineCoefficients_(degree, order) = itY->second.first->interpolateWithState(interpolationState);
                            if ( order > 0 )
                            {
                                cachedYSineCoefficients_(degree, order) = itY->second.second->interpolateWithState(interpolationState);
                            }
                        }
                        else
                        {
                            createFallback2DInterpolator( yStokesDataset_, yFallbackStokesInterpolators_, fileIndex, degree, order,
                                                          cachedYCosineCoefficients_(degree, order),
                                                          cachedYSineCoefficients_(degree, order),
                                                          effectiveRadius, solarLongitude );
                        }

                        // Z-component
                        auto itZ = zStokesInterpolators_.find(degreeOrderPair);
                        if ( itZ != zStokesInterpolators_.end() )
                        {
                            cachedZCosineCoefficients_(degree, order) = itZ->second.first->interpolateWithState(interpolationState);
                            if ( order > 0 )
                            {
                                cachedZSineCoefficients_(degree, order) = itZ->second.second->interpolateWithState(interpolationState);
                            }
                        }
                        else
                        {
                            createFallback2DInterpolator( zStokesDataset_, zFallbackStokesInterpolators_, fileIndex, degree, order,
                                                          cachedZCosineCoefficients_(degree, order),
                                                          cachedZSineCoefficients_(degree, order),
                                                          effectiveRadius, solarLongitude );
                        }
                    }
                }
            }
            else
            {
                // No interpolators available - use fallback for all coefficients and components
                for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
                {
                    for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                    {
                        createFallback2DInterpolator( xStokesDataset_, xFallbackStokesInterpolators_, fileIndex, degree, order,
                                                      cachedXCosineCoefficients_(degree, order),
                                                      cachedXSineCoefficients_(degree, order),
                                                      effectiveRadius, solarLongitude );
                        createFallback2DInterpolator( yStokesDataset_, yFallbackStokesInterpolators_, fileIndex, degree, order,
                                                      cachedYCosineCoefficients_(degree, order),
                                                      cachedYSineCoefficients_(degree, order),
                                                      effectiveRadius, solarLongitude );
                        createFallback2DInterpolator( zStokesDataset_, zFallbackStokesInterpolators_, fileIndex, degree, order,
                                                      cachedZCosineCoefficients_(degree, order),
                                                      cachedZSineCoefficients_(degree, order),
                                                      effectiveRadius, solarLongitude );
                    }
                }
            }
        }

        // Update interpolation cache after successful coefficient computation
        cachedEffectiveRadius_ = effectiveRadius;
        cachedInterpolationSolarLongitude_ = solarLongitude;
        cacheFlags_.interpolationValid = true;
    }
    // else: Cache hit - reuse previously computed coefficient matrices

    // Step 5: Compute wind components using spherical harmonics expansion
    // Note: This step always runs (can't cache because lat/lon may change between calls)
    auto* shCalc = getActiveSphericalHarmonicsCalculator();

    const double windX = shCalc->calculateSurfaceSphericalHarmonics(
        cachedXSineCoefficients_, cachedXCosineCoefficients_,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );

    const double windY = shCalc->calculateSurfaceSphericalHarmonics(
        cachedYSineCoefficients_, cachedYCosineCoefficients_,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );

    const double windZ = shCalc->calculateSurfaceSphericalHarmonics(
        cachedZSineCoefficients_, cachedZCosineCoefficients_,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );

    return Eigen::Vector3d( windX, windY, windZ );
}

/*!
 * \brief Compute complete wind vector from polynomial coefficients (vectorized, all 3 components).
 * This optimized method computes all wind components in a single pass, sharing caches
 * and reducing overhead. Preferred over calling single-component methods.
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute the wind vector [s]
 * \return Complete wind velocity vector [m/s]
 * \throws std::runtime_error If datasets are null or time is out of range
 */
Eigen::Vector3d ComaWindModel::computeWindVectorFromPolyCoefficients(
    double radius, double longitude, double latitude, double time ) const
{
    if ( !xPolyDataset_ || !yPolyDataset_ || !zPolyDataset_ )
    {
        throw std::runtime_error( "ComaWindModel: Polynomial datasets are null" );
    }

    // Get time-dependent data and solar longitude (shared for all components)
    const int fileIndex = findTimeIntervalIndex( time, xPolyDataset_ );  // Assume all datasets have same time structure
    const auto& fileMeta = xPolyDataset_->getFileMeta( fileIndex );
    const double solarLongitude = calculateSolarLongitude( time );

    // Only resize/zero coefficient matrices if dimensions changed
    const int maxDegreeAvailable = fileMeta.maxDegreeSH;
    const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : maxDegreeAvailable;
    const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : maxDegreeAvailable;

    // Check sizing flag first to avoid repeated size comparisons
    if ( !coefficientMatricesSized_ ||
         cachedXCosineCoefficients_.rows() != effectiveMaxDegree + 1 ||
         cachedXCosineCoefficients_.cols() != effectiveMaxOrder + 1 )
    {
        // Resize all component matrices
        cachedXCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedXSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedYSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZCosineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        cachedZSineCoefficients_ = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        coefficientMatricesSized_ = true;
    }
    else
    {
        // Same size - just zero the contents (faster than reallocation)
        cachedXCosineCoefficients_.setZero();
        cachedXSineCoefficients_.setZero();
        cachedYCosineCoefficients_.setZero();
        cachedYSineCoefficients_.setZero();
        cachedZCosineCoefficients_.setZero();
        cachedZSineCoefficients_.setZero();
    }

    // Evaluate polynomial coefficients for all three components
    const auto& xPolyCoefficients = xPolyDataset_->getPolyCoefficients( fileIndex );
    const auto& xShIndices = xPolyDataset_->getSHDegreeAndOrderIndices( fileIndex );
    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius, solarLongitude,
        xPolyCoefficients, xShIndices,
        fileMeta.powersInvRadius, fileMeta.referenceRadius,
        cachedXCosineCoefficients_, cachedXSineCoefficients_,
        maximumDegree_, maximumOrder_ );

    const auto& yPolyCoefficients = yPolyDataset_->getPolyCoefficients( fileIndex );
    const auto& yShIndices = yPolyDataset_->getSHDegreeAndOrderIndices( fileIndex );
    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius, solarLongitude,
        yPolyCoefficients, yShIndices,
        fileMeta.powersInvRadius, fileMeta.referenceRadius,
        cachedYCosineCoefficients_, cachedYSineCoefficients_,
        maximumDegree_, maximumOrder_ );

    const auto& zPolyCoefficients = zPolyDataset_->getPolyCoefficients( fileIndex );
    const auto& zShIndices = zPolyDataset_->getSHDegreeAndOrderIndices( fileIndex );
    simulation_setup::StokesCoefficientsEvaluator::evaluate2D(
        radius, solarLongitude,
        zPolyCoefficients, zShIndices,
        fileMeta.powersInvRadius, fileMeta.referenceRadius,
        cachedZCosineCoefficients_, cachedZSineCoefficients_,
        maximumDegree_, maximumOrder_ );

    // Compute wind components using spherical harmonics (shared calculator benefits from cache)
    auto* shCalc = getActiveSphericalHarmonicsCalculator();

    const double windX = shCalc->calculateSurfaceSphericalHarmonics(
        cachedXSineCoefficients_, cachedXCosineCoefficients_,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cachedXCosineCoefficients_.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cachedXCosineCoefficients_.cols() - 1 );

    const double windY = shCalc->calculateSurfaceSphericalHarmonics(
        cachedYSineCoefficients_, cachedYCosineCoefficients_,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cachedYCosineCoefficients_.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cachedYCosineCoefficients_.cols() - 1 );

    const double windZ = shCalc->calculateSurfaceSphericalHarmonics(
        cachedZSineCoefficients_, cachedZCosineCoefficients_,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cachedZCosineCoefficients_.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cachedZCosineCoefficients_.cols() - 1 );

    return Eigen::Vector3d( windX, windY, windZ );
}

/*!
 * \brief Find the time interval index for polynomial datasets with caching.
 * Searches through all time periods defined in the dataset files to find
 * which interval contains the specified time. Uses cached file index as
 * optimization hint to check last used file first.
 * \param time The time for which to find the corresponding time interval [s]
 * \param dataset Polynomial dataset to search through
 * \return Index of the time interval containing the given time
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaWindModel::findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaPolyDataset> dataset ) const
{
    // Determine which cached file index to use based on dataset pointer
    int& cachedFileIndex = (dataset == xPolyDataset_) ? lastXFileIndex_ :
                           (dataset == yPolyDataset_) ? lastYFileIndex_ : lastZFileIndex_;

    const std::size_t numFiles = dataset->getNumFiles();

    // Check last used file first (cache optimization)
    if ( cachedFileIndex >= 0 && cachedFileIndex < static_cast<int>( numFiles ) )
    {
        const auto& fileMeta = dataset->getFileMeta( cachedFileIndex );
        for ( const auto& period : fileMeta.timePeriods )
        {
            if ( time >= period.first && time <= period.second )
            {
                return cachedFileIndex;
            }
        }
    }

    // Cache miss - search all files
    for ( std::size_t fileIndex = 0; fileIndex < numFiles; ++fileIndex )
    {
        const auto& fileMeta = dataset->getFileMeta( fileIndex );
        for ( const auto& period : fileMeta.timePeriods )
        {
            if ( time >= period.first && time <= period.second )
            {
                cachedFileIndex = static_cast<int>( fileIndex );
                return cachedFileIndex;
            }
        }
    }
    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}

/*!
 * \brief Find the time interval index for Stokes datasets with caching and binary search.
 * Searches through all files in the dataset to find which time interval
 * contains the specified time. Uses cached file index as optimization hint
 * and binary search for efficient lookup.
 * \param time The time for which to find the corresponding time interval [s]
 * \param dataset Stokes dataset to search through
 * \return Index of the time interval containing the given time
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaWindModel::findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaStokesDataset> dataset ) const
{
    // Determine which cached file index to use based on dataset pointer
    int& cachedFileIndex = (dataset == xStokesDataset_) ? lastXFileIndex_ :
                           (dataset == yStokesDataset_) ? lastYFileIndex_ : lastZFileIndex_;

    const std::size_t numFiles = dataset->nFiles();

    // Check last used file first (temporal locality optimization)
    if ( cachedFileIndex >= 0 && cachedFileIndex < static_cast<int>( numFiles ) )
    {
        const auto& fileMeta = dataset->files()[cachedFileIndex];
        if ( time >= fileMeta.start_epoch && time <= fileMeta.end_epoch )
        {
            return cachedFileIndex;
        }
    }

    // Binary search for time interval (assumes files are ordered by time)
    int left = 0;
    int right = static_cast<int>( numFiles ) - 1;

    while ( left <= right )
    {
        const int mid = left + ( right - left ) / 2;
        const auto& fileMeta = dataset->files()[mid];

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
            cachedFileIndex = mid;
            return mid;
        }
    }

    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}

/*!
 * \brief Compute a single wind component from polynomial coefficients with caching.
 * Evaluates polynomial coefficients at the given radius and solar longitude,
 * then uses spherical harmonics to compute the wind component at the specified
 * surface location. Uses cached coefficient matrices to avoid repeated allocations.
 * \param dataset Polynomial dataset containing the coefficients
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute the wind component [s]
 * \return Wind velocity component [m/s]
 * \throws std::runtime_error If dataset is null or time is out of range
 */
double ComaWindModel::computeWindComponentFromPolyCoefficients(
    std::shared_ptr<simulation_setup::ComaPolyDataset> dataset,
    double radius, double longitude, double latitude, double time ) const
{
    if ( !dataset )
    {
        throw std::runtime_error( "ComaWindModel: dataset is null" );
    }

    // Determine which cached coefficient matrices to use based on dataset pointer
    Eigen::MatrixXd& cosineCoefficients = (dataset == xPolyDataset_) ? cachedXCosineCoefficients_ :
                                          (dataset == yPolyDataset_) ? cachedYCosineCoefficients_ : cachedZCosineCoefficients_;
    Eigen::MatrixXd& sineCoefficients = (dataset == xPolyDataset_) ? cachedXSineCoefficients_ :
                                        (dataset == yPolyDataset_) ? cachedYSineCoefficients_ : cachedZSineCoefficients_;

    const int fileIndex = findTimeIntervalIndex( time, dataset );
    const auto& fileMeta = dataset->getFileMeta( fileIndex );
    const auto& polyCoefficients = dataset->getPolyCoefficients( fileIndex );
    const auto& shIndices = dataset->getSHDegreeAndOrderIndices( fileIndex );

    const double solarLongitude = calculateSolarLongitude( time );

    // Only resize/zero coefficient matrices if dimensions changed
    const int maxDegreeAvailable = fileMeta.maxDegreeSH;
    const int effectiveMaxDegree = ( maximumDegree_ > 0 ) ? maximumDegree_ : maxDegreeAvailable;
    const int effectiveMaxOrder = ( maximumOrder_ > 0 ) ? maximumOrder_ : maxDegreeAvailable;

    // Check sizing flag first to avoid repeated size comparisons
    if ( !coefficientMatricesSized_ ||
         cosineCoefficients.rows() != effectiveMaxDegree + 1 ||
         cosineCoefficients.cols() != effectiveMaxOrder + 1 )
    {
        cosineCoefficients = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        sineCoefficients = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        coefficientMatricesSized_ = true;
    }
    else
    {
        // Same size - just zero the contents (faster than reallocation)
        cosineCoefficients.setZero();
        sineCoefficients.setZero();
    }

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

    return getActiveSphericalHarmonicsCalculator()->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cosineCoefficients.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cosineCoefficients.cols() - 1 );
}

/*!
 * \brief Compute a single wind component from Stokes coefficients with advanced optimizations.
 * Uses cached coefficient matrices, batch interpolation, and reduced interpolators for
 * efficient evaluation. Includes comprehensive caching to avoid redundant computations.
 * \param dataset Stokes dataset containing the coefficients
 * \param interpolators Pre-initialized 2D interpolators for efficient coefficient evaluation
 * \param radius Radial distance from comet center [m]
 * \param longitude Longitude in comet body-fixed frame [rad]
 * \param latitude Latitude in comet body-fixed frame [rad]
 * \param time Time at which to compute the wind component [s]
 * \return Wind velocity component [m/s]
 * \throws std::runtime_error If dataset is null or time is out of range
 */
double ComaWindModel::computeWindComponentFromStokesCoefficients(
    std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
    const std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                                 std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>& interpolators,
    double radius, double longitude, double latitude, double time ) const
{
    if ( !dataset )
    {
        throw std::runtime_error( "ComaWindModel: dataset is null" );
    }

    // Step 1: Get time-dependent properties
    const int fileIndex = findTimeIntervalIndex( time, dataset );
    const double solarLongitude = calculateSolarLongitude( time );

    // Step 2: Get dataset properties
    const int nmax = dataset->nmax();
    const double referenceRadius = dataset->getReferenceRadius(fileIndex);

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Clamp radius to reference radius - wind velocity is constant beyond reference radius
    const double effectiveRadius = std::min(radius, referenceRadius);

    // Determine which cached coefficient matrices and interpolators to use based on dataset pointer
    Eigen::MatrixXd& cosineCoefficients = (dataset == xStokesDataset_) ? cachedXCosineCoefficients_ :
                                          (dataset == yStokesDataset_) ? cachedYCosineCoefficients_ : cachedZCosineCoefficients_;
    Eigen::MatrixXd& sineCoefficients = (dataset == xStokesDataset_) ? cachedXSineCoefficients_ :
                                        (dataset == yStokesDataset_) ? cachedYSineCoefficients_ : cachedZSineCoefficients_;

    const auto& reducedInterpolators = (dataset == xStokesDataset_) ? xReducedStokesInterpolators_ :
                                       (dataset == yStokesDataset_) ? yReducedStokesInterpolators_ : zReducedStokesInterpolators_;

    auto& fallbackInterpolators = (dataset == xStokesDataset_) ? xFallbackStokesInterpolators_ :
                                  (dataset == yStokesDataset_) ? yFallbackStokesInterpolators_ : zFallbackStokesInterpolators_;

    auto& fallbackReducedInterpolators = (dataset == xStokesDataset_) ? xFallbackReducedStokesInterpolators_ :
                                         (dataset == yStokesDataset_) ? yFallbackReducedStokesInterpolators_ : zFallbackReducedStokesInterpolators_;

    // Check if coefficient interpolation can be skipped (same radius/solar longitude)
    constexpr double radiusTolerance = 1e-10;
    constexpr double solarLongitudeTolerance = 1e-10;
    constexpr double radiusToleranceSq = radiusTolerance * radiusTolerance;
    constexpr double solarLongitudeToleranceSq = solarLongitudeTolerance * solarLongitudeTolerance;

    const double radiusDiff = effectiveRadius - cachedEffectiveRadius_;
    const double solarLongitudeDiff = solarLongitude - cachedInterpolationSolarLongitude_;

    const bool sameRadius = cacheFlags_.interpolationValid &&
                           radiusDiff * radiusDiff < radiusToleranceSq;
    const bool sameSolarLongitude = cacheFlags_.interpolationValid &&
                                   solarLongitudeDiff * solarLongitudeDiff < solarLongitudeToleranceSq;

    if ( !sameRadius || !sameSolarLongitude )
    {
        // Cache miss - need to recompute coefficients
        // Step 3: Only zero coefficient matrices if dimensions match
        if ( cosineCoefficients.rows() == effectiveMaxDegree + 1 &&
             cosineCoefficients.cols() == effectiveMaxOrder + 1 )
        {
            cosineCoefficients.setZero();
            sineCoefficients.setZero();
        }
        else
        {
            cosineCoefficients = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
            sineCoefficients = Eigen::MatrixXd::Zero( effectiveMaxDegree + 1, effectiveMaxOrder + 1 );
        }

        // Step 4: Interpolate coefficients using 2D interpolation (radius is clamped to referenceRadius)
        {
            interpolationPoint2D_[0] = effectiveRadius;
            interpolationPoint2D_[1] = solarLongitude;

            // Prepare interpolation state once for all coefficients (batch mode)
            if ( !interpolators.empty() )
            {
                // Use the first available interpolator to prepare the shared interpolation state
                auto firstInterpolator = interpolators.begin()->second.first.get();
                auto interpolationState = firstInterpolator->prepareInterpolationState(interpolationPoint2D_);

                // Now use the pre-computed state for all coefficient interpolations
                for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
                {
                    for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
                    {
                        std::pair<int,int> degreeOrderPair = {degree, order};
                        auto it = interpolators.find(degreeOrderPair);

                        if ( it != interpolators.end() )
                        {
                            // Use batch interpolation with pre-computed state (much faster!)
                            cosineCoefficients(degree, order) = it->second.first->interpolateWithState(interpolationState);
                            if ( order > 0 )  // Sine coefficients only exist for m > 0
                            {
                                sineCoefficients(degree, order) = it->second.second->interpolateWithState(interpolationState);
                            }
                        }
                        else
                        {
                            // Fallback: create 2D interpolator on-the-fly using helper method
                            createFallback2DInterpolator( dataset, fallbackInterpolators, fileIndex, degree, order,
                                                          cosineCoefficients(degree, order),
                                                          sineCoefficients(degree, order),
                                                          effectiveRadius, solarLongitude );
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
                        createFallback2DInterpolator( dataset, fallbackInterpolators, fileIndex, degree, order,
                                                      cosineCoefficients(degree, order),
                                                      sineCoefficients(degree, order),
                                                      effectiveRadius, solarLongitude );
                    }
                }
            }
        }

        // Update interpolation cache after successful coefficient computation
        cachedEffectiveRadius_ = effectiveRadius;
        cachedInterpolationSolarLongitude_ = solarLongitude;
        cacheFlags_.interpolationValid = true;
    }
    // else: Cache hit - reuse previously computed coefficient matrices

    // Step 5: Compute wind component using spherical harmonics expansion
    // Note: This step always runs (can't cache because lat/lon may change)
    return getActiveSphericalHarmonicsCalculator()->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );
}

/*!
 * \brief Calculate the solar longitude in the comet body-fixed frame with caching.
 * Computes the angle between the Sun direction (projected onto the XY plane)
 * and the X-axis in the comet body-fixed coordinate system. This angle is
 * used to determine the solar heating pattern on the comet surface.
 * Results are cached to avoid redundant state function calls when time hasn't changed.
 * \param time Time at which to compute solar longitude [s]
 * \return Solar longitude angle from X-axis in XY plane [rad]
 * \throws std::runtime_error If state functions are not initialized
 */
double ComaWindModel::calculateSolarLongitude( const double time ) const
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
        throw std::runtime_error( "ComaWindModel: State functions must be initialized" );
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

/*!
 * \brief Get the active spherical harmonics calculator.
 * Returns the shared calculator from ComaModel if available, otherwise
 * returns the owned calculator. This allows efficient sharing of computation
 * resources when both wind and density models use spherical harmonics.
 * \return Pointer to the active spherical harmonics calculator (non-owning)
 */
SphericalHarmonicsCalculator* ComaWindModel::getActiveSphericalHarmonicsCalculator() const
{
    // Return shared calculator if available, otherwise use our own
    return sharedSphericalHarmonicsCalculator_ != nullptr ? sharedSphericalHarmonicsCalculator_ : sphericalHarmonicsCalculator_.get();
}

/*!
 * \brief Initialize interpolators for Stokes coefficients.
 * Pre-computes multilinear interpolators for all spherical harmonic coefficients
 * up to the specified maximum degree and order. This optimization significantly
 * improves performance during repeated wind velocity evaluations.
 */
void ComaWindModel::initializeStokesInterpolators()
{
    // This function is only called from the Stokes coefficients constructor,
    // so dataType_ is guaranteed to be STOKES_COEFFICIENTS and all Stokes datasets are guaranteed to be non-null

    // Get dataset properties (assume same for all datasets)
    const auto& radiiGrid = xStokesDataset_->radii();
    const auto& longitudeGrid = xStokesDataset_->lons();
    const int nmax = xStokesDataset_->nmax();

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Set up interpolation grids (shared for all interpolators)
    std::vector<std::vector<double>> independentGrids(2);
    independentGrids[0] = radiiGrid;     // Radius grid
    independentGrids[1] = longitudeGrid; // Solar longitude grid

    // Initialize interpolators for each (n,m) pair and each component (x, y, z)
    for ( int n = 0; n <= effectiveMaxDegree; ++n )
    {
        for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
        {
            // Create 2D grids for this coefficient (across all files for now - using first file)
            const std::size_t nRadii = radiiGrid.size();
            const std::size_t nLons = longitudeGrid.size();
            const int fileIndex = 0; // For now, use first file - this could be extended

            // X-component interpolators
            {
                boost::multi_array<double, 2> cosineGrid(boost::extents[nRadii][nLons]);
                boost::multi_array<double, 2> sineGrid(boost::extents[nRadii][nLons]);

                // Fill grids with coefficient values from x-component dataset
                for ( std::size_t r = 0; r < nRadii; ++r )
                {
                    for ( std::size_t l = 0; l < nLons; ++l )
                    {
                        auto coeffs = xStokesDataset_->getCoeff(fileIndex, r, l, n, m);
                        cosineGrid[r][l] = coeffs.first;  // Cosine coefficient
                        sineGrid[r][l] = coeffs.second;   // Sine coefficient
                    }
                }

                // Create interpolators for x-component
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

                // Store interpolators in x-component map
                std::pair<int,int> nmPair = {n, m};
                xStokesInterpolators_[nmPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
            }

            // Y-component interpolators
            {
                boost::multi_array<double, 2> cosineGrid(boost::extents[nRadii][nLons]);
                boost::multi_array<double, 2> sineGrid(boost::extents[nRadii][nLons]);

                // Fill grids with coefficient values from y-component dataset
                for ( std::size_t r = 0; r < nRadii; ++r )
                {
                    for ( std::size_t l = 0; l < nLons; ++l )
                    {
                        auto coeffs = yStokesDataset_->getCoeff(fileIndex, r, l, n, m);
                        cosineGrid[r][l] = coeffs.first;  // Cosine coefficient
                        sineGrid[r][l] = coeffs.second;   // Sine coefficient
                    }
                }

                // Create interpolators for y-component
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

                // Store interpolators in y-component map
                std::pair<int,int> nmPair = {n, m};
                yStokesInterpolators_[nmPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
            }

            // Z-component interpolators
            {
                boost::multi_array<double, 2> cosineGrid(boost::extents[nRadii][nLons]);
                boost::multi_array<double, 2> sineGrid(boost::extents[nRadii][nLons]);

                // Fill grids with coefficient values from z-component dataset
                for ( std::size_t r = 0; r < nRadii; ++r )
                {
                    for ( std::size_t l = 0; l < nLons; ++l )
                    {
                        auto coeffs = zStokesDataset_->getCoeff(fileIndex, r, l, n, m);
                        cosineGrid[r][l] = coeffs.first;  // Cosine coefficient
                        sineGrid[r][l] = coeffs.second;   // Sine coefficient
                    }
                }

                // Create interpolators for z-component
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

                // Store interpolators in z-component map
                std::pair<int,int> nmPair = {n, m};
                zStokesInterpolators_[nmPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
            }
        }
    }

    // Initialize reduced interpolators for coefficients beyond reference radius
    // These are 1D interpolators (solar longitude only)
    const std::size_t nFiles = xStokesDataset_->nFiles();

    // Resize reduced interpolator vectors to accommodate all files
    xReducedStokesInterpolators_.resize( nFiles );
    yReducedStokesInterpolators_.resize( nFiles );
    zReducedStokesInterpolators_.resize( nFiles );

    // Resize fallback cache vectors to accommodate all files
    xFallbackStokesInterpolators_.resize( nFiles );
    yFallbackStokesInterpolators_.resize( nFiles );
    zFallbackStokesInterpolators_.resize( nFiles );
    xFallbackReducedStokesInterpolators_.resize( nFiles );
    yFallbackReducedStokesInterpolators_.resize( nFiles );
    zFallbackReducedStokesInterpolators_.resize( nFiles );

    std::vector<std::vector<double>> reducedIndependentGrids(1);
    reducedIndependentGrids[0] = longitudeGrid; // Solar longitude grid only

    for ( std::size_t fileIndex = 0; fileIndex < nFiles; ++fileIndex )
    {
        for ( int degree = 0; degree <= effectiveMaxDegree; ++degree )
        {
            for ( int order = 0; order <= std::min(degree, effectiveMaxOrder); ++order )
            {
                const std::size_t numLongitudes = longitudeGrid.size();

                // X-component reduced interpolators
                {
                    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
                    boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

                    // Fill grids with reduced coefficient values for x-component
                    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
                    {
                        auto reducedCoeffs = xStokesDataset_->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
                        reducedCosineGrid[longitudeIndex] = reducedCoeffs.first;  // Cosine coefficient
                        reducedSineGrid[longitudeIndex] = reducedCoeffs.second;   // Sine coefficient
                    }

                    // Create and store reduced interpolators for x-component
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
                    xReducedStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
                }

                // Y-component reduced interpolators
                {
                    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
                    boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

                    // Fill grids with reduced coefficient values for y-component
                    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
                    {
                        auto reducedCoeffs = yStokesDataset_->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
                        reducedCosineGrid[longitudeIndex] = reducedCoeffs.first;  // Cosine coefficient
                        reducedSineGrid[longitudeIndex] = reducedCoeffs.second;   // Sine coefficient
                    }

                    // Create and store reduced interpolators for y-component
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
                    yReducedStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
                }

                // Z-component reduced interpolators
                {
                    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
                    boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

                    // Fill grids with reduced coefficient values for z-component
                    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
                    {
                        auto reducedCoeffs = zStokesDataset_->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
                        reducedCosineGrid[longitudeIndex] = reducedCoeffs.first;  // Cosine coefficient
                        reducedSineGrid[longitudeIndex] = reducedCoeffs.second;   // Sine coefficient
                    }

                    // Create and store reduced interpolators for z-component
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
                    zReducedStokesInterpolators_[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
                }
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
void ComaWindModel::createFallback2DInterpolator(
    std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>>& fallbackCache,
    const int fileIndex, const int degree, const int order,
    double& cosineCoeff, double& sineCoeff,
    const double radius, const double solarLongitude ) const
{
    std::pair<int,int> degreeOrderPair = {degree, order};

    // Check if we already created this interpolator
    auto it = fallbackCache[fileIndex].find(degreeOrderPair);
    if ( it != fallbackCache[fileIndex].end() )
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
    const auto& radiiGrid = dataset->radii();
    const auto& longitudeGrid = dataset->lons();
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
            auto coeffs = dataset->getCoeff(fileIndex, radiusIndex, longitudeIndex, degree, order);
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
    fallbackCache[fileIndex][degreeOrderPair] = {std::move(cosineInterpolator), std::move(sineInterpolator)};
}

/*!
 * \brief Helper to create 1D reduced interpolator for Stokes coefficients on-the-fly.
 * This is a fallback method used when pre-initialized reduced interpolators are not available.
 * Created interpolators are cached for reuse to avoid redundant construction.
 */
void ComaWindModel::createFallback1DInterpolator(
    std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>>& fallbackCache,
    const int fileIndex, const int degree, const int order,
    double& cosineCoeff, double& sineCoeff,
    const double solarLongitude ) const
{
    std::pair<int,int> degreeOrderPair = {degree, order};

    // Check if we already created this interpolator
    auto it = fallbackCache[fileIndex].find(degreeOrderPair);
    if ( it != fallbackCache[fileIndex].end() )
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
    const auto& longitudeGrid = dataset->lons();
    const std::size_t numLongitudes = longitudeGrid.size();

    std::vector<std::vector<double>> reducedIndependentGrids(1);
    reducedIndependentGrids[0] = longitudeGrid;

    boost::multi_array<double, 1> reducedCosineGrid(boost::extents[numLongitudes]);
    boost::multi_array<double, 1> reducedSineGrid(boost::extents[numLongitudes]);

    for ( std::size_t longitudeIndex = 0; longitudeIndex < numLongitudes; ++longitudeIndex )
    {
        auto reducedCoeffs = dataset->getReducedCoeff(fileIndex, longitudeIndex, degree, order);
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
    fallbackCache[fileIndex][degreeOrderPair] = {std::move(reducedCosineInterpolator), std::move(reducedSineInterpolator)};
}

} // namespace aerodynamics
} // namespace tudat