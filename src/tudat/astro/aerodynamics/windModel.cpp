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
 */
ComaWindModel::ComaWindModel( const simulation_setup::ComaPolyDataset& xPolyDataset,
                   const simulation_setup::ComaPolyDataset& yPolyDataset,
                   const simulation_setup::ComaPolyDataset& zPolyDataset,
                   std::shared_ptr<ComaModel> comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree,
                   const int& maximumOrder,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame ) :
        WindModel( associatedFrame ),
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
        sharedSphericalHarmonicsCalculator_( nullptr )
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
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame ) :
        WindModel( associatedFrame ),
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
        sharedSphericalHarmonicsCalculator_( nullptr )
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

        // Initialize interpolators for Stokes coefficients
        initializeStokesInterpolators();
    }

/*!
 * \brief Compute the current wind velocity vector in body-fixed coordinates.
 * This method evaluates the wind velocity at the specified location and time
 * using either polynomial or Stokes coefficient data depending on the model type.
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

        // Compute wind velocity components
        double windX, windY, windZ;

        switch ( dataType_ )
        {
            case 0: // POLYNOMIAL_COEFFICIENTS
                windX = computeWindComponentFromPolyCoefficients( xPolyDataset_, radius, currentLongitude, currentLatitude, currentTime );
                windY = computeWindComponentFromPolyCoefficients( yPolyDataset_, radius, currentLongitude, currentLatitude, currentTime );
                windZ = computeWindComponentFromPolyCoefficients( zPolyDataset_, radius, currentLongitude, currentLatitude, currentTime );
                break;

            case 1: // STOKES_COEFFICIENTS
                windX = computeWindComponentFromStokesCoefficients( xStokesDataset_, xStokesInterpolators_, radius, currentLongitude, currentLatitude, currentTime );
                windY = computeWindComponentFromStokesCoefficients( yStokesDataset_, yStokesInterpolators_, radius, currentLongitude, currentLatitude, currentTime );
                windZ = computeWindComponentFromStokesCoefficients( zStokesDataset_, zStokesInterpolators_, radius, currentLongitude, currentLatitude, currentTime );
                break;

            default:
                throw std::runtime_error( "ComaWindModel: Unknown data type" );
        }

        return Eigen::Vector3d( windX, windY, windZ );
    }

/*!
 * \brief Find the time interval index for polynomial datasets.
 * Searches through all time periods defined in the dataset files to find
 * which interval contains the specified time.
 * \param time The time for which to find the corresponding time interval [s]
 * \param dataset Polynomial dataset to search through
 * \return Index of the time interval containing the given time
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaWindModel::findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaPolyDataset> dataset ) const
{
    for ( std::size_t f = 0; f < dataset->getNumFiles(); ++f )
    {
        const auto& fileMeta = dataset->getFileMeta( f );
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

/*!
 * \brief Find the time interval index for Stokes datasets.
 * Searches through all files in the dataset to find which time interval
 * contains the specified time.
 * \param time The time for which to find the corresponding time interval [s]
 * \param dataset Stokes dataset to search through
 * \return Index of the time interval containing the given time
 * \throws std::runtime_error If no matching time interval is found
 */
int ComaWindModel::findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaStokesDataset> dataset ) const
{
    for ( std::size_t f = 0; f < dataset->nFiles(); ++f )
    {
        const auto& fileMeta = dataset->files()[f];
        if ( time >= fileMeta.start_epoch && time <= fileMeta.end_epoch )
        {
            return static_cast<int>( f );
        }
    }
    throw std::runtime_error( "Time " + std::to_string( time ) + " does not fall into any defined time period" );
}

/*!
 * \brief Compute a single wind component from polynomial coefficients.
 * Evaluates polynomial coefficients at the given radius and solar longitude,
 * then uses spherical harmonics to compute the wind component at the specified
 * surface location.
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

    const int fileIndex = findTimeIntervalIndex( time, dataset );
    const auto& fileMeta = dataset->getFileMeta( fileIndex );
    const auto& polyCoefficients = dataset->getPolyCoefficients( fileIndex );
    const auto& shIndices = dataset->getSHDegreeAndOrderIndices( fileIndex );

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

    return getActiveSphericalHarmonicsCalculator()->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        maximumDegree_ > 0 ? maximumDegree_ : cosineCoefficients.rows() - 1,
        maximumOrder_ > 0 ? maximumOrder_ : cosineCoefficients.cols() - 1 );
}

/*!
 * \brief Compute a single wind component from Stokes coefficients.
 * Uses pre-initialized interpolators to efficiently evaluate Stokes coefficients
 * at the given radius and solar longitude, then applies spherical harmonics
 * to compute the wind component. Includes distance-dependent decay for radii
 * beyond the reference radius.
 * \param dataset Stokes dataset containing the coefficients
 * \param interpolators Pre-initialized interpolators for efficient coefficient evaluation
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

    // Step 1: Find time interval index
    const int fileIndex = findTimeIntervalIndex( time, dataset );

    // Step 2: Calculate solar longitude
    const double solarLongitude = calculateSolarLongitude();

    // Step 3: Get dataset properties
    const int nmax = dataset->nmax();
    const double referenceRadius = dataset->getReferenceRadius(fileIndex);
    const auto& radiiGrid = dataset->radii();

    // Determine effective maximum degree and order
    const int effectiveMaxDegree = maximumDegree_ > 0 ? maximumDegree_ : nmax;
    const int effectiveMaxOrder = maximumOrder_ > 0 ? maximumOrder_ : nmax;

    // Initialize coefficient matrices
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero(effectiveMaxDegree + 1, effectiveMaxOrder + 1);
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero(effectiveMaxDegree + 1, effectiveMaxOrder + 1);

    // Step 4: Handle radius beyond reference radius (apply 1/r² decay)
    bool applyDecay = false;
    double interpolationRadius = radius;

    if (radius > referenceRadius)
    {
        applyDecay = true;
        interpolationRadius = radiiGrid.back();
    }

    std::vector<double> interpolationPoint = {interpolationRadius, solarLongitude};

    // Step 5: For each spherical harmonic coefficient (n,m), use pre-initialized interpolators
    for ( int n = 0; n <= effectiveMaxDegree; ++n )
    {
        for ( int m = 0; m <= std::min(n, effectiveMaxOrder); ++m )
        {
            std::pair<int,int> nmPair = {n, m};

            // Find the pre-initialized interpolators for this (n,m) pair
            auto it = interpolators.find(nmPair);
            if ( it != interpolators.end() )
            {
                // Use pre-initialized interpolators
                cosineCoefficients(n, m) = it->second.first->interpolate(interpolationPoint);

                if ( m > 0 )  // Sine coefficients only exist for m > 0
                {
                    sineCoefficients(n, m) = it->second.second->interpolate(interpolationPoint);
                }
            }
            else
            {
                // Fallback: create interpolator on-the-fly (should not happen if initialization was correct)
                const auto& radiiGrid = dataset->radii();
                const auto& longitudeGrid = dataset->lons();
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
                        auto coeffs = dataset->getCoeff(fileIndex, r, l, n, m);
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

    // Step 6: Apply logarithmic decay term for 1/r² behavior if radius > reference radius
    if (applyDecay)
    {
        simulation_setup::StokesCoefficientsEvaluator::applyDecayTerm(cosineCoefficients, radius, referenceRadius);
    }

    // Step 7: Calculate wind component using spherical harmonics
    return getActiveSphericalHarmonicsCalculator()->calculateSurfaceSphericalHarmonics(
        sineCoefficients, cosineCoefficients,
        latitude, longitude,
        effectiveMaxDegree, effectiveMaxOrder
    );
}

/*!
 * \brief Calculate the solar longitude in the comet body-fixed frame.
 * Computes the angle between the Sun direction (projected onto the XY plane)
 * and the X-axis in the comet body-fixed coordinate system. This angle is
 * used to determine the solar heating pattern on the comet surface.
 * \return Solar longitude angle from X-axis in XY plane [rad]
 * \throws std::runtime_error If state functions are not initialized
 */
double ComaWindModel::calculateSolarLongitude() const
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
}

} // namespace aerodynamics
} // namespace tudat