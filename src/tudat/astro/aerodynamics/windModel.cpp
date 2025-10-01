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

ComaWindModel::ComaWindModel( const simulation_setup::ComaPolyDataset& xPolyDataset,
                   const simulation_setup::ComaPolyDataset& yPolyDataset,
                   const simulation_setup::ComaPolyDataset& zPolyDataset,
                   std::shared_ptr<ComaModel> comaModel,
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
        sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
    {
        if ( !comaModel_ )
        {
            throw std::invalid_argument( "ComaWindModel: ComaModel must be provided" );
        }

        if ( xPolyDataset_->getNumFiles() == 0 || yPolyDataset_->getNumFiles() == 0 || zPolyDataset_->getNumFiles() == 0 )
        {
            throw std::invalid_argument( "ComaWindModel: All PolyDatasets must contain at least one file" );
        }

        if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
        {
            throw std::invalid_argument( "ComaWindModel: Maximum degree and order must be >= -1" );
        }
    }

ComaWindModel::ComaWindModel( const simulation_setup::ComaStokesDataset& xStokesDataset,
                   const simulation_setup::ComaStokesDataset& yStokesDataset,
                   const simulation_setup::ComaStokesDataset& zStokesDataset,
                   std::shared_ptr<ComaModel> comaModel,
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
        sphericalHarmonicsCalculator_( std::make_unique<SphericalHarmonicsCalculator>() )
    {
        if ( !comaModel_ )
        {
            throw std::invalid_argument( "ComaWindModel: ComaModel must be provided" );
        }

        if ( xStokesDataset_->nFiles() == 0 || yStokesDataset_->nFiles() == 0 || zStokesDataset_->nFiles() == 0 )
        {
            throw std::invalid_argument( "ComaWindModel: All StokesDatasets must contain at least one file" );
        }

        if ( maximumDegree_ < -1 || maximumOrder_ < -1 )
        {
            throw std::invalid_argument( "ComaWindModel: Maximum degree and order must be >= -1" );
        }

        // Initialize interpolators for Stokes coefficients
        initializeStokesInterpolators();
    }

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