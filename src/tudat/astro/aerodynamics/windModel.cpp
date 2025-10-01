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

} // namespace aerodynamics
} // namespace tudat