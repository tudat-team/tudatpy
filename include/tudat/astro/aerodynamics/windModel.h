/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_WIND_MODEL_H
#define TUDAT_WIND_MODEL_H

#include <Eigen/Core>

#include <memory>
#include <functional>
#include <map>

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include "tudat/math/interpolators/multiLinearInterpolator.h"

namespace tudat
{

namespace simulation_setup
{
    class ComaPolyDataset;
    class ComaStokesDataset;
}


namespace aerodynamics
{

// Forward declarations
class ComaModel;
class SphericalHarmonicsCalculator;

//! Base class for a wind model.
/*!
 *  Base class for a wind model.  The wind vector is defined as the local velocity of the atmosphere expressed in the
 *  frame corotating with the body.
 */
class WindModel
{
public:
    //! Constructor.
    WindModel( const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
               const bool includeCorotation = true ):
        associatedFrame_( associatedFrame ), includeCorotation_( includeCorotation )
    {
        if( !( associatedFrame == reference_frames::inertial_frame || associatedFrame == reference_frames::corotating_frame ||
               associatedFrame == reference_frames::vertical_frame ) )
        {
            throw std::runtime_error( "Error when creating wind model, definition must be in inertial, corotating or vertical frame" );
        }
    }

    //! Destructor.
    virtual ~WindModel( ) { }

    //! Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    virtual Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                                      const double currentLongitude,
                                                                      const double currentLatitude,
                                                                      const double currentTime ) = 0;

    /*!
     * \brief Get the reference frame associated with this wind model.
     * \return Reference frame in which the wind model is defined
     */
    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

    /*!
     * \brief Get whether atmospheric co-rotation is included.
     * \return Boolean indicating if atmospheric co-rotation is included in aerodynamic computations
     */
    bool getIncludeCorotation( ) const
    {
        return includeCorotation_;
    }

protected:
    reference_frames::AerodynamicsReferenceFrames associatedFrame_;

    //! Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations
    bool includeCorotation_;
};

/*!
 * \brief Simple wind model with constant wind velocity.
 */
class ConstantWindModel : public WindModel
{
public:
    /*!
     * \brief Constructor for constant wind model.
     * \param constantWindVelocity Constant wind velocity vector [m/s]
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    ConstantWindModel( const Eigen::Vector3d constantWindVelocity,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                       const bool includeCorotation = true ):
        WindModel( associatedFrame, includeCorotation ), constantWindVelocity_( constantWindVelocity )
    { }

    /*!
     * \brief Get the constant wind velocity vector.
     * \param currentAltitude Altitude (unused for constant wind)
     * \param currentLongitude Longitude (unused for constant wind)
     * \param currentLatitude Latitude (unused for constant wind)
     * \param currentTime Time (unused for constant wind)
     * \return Constant wind velocity vector [m/s]
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime )
    {
        return constantWindVelocity_;
    }

private:
    Eigen::Vector3d constantWindVelocity_;
};

//! Class for computing the wind velocity vector from a custom, user-defined function.
class CustomWindModel : public WindModel
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    CustomWindModel( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                     const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                     const bool includeCorotation = true ):
        WindModel( associatedFrame, includeCorotation ), windFunction_( windFunction )
    { }

    //! Destructor
    ~CustomWindModel( ) { }

    //! Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime )
    {
        return windFunction_( currentAltitude, currentLongitude, currentLatitude, currentTime );
    }

private:
    //! Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

//! Class for empty wind model (no wind, only controls co-rotation behavior)
/*!
 * Class for empty wind model that always returns zero wind velocity.
 * This is used when no physical wind model is required but co-rotation behavior needs to be specified.
 */
class EmptyWindModel : public WindModel
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    EmptyWindModel( const bool includeCorotation = true ):
        WindModel( reference_frames::vertical_frame, includeCorotation )
    { }

    //! Destructor
    ~EmptyWindModel( ) { }

    //! Function to retrieve wind velocity vector (always zero)
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere.
     * Always returns zero vector as this is an empty wind model.
     * \param currentAltitude Altitude at which wind vector is to be retrieved (unused).
     * \param currentLongitude Longitude at which wind vector is to be retrieved (unused).
     * \param currentLatitude Latitude at which wind vector is to be retrieved (unused).
     * \param currentTime Time at which wind vector is to be retrieved (unused).
     * \return Zero wind velocity vector in body-fixed, body-centered frame
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime )
    {
        return Eigen::Vector3d::Zero( );
    }
};

//! Class for computing the wind velocity vector from coma models
/*!
 * Class for computing the wind velocity vector from coma models. This class uses three separate datasets
 * for the wind velocity components in the modified vertical frame, each of which can be either polynomial or Stokes coefficients.
 *
 * The wind velocity components are defined in a modified vertical frame:
 *   - X-component: Meridional direction (North, in meridian plane)
 *   - Y-component: Zonal direction (West, completing the right-handed frame)
 *   - Z-component: Radial direction pointing OUTWARD from the nucleus (opposite to standard vertical frame)
 */
class ComaWindModel : public WindModel
{
public:
    /*!
     * \brief Constructor for polynomial coefficient datasets.
     * \param xPolyDataset Dataset for X-component wind (meridional/North, polynomial coefficients)
     * \param yPolyDataset Dataset for Y-component wind (zonal/West, polynomial coefficients)
     * \param zPolyDataset Dataset for Z-component wind (radial outward, polynomial coefficients)
     * \param comaModel Shared pointer to the ComaModel for accessing state functions
     * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     * \param cometRotationFunction Function returning comet body-fixed rotation matrix
     * \param maximumDegree Maximum degree used for computation (-1 for auto)
     * \param maximumOrder Maximum order used for computation (-1 for auto)
     * \param associatedFrame Reference frame for the wind model
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    ComaWindModel( const simulation_setup::ComaPolyDataset& xPolyDataset,
                   const simulation_setup::ComaPolyDataset& yPolyDataset,
                   const simulation_setup::ComaPolyDataset& zPolyDataset,
                   const std::shared_ptr<ComaModel>& comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree = -1,
                   const int& maximumOrder = -1,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                   const bool includeCorotation = true );

    /*!
     * \brief Constructor for Stokes coefficient datasets.
     * \param xStokesDataset Dataset for X-component wind (meridional/North, Stokes coefficients)
     * \param yStokesDataset Dataset for Y-component wind (zonal/West, Stokes coefficients)
     * \param zStokesDataset Dataset for Z-component wind (radial outward, Stokes coefficients)
     * \param comaModel Shared pointer to the ComaModel for accessing state functions
     * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     * \param cometRotationFunction Function returning comet body-fixed rotation matrix
     * \param maximumDegree Maximum degree used for computation (-1 for auto)
     * \param maximumOrder Maximum order used for computation (-1 for auto)
     * \param associatedFrame Reference frame for the wind model
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    ComaWindModel( const simulation_setup::ComaStokesDataset& xStokesDataset,
                   const simulation_setup::ComaStokesDataset& yStokesDataset,
                   const simulation_setup::ComaStokesDataset& zStokesDataset,
                   std::shared_ptr<ComaModel> comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree = -1,
                   const int& maximumOrder = -1,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                   const bool includeCorotation = true );

    //! Destructor
    ~ComaWindModel( ) = default;

    //! Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime ) override;

private:
    //! Type of data used (0: polynomial coefficients, 1: Stokes coefficients)
    int dataType_;

    //! Maximum spherical harmonic degree used for computation (-1 for auto-detect)
    int maximumDegree_;

    //! Maximum spherical harmonic order used for computation (-1 for auto-detect)
    int maximumOrder_;

    //! Polynomial coefficient dataset for X-component wind (meridional/North, used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> xPolyDataset_;
    //! Polynomial coefficient dataset for Y-component wind (zonal/West, used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> yPolyDataset_;
    //! Polynomial coefficient dataset for Z-component wind (radial outward, used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> zPolyDataset_;

    //! Stokes coefficient dataset for X-component wind (meridional/North, used when dataType_ == 1)
    std::shared_ptr<simulation_setup::ComaStokesDataset> xStokesDataset_;
    //! Stokes coefficient dataset for Y-component wind (zonal/West, used when dataType_ == 1)
    std::shared_ptr<simulation_setup::ComaStokesDataset> yStokesDataset_;
    //! Stokes coefficient dataset for Z-component wind (radial outward, used when dataType_ == 1)
    std::shared_ptr<simulation_setup::ComaStokesDataset> zStokesDataset_;

    //! Reference to the ComaModel for accessing shared spherical harmonics calculator
    std::shared_ptr<ComaModel> comaModel_;

    //! Function returning Sun state vector (position [m], velocity [m/s])
    std::function<Eigen::Vector6d()> sunStateFunction_;

    //! Function returning Comet state vector (position [m], velocity [m/s])
    std::function<Eigen::Vector6d()> cometStateFunction_;

    //! Function returning comet body-fixed to inertial rotation matrix
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

    //! Owned spherical harmonics calculator (used when not sharing with ComaModel)
    std::unique_ptr<SphericalHarmonicsCalculator> sphericalHarmonicsCalculator_;

    //! Non-owning pointer to shared spherical harmonics calculator (used when sharing with ComaModel)
    SphericalHarmonicsCalculator* sharedSphericalHarmonicsCalculator_;

    //! Pre-initialized interpolators for X-component (meridional/North) Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> xStokesInterpolators_;
    //! Pre-initialized interpolators for Y-component (zonal/West) Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> yStokesInterpolators_;
    //! Pre-initialized interpolators for Z-component (radial outward) Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> zStokesInterpolators_;

    //! Pre-initialized reduced interpolators for X-component Stokes coefficients (1D: solar longitude only, for radius > reference radius)
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> xReducedStokesInterpolators_;
    //! Pre-initialized reduced interpolators for Y-component Stokes coefficients (1D: solar longitude only, for radius > reference radius)
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> yReducedStokesInterpolators_;
    //! Pre-initialized reduced interpolators for Z-component Stokes coefficients (1D: solar longitude only, for radius > reference radius)
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> zReducedStokesInterpolators_;

    //! Cache for fallback X-component interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>> xFallbackStokesInterpolators_;
    //! Cache for fallback Y-component interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>> yFallbackStokesInterpolators_;
    //! Cache for fallback Z-component interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>> zFallbackStokesInterpolators_;

    //! Cache for fallback X-component reduced interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> xFallbackReducedStokesInterpolators_;
    //! Cache for fallback Y-component reduced interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> yFallbackReducedStokesInterpolators_;
    //! Cache for fallback Z-component reduced interpolators (created on-demand, then cached for reuse)
    mutable std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> zFallbackReducedStokesInterpolators_;

    // ========== Hot path: Frequently accessed cached values (grouped for cache locality) ==========

    //! Cached solar longitude value to avoid redundant state function calls [rad]
    mutable double cachedSolarLongitude_;

    //! Time at which solar longitude was last cached [s]
    mutable double cachedTime_;

    //! Pre-allocated coefficient matrices to avoid repeated allocations for x/y/z components
    mutable Eigen::MatrixXd cachedXCosineCoefficients_;
    mutable Eigen::MatrixXd cachedXSineCoefficients_;
    mutable Eigen::MatrixXd cachedYCosineCoefficients_;
    mutable Eigen::MatrixXd cachedYSineCoefficients_;
    mutable Eigen::MatrixXd cachedZCosineCoefficients_;
    mutable Eigen::MatrixXd cachedZSineCoefficients_;

    //! Cached file indices from last time interval search (optimization hint for each component dataset)
    mutable int lastXFileIndex_;
    mutable int lastYFileIndex_;
    mutable int lastZFileIndex_;

    //! Cached interpolation point (radius, solar longitude) to detect when cache is valid
    mutable double cachedRadius_;
    mutable double cachedEffectiveRadius_;  // Cached effective radius (clamped to ref radius) for interpolation cache
    mutable double cachedInterpolationSolarLongitude_;

    //! Flag to track if coefficient matrices have been sized (optimization to avoid repeated size checks)
    mutable bool coefficientMatricesSized_;

    //! Cached latitude/longitude
    mutable double cachedLatitude_;
    mutable double cachedLongitude_;

    //! Cached final wind vector result to avoid repeated computations
    mutable Eigen::Vector3d cachedFinalWindVector_;

    //! Cached state function results
    mutable Eigen::Vector6d cachedSunState_;
    mutable Eigen::Vector6d cachedCometState_;
    mutable Eigen::Matrix3d cachedRotationMatrix_;

    //! Pre-allocated interpolation point vectors to avoid repeated allocations
    mutable std::vector<double> interpolationPoint2D_;
    mutable std::vector<double> interpolationPoint1D_;

    //! Cache validity flags packed into a bitfield for memory efficiency
    struct CacheFlags {
        bool solarLongitudeValid : 1;
        bool interpolationValid : 1;
        bool windVectorValid : 1;
        bool stateValid : 1;

        CacheFlags() : solarLongitudeValid(false), interpolationValid(false),
                       windVectorValid(false), stateValid(false) {}
    };
    mutable CacheFlags cacheFlags_;

    // Helper methods
    /*!
     * \brief Find the index of the time interval containing a given time for polynomial datasets.
     * \param time The time for which to find the corresponding time interval
     * \param dataset Polynomial dataset to search through
     * \return Index of the time interval containing the given time
     * \throws std::runtime_error If no matching time interval is found
     */
    int findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaPolyDataset> dataset ) const;

    /*!
     * \brief Find the index of the time interval containing a given time for Stokes datasets.
     * \param time The time for which to find the corresponding time interval
     * \param dataset Stokes dataset to search through
     * \return Index of the time interval containing the given time
     * \throws std::runtime_error If no matching time interval is found
     */
    int findTimeIntervalIndex( double time, std::shared_ptr<simulation_setup::ComaStokesDataset> dataset ) const;

    /*!
     * \brief Compute complete wind vector from polynomial coefficients (vectorized, all 3 components).
     * This is the optimized method that computes all wind components (x, y, z) in a single pass,
     * sharing caches and reducing overhead. Preferred over calling single-component methods.
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute the wind vector [s]
     * \return Complete wind velocity vector [m/s]
     * \throws std::runtime_error If datasets are null or time is out of range
     */
    Eigen::Vector3d computeWindVectorFromPolyCoefficients(
        double radius, double longitude, double latitude, double time ) const;

    /*!
     * \brief Compute complete wind vector from Stokes coefficients (vectorized, all 3 components).
     * This is the optimized method that computes all wind components (x, y, z) in a single pass,
     * using batch interpolation and sharing caches for maximum performance. Preferred over calling
     * single-component methods.
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute the wind vector [s]
     * \return Complete wind velocity vector [m/s]
     * \throws std::runtime_error If datasets are null or time is out of range
     */
    Eigen::Vector3d computeWindVectorFromStokesCoefficients(
        double radius, double longitude, double latitude, double time ) const;

    /*!
     * \brief Compute wind component from polynomial coefficients (single component - internal use).
     * \param dataset Polynomial dataset containing the coefficients
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute the wind component [s]
     * \return Wind velocity component [m/s]
     * \throws std::runtime_error If dataset is null or time is out of range
     */
    double computeWindComponentFromPolyCoefficients(
        std::shared_ptr<simulation_setup::ComaPolyDataset> dataset,
        double radius, double longitude, double latitude, double time ) const;

    /*!
     * \brief Compute wind component from Stokes coefficients using pre-initialized interpolators (single component - internal use).
     * \param dataset Stokes dataset containing the coefficients
     * \param interpolators Pre-initialized interpolators for efficient coefficient evaluation
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute the wind component [s]
     * \return Wind velocity component [m/s]
     * \throws std::runtime_error If dataset is null or time is out of range
     */
    double computeWindComponentFromStokesCoefficients(
        std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
        const std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                                     std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>& interpolators,
        double radius, double longitude, double latitude, double time ) const;

    /*!
     * \brief Calculate solar longitude in comet body-fixed frame with caching.
     * \param time Time at which to compute solar longitude [s]
     * \return Solar longitude angle from X-axis in XY plane [rad]
     * \throws std::runtime_error If state functions are not initialized
     */
    double calculateSolarLongitude( double time ) const;

    /*!
     * \brief Initialize interpolators for Stokes coefficients.
     * Called only for Stokes coefficient datasets to pre-compute interpolators
     * for efficient evaluation during wind velocity calculations.
     */
    void initializeStokesInterpolators();

    /*!
     * \brief Get the active spherical harmonics calculator.
     * Returns the shared calculator if available, otherwise returns the owned calculator.
     * \return Pointer to the active spherical harmonics calculator
     */
    SphericalHarmonicsCalculator* getActiveSphericalHarmonicsCalculator() const;

    /*!
     * @brief Helper to create 2D interpolator for Stokes coefficients on-the-fly (fallback)
     * @param dataset Dataset to use for coefficient extraction
     * @param fallbackCache Fallback cache to store the created interpolator
     * @param fileIndex File index in dataset
     * @param n Spherical harmonic degree
     * @param m Spherical harmonic order
     * @param cosineCoeff Output: interpolated cosine coefficient
     * @param sineCoeff Output: interpolated sine coefficient
     * @param radius Radius for interpolation [m]
     * @param solarLongitude Solar longitude for interpolation [rad]
     */
    void createFallback2DInterpolator(
        std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
        std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>>& fallbackCache,
        int fileIndex, int n, int m,
        double& cosineCoeff, double& sineCoeff,
        double radius, double solarLongitude ) const;

    /*!
     * @brief Helper to create 1D reduced interpolator for Stokes coefficients on-the-fly (fallback)
     * @param dataset Dataset to use for coefficient extraction
     * @param fallbackCache Fallback cache to store the created interpolator
     * @param fileIndex File index in dataset
     * @param n Spherical harmonic degree
     * @param m Spherical harmonic order
     * @param cosineCoeff Output: interpolated cosine coefficient
     * @param sineCoeff Output: interpolated sine coefficient
     * @param solarLongitude Solar longitude for interpolation [rad]
     */
    void createFallback1DInterpolator(
        std::shared_ptr<simulation_setup::ComaStokesDataset> dataset,
        std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>>& fallbackCache,
        int fileIndex, int n, int m,
        double& cosineCoeff, double& sineCoeff,
        double solarLongitude ) const;
};

}  // namespace aerodynamics

}  // namespace tudat

#endif  // TUDAT_WIND_MODEL_H
