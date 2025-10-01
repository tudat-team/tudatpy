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
    WindModel( const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        associatedFrame_( associatedFrame )
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

protected:
    reference_frames::AerodynamicsReferenceFrames associatedFrame_;
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
     */
    ConstantWindModel( const Eigen::Vector3d constantWindVelocity,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModel( associatedFrame ), constantWindVelocity_( constantWindVelocity )
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
     */
    CustomWindModel( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                     const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModel( associatedFrame ), windFunction_( windFunction )
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

//! Class for computing the wind velocity vector from coma models
/*!
 * Class for computing the wind velocity vector from coma models. This class uses three separate datasets
 * for the x, y, and z wind components, each of which can be either polynomial or Stokes coefficients.
 */
class ComaWindModel : public WindModel
{
public:
    /*!
     * \brief Constructor for polynomial coefficient datasets.
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
    ComaWindModel( const simulation_setup::ComaPolyDataset& xPolyDataset,
                   const simulation_setup::ComaPolyDataset& yPolyDataset,
                   const simulation_setup::ComaPolyDataset& zPolyDataset,
                   std::shared_ptr<ComaModel> comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree = -1,
                   const int& maximumOrder = -1,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame );

    /*!
     * \brief Constructor for Stokes coefficient datasets.
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
    ComaWindModel( const simulation_setup::ComaStokesDataset& xStokesDataset,
                   const simulation_setup::ComaStokesDataset& yStokesDataset,
                   const simulation_setup::ComaStokesDataset& zStokesDataset,
                   std::shared_ptr<ComaModel> comaModel,
                   std::function<Eigen::Vector6d()> sunStateFunction,
                   std::function<Eigen::Vector6d()> cometStateFunction,
                   std::function<Eigen::Matrix3d()> cometRotationFunction,
                   const int& maximumDegree = -1,
                   const int& maximumOrder = -1,
                   const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame );

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

    //! Polynomial coefficient dataset for x-component wind (used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> xPolyDataset_;
    //! Polynomial coefficient dataset for y-component wind (used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> yPolyDataset_;
    //! Polynomial coefficient dataset for z-component wind (used when dataType_ == 0)
    std::shared_ptr<simulation_setup::ComaPolyDataset> zPolyDataset_;

    //! Stokes coefficient dataset for x-component wind (used when dataType_ == 1)
    std::shared_ptr<simulation_setup::ComaStokesDataset> xStokesDataset_;
    //! Stokes coefficient dataset for y-component wind (used when dataType_ == 1)
    std::shared_ptr<simulation_setup::ComaStokesDataset> yStokesDataset_;
    //! Stokes coefficient dataset for z-component wind (used when dataType_ == 1)
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

    //! Pre-initialized interpolators for x-component Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> xStokesInterpolators_;
    //! Pre-initialized interpolators for y-component Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> yStokesInterpolators_;
    //! Pre-initialized interpolators for z-component Stokes coefficients: map from (n,m) to (cosine, sine) interpolators
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> zStokesInterpolators_;

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
     * \brief Compute wind component from polynomial coefficients.
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
     * \brief Compute wind component from Stokes coefficients using pre-initialized interpolators.
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
     * \brief Calculate solar longitude in comet body-fixed frame.
     * \return Solar longitude angle from X-axis in XY plane [rad]
     * \throws std::runtime_error If state functions are not initialized
     */
    double calculateSolarLongitude() const;

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
};

}  // namespace aerodynamics

}  // namespace tudat

#endif  // TUDAT_WIND_MODEL_H
