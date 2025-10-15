#ifndef COMAMODEL_H
#define COMAMODEL_H


#include <memory>


#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"

#include <functional>
#include <map>


namespace tudat
{
namespace aerodynamics
{
class SphericalHarmonicsCalculator;

/*!
 * \brief Enumeration for the type of coefficient data used in coma models.
 */
enum class ComaDataType
{
    POLYNOMIAL_COEFFICIENTS,  //!< Use polynomial coefficients for density computation
    STOKES_COEFFICIENTS      //!< Use Stokes coefficients for density computation
};

/*!
 * \brief Model for computing coma density, pressure, temperature, and sound speed.
 *
 * This class implements an atmosphere model specifically designed for comet comas.
 * It can use either polynomial coefficients or Stokes coefficients to represent
 * the spatial and temporal variation of coma properties. The model uses spherical
 * harmonics to efficiently compute atmospheric properties at any location.
 */
class ComaModel final : public AtmosphereModel
{
public:
    /*!
     * \brief Constructor for polynomial coefficient data.
     * \param polyDataset Structured polynomial coefficient dataset
     * \param molecularWeight Molecular weight of the gas species [kg/mol]
     * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     * \param cometRotationFunction Function returning comet body-fixed rotation matrix
     * \param maximumDegree Maximum degree used to compute the coma density with SH (-1 for auto)
     * \param maximumOrder Maximum Order used to compute the coma density with SH (-1 for auto)
     */
    ComaModel( const simulation_setup::ComaPolyDataset& polyDataset,
               const double molecularWeight,
               std::function<Eigen::Vector6d()> sunStateFunction,
               std::function<Eigen::Vector6d()> cometStateFunction,
               std::function<Eigen::Matrix3d()> cometRotationFunction,
               const int& maximumDegree = -1,
               const int& maximumOrder = -1 );

    /*!
     * \brief Constructor for Stokes coefficient data.
     * \param stokesDataset Structured Stokes coefficient dataset
     * \param molecularWeight Molecular weight of the gas species [kg/mol]
     * \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     * \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     * \param cometRotationFunction Function returning comet body-fixed rotation matrix
     * \param maximumDegree Maximum degree used to compute the coma density with SH (-1 for auto)
     * \param maximumOrder Maximum Order used to compute the coma density with SH (-1 for auto)
     */
    ComaModel( const simulation_setup::ComaStokesDataset& stokesDataset,
               const double molecularWeight,
               std::function<Eigen::Vector6d()> sunStateFunction,
               std::function<Eigen::Vector6d()> cometStateFunction,
               std::function<Eigen::Matrix3d()> cometRotationFunction,
               const int& maximumDegree = -1,
               const int& maximumOrder = -1 );


    /*!
     * \brief Returns the local density of the coma in kg per meter^3.
     * \param radius Radius from comet center at which density is to be computed [m]
     * \param longitude Longitude in comet body-fixed frame at which density is to be computed [rad]
     * \param latitude Latitude in comet body-fixed frame at which latitude is to be computed [rad]
     * \param time Time at which density is to be computed [s]
     * \return Coma density at specified location and time [kg/m³]
     */
    double getDensity( double radius,
                       double longitude,
                       double latitude,
                       double time ) override;

    /*!
     * \brief Returns the local number density of the coma in particles per meter^3.
     * \param radius Radius from comet center at which number density is to be computed [m]
     * \param longitude Longitude in comet body-fixed frame at which number density is to be computed [rad]
     * \param latitude Latitude in comet body-fixed frame at which number density is to be computed [rad]
     * \param time Time at which number density is to be computed [s]
     * \return Coma number density at specified location and time [m^-3]
     */
    double getNumberDensity( double radius,
                             double longitude,
                             double latitude,
                             double time );

    /*!
     * \brief Returns the local pressure of the coma in Newton per meter^2.
     * \param radius Radius from comet center at which pressure is to be computed [m]
     * \param longitude Longitude in comet body-fixed frame at which pressure is to be computed [rad]
     * \param latitude Latitude in comet body-fixed frame at which pressure is to be computed [rad]
     * \param time Time at which pressure is to be computed [s]
     * \return Coma pressure at specified location and time [N/m²]
     */
    double getPressure( double radius,
                         double longitude,
                        double latitude,
                        double time ) override;

    /*!
     * \brief Returns the local temperature of the coma in Kelvin.
     * \param radius Radius from comet center at which temperature is to be computed [m]
     * \param longitude Longitude in comet body-fixed frame at which temperature is to be computed [rad]
     * \param latitude Latitude in comet body-fixed frame at which temperature is to be computed [rad]
     * \param time Time at which temperature is to be computed [s]
     * \return Coma temperature at specified location and time [K]
     */
    double getTemperature( double radius,
                               double longitude,
                           double latitude,
                           double time ) override;

    /*!
     * \brief Returns the speed of sound in the coma in m/s.
     * \param radius Radius from comet center at which speed of sound is to be computed [m]
     * \param longitude Longitude in comet body-fixed frame at which speed of sound is to be computed [rad]
     * \param latitude Latitude in comet body-fixed frame at which speed of sound is to be computed [rad]
     * \param time Time at which speed of sound is to be computed [s]
     * \return Coma speed of sound at specified location and time [m/s]
     */
    double getSpeedOfSound( double radius,
                                 double longitude,
                            double latitude,
                            double time ) override;

    /*!
     * \brief Get the data type used by this ComaModel instance.
     * \return ComaDataType indicating whether polynomial or Stokes coefficients are used
     */
    inline ComaDataType getDataType() const { return dataType_; }

    /*!
     * \brief Get a pointer to the spherical harmonics calculator for sharing.
     * \return Pointer to the spherical harmonics calculator (non-owning)
     */
    inline SphericalHarmonicsCalculator* getSphericalHarmonicsCalculator() const { return sphericalHarmonicsCalculator_.get(); }

    /*!
     * \brief Get the current solar longitude value.
     * \return Solar longitude in comet body-fixed frame [rad]
     * \note Returns the cached value from the most recent density computation. If no density has been computed yet, returns 0.0.
     */
    inline double getSolarLongitude() const { return cachedSolarLongitude_; }

private:
    // ========== Hot path: Frequently accessed cached values (grouped for cache locality) ==========

    //! Cached solar longitude value to avoid redundant state function calls [rad]
    mutable double cachedSolarLongitude_;

    //! Time at which solar longitude was last cached [s]
    mutable double cachedTime_;

    //! Pre-allocated coefficient matrices to avoid repeated allocations
    mutable Eigen::MatrixXd cachedCosineCoefficients_;
    mutable Eigen::MatrixXd cachedSineCoefficients_;

    //! Cached file index from last time interval search (optimization hint)
    mutable int lastFileIndex_;

    //! Flag indicating whether cached solar longitude is valid
    mutable bool solarLongitudeCacheValid_;

    // ========== Configuration: Model parameters (small, accessed during initialization and queries) ==========

    //! Type of data used to determine computation method (polynomial or Stokes coefficients)
    ComaDataType dataType_;

    //! Molecular weight of gas species [kg/mol]
    double molecularWeight_;

    //! Maximum spherical harmonic degree used for density computation (-1 for auto-detect)
    int maximumDegree_;

    //! Maximum spherical harmonic order used for density computation (-1 for auto-detect)
    int maximumOrder_;

    // ========== Data and computation infrastructure (large objects, less frequently accessed) ==========

    //! Polynomial coefficient dataset containing coma density data (used when dataType_ == POLYNOMIAL_COEFFICIENTS)
    std::shared_ptr<simulation_setup::ComaPolyDataset> polyDataset_;

    //! Stokes coefficient dataset containing coma density data (used when dataType_ == STOKES_COEFFICIENTS)
    std::shared_ptr<simulation_setup::ComaStokesDataset> stokesDataset_;

    //! Function returning Sun state vector (position [m], velocity [m/s]) in inertial frame
    std::function<Eigen::Vector6d()> sunStateFunction_;

    //! Function returning Comet state vector (position [m], velocity [m/s]) in inertial frame
    std::function<Eigen::Vector6d()> cometStateFunction_;

    //! Function returning comet body-fixed to inertial frame rotation matrix
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

    //! Spherical harmonics calculator with cached computations for efficient evaluation
    std::unique_ptr<SphericalHarmonicsCalculator> sphericalHarmonicsCalculator_;

    //! Pre-initialized interpolators for Stokes coefficients (only used for STOKES_COEFFICIENTS data type)
    //! Vector indexed by file, each containing a map from (n,m) pairs to cosine and sine coefficient interpolators
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>>> stokesInterpolators_;

    //! Pre-initialized interpolators for reduced Stokes coefficients (for radius > reference radius)
    //! Vector indexed by file, each containing a map from (n,m) pairs to cosine and sine coefficient interpolators (1D: solar longitude only)
    std::vector<std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 1>>>>> reducedStokesInterpolators_;

    /*!
     * @brief Find the index of the time interval that contains a given time.
     *
     * This function searches through the defined time periods to find which one contains the input time.
     *
     * @param time The time for which the corresponding time interval index is to be found.
     * @return The index of the time period in which the given time falls.
     *
     * @throws std::runtime_error If no matching time interval is found.
     */
    int findTimeIntervalIndex( double time ) const;

    /*!
     * \brief Compute number density from polynomial coefficients.
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute number density [s]
     * \return Coma number density [m^-3]
     * \throws std::runtime_error If dataset is null or time is out of range
     */
    double computeNumberDensityFromPolyCoefficients( double radius, double longitude, double latitude, double time ) const;

    /*!
     * \brief Compute number density from Stokes coefficients using pre-initialized interpolators.
     * \param radius Radial distance from comet center [m]
     * \param longitude Longitude in comet body-fixed frame [rad]
     * \param latitude Latitude in comet body-fixed frame [rad]
     * \param time Time at which to compute number density [s]
     * \return Coma number density [m^-3]
     * \throws std::runtime_error If dataset is null or time is out of range
     */
    double computeNumberDensityFromStokesCoefficients( double radius, double longitude, double latitude, double time ) const;

    /*!
     * @brief Calculate solar longitude in comet body-fixed frame with caching.
     * @param time Time at which to compute solar longitude [s]
     * @return Solar longitude [rad]
     */
    double calculateSolarLongitude( double time ) const;

    /*!
     * @brief Initialize interpolators for Stokes coefficients (called only for STOKES_COEFFICIENTS data type)
     */
    void initializeStokesInterpolators();

    /*!
     * @brief Helper to create 2D interpolator for Stokes coefficients on-the-fly (fallback)
     * @param fileIndex File index in dataset
     * @param n Spherical harmonic degree
     * @param m Spherical harmonic order
     * @param cosineCoeff Output: interpolated cosine coefficient
     * @param sineCoeff Output: interpolated sine coefficient
     * @param radius Radius for interpolation [m]
     * @param solarLongitude Solar longitude for interpolation [rad]
     */
    void createFallback2DInterpolator( int fileIndex, int n, int m,
                                       double& cosineCoeff, double& sineCoeff,
                                       double radius, double solarLongitude ) const;

    /*!
     * @brief Helper to create 1D reduced interpolator for Stokes coefficients on-the-fly (fallback)
     * @param fileIndex File index in dataset
     * @param n Spherical harmonic degree
     * @param m Spherical harmonic order
     * @param cosineCoeff Output: interpolated cosine coefficient
     * @param sineCoeff Output: interpolated sine coefficient
     * @param solarLongitude Solar longitude for interpolation [rad]
     */
    void createFallback1DInterpolator( int fileIndex, int n, int m,
                                       double& cosineCoeff, double& sineCoeff,
                                       double solarLongitude ) const;

};


/*!
 * \brief Calculator for efficient spherical harmonics evaluation.
 *
 * This class provides optimized computation of spherical harmonic expansions
 * with internal caching to improve performance for repeated evaluations.
 * It is designed to be shared between multiple models (e.g., ComaModel and
 * ComaWindModel) to avoid redundant computations.
 */
class SphericalHarmonicsCalculator
{
public:
    /*!
     * \brief Constructor for spherical harmonics calculator.
     * \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional)
     */
    explicit SphericalHarmonicsCalculator( std::string fixedReferenceFrame = "" );

    //! Destructor
    ~SphericalHarmonicsCalculator( ) = default;


    /*!
     * \brief Evaluate the surface spherical harmonics field at a given latitude and longitude.
     *
     * This function evaluates a scalar field expressed as a spherical harmonic expansion
     * (e.g., density, gas velocity, etc.) on the surface of a body, using
     * geodesy-normalised Legendre polynomials and stored spherical harmonic coefficients.
     * \param sineCoefficients Sine coefficients
     * \param cosineCoefficients Cosine coefficients
     * \param latitude Geocentric latitude [rad].
     * \param longitude Longitude [rad].
     * \param highestDegree Highest spherical harmonic degree to include in the evaluation.
     * \param highestOrder Highest spherical harmonic order to include in the evaluation.
     * \return Evaluated surface field value at the given latitude and longitude.
     */
    double calculateSurfaceSphericalHarmonics( const Eigen::MatrixXd& sineCoefficients,
                                               const Eigen::MatrixXd& cosineCoefficients,
                                               double latitude,
                                               double longitude,
                                               int highestDegree,
                                               int highestOrder );

private:
    //! Identifier for the body-fixed reference frame used for coordinate transformations
    std::string fixedReferenceFrame_;

    //! Cache object for efficient spherical harmonics calculations (stores Legendre polynomials and trigonometric values)
    basic_mathematics::SphericalHarmonicsCache sphericalHarmonicsCache_;

    //! Last used maximum degree (for cache optimization)
    int lastMaxDegree_ = -1;

    //! Last used maximum order (for cache optimization)
    int lastMaxOrder_ = -1;
};
} // end namespace aerodynamics
} // end namespace tudat


#endif //COMAMODEL_H