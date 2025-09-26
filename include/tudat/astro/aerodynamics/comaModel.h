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

enum class ComaDataType
{
    POLYNOMIAL_COEFFICIENTS,
    STOKES_COEFFICIENTS
};

class ComaModel final : public AtmosphereModel
{
public:
    /*!
     *  Constructor for polynomial coefficient data.
     *  \param polyDataset Structured polynomial coefficient dataset
     *  \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     *  \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     *  \param cometRotationFunction Function returning comet body-fixed rotation matrix
     *  \param maximumDegree Maximum degree used to compute the coma density with SH
     *  \param maximumOrder Maximum Order used to compute the coma density with SH
     */
    ComaModel( const simulation_setup::ComaPolyDataset& polyDataset,
               std::function<Eigen::Vector6d()> sunStateFunction,
               std::function<Eigen::Vector6d()> cometStateFunction,
               std::function<Eigen::Matrix3d()> cometRotationFunction,
               const int& maximumDegree = -1,
               const int& maximumOrder = -1 );

    /*!
     *  Constructor for Stokes coefficient data.
     *  \param stokesDataset Structured Stokes coefficient dataset
     *  \param sunStateFunction Function returning Sun state vector (position, velocity) [m, m/s]
     *  \param cometStateFunction Function returning Comet state vector (position, velocity) [m, m/s]
     *  \param cometRotationFunction Function returning comet body-fixed rotation matrix
     *  \param maximumDegree Maximum degree used to compute the coma density with SH
     *  \param maximumOrder Maximum Order used to compute the coma density with SH
     */
    ComaModel( const simulation_setup::ComaStokesDataset& stokesDataset,
               std::function<Eigen::Vector6d()> sunStateFunction,
               std::function<Eigen::Vector6d()> cometStateFunction,
               std::function<Eigen::Matrix3d()> cometRotationFunction,
               const int& maximumDegree = -1,
               const int& maximumOrder = -1 );


    /*!
     * Returns the local density of the coma in kg per meter^3.
     * \param radius Radius from comet center at which density is to be computed.
     * \param longitude Longitude in comet body-fixed frame at which density is to be computed.
     * \param latitude Latitude in comet body-fixed frame at which density is to be computed.
     * \param time Time at which density is to be computed.
     * \return Coma density at specified location and time.
     */
    double getDensity( double radius,
                       double longitude,
                       double latitude,
                       double time ) override;

    /*!
     * Returns the local pressure of the coma in Newton per meter^2.
     * \param radius Radius from comet center at which pressure is to be computed.
     * \param longitude Longitude in comet body-fixed frame at which pressure is to be computed.
     * \param latitude Latitude in comet body-fixed frame at which pressure is to be computed.
     * \param time Time at which pressure is to be computed.
     * \return Coma pressure at specified location and time.
     */
    double getPressure( double radius,
                         double longitude,
                        double latitude,
                        double time ) override;

    /*!
     * Returns the local temperature of the coma in Kelvin.
     * \param radius Radius from comet center at which temperature is to be computed.
     * \param longitude Longitude in comet body-fixed frame at which temperature is to be computed.
     * \param latitude Latitude in comet body-fixed frame at which temperature is to be computed.
     * \param time Time at which temperature is to be computed.
     * \return Coma temperature at specified location and time.
     */
    double getTemperature( double radius,
                               double longitude,
                           double latitude,
                           double time ) override;

    /*!
     * Returns the speed of sound in the coma in m/s.
     * @param radius Radius from comet center at which speed of sound is to be computed.
     * @param longitude Longitude in comet body-fixed frame at which speed of sound is to be computed.
     * @param latitude Latitude in comet body-fixed frame at which speed of sound is to be computed.
     * @param time Time at which speed of sound is to be computed.
     * @return Coma speed of sound at specified location and time.
     */
    double getSpeedOfSound( double radius,
                                 double longitude,
                            double latitude,
                            double time ) override;

    /*!
     * \brief Get the data type used by this ComaModel instance.
     * \return ComaDataType indicating whether polynomial or Stokes coefficients are used
     */
    ComaDataType getDataType() const { return dataType_; }

private:
    //! Type of data used (polynomial or Stokes coefficients)
    ComaDataType dataType_;

    //! Maximum degree used for computation of density
    int maximumDegree_;

    //! Maximum order used for computation of density
    int maximumOrder_;

    //! Polynomial coefficient dataset (used when dataType_ == POLYNOMIAL_COEFFICIENTS)
    std::shared_ptr<simulation_setup::ComaPolyDataset> polyDataset_;

    //! Stokes coefficient dataset (used when dataType_ == STOKES_COEFFICIENTS)
    std::shared_ptr<simulation_setup::ComaStokesDataset> stokesDataset_;

    //! Function to get Sun state [m, m/s]
    std::function<Eigen::Vector6d()> sunStateFunction_;

    //! Function to get Comet state [m, m/s]
    std::function<Eigen::Vector6d()> cometStateFunction_;

    //! Function to get comet body-fixed rotation matrix
    std::function<Eigen::Matrix3d()> cometRotationFunction_;

    //! Spherical harmonics calculator with shared cache
    std::unique_ptr<SphericalHarmonicsCalculator> sphericalHarmonicsCalculator_;

    //! Pre-initialized interpolators for Stokes coefficients (only for STOKES_COEFFICIENTS data type)
    //! Map from (n,m) pair to (cosineInterpolator, sineInterpolator) pair
    std::map<std::pair<int,int>, std::pair<std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>,
                                           std::unique_ptr<interpolators::MultiLinearInterpolator<double, double, 2>>>> stokesInterpolators_;

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

    double computeDensityFromPolyCoefficients( double radius, double longitude, double latitude, double time ) const;

    double computeDensityFromStokesCoefficients( double radius, double longitude, double latitude, double time ) const;

    /*!
     * @brief Calculate solar longitude in comet body-fixed frame.
     * @return Solar longitude [rad]
     */
    double calculateSolarLongitude() const;

    /*!
     * @brief Initialize interpolators for Stokes coefficients (called only for STOKES_COEFFICIENTS data type)
     */
    void initializeStokesInterpolators();

};


class SphericalHarmonicsCalculator
{
public:
    /*!
     *  Default constructor.
     *  \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional).
     */
    explicit SphericalHarmonicsCalculator( std::string fixedReferenceFrame = "" );

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
    //! Identifier for body-fixed reference frame
    std::string fixedReferenceFrame_;

    //! Cache object for spherical harmonics calculations
    basic_mathematics::SphericalHarmonicsCache sphericalHarmonicsCache_;
};
} // end namespace aerodynamics
} // end namespace tudat


#endif //COMAMODEL_H