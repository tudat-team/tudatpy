#ifndef COMAMODEL_H
#define COMAMODEL_H


#include <memory>

#include <cmath>
#include <utility>

#include "tudat/basics/utilityMacros.h"

#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/math/basic/sphericalHarmonics.h"

#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace aerodynamics
{
class ComaModel final : public AtmosphereModel
{
public:
    /*!
     *  Default constructor.
     *  \param polyCoefficients Coefficients from which the SH coefficients are computed
     *  \param powersInvRadius
     *  \param referenceRadius Reference radius of the polyCoefficient Matrix. Extrapolation is used
     *  to compute the density outside the reference radius
     *  \param TimePeriods  Time Periods for which each polyCoefficient Matrix is valid for
     *  \param SHDegreeAndOrder Maximum Degree/Order of SH coefficients, which can be computed from
     *  poly Coefficients
     *  \param maximumDegree Maximum degree used to compute the coma density with SH
     *  \param maximumOrder Maximum Order used to compute the coma density with SH
     *
     */
    ComaModel( const std::vector< Eigen::MatrixXd >& polyCoefficients,
               const std::vector< Eigen::ArrayXXi >& SHDegreeAndOrder,
               const std::vector< Eigen::VectorXd >& powersInvRadius,
               const std::vector< double >& referenceRadius,
               const std::vector< std::vector< double > >& TimePeriods,
               const int& maximumDegree,
               const int& maximumOrder );

    /*!
     * Returns the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude at which density is to be computed.
     * \param longitude Longitude at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which density is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric density at specified altitude.
     */
    double getDensity( double altitude,
                       // TODO: update information text
                       double longitude,
                       double latitude,
                       double time ) override;

    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude  at which pressure is to be computed.
     * \param longitude Longitude at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which pressure is to be computed (not used but included for
     * consistency with base class interface).
     * \return Atmospheric pressure at specified altitude.
     */
    double getPressure( double altitude,
                        // TODO: update information text
                        double longitude,
                        double latitude,
                        double time ) override;

    /*!
     * Returns the local temperature of the atmosphere in Kelvin.
     * \param altitude Altitude at which temperature is to be computed (not used since
     * temperature is assumed to be constant).
     * \param longitude Longitude at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \param latitude Latitude at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \param time Time at which temperature is to be computed (not used but included for
     * consistency with base class interface).
     * \return constantTemperature Atmospheric temperature at specified altitude.
     */
    double getTemperature( double altitude,
                           // TODO: update information text
                           double longitude,
                           double latitude,
                           double time ) override;

    /*!
     * Returns the speed of sound in the atmosphere in m/s.
     * @param altitude Altitude at which speed of sounds is to be computed.
     * @param longitude Longitude at which speed of sounds is to be computed (not used but included
     * for consistency with base class interface).
     * @param latitude Latitude at which speed of sounds is to be computed (not used but included
     * for consistency with base class interface).
     * @param time Time at which speed of sounds is to be computed (not used but included for
     * consistency with base class interface).
     * @return Atmospheric speed of sounds at specified altitude.
     */
    double getSpeedOfSound( double altitude,
                            // TODO: update information text
                            double longitude,
                            double latitude,
                            double time ) override;

    //! Static test wrapper for evaluating Stokes coefficients (for unit testing only)
    static void testEvaluateStokesCoefficients2D(
            double distanceToCometCentre,
            double solarLongitude,
            const Eigen::ArrayXXd& polyCoefficients,
            const Eigen::ArrayXXi& atSHDegreeAndOrder,
            const Eigen::VectorXd& atPowersInvRadius,
            double refRadius,
            Eigen::MatrixXd& cosineCoefficients,
            Eigen::MatrixXd& sineCoefficients,
            int maxDegree = -1,
            int maxOrder = -1 );

private:
    //! Maximum degree used for computation of density
    int maximumDegree_;

    //! Maximum order used for computation of density
    int maximumOrder_;

    //! Polynomial Coefficients of input files
    std::vector< Eigen::MatrixXd > polyCoefficients_;

    //! Spherical Harmonics Degree and Order of input tables
    std::vector< Eigen::ArrayXXi > SHDegreeAndOrder_;

    //! Reference radius of input table
    std::vector< double > referenceRadius_;

    //! Power Inverse Radius of input tables
    std::vector< Eigen::VectorXd > powersInvRadius_;

    //! Time periods where input tables are valid
    std::vector< std::vector< double > > TimePeriods_;

    //TODO: figure out Temperature dependency
    double constantTemperature_;

    // TODO: figure out specific gas constant dependency
    double specificGasConstant_;

    // TODO: figure out about ratio of specific heat
    double ratioOfSpecificHeats_;

    /*!
     * @brief Find the index of the time interval that contains a given time.
     *
     * This function searches through the defined time periods to find which one contains the input time.
     * Each time period is expected to be a vector of exactly two elements: the start and end times (in Julian days).
     *
     * @param time The time (in Julian days) for which the corresponding time interval index is to be found.
     * @return The index of the time period in which the given time falls.
     *
     * @throws std::runtime_error If no matching time interval is found, or if any time period does not consist of exactly two values.
     */
    int findTimeIntervalIndex( double time ) const;


    /*!
     * @brief Evaluate 2D Stokes coefficients from polynomial input.
     *
     * This function computes the Stokes coefficients \( C_{nm} \) and \( S_{nm} \) of the spherical harmonic
     * coma density expansion based on a set of temporally and radially dependent polynomial coefficients.
     * The evaluation depends on whether the evaluation radius is less than or greater than a reference radius.
     *
     * @param distanceToCometCentre Current radius (distance from nucleus centre) at which the density is evaluated [m].
     * @param solarLongitude Solar longitude (angle in radians) of the Sun direction in the body-fixed frame of the comet [rad].
     * @param polyCoefficients Matrix of stacked polynomial coefficients: (numRadialTerms × numTimeIntervals) × numCoeffs.
     * @param atDegreeAndOrder Integer matrix (2 × numCoefs) containing spherical harmonic degree (row 0)
     *        and order (row 1)
     * @param atPowersInvRadius Vector of powers of inverse radius (1/r)^k for polynomial evaluation.
     * @param refRadius Reference radius [m] used for determining which computation branch to apply.
     * @param maxDegree
     * @param maxOrder
     * @param cosineCoefficients Output matrix of cosine (C) spherical harmonic coefficients.
     * @param sineCoefficients  Output matrix of sine (S) spherical harmonic coefficients.

     */
    static void evaluateStokesCoefficients2D(
            double distanceToCometCentre,
            double solarLongitude,
            const Eigen::ArrayXXd& polyCoefficients,
            const Eigen::ArrayXXi& atDegreeAndOrder,
            const Eigen::VectorXd& atPowersInvRadius,
            double refRadius,
            Eigen::MatrixXd& cosineCoefficients,
            Eigen::MatrixXd& sineCoefficients,
            int maxDegree = -1,
            int maxOrder = -1 );

    /*!
     * @brief Evaluate a Stokes coefficient (Cnm or Snm) using a radial polynomial and temporal inverse FFT-like form.
     *
     * This function evaluates a single spherical harmonic coefficient based on a set of polynomial coefficients
     * defined over radial and temporal (Fourier) components. It performs a projection using cosine and sine
     * weights at the given solar longitude and applies a power-law scaling in radius.
     *
     * @param nucleusToSpacecraftDistance Distance from the comet nucleus to the spacecraft [m].
     * @param solarLongitude Solar longitude (in radians) of the Sun direction in the body-fixed frame [rad].
     * @param polyCoefs Matrix of shape (N × K) containing the polynomial coefficients.
     * @param atPowersInvRadius Array of exponent values used to compute powers of the inverse distance (1/r)^k.
     * @return The evaluated Stokes coefficient.
     */
    static double radialPolyvalAndTemporalIFFT(
            double nucleusToSpacecraftDistance,
            double solarLongitude,
            const Eigen::MatrixXd& polyCoefs,
            const Eigen::ArrayXd& atPowersInvRadius );


    /*!
     * @brief Evaluate a Stokes coefficient (Cnm or Snm) from a 1D radial polynomial using an inverse FFT-like form.
     *
     * This reduced version of the evaluator assumes that temporal components have already been collapsed, leaving
     * only a 1D radial polynomial. The method applies angular cosine/sine weighting at the given solar longitude
     * and a square-distance scaling.
     *
     * @param nucleusToSpacecraftDistance Distance from the comet nucleus to the spacecraft [m].
     * @param solarLongitude Solar longitude (in radians) in the body-fixed frame [rad].
     * @param polyCoefs Vector of polynomial coefficients corresponding to angular basis functions.
     * @return The evaluated Stokes coefficient.
     */
    static double radialSquaredAndTemporalIFFT(
            double nucleusToSpacecraftDistance,
            double solarLongitude,
            const Eigen::VectorXd& polyCoefs );
};


class SphericalHarmonicsDensity
{
public:
    /*!
     *  Default constructor.
     *  \param sineCoefficients Sine spherical harmonic coefficients
     *  \param cosineCoefficients Cosine spherical harmonic coefficients
     *  \param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional).
     */
    explicit SphericalHarmonicsDensity( Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Identity( 1, 1 ),
                                        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 1, 1 ),
                                        std::string fixedReferenceFrame = "" );

    ~SphericalHarmonicsDensity( ) = default;


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
    //! Cosine spherical harmonic coefficients
    Eigen::MatrixXd cosineCoefficients_;

    //! Sine spherical harmonic coefficients
    Eigen::MatrixXd sineCoefficients_;

    //! Identifier for body-fixed reference frame
    std::string fixedReferenceFrame_;

    //! Maximum degree used for computation of density
    int maximumDegree_;

    //! Maximum order used for computation of density
    int maximumOrder_;

    //! Cache object for potential calculations.
    basic_mathematics::SphericalHarmonicsCache sphericalHarmonicsCache_;
};
} // end namespace aerodynamics
} // end namespace tudat


#endif //COMAMODEL_H