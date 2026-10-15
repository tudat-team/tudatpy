#ifndef TUDAT_ATMOSPHERIC_TIDE_LOADING_H
#define TUDAT_ATMOSPHERIC_TIDE_LOADING_H
#include "tudat/io/readAtmosphericTideData.h"

namespace tudat
{
namespace basic_astrodynamics
{
//! Pure IERS Eqs.19a/b: METRE coefficients in radial/north/east rows; UT1 seconds since J2000 noon.
//! The accepted RP03 convention has midnight phase origin; no further 20-minute or S1 180-degree shift.
Eigen::Vector3d synthesizeAtmosphericTideLoading( const Eigen::Matrix< double, 3, 4 >& coefficientsMetres,
                                                  const double ut1SecondsSinceJ2000 );
//! Exact analytic derivative with respect to UT1 seconds, in metres/second.
//! A derivative with respect to an environment epoch must additionally multiply by dUT1/d(epoch).
Eigen::Vector3d differentiateAtmosphericTideLoadingWrtUt1( const Eigen::Matrix< double, 3, 4 >& coefficientsMetres,
                                                           const double ut1SecondsSinceJ2000 );

//! Pure Cartesian dX/dY/dZ of IERS Eq.20/Table7.6, in metres; UT1 seconds from J2000 noon.
//! This helper does not change any displacement's reference origin. Its RP03 applicability must be checked by a caller.
Eigen::Vector3d calculateIers2010Rp03AtmosphericGeocentreTranslation( const double ut1SecondsSinceJ2000 );

//! Cartesian derivative of the same table translation, in metres per UT1 second.
Eigen::Vector3d differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( const double ut1SecondsSinceJ2000 );

//! Validated interpolation object; no file I/O, time conversion, station latitude inference or origin conversion.
//! Duplicated seam columns are averaged, and polar rows are replaced by their Cartesian mean, on a model copy only.
//! Provider rounding makes these duplicates differ slightly (below 1 micrometre for RP03); a grid whose duplicates
//! differ by more than maximumEndpointDiscrepancyMetres is rejected as inconsistent.
class AtmosphericTideGridInterpolator
{
public:
    static constexpr double maximumEndpointDiscrepancyMetres = 1.0e-5;

    AtmosphericTideGridInterpolator( const input_output::AtmosphericTideGrid& grid );
    //! Interpolate metre coefficients at degree coordinates; missing spatial coverage raises std::out_of_range.
    Eigen::Matrix< double, 3, 4 > interpolate( const double providerLatitudeDegrees, const double providerLongitudeDegrees ) const;

    //! The original neutral grid, including both conflicting endpoints, is never modified.
    const input_output::AtmosphericTideGrid& getProviderGrid( ) const
    {
        return providerGrid_;
    }
    double getMaximumSeamDiscrepancyMetres( ) const
    {
        return seamDiscrepancy_;
    }
    double getMaximumPolarCartesianDiscrepancyMetres( ) const
    {
        return polarDiscrepancy_;
    }

private:
    input_output::AtmosphericTideGrid providerGrid_;
    std::vector< Eigen::Matrix< double, 3, 4 > > modelCoefficients_;
    bool periodic_ = false;
    bool duplicatedSeam_ = false;
    double seamDiscrepancy_ = 0;
    double polarDiscrepancy_ = 0;
};
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
