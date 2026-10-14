/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_OCEAN_POLE_TIDE_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_OCEAN_POLE_TIDE_GRAVITY_FIELD_VARIATIONS_H

#include <complex>
#include <vector>
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/io/readOceanPoleGravityData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace gravitation
{
//! Explicit IERS 6.23 forcing constants. No reference constants inferred from the A/B table.
struct OceanPoleGravityConstants {
    double angularVelocity = TUDAT_NAN;
    double gravitationalConstant = TUDAT_NAN;
    double oceanDensity = TUDAT_NAN;
    double equatorialGravity = TUDAT_NAN;
    double gravitationalParameter = TUDAT_NAN;
    double referenceRadius = TUDAT_NAN;
    std::complex< double > gamma2{ TUDAT_NAN, TUDAT_NAN };
};

//! Full normalized ocean-pole response. Wobble callback returns (m1,m2) in radians at TDB seconds.
class OceanPoleTideGravityFieldVariations : public GravityFieldVariations
{
public:
    //! loadLoveNumbers[n] is k'_n; require every selected degree, with no hidden response defaults.
    //! Degrees >=2 describe body-centered gravity. Published degree-one integrals remain in raw data.
    //! Output GM/R conversion is applied once; source constants describe the forcing potential.
    OceanPoleTideGravityFieldVariations( const input_output::OceanPoleGravityData& data,
                                         const std::function< Eigen::Vector2d( double ) >& wobbleRadians,
                                         const std::vector< double >& loadLoveNumbers,
                                         const OceanPoleGravityConstants& constants,
                                         double outputGravitationalParameter,
                                         double outputReferenceRadius,
                                         int minimumDegree,
                                         int minimumOrder,
                                         int maximumDegree,
                                         int maximumOrder );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( double time ) override;

private:
    std::function< Eigen::Vector2d( double ) > wobbleRadians_;
    Eigen::MatrixXd aReal_, aImaginary_, bReal_, bImaginary_;
    std::complex< double > gamma2_;
};

}  // namespace gravitation
}  // namespace tudat
#endif
