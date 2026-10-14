/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_EARTHSOLIDBODYTIDEGRAVITYFIELDVARIATIONS_H
#define TUDAT_EARTHSOLIDBODYTIDEGRAVITYFIELDVARIATIONS_H

#include "tudat/astro/basic_astro/tideSystem.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"

namespace tudat
{
namespace gravitation
{

//! IERS 2010 anelastic solid Earth gravity tide, equations 6.6--6.8 and tables 6.3/6.5.
/*! State callbacks share inertial origin/axes; orientation maps those axes to the
 * coefficient frame. The argument callback returns (tau,s,h,p,N',ps) in radians
 * using explicitly routed TT/UT1. Scalar epochs follow the existing gravity API.
 * Source GM/radius define the response; coefficients are rescaled once to the
 * consuming field. Only known tide-free and zero-tide nominal fields are supported.
 * This model includes no secular gravity drift or additional degree-four forcing.
 */
class EarthSolidBodyTideGravityFieldVariations : public SolidBodyTideGravityFieldVariations
{
public:
    EarthSolidBodyTideGravityFieldVariations( const std::function< Eigen::Vector6d( double ) >& deformedBodyState,
                                              const std::function< Eigen::Quaterniond( double ) >& toBodyFixedFrame,
                                              const std::vector< std::function< Eigen::Vector6d( double ) > >& deformingBodyStates,
                                              double sourceRadius,
                                              const std::function< double( ) >& sourceGravitationalParameter,
                                              const std::vector< std::function< double( ) > >& deformingBodyGravitationalParameters,
                                              const std::vector< std::string >& deformingBodies,
                                              const std::function< Eigen::Vector6d( double ) >& doodsonArguments,
                                              double fieldRadius,
                                              double fieldGravitationalParameter,
                                              basic_astrodynamics::TideSystem tideSystem );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( double time ) override;

protected:
    void addBasicSolidBodyTideCorrections( Eigen::MatrixXd& cosine, Eigen::MatrixXd& sine ) override;

private:
    std::function< Eigen::Vector6d( double ) > doodsonArguments_;
    double fieldRadius_, fieldGravitationalParameter_;
    basic_astrodynamics::TideSystem tideSystem_;
};

//! IERS 2010 section 6.4 solid pole tide, with explicitly supplied m1/m2 in arcseconds.
/*! The production adapter obtains wobble from the selected pole and shared secular
 * pole, m1=xp-xbar and m2=-(yp-ybar). Published benchmark wobble can be injected
 * without changing that production convention. Coefficient output is a 1x1 block
 * at degree 2, order 1, rescaled once from source GM/radius to the consuming field.
 */
class SolidPoleTideGravityFieldVariations : public GravityFieldVariations
{
public:
    SolidPoleTideGravityFieldVariations( const std::function< Eigen::Vector2d( double ) >& wobbleArcseconds,
                                         double sourceGravitationalParameter,
                                         double sourceRadius,
                                         double fieldGravitationalParameter,
                                         double fieldRadius );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( double time ) override;

private:
    std::function< Eigen::Vector2d( double ) > wobbleArcseconds_;
    double coefficientScale_;
};

}  // namespace gravitation
}  // namespace tudat
#endif
