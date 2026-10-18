/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_IERS202SOLIDTIDALBODYDEFORMATION_H
#define TUDAT_IERS202SOLIDTIDALBODYDEFORMATION_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>
#include <memory>

#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/interface/sofa/fundamentalArguments.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Doodson arguments for the corrected IERS DEHANTTIDEINEL solid Earth tide model.
/*!
 * Uses the STEP2DIU/STEP2LON polynomials, including their lunar-longitude
 * precession correction, and UTC fractional hours as specified by the reference
 * routine (2016-12-19 revision, including the 2013 fractional-hour correction).
 * TT centuries are formed from the input and UTC is obtained with SOFA's
 * TT-to-UTC conversion. Native model calls convert their precise TDB input to TT
 * with the model-owned converter before evaluating these scalar polynomials.
 * These arguments are specific to this station displacement convention.
 * \param secondsSinceJ2000 TT seconds since J2000.
 * \return Doodson arguments in radians: tau, s, h, p, N-prime, p-s.
 */
Eigen::Vector6d calculateIers2010SolidEarthTideArguments( const double secondsSinceJ2000 );

namespace iers_2010_parameters
{

const static double PRINCIPAL_DEGREE_TWO_LOVE_NUMBER = 0.6078;

const static double PRINCIPAL_DEGREE_THREE_LOVE_NUMBER = 0.292;

const static double PRINCIPAL_DEGREE_TWO_SHIDA_NUMBER = 0.0847;

const static double PRINCIPAL_DEGREE_THREE_SHIDA_NUMBER = 0.015;

const static double DEGREE_TWO_LATITUDE_LOVE_NUMBER = -0.0006;

const static double DEGREE_TWO_LATITUDE_SHIDA_NUMBER = 0.0002;

const static double IMAGINARY_DEGREE_TWO_DIURNAL_LOVE_NUMBER = -0.0025;

const static double IMAGINARY_DEGREE_TWO_DIURNAL_SHIDA_NUMBER = -0.0007;

const static double IMAGINARY_DEGREE_TWO_SEMIDIURNAL_LOVE_NUMBER = -0.0022;

const static double IMAGINARY_DEGREE_TWO_SEMIDIURNAL_SHIDA_NUMBER = -0.0007;

const static double DEGREE_TWO_DIURNAL_TOROIDAL_LOVE_NUMBER = 0.0012;

const static double DEGREE_TWO_SEMIDIURNAL_TOROIDAL_LOVE_NUMBER = 0.0024;

}  // namespace iers_2010_parameters

//! Calculate trigonometric expressions for correction fo site displacements, step 1 of IERS Conventions 2010, Chapter 7
std::vector< Eigen::Vector3d > calculateDisplacementCorrectionTrigonometricPart( const Eigen::Vector3d& relativeBodyState,
                                                                                 const Eigen::Vector3d& stationSphericalCoordinates,
                                                                                 const std::vector< bool >& calculateTerms );

Eigen::Vector3d calculateFirstCorrectionStepDisplacements( const Eigen::Vector3d& relativeBodyState,
                                                           const Eigen::Vector3d& stationSphericalCoordinates,
                                                           const double gravitationalParameterRatio,
                                                           const double bodyEquatorialRadius,
                                                           const std::vector< double >& correctionLoveAndShidaNumbers,
                                                           const std::vector< bool >& areFirstStepCorrectionsCalculated );

//! Function to calculate diurnal tide-frequency dependent site displacement correction to solid Earth tide displacements.
/*!
 *  Function to calculate diurnal tide-frequency dependent site displacement correction to solid Earth tide displacements,
 *  per step II as described in Section 7.1.1 of the IERS 2010 Conventions.
 *  \param tideArgument Current phase of tide for which correction is to be applied.
 *  \param stationLongitude Longitude of station at which correction is to be calculated.
 *  \param stationLatitude Latitude of station at which correction is to be calculated.
 *  \param Amplitudes of tide corrections. They should be provided as,
 *  in order, in-phase radial (i.e locally up), out-of-phase radial, in-phase transverse (i.e. locally horizontal),
 *  out-of-phase transverse.
 */
Eigen::Vector3d calculateDiurnalFrequencyDependentDisplacementCorrection( const double tideArgument,
                                                                          const double stationLongitude,
                                                                          const double stationLatitude,
                                                                          const Eigen::Vector4d& amplitudes );

//! Function to calculate long period tide-frequency dependent site displacement correction to solid Earth tide displacements.
/*!
 *  Function to calculate long period tide-frequency dependent site displacement correction to solid Earth tide displacements,
 *  per step II as described in Section 7.1.1 of the IERS 2010 Conventions.
 *  \param tideArgument Current phase of tide for which correction is to be applied.
 *  \param stationLatitude Latitude of station at which correction is to be calculated.
 *  \param Amplitudes of tide corrections. They should be provided as,
 *  in order, in-phase radial (i.e locally up), out-of-phase radial, in-phase transverse (i.e. locally horizontal),
 *  out-of-phase transverse.
 */
Eigen::Vector3d calculateLongPeriodFrequencyDependentDisplacementCorrection( const double tideArgument,
                                                                             const double stationLatitude,
                                                                             const Eigen::Vector4d& amplitudes );

//! Function to calculate tide-frequency dependent site displacement corrections to solid Earth tide displacements.
/*!
 *  Function to calculate tide-frequency dependent site displacement corrections to solid Earth tide displacements,
 *  per step II as described in Section 7.1.1 of the IERS 2010 Conventions.
 *  \param doodsonMultipliers. Doodson multipliers of the tides for which corrections are tobe applied, the 6
 *  multipliers of each tide are given on subsequent rows of the Matrix.
 *  \param doodsonArguments. Doodson arguments at time at which displacements are to be calculated.
 *  \param amplitudes Amplitudes of corrections. Four amplitudea are provided per tide (row). They should be provided as,
 *  in order, in-phase radial (i.e locally up), out-of-phase radial, in-phase transverse (i.e. locally horizontal),
 *  out-of-phase transverse.
 *  \param stationSphericalCoordinates Spherical coordinates of site where displacement corrections are to be calculated.
 */

Eigen::Vector3d calculateFrequencyDependentDisplacementCorrections( const Eigen::MatrixXd& doodsonMultipliers,
                                                                    const Eigen::Matrix< double, 6, 1 >& doodsonArguments,
                                                                    const Eigen::MatrixXd& amplitudes,
                                                                    const Eigen::Vector3d& stationSphericalCoordinates );

class Iers2010EarthDeformation : public BasicTidalBodyDeformation
{
public:
    // An empty argument callback selects the corrected native UTC-based phases. Their
    // numerical rates reject UTC offset steps; explicit custom callbacks keep their own time basis.
    Iers2010EarthDeformation( const std::function< Eigen::Vector6d( const Time& ) > deformedBodyEphemeris,
                              const std::vector< std::function< Eigen::Vector6d( const Time& ) > >& deformingBodyEphemerides,
                              const std::function< Eigen::Quaterniond( const Time& ) > deformedBodyRotation,
                              const std::function< double( ) > gravitionalParameterOfDeformedBody,
                              const std::vector< std::function< double( ) > >& gravitionalParametersOfDeformingBodies,
                              const double deformedBodyEquatorialRadius,
                              const std::map< int, std::pair< double, double > >& nominalDisplacementLoveNumbers,
                              const std::vector< double >& firstStepLatitudeDependenceTerms,
                              const std::vector< bool >& areFirstStepCorrectionsCalculated,
                              const std::vector< double >& correctionLoveAndShidaNumbers,
                              const std::string& diurnalFile = "",
                              const std::string& longPeriodFile = "",
                              const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction = {},
                              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > contextClock = nullptr );

    ~Iers2010EarthDeformation( ) {}

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    Eigen::Vector3d displacementFromArguments( const std::shared_ptr< ground_stations::GroundStationState > station,
                                               const Eigen::Vector6d& arguments );

    //! Calculate site displacements due to solid Earth tide deformation according to first step of Section 7.1.1 of IERS 2010 Conventions.
    /*!
     *  Calculate site displacements due to solid Earth tide deformation according to first step of Section 7.1.1 of IERS 2010 Conventions.
     *  The first part of step 1 (Eqs. 7.5 and 7.6) give results in planetcentered and fixed frame (ITRS), whereas the corrections
     *  from Eqs 7.8, 7.9, 7.10 and 7.11 are calculated in a local (east, north, up) frame. As such, two Vector3ds
     *  are returned in a pair and are to be combined subsequently by a different function.
     *  \param tideRaisingBodiesPositions Positions of bodies raising tides in earth-fixed, centered frame
     *  \param nominalSiteState Nominal state of site at which displacements are to be calculated.
     *  \return Pair of Vector3d containing displacemnts in, first: planet-centered frame, second: local frame (see above)
     */
    std::pair< Eigen::Vector3d, Eigen::Vector3d > calculateFirstStepDisplacements(
            const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState );

    std::vector< double > firstStepLatitudeDependenceTerms_;

    std::vector< bool > areFirstStepCorrectionsCalculated_;

    std::vector< double > correctionLoveAndShidaNumbers_;

    bool areFrequencyDependentTermsCalculated_;

    //! Matrix containing doodson multipliers of tide constitudent for which corrections are to be applied.
    /*!
     *  Matrix containing doodson multipliers of tide constitudent for which corrections are to be applied.
     *  Typical values can be found in Tables 7.3a and 7.3b of IERS 2010 Conventions.
     */
    Eigen::MatrixXd doodsonMultipliers_;

    //! Matrix containing amplitudes of tide constitudent for which corrections are to be applied.
    /*!
     *  Matrix containing amplitudes of tide constitudent for which corrections are to be applied.
     *  Typical values can be found in Tables 7.3a and 7.3b of IERS 2010 Conventions.
     */
    Eigen::MatrixXd tideCorrectionAmplitudes_;

    std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction_;

    // Empty/default callback selects native UTC phases; explicit callbacks retain their own time basis.
    bool usesNativeUtcPhase_;
};

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_IERS202SOLIDTIDALBODYDEFORMATION_H
