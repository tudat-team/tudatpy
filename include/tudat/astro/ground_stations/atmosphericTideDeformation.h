#ifndef TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H
#define TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H

#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! S1/S2 atmospheric tide loading displacements (IERS 2010 Conventions, Section 7.1.3).
/*!
 * The coefficients are taken from a grid (interpolated at the geocentric latitude and longitude of the station) or per
 * station (selected by the station name). Their radial, north and east components are geocentric (as for the Ray and Ponte
 * (2003) model). The phase is computed from UT1. The displacements are with respect to the centre of the solid Earth, unless
 * the IERS 2010 (Table 7.6) translation to the centre of mass of the whole Earth is added.
 */
class AtmosphericTideDeformation : public BodyDeformationModel
{
public:
    AtmosphericTideDeformation( const input_output::AtmosphericTideGrid& grid,
                                const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
                                const bool convertToCentreOfMass = false ):
        BodyDeformationModel( converter ), grid_( grid ), convertToCentreOfMass_( convertToCentreOfMass )
    {}

    AtmosphericTideDeformation( const std::map< std::string, input_output::AtmosphericTideCoefficients >& stationCoefficients,
                                const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
                                const bool convertToCentreOfMass = false ):
        BodyDeformationModel( converter ), stationCoefficients_( stationCoefficients ), convertToCentreOfMass_( convertToCentreOfMass )
    {}

protected:
    //! Cartesian displacement in metres.
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;

    //! Cartesian displacement and its rate (per TDB second).
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    //! Displacement and (if requested) its rate for a station.
    Eigen::Vector6d evaluate( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState >& station, bool withRate );

    //! Grid of coefficients (empty if station coefficients are used), and coefficients per station name.
    input_output::AtmosphericTideGrid grid_;
    std::map< std::string, input_output::AtmosphericTideCoefficients > stationCoefficients_;

    //! Whether the translation from the centre of the solid Earth to the centre of mass of the whole Earth is added.
    bool convertToCentreOfMass_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
