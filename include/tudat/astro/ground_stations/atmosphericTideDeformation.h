#ifndef TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H
#define TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H

#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace basic_astrodynamics
{
//! Unspecified conventions are taken from the source metadata. For a grid identified as Ray-Ponte 2003 with Farrell
//! Green functions (RP03), unspecified conventions default to geocentric lookup latitude and geocentric radial/north/east
//! components. The provider does not document this choice; geodetic and geocentric latitude differ by at most 0.19 degrees
//! on the 1-degree grid, which changes the displacement by micrometres.
//! convertToCentreOfMass adds the IERS 2010 Table 7.6 CE-to-CM translation; it requires the identified RP03 product.
struct AtmosphericTideEvaluationSettings {
    //! Latitude used to look up a station in the grid, and directions of the displacement components.
    input_output::AtmosphericTideLatitudeConvention lookupLatitudeConvention = input_output::AtmosphericTideLatitudeConvention::unspecified;
    input_output::AtmosphericTideDisplacementBasis displacementBasis = input_output::AtmosphericTideDisplacementBasis::unspecified;

    //! Whether to convert the displacement from the centre of the solid Earth to the centre of mass of the whole Earth.
    bool convertToCentreOfMass = false;
};

//! Native conventional atmospheric loading with explicit geometry and source-origin conventions.
//! Native input is precise environment TDB Time. Phase and analytic state rates use the model-owned shared clock.
//! No source file is read during evaluation; provider coefficients and source metadata remain unchanged.
class AtmosphericTideDeformation : public BodyDeformationModel
{
public:
    AtmosphericTideDeformation( const input_output::AtmosphericTideGrid& grid,
                                const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
                                const AtmosphericTideEvaluationSettings& settings = {},
                                const std::shared_ptr< BodyShapeModel >& bodyShape = nullptr );
    AtmosphericTideDeformation( const input_output::AtmosphericTideStationData& stations,
                                const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
                                const AtmosphericTideEvaluationSettings& settings = {},
                                const std::shared_ptr< BodyShapeModel >& bodyShape = nullptr );

    bool usesStationCoefficients( ) const
    {
        return gridInterpolator_ == nullptr;
    }
    const std::shared_ptr< AtmosphericTideGridInterpolator >& getGridInterpolator( ) const
    {
        return gridInterpolator_;
    }
    const input_output::AtmosphericTideStationData& getStationData( ) const
    {
        return stationData_;
    }
    const AtmosphericTideEvaluationSettings& getEvaluationSettings( ) const
    {
        return settings_;
    }
    input_output::AtmosphericLoadingReferenceOrigin getOutputReferenceOrigin( const std::string& identifier = "" ) const;

protected:
    //! Return Cartesian atmospheric displacement in metres without requesting a clock rate.
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;

    //! Return displacement and its analytic derivative in metres per TDB second.
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    //! Select coefficients/geometry once and optionally evaluate the rate with the same UT1 phase and basis.
    Eigen::Vector6d evaluateCartesianLoading( const Time& tdb,
                                              const std::shared_ptr< ground_stations::GroundStationState >& station,
                                              bool withRate );

    void validateCommonSettings( ) const;
    void validateShapeForBasis( input_output::AtmosphericTideDisplacementBasis basis, bool geodeticLookup = false ) const;
    void validateOriginConversion( const input_output::AtmosphericTideMetadata& metadata ) const;

    //! Interpolator of the grid coefficients (nullptr if the model uses station coefficients), and the station coefficients.
    std::shared_ptr< AtmosphericTideGridInterpolator > gridInterpolator_;
    input_output::AtmosphericTideStationData stationData_;

    //! Metadata of the grid and of each station, with unspecified conventions completed from the settings.
    input_output::AtmosphericTideMetadata resolvedGridMetadata_;
    std::map< std::string, input_output::AtmosphericTideMetadata > resolvedStationMetadata_;

    //! Settings for the evaluation of the displacements.
    AtmosphericTideEvaluationSettings settings_;

    //! Shape of the Earth, used to compute geodetic latitudes and directions (may be nullptr if these are not needed).
    std::shared_ptr< BodyShapeModel > bodyShape_;
};
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
