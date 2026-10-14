#ifndef TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H
#define TUDAT_ATMOSPHERIC_TIDE_DEFORMATION_H

#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace basic_astrodynamics
{
enum class AtmosphericLoadingOriginPolicy { preserveSourceOrigin, convertIdentifiedRp03CeToCm };
struct AtmosphericTideEvaluationSettings {
    input_output::AtmosphericTideLatitudeConvention lookupLatitudeConvention = input_output::AtmosphericTideLatitudeConvention::unspecified;
    input_output::AtmosphericTideDisplacementBasis displacementBasis = input_output::AtmosphericTideDisplacementBasis::unspecified;
    AtmosphericTideEndpointPolicy endpointPolicy = AtmosphericTideEndpointPolicy::requireConsistency;
    AtmosphericLoadingOriginPolicy originPolicy = AtmosphericLoadingOriginPolicy::preserveSourceOrigin;
};

//! Native conventional atmospheric loading with explicit geometry and source-origin conventions.
//! Scalar epochs are environment TDB seconds since J2000. Phase uses continuous UT1 from this model's supplied converter.
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

    Eigen::Vector3d calculateDisplacement( double tdb, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacement( double tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacement( const EarthDeformationEpoch& epoch,
                                           const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    //! Actual station topocentric ENU, independently of the provider lookup/basis conventions.
    Eigen::Vector3d calculateEnuDisplacement( const EarthDeformationEpoch& epoch,
                                              const std::shared_ptr< ground_stations::GroundStationState > station );
    Eigen::Vector3d calculateEnuDisplacement( double tdb, const std::shared_ptr< ground_stations::GroundStationState > station );
    //! Exact analytic loading velocity/state with respect to TDB seconds, using this model's continuous UT1 rate.
    //! Position is evaluated independently; rate requests reject unavailable EOP rate stencils.
    Eigen::Vector3d calculateDisplacementVelocity( double tdb, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacementVelocity( double tdb,
                                                   const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector3d calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                   const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d calculateDisplacementState( double tdb, const Eigen::Vector3d& position ) override;
    Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position ) override;
    Eigen::Vector6d calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d calculateDisplacementState( double tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;

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

private:
    Eigen::Vector3d evaluateCartesianLoading( const EarthDeformationEpoch& epoch,
                                              const std::shared_ptr< ground_stations::GroundStationState >& station,
                                              bool derivativeWrtTdb );
    void validateCommonSettings( ) const;
    void validateShapeForBasis( input_output::AtmosphericTideDisplacementBasis basis, bool geodeticLookup = false ) const;
    void validateOriginConversion( const input_output::AtmosphericTideMetadata& metadata ) const;
    std::shared_ptr< AtmosphericTideGridInterpolator > gridInterpolator_;
    input_output::AtmosphericTideStationData stationData_;
    input_output::AtmosphericTideMetadata resolvedGridMetadata_;
    std::map< std::string, input_output::AtmosphericTideMetadata > resolvedStationMetadata_;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
    AtmosphericTideEvaluationSettings settings_;
    std::shared_ptr< BodyShapeModel > bodyShape_;
};
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
