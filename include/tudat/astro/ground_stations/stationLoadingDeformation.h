#ifndef TUDAT_STATION_LOADING_DEFORMATION_H
#define TUDAT_STATION_LOADING_DEFORMATION_H

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/io/readStationLoadingData.h"

namespace tudat
{
namespace basic_astrodynamics
{
enum class StationLoadingSpatialSelection { explicitBindings, withinProviderRadius };

//! Cartesian correction r_sourceOrigin-r_targetOrigin; its inventory declares its physical component/band scope.
struct StationLoadingOriginTranslation {
    std::string identifier;
    std::string sourceOriginDefinition;
    std::string targetOriginDefinition;
    input_output::StationLoadingOrigin targetOrigin = input_output::StationLoadingOrigin::unspecified;
    input_output::StationLoadingData series;
};

struct StationLoadingEvaluationSettings {
    StationLoadingSpatialSelection spatialSelection = StationLoadingSpatialSelection::explicitBindings;

    //! Environment station ID -> arbitrary provider label. The empty ID explicitly binds an unidentified nominal station.
    std::map< std::string, std::string > stationBindings;

    //! Linear interpolation only; a gap larger than this is unavailable, even if present in an irregular source grid.
    double maximumGapSeconds = 86400.0;
    std::string outputOriginDefinition;
    std::shared_ptr< const StationLoadingOriginTranslation > originTranslation;
};

//! Precomputed loading, evaluated on continuous TAI and returned in provider Cartesian axes.
//! Native input is precise TDB Time. Rates differentiate the same linear interpolant,
//! taking the interval to the right of a knot, and to the left of the final endpoint.
class StationLoadingDeformation : public BodyDeformationModel
{
public:
    StationLoadingDeformation( const input_output::StationLoadingData& data,
                               const StationLoadingEvaluationSettings& settings,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr );
    const input_output::StationLoadingData& getData( ) const
    {
        return data_;
    }
    const StationLoadingEvaluationSettings& getEvaluationSettings( ) const
    {
        return settings_;
    }

protected:
    //! Return Cartesian displacement in metres, including valid singleton source positions.
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;

    //! Return displacement and the derivative of its linear TAI interpolant per TDB second.
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    struct PreparedGrid {
        std::vector< Time > tai, tdb;
    };
    PreparedGrid prepareGrid( const input_output::StationLoadingData& data ) const;
    const input_output::StationLoadingSite& selectSite( const Eigen::Vector3d& position, const std::string& stationId ) const;

    //! Select one provider point and calculate displacement plus optional rate, applying its basis/origin once.
    Eigen::Vector6d evaluate( const Time& tdb, const Eigen::Vector3d& position, const std::string& stationId, bool withRate ) const;

    //! Evaluate one physical TAI interval; the optional velocity is returned per TAI second.
    Eigen::Vector6d interpolate( const input_output::StationLoadingSite& site,
                                 const PreparedGrid& grid,
                                 const Time& tdb,
                                 bool withRate ) const;

    input_output::StationLoadingData data_;
    StationLoadingEvaluationSettings settings_;
    PreparedGrid grid_, translationGrid_;
    std::map< std::string, Eigen::Matrix3d > sourceBases_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
