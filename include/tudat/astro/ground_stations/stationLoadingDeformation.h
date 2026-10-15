#ifndef TUDAT_STATION_LOADING_DEFORMATION_H
#define TUDAT_STATION_LOADING_DEFORMATION_H

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/io/readStationLoadingData.h"

namespace tudat
{
namespace basic_astrodynamics
{
//! How a ground station is matched to a site of a station loading series: by an explicit list of pairs, or by finding the
//! site within whose validity radius the station lies.
enum class StationLoadingSpatialSelection { explicitBindings, withinProviderRadius };

//! Time series of the translation between two origins (r_sourceOrigin - r_targetOrigin), used to convert station loading
//! displacements to another origin.
struct StationLoadingOriginTranslation {
    //! Name of the translation series.
    std::string identifier;

    //! Descriptions of the origin of the displacements and of the origin to which they are converted, and the latter origin.
    std::string sourceOriginDefinition;
    std::string targetOriginDefinition;
    input_output::StationLoadingOrigin targetOrigin = input_output::StationLoadingOrigin::unspecified;

    //! Translation, as a station loading series with a single site.
    input_output::StationLoadingData series;
};

//! Settings for computing station displacements from a station loading series.
struct StationLoadingEvaluationSettings {
    //! How ground stations are matched to sites of the series.
    StationLoadingSpatialSelection spatialSelection = StationLoadingSpatialSelection::explicitBindings;

    //! Environment station ID -> arbitrary provider label. The empty ID explicitly binds an unidentified nominal station.
    std::map< std::string, std::string > stationBindings;

    //! Linear interpolation only; a gap larger than this is unavailable, even if present in an irregular source grid.
    double maximumGapSeconds = 86400.0;

    //! Description of the origin of the output displacements, and the translation used to convert them to that origin
    //! (nullptr if the origin of the series is used).
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
    //! Epochs of a station loading series, converted to TAI and TDB.
    struct PreparedGrid {
        //! TAI and TDB of each epoch.
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

    //! Station loading series, and the settings for evaluating it.
    input_output::StationLoadingData data_;
    StationLoadingEvaluationSettings settings_;

    //! Epochs of the series and of the origin translation, in TAI and TDB.
    PreparedGrid grid_, translationGrid_;

    //! Rotation from the component directions of each site to body-fixed Cartesian axes, with the site identifier as key.
    std::map< std::string, Eigen::Matrix3d > sourceBases_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
