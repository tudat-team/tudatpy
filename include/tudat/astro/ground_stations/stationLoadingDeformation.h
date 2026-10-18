#ifndef TUDAT_STATION_LOADING_DEFORMATION_H
#define TUDAT_STATION_LOADING_DEFORMATION_H

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/io/readStationLoadingData.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Displacements of ground stations from a station loading series (for instance non-tidal atmospheric loading).
/*!
 * Each ground station is associated with a site of the series: the site given for its name in the station bindings or, if
 * no bindings are given, the nearest site within whose validity radius it lies. The radial, east and north displacements
 * of the site are interpolated linearly in TAI, and rotated to body-fixed Cartesian axes using the position of the site.
 * The rate is that of the linear interpolant (at a sample, that of the interval to its right). Epochs outside the series,
 * or in an interval longer than the maximum gap, raise an exception.
 */
class StationLoadingDeformation : public BodyDeformationModel
{
public:
    //! Constructor
    /*!
     * \param data Station loading series
     * \param stationBindings Site identifier of the series per ground station name (empty to select sites by their radius)
     * \param maximumGapSeconds Maximum interval between two samples over which the series is interpolated
     * \param converter Time scale converter (if nullptr, a converter without EOP data is used)
     */
    StationLoadingDeformation( const input_output::StationLoadingData& data,
                               const std::map< std::string, std::string >& stationBindings = {},
                               double maximumGapSeconds = 86400.0,
                               std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr );

protected:
    //! Cartesian displacement in metres.
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;

    //! Cartesian displacement and its rate (per TDB second).
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    //! Displacement and (if requested) its rate for a station.
    Eigen::Vector6d evaluate( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState >& station, bool withRate ) const;

    //! Station loading series, and the site identifier per ground station name.
    input_output::StationLoadingData data_;
    std::map< std::string, std::string > stationBindings_;

    //! Maximum interval between two samples over which the series is interpolated.
    double maximumGapSeconds_;

    //! Epochs of the series in TDB.
    std::vector< Time > tdbEpochs_;

    //! Rotation from the radial, east and north directions of each site to body-fixed Cartesian axes.
    std::map< std::string, Eigen::Matrix3d > siteBases_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
