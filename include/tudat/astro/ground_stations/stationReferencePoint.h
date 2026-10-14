/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_REFERENCE_POINT_H
#define TUDAT_STATION_REFERENCE_POINT_H
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace ground_stations
{
struct StationReferencePointIdentity {
    input_output::SinexSiteKey site;
    std::string domesId;
    StationReferencePointKind kind = StationReferencePointKind::marker;
    std::string name;
};
enum class StationLocalVertical { geocentricRadial, geodeticEllipsoidNormal, explicitFixedBasis };
struct StationReferencePointGeometry {
    std::string frame;
    StationLocalVertical vertical = StationLocalVertical::explicitFixedBasis;
    Eigen::Matrix3d enuToTerrestrial = Eigen::Matrix3d::Identity( );
    Eigen::Vector3d anchorMetres = Eigen::Vector3d::Zero( );
    std::optional< double > equatorialRadiusMetres;
    std::optional< double > flattening;
};
StationReferencePointGeometry createGeocentricStationReferencePointGeometry( const Eigen::Vector3d& anchor, const std::string& frame );
StationReferencePointGeometry createGeodeticStationReferencePointGeometry( const Eigen::Vector3d& anchor,
                                                                           const std::string& frame,
                                                                           double equatorialRadiusMetres,
                                                                           double flattening );
//! Explicit bounded context for missing source header limits, with its own provenance.
struct SinexEquipmentHeaderTimeContext {
    std::optional< Time > startTai;
    std::optional< Time > endTai;
    std::string source;
};
enum class StationReferencePointOffsetBasis { une, xyz };
enum class StationReferencePointRatePolicy { rejectAtSwitch, rightHandAtSwitch };
struct StationReferencePointOffset {
    StationReferencePointIdentity from, to;
    StationHistoryInterval validity;
    std::string frame;
    StationReferencePointOffsetBasis basis = StationReferencePointOffsetBasis::xyz;
    Eigen::Vector3d componentsMetres = Eigen::Vector3d::Zero( );
    std::string equipmentType, equipmentSerialNumber;
    std::string calibrationId, signal;
    std::string source;
    input_output::SinexSourceRecord provenance;
    input_output::SinexSourceRecord equipmentProvenance;
    std::optional< SinexEquipmentHeaderTimeContext > headerTimeContext;
    //! Optional time-dependent reference vectors in the declared basis; state rates are d/dTDB.
    std::function< Eigen::Vector3d( const StationHistoryTime& ) > positionFunction;
    std::function< Eigen::Vector6d( const StationHistoryTime& ) > stateFunction;
    bool directionDependent = false;
};
class StationReferencePointOffsetHistory
{
public:
    StationReferencePointOffsetHistory( std::vector< StationReferencePointOffset > offsets,
                                        std::optional< StationReferencePointGeometry > geometry = std::nullopt,
                                        StationReferencePointRatePolicy ratePolicy = StationReferencePointRatePolicy::rejectAtSwitch );
    Eigen::Vector3d getPositionCorrection( const StationHistoryTime& time ) const;
    Eigen::Vector6d getStateCorrection( const StationHistoryTime& time ) const;
    const StationReferencePointOffset& getActiveOffset( const Time& tai ) const;
    const std::vector< StationReferencePointOffset >& getOffsets( ) const
    {
        return offsets_;
    }
    const std::optional< StationReferencePointGeometry >& getGeometry( ) const
    {
        return geometry_;
    }
    StationReferencePointRatePolicy getRatePolicy( ) const
    {
        return ratePolicy_;
    }

private:
    Eigen::Vector3d toTerrestrial( const Eigen::Vector3d& components, StationReferencePointOffsetBasis basis ) const;
    std::vector< StationReferencePointOffset > offsets_;
    std::optional< StationReferencePointGeometry > geometry_;
    StationReferencePointRatePolicy ratePolicy_;
};

//! Compose one declared reference-point vector once with its source trajectory.
class StationReferencePointTrajectory : public StationHistoryTrajectory
{
public:
    StationReferencePointTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                     std::shared_ptr< StationReferencePointOffsetHistory > offsets );
    Eigen::Vector3d getPosition( const StationHistoryTime& time ) const override;
    Eigen::Vector6d getState( const StationHistoryTime& time ) const override;
    const StationHistoryMetadata& getMetadata( ) const override
    {
        return outputMetadata_;
    }
    bool hasFittedSeasonalTerms( ) const override
    {
        return source_->hasFittedSeasonalTerms( );
    }
    const StationReferencePointIdentity& getReferencePointIdentity( ) const
    {
        return offsets_->getOffsets( ).front( ).to;
    }
    const std::shared_ptr< StationHistoryTrajectory >& getSourceTrajectory( ) const
    {
        return source_;
    }
    const std::shared_ptr< StationReferencePointOffsetHistory >& getOffsetHistory( ) const
    {
        return offsets_;
    }

private:
    std::shared_ptr< StationHistoryTrajectory > source_;
    std::shared_ptr< StationReferencePointOffsetHistory > offsets_;
    StationHistoryMetadata outputMetadata_;
};
StationReferencePointIdentity getStationReferencePointIdentity( const std::shared_ptr< StationHistoryTrajectory >& trajectory );
std::shared_ptr< StationReferencePointTrajectory > createStationReferencePointTrajectory(
        std::shared_ptr< StationHistoryTrajectory > source,
        std::shared_ptr< StationReferencePointOffsetHistory > offsets );
//! GNSS marker-to-ARP interpretation of SINEX equipment/eccentricity records, never a phase centre.
std::shared_ptr< StationReferencePointTrajectory > createSinexAntennaReferencePointTrajectory(
        std::shared_ptr< StationHistoryTrajectory > marker,
        const input_output::SinexStationHistory& equipment,
        std::string referencePointName,
        std::optional< StationReferencePointGeometry > geometry = std::nullopt,
        std::optional< SinexEquipmentHeaderTimeContext > headerContext = std::nullopt,
        StationReferencePointRatePolicy ratePolicy = StationReferencePointRatePolicy::rejectAtSwitch,
        const std::optional< std::string >& domesSelector = std::nullopt,
        const std::optional< input_output::SinexSiteKey >& equipmentSite = std::nullopt );
}  // namespace ground_stations
}  // namespace tudat
#endif
