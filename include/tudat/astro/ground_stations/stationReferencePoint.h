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

//! At an offset switch, position and state use the new (right-hand) record, as elsewhere at station-history discontinuities.
class StationReferencePointOffsetHistory
{
public:
    StationReferencePointOffsetHistory( std::vector< StationReferencePointOffset > offsets,
                                        std::optional< StationReferencePointGeometry > geometry = std::nullopt );
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

private:
    Eigen::Vector3d toTerrestrial( const Eigen::Vector3d& components, StationReferencePointOffsetBasis basis ) const;
    std::vector< StationReferencePointOffset > offsets_;
    std::optional< StationReferencePointGeometry > geometry_;
};

//! Compose one declared reference-point vector once with its source trajectory.
class StationReferencePointTrajectory : public StationHistoryTrajectory
{
public:
    //! Retain the source and bind one validated reference-point offset history.
    StationReferencePointTrajectory( std::shared_ptr< StationHistoryTrajectory > source,
                                     std::shared_ptr< StationReferencePointOffsetHistory > offsets );

    //! Add the active reference vector to the source position.
    Eigen::Vector3d getPosition( const StationHistoryTime& time ) const override;

    //! Add the active reference vector and its TDB velocity to the source state.
    Eigen::Vector6d getState( const StationHistoryTime& time ) const override;

    //! Describe the evaluated reference point while retaining source frame and provenance.
    const StationHistoryMetadata& getMetadata( ) const override
    {
        return outputMetadata_;
    }

    //! Return the directed offset's output point identity.
    const StationReferencePointIdentity& getReferencePointIdentity( ) const
    {
        return offsets_->getOffsets( ).front( ).to;
    }

    //! Inspect the source retained by this ordered reference-point composition.
    const std::shared_ptr< StationHistoryTrajectory >& getSourceTrajectory( ) const
    {
        return source_;
    }

    //! Inspect the validated offsets and their equipment association.
    const std::shared_ptr< StationReferencePointOffsetHistory >& getOffsetHistory( ) const
    {
        return offsets_;
    }

private:
    //! Immutable offset history and metadata describing the evaluated output point.
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
        const std::optional< std::string >& domesSelector = std::nullopt,
        const std::optional< input_output::SinexSiteKey >& equipmentSite = std::nullopt );
}  // namespace ground_stations
}  // namespace tudat
#endif
