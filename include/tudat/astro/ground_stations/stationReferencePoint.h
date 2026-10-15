/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_STATION_REFERENCE_POINT_H
#define TUDAT_STATION_REFERENCE_POINT_H
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace ground_stations
{

//! Identifier of a point of a station (for instance its marker or its antenna reference point).
struct StationReferencePointIdentity {
    //! Station identifier and DOMES number.
    input_output::SinexSiteKey site;
    std::string domesId;

    //! Kind of point, and its name.
    StationReferencePointKind kind = StationReferencePointKind::marker;
    std::string name;
};

//! Definition of the local vertical used for up/north/east offsets at a station.
enum class StationLocalVertical { geocentricRadial, geodeticEllipsoidNormal, explicitFixedBasis };

//! Local directions (east, north, up) at a station, used to convert up/north/east offsets to terrestrial coordinates.
struct StationReferencePointGeometry {
    //! Terrestrial reference frame, and the definition of the local vertical.
    std::string frame;
    StationLocalVertical vertical = StationLocalVertical::explicitFixedBasis;

    //! Rotation from east/north/up to terrestrial Cartesian axes, and the position (metres) at which it was computed.
    Eigen::Matrix3d enuToTerrestrial = Eigen::Matrix3d::Identity( );
    Eigen::Vector3d anchorMetres = Eigen::Vector3d::Zero( );

    //! Equatorial radius (metres) and flattening of the ellipsoid, for a geodetic vertical.
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
    //! Start and end (TAI) used for equipment records whose start or end refers to the missing header limits.
    std::optional< Time > startTai;
    std::optional< Time > endTai;

    //! Description of where these limits come from.
    std::string source;
};

//! Component order of a reference point offset: up, north, east (UNE) or Cartesian X, Y, Z.
enum class StationReferencePointOffsetBasis { une, xyz };

//! Offset from one point of a station to another (for instance from the marker to the antenna reference point), valid
//! in a given time interval.
struct StationReferencePointOffset {
    //! Points from which and to which the offset is defined, and the interval in which it is valid.
    StationReferencePointIdentity from, to;
    StationHistoryInterval validity;

    //! Terrestrial reference frame, component order and components (metres) of the offset.
    std::string frame;
    StationReferencePointOffsetBasis basis = StationReferencePointOffsetBasis::xyz;
    Eigen::Vector3d componentsMetres = Eigen::Vector3d::Zero( );

    //! Type and serial number of the equipment (for instance the antenna) to which the offset belongs.
    std::string equipmentType, equipmentSerialNumber;

    //! Calibration and signal to which the offset applies (required for an offset to a calibrated instrument point).
    std::string calibrationId, signal;

    //! Name of the source of the offset, the lines from which the offset and the equipment were read, and the time limits
    //! used when the source header gives none.
    std::string source;
    input_output::SinexSourceRecord provenance;
    input_output::SinexSourceRecord equipmentProvenance;
    std::optional< SinexEquipmentHeaderTimeContext > headerTimeContext;

    //! Optional time-dependent reference vectors in the declared basis; state rates are d/dTDB.
    std::function< Eigen::Vector3d( const StationHistoryTime& ) > positionFunction;
    std::function< Eigen::Vector6d( const StationHistoryTime& ) > stateFunction;

    //! Whether the offset depends on the direction of the signal; such offsets are not supported.
    bool directionDependent = false;
};

//! Offsets between two points of a station over time (for instance with a new offset after each antenna change).
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

    //! Offsets, ordered in time, and the local directions used for up/north/east offsets (if any).
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
