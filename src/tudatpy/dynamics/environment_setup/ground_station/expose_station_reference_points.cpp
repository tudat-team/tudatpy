/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_ground_station.h"
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "tudat/simulation/environment_setup/createStationReferencePoint.h"
namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace ground_station
{
namespace gs = tudat::ground_stations;
void expose_station_reference_points( py::module& m )
{
    py::class_< gs::StationReferencePointIdentity >( m, "StationReferencePointIdentity" )
            .def( py::init<>( ) )
            .def_readwrite( "site", &gs::StationReferencePointIdentity::site )
            .def_readwrite( "domes_id", &gs::StationReferencePointIdentity::domesId )
            .def_readwrite( "kind", &gs::StationReferencePointIdentity::kind )
            .def_readwrite( "name", &gs::StationReferencePointIdentity::name );
    py::enum_< gs::StationLocalVertical >( m, "StationLocalVertical" )
            .value( "geocentric_radial", gs::StationLocalVertical::geocentricRadial )
            .value( "geodetic_ellipsoid_normal", gs::StationLocalVertical::geodeticEllipsoidNormal )
            .value( "explicit_fixed_basis", gs::StationLocalVertical::explicitFixedBasis );
    py::class_< gs::StationReferencePointGeometry >( m, "StationReferencePointGeometry" )
            .def( py::init<>( ) )
            .def_readwrite( "frame", &gs::StationReferencePointGeometry::frame )
            .def_readwrite( "vertical", &gs::StationReferencePointGeometry::vertical )
            .def_readwrite( "enu_to_terrestrial", &gs::StationReferencePointGeometry::enuToTerrestrial )
            .def_readwrite( "anchor_metres", &gs::StationReferencePointGeometry::anchorMetres )
            .def_readwrite( "equatorial_radius_metres", &gs::StationReferencePointGeometry::equatorialRadiusMetres )
            .def_readwrite( "flattening", &gs::StationReferencePointGeometry::flattening );
    py::class_< gs::SinexEquipmentHeaderTimeContext >( m, "SinexEquipmentHeaderTimeContext" )
            .def( py::init<>( ) )
            .def_readwrite( "start_tai", &gs::SinexEquipmentHeaderTimeContext::startTai )
            .def_readwrite( "end_tai", &gs::SinexEquipmentHeaderTimeContext::endTai )
            .def_readwrite( "source", &gs::SinexEquipmentHeaderTimeContext::source );
    py::enum_< gs::StationReferencePointOffsetBasis >( m, "StationReferencePointOffsetBasis" )
            .value( "une", gs::StationReferencePointOffsetBasis::une )
            .value( "xyz", gs::StationReferencePointOffsetBasis::xyz );
    py::enum_< gs::StationReferencePointRatePolicy >( m, "StationReferencePointRatePolicy" )
            .value( "reject_at_switch", gs::StationReferencePointRatePolicy::rejectAtSwitch )
            .value( "right_hand_at_switch", gs::StationReferencePointRatePolicy::rightHandAtSwitch );
    py::class_< gs::StationReferencePointOffset >( m, "StationReferencePointOffset" )
            .def( py::init<>( ) )
            .def_readwrite( "from_point", &gs::StationReferencePointOffset::from )
            .def_readwrite( "to_point", &gs::StationReferencePointOffset::to )
            .def_readwrite( "validity", &gs::StationReferencePointOffset::validity )
            .def_readwrite( "frame", &gs::StationReferencePointOffset::frame )
            .def_readwrite( "basis", &gs::StationReferencePointOffset::basis )
            .def_readwrite( "components_metres", &gs::StationReferencePointOffset::componentsMetres )
            .def_readwrite( "equipment_type", &gs::StationReferencePointOffset::equipmentType )
            .def_readwrite( "equipment_serial_number", &gs::StationReferencePointOffset::equipmentSerialNumber )
            .def_readwrite( "calibration_id", &gs::StationReferencePointOffset::calibrationId )
            .def_readwrite( "signal", &gs::StationReferencePointOffset::signal )
            .def_readwrite( "source", &gs::StationReferencePointOffset::source )
            .def_readwrite( "position_function", &gs::StationReferencePointOffset::positionFunction )
            .def_readwrite( "state_function", &gs::StationReferencePointOffset::stateFunction )
            .def_readwrite( "direction_dependent", &gs::StationReferencePointOffset::directionDependent )
            .def_readonly( "provenance", &gs::StationReferencePointOffset::provenance )
            .def_readonly( "equipment_provenance", &gs::StationReferencePointOffset::equipmentProvenance )
            .def_readonly( "header_time_context", &gs::StationReferencePointOffset::headerTimeContext );
    py::class_< gs::StationReferencePointOffsetHistory, std::shared_ptr< gs::StationReferencePointOffsetHistory > >(
            m, "StationReferencePointOffsetHistory" )
            .def( py::init< std::vector< gs::StationReferencePointOffset >,
                            std::optional< gs::StationReferencePointGeometry >,
                            gs::StationReferencePointRatePolicy >( ),
                  py::arg( "offsets" ),
                  py::arg( "geometry" ) = py::none( ),
                  py::arg( "rate_policy" ) = gs::StationReferencePointRatePolicy::rejectAtSwitch )
            .def( "position_correction", &gs::StationReferencePointOffsetHistory::getPositionCorrection, py::arg( "source_time" ) )
            .def( "state_correction", &gs::StationReferencePointOffsetHistory::getStateCorrection, py::arg( "source_time" ) )
            .def_property_readonly( "offsets", &gs::StationReferencePointOffsetHistory::getOffsets, py::return_value_policy::copy )
            .def_property_readonly( "geometry", []( const gs::StationReferencePointOffsetHistory& h ) { return h.getGeometry( ); } )
            .def_property_readonly( "rate_policy", &gs::StationReferencePointOffsetHistory::getRatePolicy );
    py::class_< gs::StationReferencePointTrajectory, std::shared_ptr< gs::StationReferencePointTrajectory >, gs::StationHistoryTrajectory >(
            m, "StationReferencePointTrajectory" )
            .def_property_readonly( "reference_point_identity",
                                    []( const gs::StationReferencePointTrajectory& p ) { return p.getReferencePointIdentity( ); } )
            .def_property_readonly( "source_trajectory", &gs::StationReferencePointTrajectory::getSourceTrajectory )
            .def_property_readonly( "offset_history", &gs::StationReferencePointTrajectory::getOffsetHistory );
    m.def( "geocentric_reference_point_geometry",
           &gs::createGeocentricStationReferencePointGeometry,
           py::arg( "anchor_metres" ),
           py::arg( "frame" ) );
    m.def( "geodetic_reference_point_geometry",
           &gs::createGeodeticStationReferencePointGeometry,
           py::arg( "anchor_metres" ),
           py::arg( "frame" ),
           py::arg( "equatorial_radius_metres" ),
           py::arg( "flattening" ) );
    m.def( "station_reference_point_identity", &gs::getStationReferencePointIdentity, py::arg( "trajectory" ) );
    m.def( "station_reference_point_trajectory",
           &gs::createStationReferencePointTrajectory,
           py::arg( "source" ),
           py::arg( "offset_history" ) );
    m.def( "sinex_antenna_reference_point",
           &gs::createSinexAntennaReferencePointTrajectory,
           py::arg( "marker" ),
           py::arg( "equipment" ),
           py::arg( "reference_point_name" ),
           py::arg( "geometry" ) = py::none( ),
           py::arg( "header_context" ) = py::none( ),
           py::arg( "rate_policy" ) = gs::StationReferencePointRatePolicy::rejectAtSwitch,
           py::arg( "domes_selector" ) = py::none( ),
           py::arg( "equipment_site" ) = py::none( ) );
    m.def( "reference_point_station",
           &tudat::simulation_setup::stationReferencePointSettings,
           py::arg( "trajectory" ),
           py::arg( "nominal_position" ),
           py::arg( "source_time_function" ) = py::none( ),
           "Create a named station for the explicitly selected instrument point, including complete absolute trajectory motion." );
}
}  // namespace ground_station
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
