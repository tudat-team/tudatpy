/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_sp3.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/io/readSp3File.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace sp3
{

//! Expose the SP3 file reader and parsed-file container in the data-input SP3 submodule.
void expose_sp3( py::module& m )
{
    py::class_< tudat::input_output::Sp3FileContents, std::shared_ptr< tudat::input_output::Sp3FileContents > >(
            m, "Sp3FileContents", R"doc(
Parsed contents and metadata of an SP3 precise-orbit file.

Object that stores the contents of an SP3 file, as created by
:func:`~tudatpy.data_input.environment_data.sp3.read_sp3_file`.
)doc" )
            .def_property_readonly(
                    "format_version",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return std::string( 1, contents.formatVersion ); },
                    R"doc(
**read-only**

SP3 format version declared in the file header (``a``, ``b``, ``c`` or ``d``).

:type: str
)doc" )
            .def_property_readonly(
                    "has_velocity_records",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.hasVelocityRecords; },
                    R"doc(
**read-only**

Whether the file contains velocity records.

:type: bool
)doc" )
            .def_property_readonly(
                    "velocities_were_derived",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.velocitiesWereDerived; },
                    R"doc(
**read-only**

Whether the velocities in :attr:`satellite_states` were reconstructed from the positions, which is the case for a
position-only file.

:type: bool
)doc" )
            .def_property_readonly(
                    "start_epoch",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.startEpoch; },
                    R"doc(
**read-only**

Start epoch declared in the file header, in seconds since :attr:`reference_julian_day`, in the time scale
:attr:`time_scale`.

:type: float
)doc" )
            .def_property_readonly(
                    "reference_julian_day",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.referenceJulianDay; },
                    R"doc(
**read-only**

Julian day relative to which all epochs in this object are expressed.

:type: float
)doc" )
            .def_property_readonly(
                    "declared_number_of_epochs",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredNumberOfEpochs; },
                    R"doc(
**read-only**

Number of epochs declared in the file header.

:type: int
)doc" )
            .def_property_readonly(
                    "declared_number_of_satellites",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredNumberOfSatellites; },
                    R"doc(
**read-only**

Number of satellites declared in the file header.

:type: int
)doc" )
            .def_property_readonly(
                    "declared_epoch_interval",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredEpochInterval; },
                    R"doc(
**read-only**

Interval between epochs declared in the file header, in seconds.

:type: float
)doc" )
            .def_property_readonly(
                    "frame_name",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.frameName; },
                    R"doc(
**read-only**

Reference frame tag declared in the file header, as written in the file (for instance ``IGS14``). Empty if the file
declares no reference frame.

:type: str
)doc" )
            .def_property_readonly(
                    "analysis_center",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.analysisCenter; },
                    R"doc(
**read-only**

Tag of the producer or analysis center declared in the file header. Empty if the file declares none.

:type: str
)doc" )
            .def_property_readonly(
                    "time_scale",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.timeScale; },
                    R"doc(
**read-only**

Time-system tag declared in the file header (``GPS``, ``GLO``, ``GAL``, ``BDT``, ``TAI``, ``UTC``, ``IRN`` or ``QZS``). The
legacy ``ccc`` tag is stored as ``GPS``.

:type: str
)doc" )
            .def_property_readonly(
                    "satellite_states",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.satelliteStates; },
                    R"doc(
**read-only**

State history of each satellite, with the satellite ID (for instance ``G01``) as the first key, and the epoch as the
second key. Epochs are in seconds since :attr:`reference_julian_day`, in the time scale :attr:`time_scale`. Each state is a
Cartesian state, with position in meters and velocity in meters per second, in the frame :attr:`frame_name`. Missing values
are NaN.

:type: dict[str, dict[float, numpy.ndarray[numpy.float64[6, 1]]]]
)doc" )
            .def( "get_satellite_state_history",
                  &tudat::input_output::Sp3FileContents::getSatelliteStateHistory,
                  py::arg( "satellite_id" ),
                  R"doc(
Return the state history for one satellite.

States are returned in the reference frame and time system declared by the SP3
file. No frame or time conversion is performed by the data-input module. An error
is raised if the satellite is not in the file, or if one of its states contains a
missing value.

Parameters
----------
satellite_id : str
    ID of the satellite, as written in the file without surrounding whitespace (for instance ``G01``).

Returns
-------
dict[float, numpy.ndarray[numpy.float64[6, 1]]]
    State history of the satellite. Keys are epochs in seconds since :attr:`reference_julian_day`, in the time scale
    :attr:`time_scale`. Values are Cartesian states, with position in meters and velocity in meters per second, in the
    frame :attr:`frame_name`.
)doc" );

    m.def( "read_sp3_file",
           &tudat::input_output::readSp3File,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read an SP3 precise-orbit file.

Epochs are expressed in seconds since ``reference_julian_day`` in the time
scale declared by the file. Positions and velocities are converted to meters
and meters per second, respectively. For position-only files, velocities are
reconstructed with second-order three-point finite differences, which requires
at least three epochs. SP3-a through SP3-d and the GPS, GLO, GAL, BDT, TAI, UTC,
IRN, and QZS time-system tags are supported. Reference-frame and time-system tags
are retained as metadata; the loader does not transform either one. Clock and
covariance records are not read.

Parameters
----------
file_name : str
    Path to the SP3 file.
reference_julian_day : float, default=2451545.0
    Julian day relative to which the epochs are expressed. The default value is J2000.

Returns
-------
:class:`~tudatpy.data_input.environment_data.sp3.Sp3FileContents`
    Contents and metadata of the SP3 file.
)doc" );
}

}  // namespace sp3
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
