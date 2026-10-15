#include "expose_earth_orientation.h"
#include "tudat/io/readEarthOrientationData.h"
#include <pybind11/stl.h>
namespace py = pybind11;
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace earth_orientation
{
void expose_earth_orientation( py::module& module )
{
    using namespace tudat::input_output;
    py::enum_< EarthOrientationFormat >( module, "EarthOrientationFormat", R"doc(

Enumeration of supported Earth orientation parameter (EOP) file formats.

The ``c04`` value covers both the legacy 16-column and the modern 21-column IERS C04 formats (the column count is
detected automatically). The ``finals2000a`` value denotes the fixed-width IERS ``finals2000A`` format, of which only
the Bulletin A values are used.

      )doc" )
            .value( "c04", EarthOrientationFormat::c04, "IERS C04 series (legacy 16-column or modern 21-column)." )
            .value( "finals2000a", EarthOrientationFormat::finals2000a, "IERS finals2000A (Bulletin A values)." );

    py::class_< EarthOrientationRecord >( module, "EarthOrientationRecord", R"doc(

Earth orientation parameters at a single UTC epoch, as read from file.

All quantities are converted to SI units: angles in radians, angle rates in radians per second, and UT1-UTC and
excess length of day in seconds. Optional quantities that are not given in the file are ``None``; a zero in the file is
retained as a measured value.

      )doc" )
            .def_readonly( "utc_modified_julian_date", &EarthOrientationRecord::utcModifiedJulianDate, R"doc(

         **read-only**

         UTC epoch of the record, as a modified Julian date (including the hour field of the modern C04 format).

         :type: float
      )doc" )
            .def_readonly( "polar_motion", &EarthOrientationRecord::polarMotion, R"doc(

         **read-only**

         Polar motion coordinates :math:`(x_{p}, y_{p})` in radians.

         :type: list[float]
      )doc" )
            .def_readonly( "ut1_minus_utc", &EarthOrientationRecord::ut1MinusUtc, R"doc(

         **read-only**

         UT1-UTC in seconds.

         :type: float
      )doc" )
            .def_readonly( "length_of_day", &EarthOrientationRecord::lengthOfDay, R"doc(

         **read-only**

         Excess length of day in seconds, or ``None`` if not given.

         :type: float | None
      )doc" )
            .def_readonly( "celestial_pole_offsets", &EarthOrientationRecord::celestialPoleOffsets, R"doc(

         **read-only**

         Celestial pole offsets :math:`(dX, dY)` in radians, or ``None`` if not given. Files that give the older
         :math:`(d\psi, d\epsilon)` offsets are rejected.

         :type: list[float] | None
      )doc" )
            .def_readonly( "polar_motion_rates", &EarthOrientationRecord::polarMotionRates, R"doc(

         **read-only**

         Polar motion rates in radians per second (modern C04 only), or ``None`` if not given.

         :type: list[float] | None
      )doc" )
            .def_readonly( "polar_motion_uncertainty", &EarthOrientationRecord::polarMotionUncertainty, R"doc(

         **read-only**

         Uncertainty of the polar motion coordinates in radians, or ``None`` if not given.

         :type: list[float] | None
      )doc" )
            .def_readonly( "ut1_minus_utc_uncertainty", &EarthOrientationRecord::ut1MinusUtcUncertainty, R"doc(

         **read-only**

         Uncertainty of UT1-UTC in seconds, or ``None`` if not given.

         :type: float | None
      )doc" )
            .def_readonly( "length_of_day_uncertainty", &EarthOrientationRecord::lengthOfDayUncertainty, R"doc(

         **read-only**

         Uncertainty of the excess length of day in seconds, or ``None`` if not given.

         :type: float | None
      )doc" )
            .def_readonly( "celestial_pole_offset_uncertainty", &EarthOrientationRecord::celestialPoleOffsetUncertainty, R"doc(

         **read-only**

         Uncertainty of the celestial pole offsets in radians, or ``None`` if not given.

         :type: list[float] | None
      )doc" )
            .def_readonly( "polar_motion_rate_uncertainty", &EarthOrientationRecord::polarMotionRateUncertainty, R"doc(

         **read-only**

         Uncertainty of the polar motion rates in radians per second, or ``None`` if not given.

         :type: list[float] | None
      )doc" )
            .def_readonly( "polar_motion_status", &EarthOrientationRecord::polarMotionStatus, R"doc(

         **read-only**

         Bulletin A flag of the polar motion values (``I`` for IERS/interpolated, ``P`` for predicted). Empty for C04.

         :type: str
      )doc" )
            .def_readonly( "ut1_status", &EarthOrientationRecord::ut1Status, R"doc(

         **read-only**

         Bulletin A flag of the UT1-UTC value (``I`` or ``P``). Empty for C04.

         :type: str
      )doc" )
            .def_readonly( "celestial_pole_status", &EarthOrientationRecord::celestialPoleStatus, R"doc(

         **read-only**

         Bulletin A flag of the celestial pole offsets (``I`` or ``P``). Empty for C04.

         :type: str
      )doc" );

    py::class_< EarthOrientationData >( module, "EarthOrientationData", R"doc(

Earth orientation parameters read from a file, with their source information.

The records are stored as read: no interpolation, extrapolation or padding is applied, and the coverage is exactly
that of the file. Objects of this class are created by :func:`read_earth_orientation_data`, and can be supplied to the
Earth rotation model settings (see :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`).

      )doc" )
            .def_readonly( "format", &EarthOrientationData::format, R"doc(

         **read-only**

         Format of the file from which the data was read.

         :type: EarthOrientationFormat
      )doc" )
            .def_readonly( "source", &EarthOrientationData::source, R"doc(

         **read-only**

         Name of the file (or other source) from which the data was read.

         :type: str
      )doc" )
            .def_readonly( "format_description", &EarthOrientationData::formatDescription, R"doc(

         **read-only**

         Description of the detected file format (for instance, the C04 column layout).

         :type: str
      )doc" )
            .def_readonly( "header", &EarthOrientationData::header, R"doc(

         **read-only**

         Header lines of the file, unmodified.

         :type: list[str]
      )doc" )
            .def_readonly( "terrestrial_frame", &EarthOrientationData::terrestrialFrame, R"doc(

         **read-only**

         Terrestrial reference frame of the data, taken from the C04 header when given (default ``ITRS``). A frame
         realization is never inferred from the file name.

         :type: str
      )doc" )
            .def_readonly( "celestial_frame", &EarthOrientationData::celestialFrame, R"doc(

         **read-only**

         Celestial reference frame of the data (``GCRS``).

         :type: str
      )doc" )
            .def_readonly( "celestial_pole_convention", &EarthOrientationData::celestialPoleConvention, R"doc(

         **read-only**

         Convention of the celestial pole offsets (``dX,dY``).

         :type: str
      )doc" )
            .def_readonly( "records", &EarthOrientationData::records, R"doc(

         **read-only**

         Records, with the UTC modified Julian date of each record as key.

         :type: dict[float, EarthOrientationRecord]
      )doc" );

    module.def( "read_earth_orientation_data",
                py::overload_cast< const std::string&, EarthOrientationFormat >( &readEarthOrientationData ),
                py::arg( "file_name" ),
                py::arg( "format" ) = EarthOrientationFormat::c04,
                R"doc(

Function to read Earth orientation parameters from a file.

Function to read Earth orientation parameters (polar motion, UT1-UTC, length of day, celestial pole offsets and, where
available, their rates and uncertainties) from an IERS C04 or ``finals2000A`` file. For C04 files, the legacy 16-column and
the modern 21-column formats are detected automatically. For ``finals2000A`` files, only the Bulletin A values are used:
Bulletin B columns are never substituted, and future dates without values do not extend the coverage.

The data are stored without interpolation or extrapolation. Malformed files (inconsistent calendar dates and MJDs,
duplicate or unsorted epochs, incomplete coordinate pairs, truncated records, non-finite values or negative
uncertainties) raise an exception that gives the file and line number.


Parameters
----------
file_name : str
    Path of the file that is to be read.
format : EarthOrientationFormat, default = EarthOrientationFormat.c04
    Format of the file.

Returns
-------
EarthOrientationData
    Records read from the file, with the file's source information.


Examples
--------
.. code-block:: python

   from tudatpy.data_input.environment_data import earth_orientation

   data = earth_orientation.read_earth_orientation_data(
       "eopc04.txt", earth_orientation.EarthOrientationFormat.c04 )
   first_record = data.records[ min( data.records ) ]
   print( first_record.polar_motion, first_record.ut1_minus_utc )

      )doc" );
}
}  // namespace earth_orientation
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
