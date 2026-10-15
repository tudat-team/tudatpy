#include "expose_terrestrial_frames.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace terrestrial_frames
{
void expose_terrestrial_frames( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::enum_< io::HelmertRotationConvention >( module, "HelmertRotationConvention", R"doc(

Enumeration of sign conventions for the rotation angles of a Helmert transformation.

      )doc" )
            .value( "unspecified", io::HelmertRotationConvention::unspecified, "Convention not declared." )
            .value( "position_vector",
                    io::HelmertRotationConvention::positionVector,
                    "Rotation of the position vector (the convention of the IGN/ITRF tables)." )
            .value( "coordinate_frame", io::HelmertRotationConvention::coordinateFrame, "Rotation of the coordinate frame." );

    py::class_< io::TerrestrialFrameIdentity >( module, "TerrestrialFrameIdentity", R"doc(

Identification of a terrestrial reference frame realization (for instance ITRF2020).

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "realization", &io::TerrestrialFrameIdentity::realization, R"doc(

         Name of the realization (for instance ``ITRF2020``).

         :type: str
      )doc" )
            .def_readwrite( "release", &io::TerrestrialFrameIdentity::release, R"doc(

         Release of the catalogue in which the realization is given. A different release of the same realization does
         not denote a different frame.

         :type: str
      )doc" );

    py::class_< io::TerrestrialFrameTransformationSource >( module, "TerrestrialFrameTransformationSource", R"doc(

Description of the source of a table of terrestrial frame transformation parameters.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::TerrestrialFrameTransformationSource::source, R"doc(

         Name of the file (or other source) of the table.

         :type: str
      )doc" )
            .def_readwrite( "source_release", &io::TerrestrialFrameTransformationSource::sourceRelease, R"doc(

         Release of the table.

         :type: str
      )doc" )
            .def_readwrite( "source_sha256", &io::TerrestrialFrameTransformationSource::sourceSha256, R"doc(

         SHA-256 checksum of the source file.

         :type: str
      )doc" );

    py::class_< io::TerrestrialFrameTransformationParameters >( module, "TerrestrialFrameTransformationParameters", R"doc(

Parameters of a linearized (14-parameter) Helmert transformation between two terrestrial frame realizations.

The parameters and their rates are stored in SI units; the units of the source are retained. The transformation is
evaluated with :func:`~tudatpy.astro.frame_conversion.evaluate_linearized_helmert`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source_frame", &io::TerrestrialFrameTransformationParameters::sourceFrame, R"doc(

         Frame from which the parameters transform.

         :type: TerrestrialFrameIdentity
      )doc" )
            .def_readwrite( "target_frame", &io::TerrestrialFrameTransformationParameters::targetFrame, R"doc(

         Frame to which the parameters transform.

         :type: TerrestrialFrameIdentity
      )doc" )
            .def_readwrite( "reference_decimal_year", &io::TerrestrialFrameTransformationParameters::referenceDecimalYear, R"doc(

         Reference epoch of the parameters, as a decimal year.

         :type: float
      )doc" )
            .def_readwrite( "translation_metres", &io::TerrestrialFrameTransformationParameters::translationMetres, R"doc(

         Translation :math:`(T_{x}, T_{y}, T_{z})` in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "rotation_radians", &io::TerrestrialFrameTransformationParameters::rotationRadians, R"doc(

         Rotation angles :math:`(R_{x}, R_{y}, R_{z})` in radians.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "scale", &io::TerrestrialFrameTransformationParameters::scale, R"doc(

         Scale difference :math:`D` (dimensionless).

         :type: float
      )doc" )
            .def_readwrite( "translation_rate_metres_per_second",
                            &io::TerrestrialFrameTransformationParameters::translationRateMetresPerSecond,
                            R"doc(

         Rate of the translation in metres per second.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite(
                    "rotation_rate_radians_per_second", &io::TerrestrialFrameTransformationParameters::rotationRateRadiansPerSecond, R"doc(

         Rate of the rotation angles in radians per second.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "scale_rate_per_second", &io::TerrestrialFrameTransformationParameters::scaleRatePerSecond, R"doc(

         Rate of the scale difference per second.

         :type: float
      )doc" )
            .def_readwrite( "parameter_rate_year_seconds", &io::TerrestrialFrameTransformationParameters::parameterRateYearSeconds, R"doc(

         Length of the year in which the source gives its rates, in seconds (a Julian year of 365.25 days). Rates are converted to SI with this value; it does not define a time scale.

         :type: float
      )doc" )
            .def_readwrite( "rotation_convention", &io::TerrestrialFrameTransformationParameters::rotationConvention, R"doc(

         Sign convention of the rotation angles.

         :type: HelmertRotationConvention
      )doc" )
            .def_readwrite( "original_translation_units", &io::TerrestrialFrameTransformationParameters::originalTranslationUnits, R"doc(

         Units of the translation in the source.

         :type: str
      )doc" )
            .def_readwrite( "original_rotation_units", &io::TerrestrialFrameTransformationParameters::originalRotationUnits, R"doc(

         Units of the rotation angles in the source.

         :type: str
      )doc" )
            .def_readwrite( "original_scale_units", &io::TerrestrialFrameTransformationParameters::originalScaleUnits, R"doc(

         Units of the scale difference in the source.

         :type: str
      )doc" )
            .def_readwrite( "original_rate_time_units", &io::TerrestrialFrameTransformationParameters::originalRateTimeUnits, R"doc(

         Time unit of the rates in the source.

         :type: str
      )doc" )
            .def_readwrite( "provenance", &io::TerrestrialFrameTransformationParameters::provenance, R"doc(

         Source of the parameters.

         :type: TerrestrialFrameTransformationSource
      )doc" )
            .def_readwrite( "original_parameter_line", &io::TerrestrialFrameTransformationParameters::originalParameterLine, R"doc(

         Line of the source that gives the parameters, unmodified.

         :type: str
      )doc" )
            .def_readwrite( "original_rate_line", &io::TerrestrialFrameTransformationParameters::originalRateLine, R"doc(

         Line of the source that gives the rates, unmodified.

         :type: str
      )doc" )
            .def_readwrite( "parameter_source_line", &io::TerrestrialFrameTransformationParameters::parameterSourceLine, R"doc(

         Line number of the parameters in the source.

         :type: int
      )doc" )
            .def_readwrite( "rate_source_line", &io::TerrestrialFrameTransformationParameters::rateSourceLine, R"doc(

         Line number of the rates in the source.

         :type: int
      )doc" )
            .def( "validate", &io::validateTerrestrialFrameTransformationParameters, R"doc(

         Check that the parameters are finite and consistent, and raise an exception if not.

      )doc" );

    py::class_< io::TerrestrialFrameTransformationData >( module, "TerrestrialFrameTransformationData", R"doc(

Table of transformation parameters from a common reference frame realization to a set of other realizations.

Objects of this class are created by :func:`read_terrestrial_frame_transformation_data`,
:func:`official_itrf2014_parameters` and :func:`official_itrf2020_parameters`, and used by
:class:`~tudatpy.astro.frame_conversion.TerrestrialFrameTransformationModel`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "reference_frame", &io::TerrestrialFrameTransformationData::referenceFrame, R"doc(

         Common realization from which all parameters in the table transform.

         :type: TerrestrialFrameIdentity
      )doc" )
            .def_readwrite( "provenance", &io::TerrestrialFrameTransformationData::provenance, R"doc(

         Source of the table.

         :type: TerrestrialFrameTransformationSource
      )doc" )
            .def_readwrite( "source_headers", &io::TerrestrialFrameTransformationData::sourceHeaders, R"doc(

         Header lines and notices of the source, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "records", &io::TerrestrialFrameTransformationData::records, R"doc(

         Parameters for each target realization, with the name of the target realization as key.

         :type: dict[str, TerrestrialFrameTransformationParameters]
      )doc" )
            .def( "get_parameters",
                  &io::TerrestrialFrameTransformationData::getParameters,
                  py::arg( "target_realization" ),
                  py::return_value_policy::reference_internal,
                  R"doc(

         Function to retrieve the parameters for a single target realization.

         Parameters
         ----------
         target_realization : str
             Name of the target realization. An unknown name raises an exception.

         Returns
         -------
         TerrestrialFrameTransformationParameters
             Parameters of the transformation from the reference realization to the target realization.

      )doc" )
            .def( "validate", &io::validateTerrestrialFrameTransformationData, R"doc(

         Check that all records are valid and consistent with the reference realization, and raise an exception if not.

      )doc" );

    module.def( "read_terrestrial_frame_transformation_data",
                py::overload_cast< const std::string&, const io::TerrestrialFrameTransformationSource& >(
                        &io::readTerrestrialFrameTransformationData ),
                py::arg( "file_path" ),
                py::arg( "source" ) = io::TerrestrialFrameTransformationSource( ),
                R"doc(

Function to read a table of terrestrial frame transformation parameters.

Function to read a table of terrestrial frame transformation parameters in the format published by IGN for the ITRF
realizations (for instance ``Transfo-ITRF2020_TRFs.txt``). The units, the column order and the transformation equation
declared in the file are checked. Translations, rotations, scale and their rates are converted to SI units, using a
Julian year (365.25 days) for the rates. Malformed, non-finite, duplicated or truncated records raise an exception that
gives the file and line.


Parameters
----------
file_path : str
    Path of the file that is to be read.
source : TerrestrialFrameTransformationSource, default = TerrestrialFrameTransformationSource()
    Description of the source (release, checksum), which is stored with the parameters.

Returns
-------
TerrestrialFrameTransformationData
    Table read from the file.

      )doc" );
    module.def( "official_itrf2014_parameters", &io::readOfficialItrf2014TransformationData, R"doc(

Function to read the official IGN table of transformation parameters from ITRF2014 to earlier realizations.

The table is read from the ``station_locations`` directory of the Tudat resources.

Returns
-------
TerrestrialFrameTransformationData
    Table of transformation parameters with ITRF2014 as reference realization.

      )doc" );
    module.def( "official_itrf2020_parameters", &io::readOfficialItrf2020TransformationData, R"doc(

Function to read the official IGN table of transformation parameters from ITRF2020 to earlier realizations.

The table is read from the ``station_locations`` directory of the Tudat resources.

Returns
-------
TerrestrialFrameTransformationData
    Table of transformation parameters with ITRF2020 as reference realization.


Examples
--------
.. code-block:: python

   from tudatpy.data_input.environment_data import terrestrial_frames

   data = terrestrial_frames.official_itrf2020_parameters( )
   parameters = data.get_parameters( "ITRF2014" )
   print( parameters.reference_decimal_year )

      )doc" );
}
}  // namespace terrestrial_frames
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
