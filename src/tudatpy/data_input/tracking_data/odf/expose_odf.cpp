#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_odf.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>

#include "scalarTypes.h"
#include "tudat/io/preProcessOdfFile.h"
#include "tudat/io/readOdfFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace odf
{

void expose_odf( py::module& m )
{
    py::enum_< tio::OdfDataType >( m, "OdfDataType", R"doc(Possible data types in orbit section of ODF file)doc" )
            .value( "narrowband_spacecraft_vlbi_doppler_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_doppler_mode )
            .value( "narrowband_spacecraft_vlbi_phase_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_phase_mode )
            .value( "narrowband_quasar_vlbi_doppler_mode", tio::OdfDataType::narrowband_quasar_vlbi_doppler_mode )
            .value( "narrowband_quasar_vlbi_phase_mode", tio::OdfDataType::narrowband_quasar_vlbi_phase_mode )
            .value( "wideband_spacecraft_vlbi", tio::OdfDataType::wideband_spacecraft_vlbi )
            .value( "wideband_quasar_vlbi", tio::OdfDataType::wideband_quasar_vlbi )
            .value( "one_way_doppler", tio::OdfDataType::one_way_doppler )
            .value( "two_way_doppler", tio::OdfDataType::two_way_doppler )
            .value( "three_way_doppler", tio::OdfDataType::three_way_doppler )
            .value( "one_way_total_count_phase", tio::OdfDataType::one_way_total_count_phase )
            .value( "two_way_total_count_phase", tio::OdfDataType::two_way_total_count_phase )
            .value( "three_way_total_count_phase", tio::OdfDataType::three_way_total_count_phase )
            .value( "pra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::pra_planetary_operational_discrete_spectrum_range )
            .value( "sra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::sra_planetary_operational_discrete_spectrum_range )
            .value( "re_range", tio::OdfDataType::re_range )
            .value( "azimuth_angle", tio::OdfDataType::azimuth_angle )
            .value( "elevation_angle", tio::OdfDataType::elevation_angle )
            .value( "hour_angle", tio::OdfDataType::hour_angle )
            .value( "declination_angle", tio::OdfDataType::declination_angle )
            .value( "x_angle_east", tio::OdfDataType::x_angle_east )
            .value( "y_angle_east", tio::OdfDataType::y_angle_east )
            .value( "x_angle_south", tio::OdfDataType::x_angle_south )
            .value( "y_angle_south", tio::OdfDataType::y_angle_south );

    py::class_< tio::OdfCommonDataBlock, std::shared_ptr< tio::OdfCommonDataBlock > >( m, "OdfCommonDataBlock", R"doc(
         Base class observable-independent ODF data containers

                 The data section of an ODF is split into blocks (lines), each associated with
                 an observation epoch. The first elements of each block (e.g. observation epoch,
                 value of the observable) are common for all the observable types, but the
                 values in the remaining columns will have different meanings for different
                 types of observations. This class serves as interface to the observable-independent
                 part of an ODF data block. The different classes inheriting from OdfDataSpecificBlock
                 provide interfaces to the observable-specific part of the blocks.
      )doc" )
            .def_property_readonly( "observable_time", &tio::OdfCommonDataBlock::getObservableTime )
            .def_property_readonly( "observable_value", &tio::OdfCommonDataBlock::getObservableValue )
            .def_property_readonly( "receiving_station_downlink_delay", &tio::OdfCommonDataBlock::getReceivingStationDownlinkDelay )
            .def_readonly( "format_id", &tio::OdfCommonDataBlock::formatId_ )
            .def_readonly( "receiving_station_id", &tio::OdfCommonDataBlock::receivingStationId_ )
            .def_readonly( "transmitting_station_id", &tio::OdfCommonDataBlock::transmittingStationId_ )
            .def_readonly( "transmitting_station_network_id", &tio::OdfCommonDataBlock::transmittingStationNetworkId_ )
            .def_readonly( "data_type", &tio::OdfCommonDataBlock::dataType_ )
            .def_readonly( "downlink_band_id", &tio::OdfCommonDataBlock::downlinkBandId_ )
            .def_readonly( "uplink_band_id", &tio::OdfCommonDataBlock::uplinkBandId_ )
            .def_readonly( "reference_band_id", &tio::OdfCommonDataBlock::referenceBandId_ )
            .def_readonly( "is_invalid", &tio::OdfCommonDataBlock::validity_ )
            .def(
                    "print_data_block",
                    []( tio::OdfCommonDataBlock& block, const std::string& output_file ) {
                        std::ofstream output_stream( output_file );
                        block.printDataBlock( output_stream );
                    },
                    py::arg( "output_file" ),
                    R"doc(
         Write the contents of the data block to a text file

                           The file is created if it does not exist, and it can have, for example, txt extension

                           Parameters
                           ----------
                           output_file : str
                               Contents will be written to the file defined by this path
      )doc" );

    py::class_< tio::OdfDataSpecificBlock, std::shared_ptr< tio::OdfDataSpecificBlock > >( m, "OdfDataSpecificBlock", R"doc(
         Base class observable-dependent ODF data containers

                 The data section of an ODF is split into blocks (lines), each associated with
                 an observation epoch. The first elements of each block (e.g. observation epoch,
                 value of the observable) are common for all the observable types, but the
                 values in the remaining columns will have different meanings for different
                 types of observations. This base class serves as parent for interfaces to the
                 observable-specific part of an ODF data block. The interface to the common
                 part of the blocks is provided by the OdfCommonDataBlock class.
      )doc" );

    py::class_< tio::OdfDopplerDataBlock, std::shared_ptr< tio::OdfDopplerDataBlock >, tio::OdfDataSpecificBlock >(
            m, "OdfDopplerDataBlock", R"doc(Container for ODF Doppler-specific data)doc" )
            .def_property_readonly( "receiver_channel", &tio::OdfDopplerDataBlock::getReceiverChannel )
            .def_property_readonly( "spacecraft_id", &tio::OdfDopplerDataBlock::getSpacecraftId )
            .def_property_readonly( "receiver_exciter_flag", &tio::OdfDopplerDataBlock::getReceiverExciterFlag )
            .def_property_readonly( "reference_frequency", &tio::OdfDopplerDataBlock::getReferenceFrequency )
            .def_property_readonly( "compression_time", &tio::OdfDopplerDataBlock::getCompressionTime )
            .def_property_readonly( "transmitting_station_uplink_delay", &tio::OdfDopplerDataBlock::getTransmittingStationUplinkDelay );

    py::class_< tio::OdfDataBlock, std::shared_ptr< tio::OdfDataBlock > >(
            m, "OdfDataBlock", R"doc(Contents of a line of the data section of an ODF)doc" )
            .def_property_readonly( "observable_specific_data_block", &tio::OdfDataBlock::getObservableSpecificDataBlock )
            .def_property_readonly( "common_data_block", &tio::OdfDataBlock::getCommonDataBlock )
            .def(
                    "print_data_block",
                    []( tio::OdfDataBlock& block, const std::string& output_file ) {
                        std::ofstream output_stream( output_file );
                        block.printDataBlock( output_stream );
                    },
                    py::arg( "output_file" ),
                    R"doc(
         Write the contents of the data block to a text file

                           The file is created if it does not exist, and it can have, for example, txt extension

                           Parameters
                           ----------
                           output_file : str
                               Contents will be written to the file defined by this path
      )doc" );

    py::class_< tio::OdfRampBlock, std::shared_ptr< tio::OdfRampBlock > >(
            m, "OdfRampBlock", R"doc(Contents of a line of the ramp section of an ODF)doc" )
            .def_property_readonly( "ramp_start_frequency", &tio::OdfRampBlock::getRampStartFrequency )
            .def_property_readonly( "ramp_rate", &tio::OdfRampBlock::getRampRate )
            .def_property_readonly( "ramp_start_epoch", &tio::OdfRampBlock::getRampStartTime )
            .def_property_readonly( "ramp_end_epoch", &tio::OdfRampBlock::getRampEndTime )
            .def_property_readonly( "transmitting_station_id", &tio::OdfRampBlock::getTransmittingStationId );

    py::class_< tio::OdfClockOffsetBlock, std::shared_ptr< tio::OdfClockOffsetBlock > >( m, "OdfClockOffsetBlock" )
            .def_property_readonly( "start_time", &tio::OdfClockOffsetBlock::getStartTime )
            .def_property_readonly( "end_time", &tio::OdfClockOffsetBlock::getEndTime )
            .def_property_readonly( "clock_offset", &tio::OdfClockOffsetBlock::getClockOffset );

    py::class_< tio::OdfRawFileContents, std::shared_ptr< tio::OdfRawFileContents > >(
            m, "RawOdfFileContents", R"doc(Raw contents of a parsed ODF file.)doc" )
            .def_property_readonly( "data_blocks", &tio::OdfRawFileContents::getDataBlocks, R"doc(
         **read-only**

         ODF data blocks.

         :type: list[OdfDataBlock]
      )doc" )
            .def_property_readonly( "ramp_blocks", &tio::OdfRawFileContents::getRampBlocks, R"doc(
         **read-only**

         ODF ramp blocks.

         :type: list[OdfRampBlock]
      )doc" )
            .def_property_readonly( "clock_offset_blocks", &tio::OdfRawFileContents::getClockOffsetBlocks, R"doc(
         **read-only**

         ODF clock-offset blocks.

         :type: list[OdfClockOffsetBlock]
      )doc" )
            .def_readonly( "file_reference_date", &tio::OdfRawFileContents::fileReferenceDate_, R"doc(
         **read-only**

         ODF file reference date.

         :type: int
      )doc" )
            .def_readonly( "file_reference_time", &tio::OdfRawFileContents::fileReferenceTime_, R"doc(
         **read-only**

         ODF file reference time.

         :type: float
      )doc" )
            .def( "write_to_text_file",
                  &tio::OdfRawFileContents::writeOdfToTextFile,
                  py::arg( "output_file" ),
                  R"doc(
         Write the raw ODF contents to a text file.

         Parameters
         ----------
         output_file : str
             Path of the text file to write.

         Returns
         -------
         None
      )doc" );

    m.def( "read_raw_odf_file_contents",
           &tio::readOdfFile,
           py::arg( "file_name" ),
           R"doc(
         Read an ODF file without converting it to tracking data.

         This is a supporting reader for inspecting raw ODF contents. In the
         typical Tudat workflow, call :func:`read_odf_data` instead, so the file
         contents are converted to Tudat tracking-data objects.

         Parameters
         ----------
         file_name : str
             Path to the ODF file.

         Returns
         -------
         RawOdfFileContents
             Raw parsed ODF file contents.
      )doc" );
    m.def( "read_odf_data",
           &tio::loadOdfFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "odf_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "verbose_output" ) = true,
           R"doc(
         Read ODF files and convert them to tracking-data containers.

         ODF files are binary Deep Space Network (DSN) orbit-data files containing
         radiometric spacecraft tracking data. The format is described in ``820-013,
         TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E,
         2008, JPL/DSN``.
         The reader parses the raw ODF contents, which can also be inspected with
         :func:`read_raw_odf_file_contents`, and converts the selected tracking
         records to Tudat ``TrackingData`` and ``TrackingSupplementaryData``
         objects.

         Parameters
         ----------
         odf_file_names : list[str]
             Paths to ODF files.
         spacecraft_name : str
             Name assigned to the spacecraft link end (to be used for all files).
         earth_name : str, default="Earth"
             Name assigned to the Earth body (to be used only when it is needed
             to place the ground stations on a body not named 'Earth').
         verbose_output : bool, default=True
             Whether reader progress is printed (to be used for all files).

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );
}

}  // namespace odf

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
