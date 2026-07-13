/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_fdets.h"

#include "scalarTypes.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{
namespace data
{
namespace tracking
{
namespace fdets
{

void expose_fdets( py::module& m )
{
    py::enum_< tudat::input_output::FdetDateFormat >( m, "FdetDateFormat", R"doc(Date format used in an Fdets file.)doc" )
            .value( "datetime_string",
                    tudat::input_output::FdetDateFormat::datetime_string,
                    R"doc(Date is provided in a single UTC datetime string column.)doc" )
            .value( "pair_of_numbers",
                    tudat::input_output::FdetDateFormat::pair_of_numbers,
                    R"doc(Date is provided as a pair of numeric columns. This format is not currently supported.)doc" )
            .export_values( );

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, tio::FdetDateFormat >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "date_format" ) = tio::FdetDateFormat::datetime_string,
           R"doc(Load contents of Fdets file into object

           This function loads the contents of an Fdets file. For the currently supported `datetime_string` date format, the following structure is assumed:

           - UTC datetime string
           - Signal to noise ratio [-]
           - Normalized spectral max (Not sure what this is)
           - Doppler measured frequency [Hz]
           - Doppler noise [Hz]

           If the file contains an additional leading scan-number column, this column is detected automatically.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           date_format : FdetDateFormat
               Date format used in the Fdets file

           Returns
           -------
           fdets_contents : TrackingTxtFileContents
               Dictionary with contents of the Fdets file as lists of strings
           )doc" );

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, const std::vector< std::string >& >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           R"doc(Load contents of Fdets file into object using explicit column identifiers.

           .. deprecated::
              Passing explicit column identifiers is deprecated. Use the `date_format` argument instead.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           column_types : List[str]
               Identifiers of the columns present in the Fdets file

           Returns
           -------
           fdets_contents : TrackingTxtFileContents
               Dictionary with contents of the Fdets file as lists of strings
           )doc" );

    m.def( "read_fdets_files",
           static_cast< std::pair< std::vector< std::shared_ptr< tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >,
                                   std::vector< std::shared_ptr< tudat::data::TrackingSupplementaryData > > > ( * )(
                   const std::vector< std::string >&,
                   const std::vector< double >&,
                   const tio::FdetDateFormat,
                   const std::string&,
                   const std::vector< std::string >&,
                   const std::vector< std::string >&,
                   const std::string& ) >( &tio::readFdetsFiles< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "fdets_file_names" ),
           py::arg( "base_frequencies" ),
           py::arg( "date_format" ),
           py::arg( "spacecraft_name" ),
           py::arg( "transmitting_station_names" ),
           py::arg( "receiving_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           R"doc(Load FDETS files into tracking data and supplementary data objects.)doc" );
}

}  // namespace fdets
}  // namespace tracking
}  // namespace data
}  // namespace tudatpy
