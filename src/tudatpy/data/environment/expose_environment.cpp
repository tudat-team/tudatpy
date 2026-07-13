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
#include "expose_environment.h"

#include "tudat/io/trackingSupplementaryData.h"

#include "coma/expose_coma.h"
#include "missions/expose_missions.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data
{
namespace environment
{

void expose_environment( py::module& m )
{
    auto coma_submodule = m.def_submodule( "coma" );
    coma::expose_coma( coma_submodule );

    auto missions_submodule = m.def_submodule( "missions" );
    missions::expose_missions( missions_submodule );

    py::class_< tudat::data::FrequencySupplementaryData, std::shared_ptr< tudat::data::FrequencySupplementaryData > >(
            m, "FrequencySupplementaryData", R"doc(Base class for frequency supplementary data containers.)doc" );

    py::class_< tudat::data::RampedFrequencySupplementaryData,
                std::shared_ptr< tudat::data::RampedFrequencySupplementaryData >,
                tudat::data::FrequencySupplementaryData >( m, "RampedFrequencySupplementaryData", R"doc(
 Container for a set of ramped (piecewise-linear) frequency intervals.
    )doc" )
            .def( py::init<>( ) )
            .def( "add_frequency_ramp",
                  &tudat::data::RampedFrequencySupplementaryData::addFrequencyRamp,
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ),
                  R"doc(Adds a single frequency ramp interval to the container.)doc" );

    py::class_< tudat::data::TrackingSupplementaryData, std::shared_ptr< tudat::data::TrackingSupplementaryData > >(
            m, "TrackingSupplementaryData", R"doc(Tracking supplementary data container.)doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::string&, const std::string& >( ), py::arg( "body_name" ), py::arg( "reference_point_name" ) )
            .def_property_readonly( "body_name", &tudat::data::TrackingSupplementaryData::getBodyName )
            .def_property_readonly( "reference_point_name", &tudat::data::TrackingSupplementaryData::getReferencePointName )
            .def( "set_frequency_supplementary_data",
                  &tudat::data::TrackingSupplementaryData::setFrequencySupplementaryData,
                  py::arg( "frequency_supplementary_data" ),
                  R"doc(Sets the list of frequency supplementary data entries (e.g. ``RampedFrequencySupplementaryData``).)doc" );
}

}  // namespace environment
}  // namespace data
}  // namespace tudatpy
