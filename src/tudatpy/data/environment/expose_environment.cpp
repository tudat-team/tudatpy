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

    py::class_< tudat::data::TrackingSupplementaryData, std::shared_ptr< tudat::data::TrackingSupplementaryData > >(
            m, "TrackingSupplementaryData", R"doc(Tracking supplementary data container.)doc" )
            .def_property_readonly( "body_name", &tudat::data::TrackingSupplementaryData::getBodyName )
            .def_property_readonly( "reference_point_name", &tudat::data::TrackingSupplementaryData::getReferencePointName );
}

}  // namespace environment
}  // namespace data
}  // namespace tudatpy
