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
#include "expose_observations_wrapper_bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"

#include "tudat/simulation/estimation_setup/processSumLmkFiles.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_wrapper
{

void expose_observations_wrapper_sum_lmk_bindings( py::module& m )
{
    py::class_< tom::SumLmkObservationConversionSettings >( m, "SumLmkObservationConversionSettings", R"doc(
        Settings for converting SPC-style SUM/LMK optical-landmark files to Tudat pixel-coordinate observations.
        )doc" )
            .def( py::init< const std::string&, const std::string& >( ), py::arg( "target_body_name" ), py::arg( "receiver_body_name" ) )
            .def_readwrite( "target_body_name", &tom::SumLmkObservationConversionSettings::targetBodyName_ )
            .def_readwrite( "receiver_body_name", &tom::SumLmkObservationConversionSettings::receiverBodyName_ )
            .def_readwrite( "body_fixed_camera_position", &tom::SumLmkObservationConversionSettings::bodyFixedCameraPosition_ )
            .def_readwrite( "validate_spacecraft_object_geometry",
                            &tom::SumLmkObservationConversionSettings::validateSpacecraftObjectGeometry_ );

    py::class_< tom::SumLmkObservationConversionResult< STATE_SCALAR_TYPE, TIME_TYPE > >( m, "SumLmkObservationConversionResult", R"doc(
        Result of a SUM/LMK observation conversion: the observation collection, matching observation
        model settings, per-image pointing a-priori inverse covariance entries, and image-to-camera map.
        )doc" )
            .def_readonly( "observation_collection",
                           &tom::SumLmkObservationConversionResult< STATE_SCALAR_TYPE, TIME_TYPE >::observationCollection_,
                           R"doc(Pixel-coordinate observation collection (observed values), keyed by (image, landmark) link ends.)doc" )
            .def_readonly( "observation_model_settings",
                           &tom::SumLmkObservationConversionResult< STATE_SCALAR_TYPE, TIME_TYPE >::observationModelSettings_,
                           R"doc(Observation model settings matching the collection (one pixel_coordinates setting per link end).)doc" )
            .def_readonly(
                    "inverse_apriori_covariance_diagonal_entries",
                    &tom::SumLmkObservationConversionResult< STATE_SCALAR_TYPE, TIME_TYPE >::inverseAprioriCovarianceDiagonalEntries_,
                    R"doc(Per-image camera_pointing_correction inverse a-priori covariance diagonal entries derived from SIGMA_PTG.)doc" )
            .def_readonly( "image_id_to_camera_name",
                           &tom::SumLmkObservationConversionResult< STATE_SCALAR_TYPE, TIME_TYPE >::imageIdToCameraName_,
                           R"doc(Map from SUM image ID to the registered Camera_<imageId> reference-point name.)doc" );

    m.def( "create_sum_lmk_observation_collection",
           py::overload_cast< const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const tss::SystemOfBodies&,
                              const tom::SumLmkObservationConversionSettings& >(
                   &tom::createSumLmkObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "sum_files" ),
           py::arg( "lmk_files" ),
           py::arg( "bodies" ),
           py::arg( "conversion_settings" ),
           R"doc(
        Create a pixel-coordinate observation collection from SPC-style SUM/LMK files.

        Reads the given SUM and LMK files, registers one ``Camera_<imageId>`` per image on the receiver
        body and each landmark as a body-fixed ground station on the target body (both added to ``bodies``
        in place), and returns the converted observations together with their matching observation model
        settings and the per-image SIGMA_PTG pointing a-priori. Each (image, landmark) pair becomes one
        observation; all landmarks of an image share the same camera (and thus the same pointing parameter).
        )doc" );

    m.def( "create_sum_lmk_observation_model_settings",
           &tom::createSumLmkObservationModelSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "light_time_corrections" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           R"doc(
        Build pixel-coordinate observation model settings matching every (image, landmark) link end in a
        SUM/LMK observation collection (light-time geometric single-leg on, stellar aberration off).
        )doc" );

    m.def( "compute_sum_lmk_residuals",
           &tom::computeSumLmkResiduals< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "bodies" ),
           py::arg( "light_time_corrections" ) = std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           R"doc(
        Compute observed-minus-computed (O-C) pixel residuals for a SUM/LMK observation collection given a
        fixed environment (e.g. a SPICE spacecraft trajectory). Residuals are stored in the collection
        (accessible via ``get_concatenated_residuals`` / ``get_rms_residuals``) and returned concatenated.
        )doc" );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
