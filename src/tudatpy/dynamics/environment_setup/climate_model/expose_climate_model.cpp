/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_climate_model.h"
#include <tudat/simulation/environment_setup.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tudat/simulation/environment_setup/createClimateModel.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace climate_model
{

void expose_climate_model_setup( py::module &m )
{
    py::class_< tss::ClimateModelSettings, std::shared_ptr< tss::ClimateModelSettings > >( m,
                                                                                     "ClimateModelSettings",
                                                                                     R"doc(

         Base class for providing settings for climate model.

      )doc" );
      
#if TUDAT_BUILD_WITH_MCD_INTERFACE

    m.def( "mars_climate_database",
           &tss::marsClimateDatabaseClimateModelSettings,
           py::arg( "mcd_data_path" ) = "",
           py::arg( "dust_scenario" ) = 1,
           py::arg( "perturbation_key" ) = 0,
           py::arg( "perturbation_seed" ) = 0.0,
           py::arg( "gravity_wave_length" ) = 0.0,
           py::arg( "high_resolution_mode" ) = 0,
           R"doc(

 Function for creating oblate spherical body shape model settings.

 Function for settings object, defining oblate spherical body shape model from equatorial radius and flattening parameter.


 Parameters
 ----------
 mcd_data_path : str, default = ""
     MCD data path
 dust_scenario : int, default = 1
     Dust and solar EUV scenario (1-8 or 24-35)
 perturbation_key : int, default = 0
     Perturbation type (0-5)
 perturbation_seed : double, default = 0.0
     Random seed or scaling factor
 gravity_wave_length : double, default = 0.0
     Gravity wave wavelength in meters
high_resolution_mode : int, default = 0
    High resolution topography flag (0 or 1), requires MOLA dataset
 Returns
 -------
 ClimateModelSettings
     )doc" );

#endif

}

}  // namespace climate_model
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy