/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_propagation_bindings.h"

#include <tudat/astro/basic_astro/accelerationModel.h>
#include <tudat/astro/basic_astro/accelerationModelTypes.h>
#include <tudat/astro/basic_astro/massRateModel.h>
#include <tudat/astro/basic_astro/torqueModel.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace dynamics
{
namespace propagation
{

void expose_propagation_state_utility_types( py::module& m )
{
    py::class_< tba::TorqueModel, std::shared_ptr< tba::TorqueModel > >( m, "TorqueModel" );

    using AccelerationModel = tba::AccelerationModel< Eigen::Vector3d >;
    py::class_< AccelerationModel, std::shared_ptr< AccelerationModel > >( m, "AccelerationModel" )
            .def( "update_members", &AccelerationModel::updateMembers, py::arg( "current_time" ) )
            .def_property_readonly( "acceleration", &AccelerationModel::getAcceleration )
            .def(
                    "update_and_get_acceleration",
                    []( AccelerationModel& accelerationModel, const double currentTime ) {
                        accelerationModel.updateMembers( currentTime );
                        return accelerationModel.getAcceleration( );
                    },
                    py::arg( "current_time" ) )
            .def_property_readonly( "acceleration_type", []( const std::shared_ptr< AccelerationModel >& accelerationModel ) {
                return tba::getAccelerationModelType( accelerationModel );
            } );

    py::class_< tba::MassRateModel, std::shared_ptr< tba::MassRateModel > >( m, "MassRateModel" );
}

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy
