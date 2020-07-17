//
// Created by ggarrett on 02-05-20.
//

#include "expose_ephemerides.h"

namespace py = pybind11;
using namespace tudat::ephemerides;

namespace tudatpy {

    void expose_ephemerides(py::module &m) {

    	py::enum_<ApproximatePlanetPositionsBase::BodiesWithEphemerisData>(m, "BodiesWithEphemerisData")
    	        .value("mercury", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury)
				.value("venus", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus)
				.value("earth_moon_barycenter", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter)
				.value("mars", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars)
				.value("jupiter", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter)
				.value("saturn", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn)
				.value("uranus", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::uranus)
				.value("neptune", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune)
				.value("pluto", ApproximatePlanetPositionsBase::BodiesWithEphemerisData::pluto)
				.export_values();

        py::class_<RotationalEphemeris,
                std::shared_ptr<RotationalEphemeris>> RotationalEphemeris_(m, "RotationalEphemeris");

		py::class_<ApproximatePlanetPositions, std::shared_ptr<ApproximatePlanetPositions>>(m, "ApproximatePlanetPositionsEphemeris")
		.def(py::init<
		        ApproximatePlanetPositionsBase::BodiesWithEphemerisData,
		        const double
		        >(),
		        py::arg("body_with_ephemeris_data"),
		        py::arg("sun_gravitational_parameter") = 1.32712440018e20
			)
		.def("get_cartesian_state", &ApproximatePlanetPositions::getCartesianState,
				py::arg("seconds_since_epoch")
			)
		.def("get_keplerian_state", &ApproximatePlanetPositions::getKeplerianStateFromEphemeris,
				py::arg("seconds_since_epoch")
			);

        m.def("transform_state_to_global_frame",
              &tudat::ephemerides::transformStateToGlobalFrame<double, double>,
              py::arg("state_in_local_frame"),
              py::arg("current_time"),
              py::arg("rotational_ephemeris")
        );

    };

}
