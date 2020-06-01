//
// Created by ggarrett on 29-04-20.
//

#include "expose_orbital_element_conversions.h"

namespace tudatpy {
    void expose_orbital_element_conversions(py::module &m) {


        py::enum_<toec::KeplerianElementIndices>(m, "KeplerianElementIndices")
                .value("semi_major_axis_index", toec::KeplerianElementIndices::semiMajorAxisIndex)
                .value("eccentricity_index", toec::KeplerianElementIndices::eccentricityIndex)
                .value("inclination_index", toec::KeplerianElementIndices::inclinationIndex)
                .value("argument_of_periapsis_index", toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
                .value("longitude_of_ascending_node_index",
                       toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
                .value("true_anomaly_index", toec::KeplerianElementIndices::trueAnomalyIndex)
                .value("semi_latus_rectum_index", toec::KeplerianElementIndices::semiLatusRectumIndex)
                .export_values();

        py::enum_<toec::SphericalOrbitalStateElementIndices>(m, "SphericalOrbitalStateElementIndices")
                .value("radius_index", toec::SphericalOrbitalStateElementIndices::radiusIndex)
                .value("latitude_index", toec::SphericalOrbitalStateElementIndices::latitudeIndex)
                .value("longitude_index", toec::SphericalOrbitalStateElementIndices::longitudeIndex)
                .value("speed_index", toec::SphericalOrbitalStateElementIndices::speedIndex)
                .value("flight_path_index", toec::SphericalOrbitalStateElementIndices::flightPathIndex)
                .value("heading_angle_index", toec::SphericalOrbitalStateElementIndices::headingAngleIndex)
                .export_values();

        m.def("convert_keplerian_to_cartesian_elements",
              &toec::convertKeplerianToCartesianElements<>,
              py::arg("keplerian_elements"),
              py::arg("central_body_gravitational_parameter"));

        m.def("convert_spherical_orbital_to_cartesian_state",
              &toec::convertSphericalOrbitalToCartesianState<>,
              py::arg("spherical_orbital_state"));

    };
}
