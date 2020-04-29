//
// Created by ggarrett on 29-04-20.
//

#include "expose_orbital_element_conversions.h"

namespace tudatpy {
    void expose_orbital_element_conversions(py::module &m){

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


        m.def("convert_keplerian_to_cartesian_elements", &toec::convertKeplerianToCartesianElements<>);

    };
}
