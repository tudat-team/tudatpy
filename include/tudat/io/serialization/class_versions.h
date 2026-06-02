#ifndef TUDAT_SERIALIZATION_CLASS_VERSIONS_H
#define TUDAT_SERIALIZATION_CLASS_VERSIONS_H

#include <cstdint>

namespace tudat
{

namespace serialization
{

//! Central manifest of cereal class versions for serializable Tudat types.
struct ClassVersions {
    static constexpr std::uint32_t simulation_results_dd = 1;
    static constexpr std::uint32_t simulation_results_dt = 1;

    static constexpr std::uint32_t single_arc_dynamics_results = 1;
    static constexpr std::uint32_t single_arc_variational_results = 1;
    static constexpr std::uint32_t multi_arc_dynamics_results = 1;
    static constexpr std::uint32_t multi_arc_variational_results = 1;
    static constexpr std::uint32_t hybrid_arc_dynamics_results = 1;
    static constexpr std::uint32_t hybrid_arc_variational_results = 1;
};

}  // namespace serialization

}  // namespace tudat

#endif  // TUDAT_SERIALIZATION_CLASS_VERSIONS_H
