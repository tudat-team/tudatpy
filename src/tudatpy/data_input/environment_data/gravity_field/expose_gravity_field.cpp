#include "expose_gravity_field.h"

#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "tudat/io/readIcgemGravityField.h"
#include "tudat/io/readOceanPoleGravityData.h"
#include "tudat/io/readHarmonicPotentialData.h"

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace gravity_field
{

//! Expose published coefficients without normalizing or adding a physical response.
void expose_gravity_field( pybind11::module& module )
{
    // Use the same provider types as the native readers and consuming settings.
    namespace py = pybind11;

    namespace io = tudat::input_output;

    namespace astro = tudat::basic_astrodynamics;

    py::enum_< astro::TideSystem >( module, "TideSystem", "Permanent-tide convention declared by the provider." )
            .value( "unknown", astro::TideSystem::unknown, "No convention declared." )
            .value( "tide_free", astro::TideSystem::tideFree, "Permanent tidal deformation removed." )
            .value( "zero_tide", astro::TideSystem::zeroTide, "Permanent deformation retained, direct permanent potential removed." )
            .value( "mean_tide", astro::TideSystem::meanTide, "Permanent deformation and direct tidal potential retained." );

    py::enum_< io::IcgemNormalization >( module, "IcgemNormalization", "Normalization stated in an ICGEM header." )
            .value( "fully_normalized", io::IcgemNormalization::fullyNormalized )
            .value( "unnormalized", io::IcgemNormalization::unnormalized );

    py::enum_< io::IcgemCoefficientErrors >( module, "IcgemCoefficientErrors", "Uncertainty columns declared in an ICGEM header." )
            .value( "none", io::IcgemCoefficientErrors::none )
            .value( "calibrated", io::IcgemCoefficientErrors::calibrated )
            .value( "formal", io::IcgemCoefficientErrors::formal )
            .value( "calibrated_and_formal", io::IcgemCoefficientErrors::calibratedAndFormal );

    py::class_< io::IcgemGravityCoefficient >( module, "IcgemGravityCoefficient", "One published dimensionless static Stokes coefficient." )
            .def( py::init<>( ) )
            .def_readwrite( "cosine", &io::IcgemGravityCoefficient::cosine, "Cosine coefficient in the declared normalization." )
            .def_readwrite( "sine", &io::IcgemGravityCoefficient::sine, "Sine coefficient in the declared normalization." )
            .def_readwrite(
                    "calibrated_errors", &io::IcgemGravityCoefficient::calibratedErrors, "Calibrated cosine/sine uncertainties, or None." )
            .def_readwrite( "formal_errors", &io::IcgemGravityCoefficient::formalErrors, "Formal cosine/sine uncertainties, or None." )
            .def_readwrite( "source_line", &io::IcgemGravityCoefficient::sourceLine, "Line number in the provider file." )
            .def_readwrite( "original_record", &io::IcgemGravityCoefficient::originalRecord, "Unmodified provider record." );

    py::class_< io::IcgemGravityFieldData >(
            module, "IcgemGravityFieldData", "Sparse static ICGEM records; missing entries remain distinct from zero." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::IcgemGravityFieldData::source, "Provider filename or source description." )
            .def_readwrite( "header_lines", &io::IcgemGravityFieldData::headerLines, "Unmodified header lines." )
            .def_readwrite( "product_type", &io::IcgemGravityFieldData::productType, "Declared product type." )
            .def_readwrite( "model_name", &io::IcgemGravityFieldData::modelName, "Declared model name." )
            .def_readwrite( "format", &io::IcgemGravityFieldData::format, "Declared ICGEM format version." )
            .def_readwrite( "gravitational_parameter",
                            &io::IcgemGravityFieldData::gravitationalParameter,
                            "Source GM in cubic metres per squared second." )
            .def_readwrite( "reference_radius", &io::IcgemGravityFieldData::referenceRadius, "Source reference radius in metres." )
            .def_readwrite( "maximum_degree", &io::IcgemGravityFieldData::maximumDegree, "Maximum degree declared by the provider." )
            .def_readwrite( "errors", &io::IcgemGravityFieldData::errors, "Declared uncertainty columns." )
            .def_readwrite( "normalization", &io::IcgemGravityFieldData::normalization, "Declared coefficient normalization." )
            .def_readwrite(
                    "tide_system", &io::IcgemGravityFieldData::tideSystem, "Declared permanent-tide convention; unknown remains unknown." )
            .def_readwrite(
                    "coefficients", &io::IcgemGravityFieldData::coefficients, "Dictionary of (degree, order) to published coefficients." )
            .def( "validate", &io::validateIcgemGravityFieldData, "Check header declarations and finite, valid coefficient records." );

    py::class_< io::OceanPoleGravityCoefficient >(
            module, "OceanPoleGravityCoefficient", "Normalized complex A/B equilibrium integrals in IERS equation 6.23." )
            .def( py::init<>( ) )
            .def_readwrite( "a", &io::OceanPoleGravityCoefficient::a, "Dimensionless complex A integral." )
            .def_readwrite( "b", &io::OceanPoleGravityCoefficient::b, "Dimensionless complex B integral." )
            .def_readwrite( "source_line", &io::OceanPoleGravityCoefficient::sourceLine, "Line number in the provider file." );

    py::class_< io::OceanPoleGravityData >(
            module, "OceanPoleGravityData", "Published ocean-pole gravity integrals, without a load response or forcing constants." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::OceanPoleGravityData::source, "Provider filename or source description." )
            .def_readwrite( "header_lines", &io::OceanPoleGravityData::headerLines, "Unmodified header lines." )
            .def_readwrite(
                    "coefficients", &io::OceanPoleGravityData::coefficients, "Dictionary of (degree, order) to equilibrium integrals." )
            .def_readwrite( "maximum_degree", &io::OceanPoleGravityData::maximumDegree, "Maximum published degree." )
            .def( "validate", &io::validateOceanPoleGravityData, "Check indices and finite published integrals." );

    py::enum_< io::HarmonicPotentialProfile >(
            module, "HarmonicPotentialProfile", "Explicit provider column, inventory and unit conventions." )
            .value( "iers_potential", io::HarmonicPotentialProfile::iersPotential )
            .value( "fes2004", io::HarmonicPotentialProfile::fes2004 )
            .value( "fes2014b_v1", io::HarmonicPotentialProfile::fes2014bV1 )
            .value( "costg_eot11a", io::HarmonicPotentialProfile::costgEot11a );

    py::class_< io::HarmonicPotentialCoefficient >(
            module, "HarmonicPotentialCoefficient", "Published Cp, Sp, Cm, Sm in declared Stokes units." )
            .def( py::init<>( ) )
            .def_readwrite( "values", &io::HarmonicPotentialCoefficient::values, "Ordered Cp, Sp, Cm, Sm values, before unit conversion." )
            .def_readwrite( "source_line", &io::HarmonicPotentialCoefficient::sourceLine, "Line number in the provider file." );

    py::class_< io::HarmonicPotentialConstituent >(
            module, "HarmonicPotentialConstituent", "One provider wave and its sparse published coefficient records." )
            .def( py::init<>( ) )
            .def_readwrite( "original_doodson", &io::HarmonicPotentialConstituent::originalDoodson, "Doodson identifier as published." )
            .def_readwrite( "doodson", &io::HarmonicPotentialConstituent::doodson, "Canonical Doodson identifier." )
            .def_readwrite( "name", &io::HarmonicPotentialConstituent::name, "Provider wave name." )
            .def_readwrite( "multipliers", &io::HarmonicPotentialConstituent::multipliers, "Six integer Doodson multipliers." )
            .def_readwrite( "coefficients",
                            &io::HarmonicPotentialConstituent::coefficients,
                            "Dictionary of (degree, order) to potential amplitudes." )
            .def_readwrite(
                    "minimum_degree", &io::HarmonicPotentialConstituent::minimumDegree, "Profile lower degree, or -1 when unspecified." )
            .def_readwrite(
                    "maximum_degree", &io::HarmonicPotentialConstituent::maximumDegree, "Profile upper degree, or -1 when unspecified." )
            .def_readwrite(
                    "maximum_order", &io::HarmonicPotentialConstituent::maximumOrder, "Profile upper order, or -1 when unspecified." );

    py::class_< io::HarmonicPotentialData >(
            module, "HarmonicPotentialData", "Ordered provider waves with no phase, loading or nodal corrections added." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::HarmonicPotentialData::source, "Provider filename or source description." )
            .def_readwrite( "header_lines", &io::HarmonicPotentialData::headerLines, "Unmodified header lines." )
            .def_readwrite( "profile", &io::HarmonicPotentialData::profile, "Selected provider convention." )
            .def_readwrite( "stokes_unit_scale",
                            &io::HarmonicPotentialData::stokesUnitScale,
                            "Multiplier converting published amplitudes to dimensionless Stokes coefficients." )
            .def_readwrite( "constituents", &io::HarmonicPotentialData::constituents, "Waves in provider order." )
            .def( "validate", &io::validateHarmonicPotentialData, "Check wave identifiers, declared units and finite coefficients." );

    module.def( "read_icgem_gravity_field",
                py::overload_cast< const std::string& >( &io::readIcgemGravityField ),
                py::arg( "file_name" ),
                "Read static ICGEM gfc records, retaining normalization, tide convention and optional uncertainties. Temporal records "
                "reject." );

    module.def( "read_ocean_pole_gravity_data",
                py::overload_cast< const std::string& >( &io::readOceanPoleGravityData ),
                py::arg( "file_name" ),
                "Read normalized complex equilibrium integrals from the six-column IERS A/B table. Files must already be decompressed." );

    module.def( "read_harmonic_potential_data",
                py::overload_cast< const std::string&, io::HarmonicPotentialProfile >( &io::readHarmonicPotentialData ),
                py::arg( "file_name" ),
                py::arg( "profile" ) = io::HarmonicPotentialProfile::iersPotential,
                "Read normalized IERS potential columns using an explicit provider profile; water-height coefficients are not supported." );
}

}  // namespace gravity_field
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
