#include "expose_gravity_field.h"

#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "tudat/io/readIcgemGravityField.h"
#include "tudat/io/readOceanPoleGravityData.h"
#include "tudat/io/readHarmonicPotentialData.h"
#include "tudat/io/readAod1bGravityField.h"

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

    py::class_< io::IcgemGravityCoefficient >(
            module, "IcgemGravityCoefficient", "Published coefficient values and provenance; units depend on the containing record." )
            .def( py::init<>( ) )
            .def_readwrite(
                    "cosine", &io::IcgemGravityCoefficient::cosine, "Cosine value in the record's declared normalization and units." )
            .def_readwrite( "sine", &io::IcgemGravityCoefficient::sine, "Sine value in the record's declared normalization and units." )
            .def_readwrite(
                    "calibrated_errors", &io::IcgemGravityCoefficient::calibratedErrors, "Calibrated cosine/sine uncertainties, or None." )
            .def_readwrite( "formal_errors", &io::IcgemGravityCoefficient::formalErrors, "Formal cosine/sine uncertainties, or None." )
            .def_readwrite( "source_line", &io::IcgemGravityCoefficient::sourceLine, "Line number in the provider file." )
            .def_readwrite( "original_record", &io::IcgemGravityCoefficient::originalRecord, "Unmodified provider record." );

    py::class_< io::IcgemTemporalCoefficient >(
            module, "IcgemTemporalCoefficient", "One gfct/trnd/dot/asin/acos provider record with independent printed dates and units." )
            .def( py::init<>( ) )
            .def_readwrite( "type", &io::IcgemTemporalCoefficient::type, "Original temporal record keyword." )
            .def_readwrite( "degree", &io::IcgemTemporalCoefficient::degree, "Spherical-harmonic degree." )
            .def_readwrite( "order", &io::IcgemTemporalCoefficient::order, "Spherical-harmonic order." )
            .def_readwrite( "coefficient",
                            &io::IcgemTemporalCoefficient::coefficient,
                            "Supplied values, uncertainties and original source record." )
            .def_readwrite( "reference_date",
                            &io::IcgemTemporalCoefficient::referenceDate,
                            "Printed reference/start date without an inferred time scale; only gfct supplies this in format 1.0." )
            .def_readwrite( "validity_end_date",
                            &io::IcgemTemporalCoefficient::validityEndDate,
                            "Independent printed format-2.0 validity end, or an empty string." )
            .def_readwrite( "period_years",
                            &io::IcgemTemporalCoefficient::periodYears,
                            "Supplied asin/acos period in years, or None; no year convention is inferred." );

    py::class_< io::IcgemGravityFieldData >(
            module, "IcgemGravityFieldData", "Sparse static and temporal ICGEM records; missing entries remain distinct from zero." )
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
            .def_readwrite( "temporal_coefficients",
                            &io::IcgemGravityFieldData::temporalCoefficients,
                            "Original ordered temporal records, without date conversion or evaluation." )
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

    // RL07 retains printed GPS epochs and distinct physical components before model construction.
    py::enum_< io::Aod1bComponent >(
            module, "Aod1bComponent", "Published non-tidal mass component; global_mass is ATM plus OCN, while OBA has a different domain." )
            .value( "atmosphere", io::Aod1bComponent::atmosphere )
            .value( "ocean", io::Aod1bComponent::ocean )
            .value( "global_mass", io::Aod1bComponent::global )
            .value( "ocean_bottom_pressure", io::Aod1bComponent::oceanBottomPressure );

    py::class_< io::Aod1bGravityCoefficient >(
            module, "Aod1bGravityCoefficient", "One supplied RL07 coefficient pair and original source line." )
            .def( py::init<>( ) )
            .def_readwrite(
                    "values", &io::Aod1bGravityCoefficient::values, "Dimensionless cosine/sine pair in the declared normalization." )
            .def_readwrite( "source_line", &io::Aod1bGravityCoefficient::sourceLine, "One-based provider line number." );

    py::class_< io::Aod1bGravitySample >(
            module, "Aod1bGravitySample", "One published component at a printed GPS epoch; missing entries remain distinct from zero." )
            .def( py::init<>( ) )
            .def_readwrite( "component", &io::Aod1bGravitySample::component, "Mass component declared in the DATA SET record." )
            .def_readwrite( "gps_date", &io::Aod1bGravitySample::gpsDate, "Original provider date/time token in GPS time." )
            .def_readwrite( "gps_seconds_since_j2000",
                            &io::Aod1bGravitySample::gpsSecondsSinceJ2000,
                            "Continuous GPS seconds from 2000-01-01 12:00:00 GPS." )
            .def_readwrite( "coefficients",
                            &io::Aod1bGravitySample::coefficients,
                            "Dictionary of supplied (degree, order) coefficient pairs, including degrees zero and one." );

    py::class_< io::Aod1bGravityFieldData >(
            module, "Aod1bGravityFieldData", "Neutral RL07 potential anomalies relative to the provider's 2007–2014 mean." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::Aod1bGravityFieldData::source, "Source filename or provider identifier." )
            .def_readwrite( "header_lines", &io::Aod1bGravityFieldData::headerLines, "Unmodified header including END OF HEADER." )
            .def_readwrite( "release", &io::Aod1bGravityFieldData::release, "Declared release number; this reader supports RL07." )
            .def_readwrite( "software_version", &io::Aod1bGravityFieldData::softwareVersion, "Provider software version." )
            .def_readwrite( "gravitational_parameter",
                            &io::Aod1bGravityFieldData::gravitationalParameter,
                            "Source GM in cubic metres per squared second." )
            .def_readwrite( "reference_radius", &io::Aod1bGravityFieldData::referenceRadius, "Source reference radius in metres." )
            .def_readwrite( "fully_normalized", &io::Aod1bGravityFieldData::fullyNormalized, "Declared coefficient normalization." )
            .def_readwrite( "maximum_degree", &io::Aod1bGravityFieldData::maximumDegree, "Declared maximum degree." )
            .def_readwrite( "samples",
                            &io::Aod1bGravityFieldData::samples,
                            "Original component/epoch records without component sums or normalization changes." )
            .def( "validate",
                  &io::validateAod1bGravityFieldData,
                  "Check RL07 metadata, GPS dates, coefficient presence and finite values." );

    module.def( "read_aod1b_gravity_field",
                py::overload_cast< const std::string& >( &io::readAod1bGravityField ),
                py::arg( "file_name" ),
                "Read one decompressed ASCII RL07 daily product, retaining its source metadata and printed GPS epochs." );

    module.def(
            "read_aod1b_gravity_fields",
            &io::readAod1bGravityFields,
            py::arg( "file_names" ),
            "Read successive local daily products without merging their raw records. Gzip decompression is the caller's responsibility." );

    module.def( "read_icgem_gravity_field",
                py::overload_cast< const std::string& >( &io::readIcgemGravityField ),
                py::arg( "file_name" ),
                "Read static and temporal ICGEM records, retaining normalization, tide convention, uncertainties and independent printed "
                "dates. Only the static field adapter rejects temporal data." );

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
