#include "expose_terrestrial_frame_transformations.h"
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace terrestrial_frame_transformations
{
void expose_terrestrial_frame_transformations( pybind11::module& module )
{
    namespace py = pybind11;
    namespace rf = tudat::reference_frames;
    py::class_< rf::TerrestrialFrameEvaluationEpoch >( module,
                                                       "TerrestrialFrameEvaluationEpoch",
                                                       "Published decimal year and its positive derivative relative to the input state's "
                                                       "seconds; no implicit TT/TDB/calendar mapping." )
            .def( py::init( []( double year, double rate ) { return rf::TerrestrialFrameEvaluationEpoch{ year, rate }; } ),
                  py::arg( "decimal_year" ),
                  py::arg( "decimal_year_rate_with_respect_to_input_seconds" ) )
            .def_readwrite( "decimal_year", &rf::TerrestrialFrameEvaluationEpoch::decimalYear )
            .def_readwrite( "decimal_year_rate_with_respect_to_input_seconds",
                            &rf::TerrestrialFrameEvaluationEpoch::decimalYearRateWithRespectToInputSeconds );
    py::class_< rf::TerrestrialFrameAffineTransformation >( module, "TerrestrialFrameAffineTransformation" )
            .def( py::init<>( ) )
            .def_readwrite( "translation", &rf::TerrestrialFrameAffineTransformation::translation )
            .def_readwrite( "matrix", &rf::TerrestrialFrameAffineTransformation::matrix )
            .def_readwrite( "translation_rate", &rf::TerrestrialFrameAffineTransformation::translationRate )
            .def_readwrite( "matrix_rate", &rf::TerrestrialFrameAffineTransformation::matrixRate )
            .def( "transform_position", &rf::TerrestrialFrameAffineTransformation::transformPosition, py::arg( "position" ) )
            .def( "transform_state", &rf::TerrestrialFrameAffineTransformation::transformState, py::arg( "state" ) )
            .def( "inverse", &rf::TerrestrialFrameAffineTransformation::inverse )
            .def( "compose_after", &rf::TerrestrialFrameAffineTransformation::composeAfter, py::arg( "before" ) );
    py::enum_< rf::ItrfTransformationCatalogue >( module, "ItrfTransformationCatalogue" )
            .value( "ign_itrf2014", rf::ItrfTransformationCatalogue::ignItrf2014 )
            .value( "ign_itrf2020", rf::ItrfTransformationCatalogue::ignItrf2020 );
    module.def( "conventional_itrf_parameters", &rf::getConventionalItrfTransformationData, py::arg( "catalogue" ) );
    module.def( "evaluate_linearized_helmert", &rf::evaluateLinearizedHelmert, py::arg( "parameters" ), py::arg( "epoch" ) );
    py::class_< rf::TerrestrialFrameTransformationModel >( module, "TerrestrialFrameTransformationModel" )
            .def( py::init< const tudat::input_output::TerrestrialFrameTransformationData&,
                            rf::TerrestrialFrameTransformationModel::Clock >( ),
                  py::arg( "data" ),
                  py::arg( "clock" ) = rf::TerrestrialFrameTransformationModel::Clock( ),
                  "Use preloaded parameters and an optional model-owned callback mapping Time to an explicit evaluation epoch. A clock is "
                  "required by at_time methods. All state velocities are with respect to the callback's input seconds; input release "
                  "provenance does not invent a different frame datum." )
            .def( "get_transformation",
                  &rf::TerrestrialFrameTransformationModel::getTransformation,
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "epoch" ) )
            .def( "get_transformation_at_time",
                  &rf::TerrestrialFrameTransformationModel::getTransformationAtTime,
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "input_time" ) )
            .def( "transform_position",
                  &rf::TerrestrialFrameTransformationModel::transformPosition,
                  py::arg( "position" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "epoch" ) )
            .def( "transform_state",
                  &rf::TerrestrialFrameTransformationModel::transformState,
                  py::arg( "state" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "epoch" ) )
            .def( "transform_position_at_time",
                  &rf::TerrestrialFrameTransformationModel::transformPositionAtTime,
                  py::arg( "position" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "input_time" ) )
            .def( "transform_state_at_time",
                  &rf::TerrestrialFrameTransformationModel::transformStateAtTime,
                  py::arg( "state" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "input_time" ) )
            .def_property_readonly( "data", []( const rf::TerrestrialFrameTransformationModel& model ) { return model.getData( ); } );
}
}  // namespace terrestrial_frame_transformations
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
