/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "expose_data.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
// #include <pybind11/native_enum.h>
#include <tudat/io/basicInputOutput.h>
#include <string>
#include <vector>

#include "tudat/io/missileDatcomData.h"
#include "tudat/io/comaModelInputOutput.h"
#include "tudat/io/readCrdFile.h"
#include "tudat/io/readHistoryFromFile.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/io/readSinexFile.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/io/readVariousPdsFiles.h"
#include "tudat/io/solarActivityData.h"

namespace py = pybind11;
namespace tio = tudat::input_output;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{

namespace data
{

void expose_data( py::module& m )
{
    py::module_::import( "tudatpy.kernel.math.interpolators" ).attr( "InterpolatorSettings" );
    // py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );
    // py::object cubic_spline_interpolation =
    //         (py::object)py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );

    m.def( "get_resource_path",
           &tudat::paths::get_resource_path,
           R"doc(

 Get the path at which tudat resources are located.

 Returns
 -------
 str
     Local path at which tudat resources are located.






     )doc" );
    m.def( "get_ephemeris_path",
           &tudat::paths::getEphemerisDataFilesPath,
           R"doc(

 Get the path at which the ephemeris used by tudat are located.

 Returns
 -------
 str
     Local path at which the tudat ephemeris resources are located.






     )doc" );
    m.def( "get_earth_orientation_path",
           &tudat::paths::getEarthOrientationDataFilesPath,
           R"doc(

 Get the path at which the Earth orientation resources used by tudat are located.

 Returns
 -------
 str
     Local path at which tudat Earth orientation resources are located.






     )doc" );
    m.def( "get_quadrature_path",
           &tudat::paths::getQuadratureDataPath,
           R"doc(

 Get the path at which the Gaussian quadrature resources are located.

 Returns
 -------
 str
     Local path at which tudat Gaussian quadrature resources are located.






     )doc" );
    m.def( "get_spice_kernel_path",
           &tudat::paths::getSpiceKernelPath,
           R"doc(

 Get the path at which the SPICE kernel used by tudat is located.

 Returns
 -------
 str
     Local path at which the SPICE kernel is located.






     )doc" );
    m.def( "get_atmosphere_tables_path",
           &tudat::paths::getAtmosphereTablesPath,
           R"doc(

 Get the path at which tudat atmosphere tables are located.

 Returns
 -------
 str
     Local path at which tudat atmosphere tables are located.






     )doc" );
    m.def( "get_gravity_models_path",
           &tudat::paths::getGravityModelsPath,
           R"doc(

 Get the path at which tudat gravity models are located.

 Returns
 -------
 str
     Local path at which tudat gravity models are located.






     )doc" );
    m.def( "get_space_weather_path",
           &tudat::paths::getSpaceWeatherDataPath,
           R"doc(

 Get the path at which tudat space weather is located.

 Returns
 -------
 str
     Local path at which tudat space weather is located.






     )doc" );



    {
        py::module_ comaModel = m.def_submodule( "coma_model" );
        py::module_& m = comaModel;

    // === Coma processing: datasets (minimal shells so Python can hold them) ===
    py::class_< tss::ComaPolyDataset >( m,
                                        "ComaPolyDataset",
                                        R"doc(

 Polynomial coefficient dataset for coma atmosphere density model.

 This class holds polynomial coefficients that describe the spherical harmonic expansion of gas
 density in a cometary coma. The coefficients represent the spatial variation of density as a
 function of position relative to the comet nucleus. The dataset can contain data from multiple
 files, each covering a different time period, enabling time-dependent coma modeling.

 Polynomial coefficients are the raw input format and need to be evaluated to produce Stokes
 coefficients for use in the coma model. This evaluation is typically performed automatically
 when creating a coma atmosphere model.

 .. note:: This class cannot be directly instantiated by the user. Create instances using
           :func:`~tudatpy.data.coma_model.coma_model_file_processor`
           and its ``create_poly_coefficient_dataset()`` method.


 Examples
 --------
 Create a polynomial dataset from coefficient files:

 .. code-block:: python

   # Define paths to polynomial coefficient files
   file_paths = [
       "coma_data/h2o_poly_epoch1.txt",
       "coma_data/h2o_poly_epoch2.txt"
   ]

   # Create file processor
   processor = data.coma_model.coma_model_file_processor(file_paths)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Use the dataset to create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
       poly_data=poly_dataset,
       molecular_weight=0.018015)  # H2O molecular weight


      )doc" );

    py::class_< tss::ComaStokesDataset >( m,
                                          "ComaStokesDataset",
                                          R"doc(

 Stokes coefficient (spherical harmonics) dataset for coma atmosphere density model.

 This class holds precomputed Stokes coefficients (spherical harmonic coefficients) that describe
 the gas density distribution in a cometary coma. The coefficients are evaluated at a specific grid
 of radii and solar longitudes, enabling efficient interpolation during simulation. The dataset can
 contain data from multiple files, each covering a different time period.

 Stokes coefficients provide a direct representation of the spherical harmonic expansion:

 .. math::
     \\rho(r, \\theta, \\phi) = \\sum_{n=0}^{N} \\sum_{m=0}^{n} [C_{nm}(r) \\cos(m\\phi) + S_{nm}(r) \\sin(m\\phi)] P_{nm}(\\cos\\theta)

 where :math:`C_{nm}` and :math:`S_{nm}` are the Stokes coefficients, :math:`P_{nm}` are associated
 Legendre polynomials, and :math:`(r, \\theta, \\phi)` are spherical coordinates.

 .. note:: This class cannot be directly instantiated by the user. Create instances using
           :func:`~tudatpy.data.coma_model.coma_model_file_processor`
           and its ``create_coma_stokes_dataset()`` method, or load from pre-existing CSV files.


 Examples
 --------
 Create a Stokes dataset by transforming polynomial coefficients:

 .. code-block:: python

   # Create file processor from polynomial files
   processor = data.coma_model.coma_model_file_processor(
       ["coma_data/h2o_poly.txt"])

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use the dataset to create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

 Alternatively, load from pre-existing Stokes CSV files:

 .. code-block:: python

   # Create file processor from existing Stokes CSV files
   processor = data.coma_model.coma_model_file_processor(
       input_dir="coma_data/stokes_csv",
       prefix="stokes")

   # Load Stokes dataset from files
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])

   # Use the dataset
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc" );


    py::class_< tss::ComaWindDatasetCollection >( m,
                                                   "ComaWindDatasetCollection",
                                                   R"doc(
        Collection of three coma datasets for wind model (x, y, z components).

        This class holds three datasets (one for each spatial component of the wind velocity)
        that are used together to construct a ComaWindModel. All three datasets must be of
        the same type (either all polynomial or all Stokes coefficients).
        )doc" )
        .def("is_poly", &tss::ComaWindDatasetCollection::isPoly,
             R"doc(

 Check if the collection contains polynomial coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.data.coma_model.ComaPolyDataset` objects,
     False if it contains :class:`~tudatpy.data.coma_model.ComaStokesDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files
   processor = data.coma_model.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_poly_coefficient_datasets()

   # Check dataset type
   if wind_datasets.is_poly():
       print("Collection contains polynomial datasets")


      )doc")
        .def("is_stokes", &tss::ComaWindDatasetCollection::isStokes,
             R"doc(

 Check if the collection contains Stokes coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.data.coma_model.ComaStokesDataset` objects,
     False if it contains :class:`~tudatpy.data.coma_model.ComaPolyDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files and convert to Stokes
   processor = data.coma_model.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Check dataset type
   if wind_datasets.is_stokes():
       print("Collection contains Stokes datasets")


      )doc");


    // === Coma processing: file processor ===
    py::class_< tss::ComaModelFileProcessor, std::shared_ptr< tss::ComaModelFileProcessor > >( m,
                                               "ComaModelFileProcessor",
                                               R"doc(

 Processor for creating coma atmosphere datasets from coefficient files.

 This class provides a high-level interface for loading and processing coma model data from files
 containing either polynomial coefficients or pre-computed Stokes (spherical harmonic) coefficients.
 It handles two distinct workflows:

 1. **Polynomial coefficient workflow**: Load polynomial coefficients from text files, optionally
    transform them to Stokes coefficients at specified radii and solar longitudes, and create datasets
    for use in coma atmosphere models.

 2. **Stokes coefficient workflow**: Load pre-computed Stokes coefficients from CSV files that were
    previously generated and saved.

 The processor automatically manages file reading, data validation, and coordinate transformations,
 simplifying the setup of complex coma atmosphere models.

 .. note:: This class cannot be directly instantiated. Create instances using the factory functions
           :func:`~tudatpy.data.coma_model.coma_model_file_processor`.


 Examples
 --------
 Polynomial coefficient workflow:

 .. code-block:: python

   # Create processor from polynomial coefficient files
   poly_files = ["h2o_epoch1.txt", "h2o_epoch2.txt"]
   processor = data.coma_model.coma_model_file_processor(poly_files)

   # Create polynomial dataset directly
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Or transform to Stokes coefficients
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

 Stokes coefficient workflow:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   processor = data.coma_model.coma_model_file_processor(
       input_dir="coma_stokes_data",
       prefix="stokes")

   # Load Stokes dataset
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc" )
        .def("create_poly_coefficient_dataset",
             &tss::ComaModelFileProcessor::createPolyCoefDataset,
             R"doc(

 Create polynomial coefficient dataset from the loaded files.

 Reads and processes polynomial coefficient files to create a :class:`~tudatpy.data.coma_model.ComaPolyDataset`.
 This method is only available when the processor was constructed from polynomial coefficient files
 using the file path variant of :func:`~tudatpy.data.coma_model.coma_model_file_processor`.

 Returns
 -------
 ComaPolyDataset
     Dataset containing polynomial coefficients for all loaded files.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes coefficient CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = data.coma_model.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
       poly_data=poly_dataset,
       molecular_weight=0.018015)


      )doc")
        .def("create_coma_stokes_dataset",
             py::overload_cast<>(&tss::ComaModelFileProcessor::createSHDataset, py::const_),
             R"doc(

 Create Stokes coefficient dataset from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.data.coma_model.coma_model_file_processor`.
 It returns the preloaded Stokes dataset directly.

 Returns
 -------
 ComaStokesDataset
     Dataset containing preloaded Stokes coefficients from CSV files.

 Raises
 ------
 RuntimeError
     If the processor was constructed from polynomial coefficient files. Use the parameterized version
     :meth:`create_coma_stokes_dataset(radii_m, sol_longitudes_deg, ...)` for polynomial files instead.

 Examples
 --------
 .. code-block:: python

   # Create processor from Stokes CSV files
   processor = data.coma_model.coma_model_file_processor(
       input_dir="stokes_data",
       prefix="stokes")

   # Load Stokes dataset (no parameters needed)
   stokes_dataset = processor.create_coma_stokes_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

      )doc")
        .def("create_coma_stokes_dataset",
             py::overload_cast<const std::vector<double>&, const std::vector<double>&, const int, const int, const bool, const bool>(
                 &tss::ComaModelFileProcessor::createSHDataset, py::const_),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             py::arg("compute_reduced_coeffs") = true,
             py::arg("is_log2") = true,
             R"doc(

 Create Stokes coefficient dataset by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.data.coma_model.coma_model_file_processor`. It transforms
 polynomial coefficients to Stokes coefficients by evaluating the spherical harmonic expansion at the
 specified grid of radii and solar longitudes.

 Parameters
 ----------
 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include. Set to -1 to use maximum available.

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include. Set to -1 to use maximum available.

 is_log2 : bool, default = True
     Whether the polynomial coefficients represent log2-transformed data.
     This affects how the 1/r² decay term is applied for radii beyond the reference radius.
     Set to False if the coefficients were fitted to non-log2-transformed (linear) data.

 Returns
 -------
 ComaStokesDataset
     Dataset containing Stokes coefficients.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files. Use the parameterless version
     :meth:`create_coma_stokes_dataset()` for Stokes CSV files instead.

 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = data.coma_model.coma_model_file_processor(
       ["h2o_poly.txt"])

   # Transform to Stokes at specific radii and solar longitudes
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0, 20000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

      )doc")
        .def("create_coma_spherical_harmonic_coefficient_files",
             &tss::ComaModelFileProcessor::createSHFiles,
             py::arg("output_dir"),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             py::arg("compute_reduced_coeffs") = true,
             py::arg("is_log2") = true,
             R"doc(

 Create and save Stokes coefficient CSV files from polynomial coefficients.

 Transforms polynomial coefficients to Stokes coefficients and saves them as CSV files in the specified
 output directory. This is useful for pre-computing Stokes coefficients to avoid repeated transformations
 during multiple simulation runs. The generated CSV files can later be loaded using a processor created
 with the directory-based variant of :func:`~tudatpy.data.coma_model.coma_model_file_processor`.

 This method is only available when the processor was constructed from polynomial coefficient files.


 Parameters
 ----------
 output_dir : str
     Directory path where Stokes coefficient CSV files will be saved. The directory will be created
     if it does not exist.

 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include. Set to -1 to automatically use the maximum
     degree available in the polynomial data.

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include. Set to -1 to automatically use the maximum
     order available in the polynomial data.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create processor from polynomial files
   processor = data.coma_model.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Pre-compute and save Stokes coefficients to CSV files
   processor.create_coma_spherical_harmonic_coefficient_files(
       output_dir="stokes_output",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Later, load from the saved CSV files
   processor_from_csv = data.coma_model.coma_model_file_processor(
       input_dir="stokes_output",
       prefix="stokes")
   stokes_dataset = processor_from_csv.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc");

    m.def(
            "coma_model_file_processor",
            py::overload_cast<const std::vector<std::string>&>(&tss::comaModelFileProcessorFromPolyFiles),
            py::arg( "file_paths" ),
            R"doc(

 Function for creating coma model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.data.coma_model.ComaModelFileProcessor` that loads
 and processes polynomial coefficient files. These files contain spherical harmonic coefficients in
 polynomial form that describe gas density distributions in a cometary coma.

 The processor can create polynomial datasets directly, or transform them to Stokes coefficients at
 specified radii and solar longitudes. Multiple files can be provided to cover different time periods,
 enabling time-dependent coma modeling.


 Parameters
 ----------
 file_paths : list[str]
     List of file paths to polynomial coefficient files. Each file should contain polynomial coefficients
     for spherical harmonic expansion of gas density. Files may cover different time periods.

 Returns
 -------
 ComaModelFileProcessor
     Processor configured for polynomial coefficient files, capable of creating both polynomial and
     Stokes datasets.

 Raises
 ------
 ValueError
     If file_paths is an empty list.
 RuntimeError
     If any of the specified files do not exist or cannot be opened.


 Examples
 --------
 Create processor and use polynomial dataset directly:

 .. code-block:: python

   # Define paths to polynomial coefficient files
   poly_files = [
       "coma_data/h2o_poly_epoch1.txt",
       "coma_data/h2o_poly_epoch2.txt"
   ]

   # Create file processor
   processor = data.coma_model.coma_model_file_processor(poly_files)

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
       poly_data=poly_dataset,
       molecular_weight=0.018015,  # H2O in kg/mol
       max_degree=10,
       max_order=10)

   # Apply to body
   body_settings.get("67P").atmosphere_settings = coma_settings

 Transform to Stokes coefficients:

 .. code-block:: python

   # Create processor
   processor = data.coma_model.coma_model_file_processor(poly_files)

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc"
            );

    m.def(
            "coma_model_file_processor",
            py::overload_cast<const std::string&, const std::string&>(&tss::comaModelFileProcessorFromSHFiles),
            py::arg( "input_dir" ),
            py::arg( "prefix" ) = "stokes",
            R"doc(

 Function for creating coma model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.data.coma_model.ComaModelFileProcessor` that loads
 pre-computed Stokes (spherical harmonic) coefficients from CSV files. These files typically contain
 Stokes coefficients that were previously generated and saved using the ``create_coma_spherical_harmonic_coefficient_files()`` method
 of a processor created from polynomial files.

 This variant is useful when you want to avoid repeated transformation of polynomial coefficients to
 Stokes coefficients across multiple simulation runs. The CSV files contain pre-evaluated coefficients
 at a fixed grid of radii and solar longitudes.


 Parameters
 ----------
 input_dir : str
     Directory path containing the Stokes coefficient CSV files to load. All CSV files in this
     directory with the specified prefix will be loaded.

 prefix : str, default = "stokes"
     File name prefix for the CSV files. Files are expected to be named as ``{prefix}_0.csv``,
     ``{prefix}_1.csv``, etc.

 Returns
 -------
 ComaModelFileProcessor
     Processor configured for Stokes coefficient CSV files, capable of loading pre-computed
     Stokes datasets.

 Raises
 ------
 RuntimeError
     If the directory does not exist, is not a directory, or contains no matching CSV files.


 Examples
 --------
 Load from previously saved Stokes CSV files:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   processor = data.coma_model.coma_model_file_processor(
       input_dir="coma_data/stokes_precomputed",
       prefix="stokes")

   # Load Stokes dataset (radii and longitudes are read from files)
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from CSV files
       sol_longitudes_deg=[])  # Ignored when loading from CSV files

   # Create coma atmosphere settings
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)  # H2O in kg/mol

   # Apply to body
   body_settings.get("67P").atmosphere_settings = coma_settings

 Complete workflow from polynomial to saved Stokes to loading:

 .. code-block:: python

   # Step 1: Create and save Stokes coefficients from polynomial files
   poly_processor = data.coma_model.coma_model_file_processor(
       ["h2o_poly.txt"])
   poly_processor.create_coma_spherical_harmonic_coefficient_files(
       output_dir="stokes_saved",
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Step 2: Later, load from saved Stokes files
   stokes_processor = data.coma_model.coma_model_file_processor(
       input_dir="stokes_saved",
       prefix="stokes")
   stokes_dataset = stokes_processor.create_coma_stokes_dataset(
       radii_m=[],
       sol_longitudes_deg=[])

   # Step 3: Use in simulation
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc"
            );

    // === Coma wind processing: file processor ===
    py::class_< tss::ComaWindModelFileProcessor, std::shared_ptr< tss::ComaWindModelFileProcessor > >( m,
                                               "ComaWindModelFileProcessor",
                                               R"doc(
        Processor for creating wind model datasets from three component file sources.

        This class manages the creation of ComaWindDatasetCollection from three sets of files
        (one for each spatial component: x, y, z). It provides a simplified interface for
        wind model setup by handling all three components together.
        )doc" )
        .def("create_poly_coefficient_datasets",
             &tss::ComaWindModelFileProcessor::createPolyCoefDatasets,
             R"doc(

 Create polynomial coefficient dataset collection for all three wind components.

 Reads and processes polynomial coefficient files for x, y, and z wind velocity components to create
 a :class:`~tudatpy.data.coma_model.ComaWindDatasetCollection`. This method is
 only available when the processor was constructed from polynomial coefficient files.

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing x, y, z polynomial datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If processor was constructed from Stokes coefficient CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Define paths to polynomial coefficient files for each component
   x_files = ["wind_x_epoch1.txt", "wind_x_epoch2.txt"]
   y_files = ["wind_y_epoch1.txt", "wind_y_epoch2.txt"]
   z_files = ["wind_z_epoch1.txt", "wind_z_epoch2.txt"]

   # Create wind file processor
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_files, y_files, z_files)

   # Create polynomial dataset collection
   poly_datasets = wind_processor.create_poly_coefficient_datasets()


      )doc")
        .def("create_coma_stokes_dataset",
             py::overload_cast<>(&tss::ComaWindModelFileProcessor::createSHDatasets, py::const_),
             R"doc(

 Create Stokes coefficient dataset collection from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.data.coma_model.coma_wind_file_processor`.
 It returns the preloaded Stokes datasets for all three wind components (x, y, z) directly.

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing preloaded x, y, z Stokes datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If the processor was constructed from polynomial coefficient files. Use the parameterized version
     :meth:`create_coma_stokes_dataset(radii_m, sol_longitudes_deg, ...)` for polynomial files instead.

 Examples
 --------
 .. code-block:: python

   # Create wind processor from Stokes CSV files
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_input_dir="wind_x_stokes",
       y_input_dir="wind_y_stokes",
       z_input_dir="wind_z_stokes",
       prefix="stokes")

   # Load Stokes datasets (no parameters needed)
   stokes_datasets = wind_processor.create_coma_stokes_dataset()

   # Use in coma wind model
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=stokes_datasets,
       requested_max_degree=10,
       requested_max_order=10)

      )doc")
        .def("create_coma_stokes_dataset",
             py::overload_cast<const std::vector<double>&, const std::vector<double>&, const int, const int>(
                 &tss::ComaWindModelFileProcessor::createSHDatasets, py::const_),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create Stokes coefficient dataset collection by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.data.coma_model.coma_wind_file_processor`. It transforms
 polynomial coefficients to Stokes coefficients for all three components (x, y, z) by evaluating
 at the specified grid of radii and solar longitudes.

 Parameters
 ----------
 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include (-1 for automatic determination).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include (-1 for automatic determination).

 Returns
 -------
 ComaWindDatasetCollection
     Collection containing x, y, z Stokes datasets for wind velocity components.

 Raises
 ------
 RuntimeError
     If the processor was constructed from Stokes CSV files. Use the parameterless version
     :meth:`create_coma_stokes_dataset()` for Stokes CSV files instead.

 Examples
 --------
 .. code-block:: python

   # Create wind processor from polynomial files
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths=["wind_x.txt"],
       y_file_paths=["wind_y.txt"],
       z_file_paths=["wind_z.txt"])

   # Transform to Stokes coefficients at specific radii and solar longitudes
   stokes_datasets = wind_processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Use in coma wind model
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=stokes_datasets,
       requested_max_degree=10,
       requested_max_order=10)

      )doc")
        .def("create_coma_spherical_harmonic_coefficient_files",
             &tss::ComaWindModelFileProcessor::createSHFiles,
             py::arg("x_output_dir"),
             py::arg("y_output_dir"),
             py::arg("z_output_dir"),
             py::arg("radii_m"),
             py::arg("sol_longitudes_deg"),
             py::arg("requested_max_degree") = -1,
             py::arg("requested_max_order") = -1,
             R"doc(

 Create and save Stokes coefficient CSV files for all three wind components.

 Transforms polynomial coefficients to Stokes coefficients for x, y, and z wind velocity components
 and saves them as CSV files in separate output directories. This is useful for pre-computing Stokes
 coefficients to avoid repeated transformations across multiple simulation runs.

 This method is only available when the processor was constructed from polynomial coefficient files.


 Parameters
 ----------
 x_output_dir : str
     Directory path where x-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 y_output_dir : str
     Directory path where y-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 z_output_dir : str
     Directory path where z-component Stokes coefficient CSV files will be saved. The directory
     will be created if it does not exist.

 radii_m : list[float]
     Vector of radii at which to evaluate Stokes coefficients [m].

 sol_longitudes_deg : list[float]
     Vector of solar longitudes at which to evaluate Stokes coefficients [degrees].

 requested_max_degree : int, default = -1
     Maximum spherical harmonic degree to include (-1 for automatic determination from data).

 requested_max_order : int, default = -1
     Maximum spherical harmonic order to include (-1 for automatic determination from data).

 Raises
 ------
 RuntimeError
     If processor was constructed from Stokes CSV files instead of polynomial files.


 Examples
 --------
 .. code-block:: python

   # Create wind processor from polynomial files
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths=["wind_x_epoch1.txt"],
       y_file_paths=["wind_y_epoch1.txt"],
       z_file_paths=["wind_z_epoch1.txt"])

   # Pre-compute and save Stokes coefficients for all components
   wind_processor.create_coma_spherical_harmonic_coefficient_files(
       x_output_dir="stokes_wind/x_component",
       y_output_dir="stokes_wind/y_component",
       z_output_dir="stokes_wind/z_component",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Later, load from the saved CSV files
   wind_processor_from_csv = data.coma_model.coma_wind_file_processor(
       x_input_dir="stokes_wind/x_component",
       y_input_dir="stokes_wind/y_component",
       z_input_dir="stokes_wind/z_component",
       prefix="stokes")


      )doc")
        .def("is_poly_type", &tss::ComaWindModelFileProcessor::isPolyType,
             R"doc(Check if this processor works with polynomial coefficient files.)doc")
        .def("is_stokes_type", &tss::ComaWindModelFileProcessor::isStokesType,
             R"doc(Check if this processor works with Stokes coefficient files.)doc");

    m.def(
            "coma_wind_file_processor",
            py::overload_cast<const std::vector<std::string>&,
                              const std::vector<std::string>&,
                              const std::vector<std::string>&>(&tss::comaWindModelFileProcessorFromPolyFiles),
            py::arg( "x_file_paths" ),
            py::arg( "y_file_paths" ),
            py::arg( "z_file_paths" ),
            R"doc(

 Function for creating coma wind model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.data.coma_model.ComaWindModelFileProcessor` that loads
 and processes polynomial coefficient files for all three wind velocity components (x, y, z). These files
 contain spherical harmonic coefficients in polynomial form that describe wind velocity distributions
 in a cometary coma.

 The processor can create polynomial dataset collections directly, or transform them to Stokes
 coefficients at specified radii and solar longitudes. Multiple files can be provided for each component
 to cover different time periods.


 Parameters
 ----------
 x_file_paths : list[str]
     List of file paths for x-component polynomial coefficients. Files may cover different time periods.

 y_file_paths : list[str]
     List of file paths for y-component polynomial coefficients. Files may cover different time periods.

 z_file_paths : list[str]
     List of file paths for z-component polynomial coefficients. Files may cover different time periods.

 Returns
 -------
 ComaWindModelFileProcessor
     Processor configured for polynomial coefficient files, capable of creating both polynomial and
     Stokes dataset collections.

 Raises
 ------
 ValueError
     If any of the file path lists is empty.
 RuntimeError
     If any of the specified files do not exist or cannot be opened.


 Examples
 --------
 Create wind model from polynomial files:

 .. code-block:: python

   # Define paths to polynomial coefficient files for each wind component
   x_files = ["wind_data/vx_epoch1.txt", "wind_data/vx_epoch2.txt"]
   y_files = ["wind_data/vy_epoch1.txt", "wind_data/vy_epoch2.txt"]
   z_files = ["wind_data/vz_epoch1.txt", "wind_data/vz_epoch2.txt"]

   # Create wind file processor
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths=x_files,
       y_file_paths=y_files,
       z_file_paths=z_files)

   # Transform to Stokes coefficients
   wind_datasets = wind_processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Create coma wind model settings
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets,
       requested_max_degree=10,
       requested_max_order=10,
       associated_reference_frame=environment.AerodynamicsReferenceFrames.vertical_frame,
       include_corotation=True)

   # Apply to atmosphere
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind

 Pre-compute and save Stokes coefficients:

 .. code-block:: python

   # Create processor
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths=x_files,
       y_file_paths=y_files,
       z_file_paths=z_files)

   # Save Stokes coefficients to CSV files for later use
   wind_processor.create_coma_spherical_harmonic_coefficient_files(
       x_output_dir="stokes_wind/x",
       y_output_dir="stokes_wind/y",
       z_output_dir="stokes_wind/z",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)


      )doc"
            );

    m.def(
            "coma_wind_file_processor",
            py::overload_cast<const std::string&,
                              const std::string&,
                              const std::string&,
                              const std::string&>(&tss::comaWindModelFileProcessorFromSHFiles),
            py::arg( "x_input_dir" ),
            py::arg( "y_input_dir" ),
            py::arg( "z_input_dir" ),
            py::arg( "prefix" ) = "stokes",
            R"doc(

 Function for creating coma wind model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.data.coma_model.ComaWindModelFileProcessor` that loads
 pre-computed Stokes (spherical harmonic) coefficients from CSV files for all three wind velocity
 components (x, y, z). These files typically contain Stokes coefficients that were previously generated
 and saved using the ``create_coma_spherical_harmonic_coefficient_files()`` method of a processor created from polynomial files.

 This variant is useful when you want to avoid repeated transformation of polynomial coefficients to
 Stokes coefficients across multiple simulation runs. The CSV files contain pre-evaluated coefficients
 at a fixed grid of radii and solar longitudes.


 Parameters
 ----------
 x_input_dir : str
     Directory path containing x-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 y_input_dir : str
     Directory path containing y-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 z_input_dir : str
     Directory path containing z-component Stokes coefficient CSV files. All CSV files in this
     directory with the specified prefix will be loaded.

 prefix : str, default = "stokes"
     File name prefix for the CSV files. Files are expected to be named as ``{prefix}_0.csv``,
     ``{prefix}_1.csv``, etc.

 Returns
 -------
 ComaWindModelFileProcessor
     Processor configured for Stokes coefficient CSV files, capable of loading pre-computed
     Stokes dataset collections.

 Raises
 ------
 RuntimeError
     If any of the directories do not exist, are not directories, or contain no matching CSV files.


 Examples
 --------
 Load from previously saved Stokes CSV files:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   wind_processor = data.coma_model.coma_wind_file_processor(
       x_input_dir="wind_stokes/x_component",
       y_input_dir="wind_stokes/y_component",
       z_input_dir="wind_stokes/z_component",
       prefix="stokes")

   # Load Stokes datasets (radii and longitudes are read from files)
   wind_datasets = wind_processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from CSV files
       sol_longitudes_deg=[])  # Ignored when loading from CSV files

   # Create coma wind model settings
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets,
       associated_reference_frame=environment.AerodynamicsReferenceFrames.vertical_frame)

   # Apply to atmosphere
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind

 Complete workflow from polynomial to saved Stokes to loading:

 .. code-block:: python

   # Step 1: Create and save Stokes coefficients from polynomial files
   poly_processor = data.coma_model.coma_wind_file_processor(
       x_file_paths=["vx.txt"],
       y_file_paths=["vy.txt"],
       z_file_paths=["vz.txt"])
   poly_processor.create_coma_spherical_harmonic_coefficient_files(
       x_output_dir="wind_stokes/x",
       y_output_dir="wind_stokes/y",
       z_output_dir="wind_stokes/z",
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Step 2: Later, load from saved Stokes files
   stokes_processor = data.coma_model.coma_wind_file_processor(
       x_input_dir="wind_stokes/x",
       y_input_dir="wind_stokes/y",
       z_input_dir="wind_stokes/z",
       prefix="stokes")
   wind_datasets = stokes_processor.create_coma_stokes_dataset(
       radii_m=[],
       sol_longitudes_deg=[])

   # Step 3: Use in simulation
   coma_wind = environment_setup.atmosphere.coma_wind_model(
       dataset_collection=wind_datasets)
   body_settings.get("67P").atmosphere_settings.wind_settings = coma_wind


      )doc"
            );




    }

    m.def( "read_vector_history_from_file",
           &tudat::input_output::readVectorHistoryFromFile< double, double >,
           py::arg( "vector_size" ),
           py::arg( "file_name" ),
           R"doc(

 Read a vector history from a file.


 Parameters
 ----------
 vector_size : int
     Size of the vector at each epoch.
 file_name : str
     Name of the file containing the vector history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the vector at the given epoch.






     )doc" );

    m.def( "read_matrix_history_from_file",
           &tudat::input_output::readMatrixHistoryFromFile< double, double >,
           py::arg( "matrix_rows" ),
           py::arg( "matrix_columns" ),
           py::arg( "file_name" ),
           R"doc(

 Read a matrix history from a file.


 Parameters
 ----------
 matrix_rows : int
     Number of rows in the matrix at each epoch.
 matrix_columns : int
     Number of columns in the matrix at each epoch.
 file_name : str
     Name of the file containing the matrix history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the matrix at the given epoch.






     )doc" );

    py::class_< tio::CrdPassConfigurationData >( m, "CrdPassConfigurationData", R"doc(
Container for CRD pass-level metadata.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_name", &tio::CrdPassConfigurationData::stationName_ )
            .def_readwrite( "cdp_pad_id", &tio::CrdPassConfigurationData::cdpPadId_ )
            .def_readwrite( "target_name", &tio::CrdPassConfigurationData::targetName_ )
            .def_readwrite( "start_year", &tio::CrdPassConfigurationData::startYear_ )
            .def_readwrite( "start_month", &tio::CrdPassConfigurationData::startMonth_ )
            .def_readwrite( "start_day", &tio::CrdPassConfigurationData::startDay_ )
            .def_readwrite( "end_year", &tio::CrdPassConfigurationData::endYear_ )
            .def_readwrite( "end_month", &tio::CrdPassConfigurationData::endMonth_ )
            .def_readwrite( "end_day", &tio::CrdPassConfigurationData::endDay_ )
            .def_readwrite( "transmit_wavelength_nm", &tio::CrdPassConfigurationData::transmitWavelengthNm_ );

    py::class_< tio::CrdNormalPointRecord >( m, "CrdNormalPointRecord", R"doc(
Container for a CRD normal-point (record ``11``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdNormalPointRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdNormalPointRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdNormalPointRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdNormalPointRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdNormalPointRecord::epochEvent_ )
            .def_readwrite( "normal_point_window_length", &tio::CrdNormalPointRecord::normalPointWindowLength_ )
            .def_readwrite( "number_of_returns", &tio::CrdNormalPointRecord::numberOfReturns_ )
            .def_readwrite( "bin_rms", &tio::CrdNormalPointRecord::binRms_ );

    py::class_< tio::CrdFullRateRecord >( m, "CrdFullRateRecord", R"doc(
Container for a CRD full-rate observation (record ``10``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdFullRateRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdFullRateRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdFullRateRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdFullRateRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdFullRateRecord::epochEvent_ )
            .def_readwrite( "filter_flag", &tio::CrdFullRateRecord::filterFlag_ )
            .def_readwrite( "detector_channel", &tio::CrdFullRateRecord::detectorChannel_ )
            .def_readwrite( "stop_number", &tio::CrdFullRateRecord::stopNumber_ );

    py::class_< tio::CrdMeteoRecord >( m, "CrdMeteoRecord", R"doc(
Container for a CRD meteorological record (record ``20``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdMeteoRecord::secondOfDay_ )
            .def_readwrite( "pressure", &tio::CrdMeteoRecord::pressure_ )
            .def_readwrite( "temperature", &tio::CrdMeteoRecord::temperature_ )
            .def_readwrite( "humidity", &tio::CrdMeteoRecord::humidity_ );

    py::class_< tio::CrdPassData >( m, "CrdPassData", R"doc(
Container for the measurement and meteorological data of a CRD pass.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "full_rate_data", &tio::CrdPassData::fullRateData_ )
            .def_readwrite( "normal_point_data", &tio::CrdPassData::normalPointData_ )
            .def_readwrite( "meteorological_data", &tio::CrdPassData::meteorologicalData_ );

    py::class_< tio::CrdPass >( m, "CrdPass", R"doc(
Container for a CRD pass, including configuration and data records.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "configuration", &tio::CrdPass::configuration_ )
            .def_readwrite( "data", &tio::CrdPass::data_ );

    m.def( "convert_crd_two_way_time_of_flight_to_slr_range",
           &tio::convertCrdTwoWayTimeOfFlightToSlrRange,
           py::arg( "two_way_time_of_flight" ),
           R"doc(
Convert a CRD two-way light time (seconds) into one-way SLR range (meters).
           )doc" );

    m.def( "read_crd_file",
           &tio::readCrdFile,
           py::arg( "file_name" ),
           R"doc(
Read a CRD file and return the parsed pass list.
           )doc" );

    m.def( "read_crd_files",
           &tio::readCrdFiles,
           py::arg( "file_names" ),
           R"doc(
Read multiple CRD files and concatenate all parsed passes.
           )doc" );

    m.def( "group_crd_data_per_target",
           &tio::groupCrdDataPerTarget,
           py::arg( "crd_passes" ),
           R"doc(
Group CRD passes by target name.
           )doc" );

    m.def( "group_crd_data_per_station",
           &tio::groupCrdDataPerStation,
           py::arg( "crd_passes" ),
           py::arg( "monument_id_to_ground_station_name_map" ),
           R"doc(
Group CRD passes by station name using a monument-ID to station-name map.
           )doc" );

    m.def( "extract_normal_point_measurements",
           py::overload_cast< const tio::CrdPass& >( &tio::extractNormalPointMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract normal-point measurements from one CRD pass as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_normal_point_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractNormalPointMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract normal-point measurements from multiple CRD passes as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_full_rate_measurements",
           py::overload_cast< const tio::CrdPass& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract full-rate measurements from one CRD pass as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_full_rate_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract full-rate measurements from multiple CRD passes as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "get_station_wavelengths",
           &tio::getStationWavelengths,
           py::arg( "grouped_data" ),
           R"doc(
Extract station transmit wavelengths from grouped CRD pass data.
           )doc" );

    py::class_< tio::SinexStationState >( m, "SinexStationState", R"doc(
Container for station state data parsed from a SINEX file.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "site_code", &tio::SinexStationState::siteCode_ )
            .def_readwrite( "domes_id", &tio::SinexStationState::domesId_ )
            .def_readwrite( "position", &tio::SinexStationState::position_ )
            .def_readwrite( "velocity", &tio::SinexStationState::velocity_ )
            .def_readwrite( "reference_epoch", &tio::SinexStationState::referenceEpoch_ );

    py::class_< tio::SinexStationEccentricity >( m, "SinexStationEccentricity", R"doc(
Container for SINEX station eccentricity data.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "domes_id", &tio::SinexStationEccentricity::domesId_ )
            .def_readwrite( "station_code", &tio::SinexStationEccentricity::stationCode_ )
            .def_readwrite( "eccentricity", &tio::SinexStationEccentricity::eccentricity_ )
            .def_readwrite( "start_epoch", &tio::SinexStationEccentricity::startEpoch_ )
            .def_readwrite( "end_epoch", &tio::SinexStationEccentricity::endEpoch_ )
            .def_readwrite( "has_open_end", &tio::SinexStationEccentricity::hasOpenEnd_ );

    py::class_< tio::IlrsStationRegistryEntry >( m, "IlrsStationRegistryEntry", R"doc(
Container for one ILRS station registry entry parsed from SINEX ``SITE/ID``.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_code", &tio::IlrsStationRegistryEntry::stationCode_ )
            .def_readwrite( "station_name", &tio::IlrsStationRegistryEntry::stationName_ )
            .def_readwrite( "domes_id", &tio::IlrsStationRegistryEntry::domesId_ )
            .def_readwrite( "approximate_longitude", &tio::IlrsStationRegistryEntry::approximateLongitude_ )
            .def_readwrite( "approximate_latitude", &tio::IlrsStationRegistryEntry::approximateLatitude_ )
            .def_readwrite( "approximate_height", &tio::IlrsStationRegistryEntry::approximateHeight_ );

    m.def( "convert_sinex_datetime_to_seconds_since_epoch",
           &tio::convertSinexDateTimeToSecondsSinceEpoch,
           py::arg( "date_time" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Convert a SINEX epoch string (``YY:DDD:SSSSS``) to seconds since a reference Julian day.
           )doc" );

    m.def( "read_sinex_station_data",
           &tio::readSinexStationData,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station positions, velocities and reference epochs from a SINEX file.
           )doc" );

    m.def( "read_sinex_station_eccentricities",
           &tio::readSinexStationEccentricities,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station eccentricity history from a SINEX file.
           )doc" );

    m.def( "read_ilrs_station_registry_from_sinex_site_id",
           &tio::readIlrsStationRegistryFromSinexSiteId,
           py::arg( "file_name" ),
           R"doc(
Read the ILRS station registry from a SINEX ``SITE/ID`` block as ``{station_code: entry}``.
           )doc" );

    m.def( "read_domes_id_numbers",
           &tio::readDomesIdNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from station name to DOMES id.
           )doc" );

    m.def( "read_monument_numbers",
           &tio::readMonumentNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from monument/station code to station name.
           )doc" );

    m.def( "read_ground_station_names",
           &tio::readGroundStationNames,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from DOMES id to station name.
           )doc" );

    py::enum_< tudat::input_output::TrackingDataType >( m, "TrackingDataType", R"doc(No documentation available.)doc" )
            .value( "year", tudat::input_output::TrackingDataType::year, R"doc(No documentation available.)doc" )
            .value( "month", tudat::input_output::TrackingDataType::month, R"doc(No documentation available.)doc" )
            .value( "day", tudat::input_output::TrackingDataType::day, R"doc(No documentation available.)doc" )
            .value( "hour", tudat::input_output::TrackingDataType::hour, R"doc(No documentation available.)doc" )
            .value( "minute", tudat::input_output::TrackingDataType::minute, R"doc(No documentation available.)doc" )
            .value( "second", tudat::input_output::TrackingDataType::second, R"doc(No documentation available.)doc" )
            .value( "observation_time_scale",
                    tudat::input_output::TrackingDataType::observation_time_scale,
                    R"doc(No documentation available.)doc" )
            .value( "file_name", tudat::input_output::TrackingDataType::file_name, R"doc(No documentation available.)doc" )
            .value( "n_way_light_time", tudat::input_output::TrackingDataType::n_way_light_time, R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_delay",
                    tudat::input_output::TrackingDataType::light_time_measurement_delay,
                    R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_accuracy",
                    tudat::input_output::TrackingDataType::light_time_measurement_accuracy,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_transmitting_station_nr",
                    tudat::input_output::TrackingDataType::dsn_transmitting_station_nr,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_receiving_station_nr",
                    tudat::input_output::TrackingDataType::dsn_receiving_station_nr,
                    R"doc(No documentation available.)doc" )
            .value( "observation_body", tudat::input_output::TrackingDataType::observation_body, R"doc(No documentation available.)doc" )
            .value( "observed_body", tudat::input_output::TrackingDataType::observed_body, R"doc(No documentation available.)doc" )
            .value( "spacecraft_id", tudat::input_output::TrackingDataType::spacecraft_id, R"doc(No documentation available.)doc" )
            .value( "spacecraft_name", tudat::input_output::TrackingDataType::spacecraft_name, R"doc(No documentation available.)doc" )
            .value( "planet_nr", tudat::input_output::TrackingDataType::planet_nr, R"doc(No documentation available.)doc" )
            .value( "tdb_reception_time_j2000",
                    tudat::input_output::TrackingDataType::tdb_reception_time_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "utc_reception_time_j2000",
                    tudat::input_output::TrackingDataType::utc_reception_time_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "utc_ramp_referencee_j2000",
                    tudat::input_output::TrackingDataType::utc_ramp_referencee_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "tdb_spacecraft_j2000",
                    tudat::input_output::TrackingDataType::tdb_spacecraft_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "x_planet_frame", tudat::input_output::TrackingDataType::x_planet_frame, R"doc(No documentation available.)doc" )
            .value( "y_planet_frame", tudat::input_output::TrackingDataType::y_planet_frame, R"doc(No documentation available.)doc" )
            .value( "z_planet_frame", tudat::input_output::TrackingDataType::z_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vx_planet_frame", tudat::input_output::TrackingDataType::vx_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vy_planet_frame", tudat::input_output::TrackingDataType::vy_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vz_planet_frame", tudat::input_output::TrackingDataType::vz_planet_frame, R"doc(No documentation available.)doc" )
            .value( "residual_de405", tudat::input_output::TrackingDataType::residual_de405, R"doc(No documentation available.)doc" )
            .value( "spacecraft_transponder_delay",
                    tudat::input_output::TrackingDataType::spacecraft_transponder_delay,
                    R"doc(No documentation available.)doc" )
            .value( "uplink_frequency", tudat::input_output::TrackingDataType::uplink_frequency, R"doc(No documentation available.)doc" )
            .value( "downlink_frequency",
                    tudat::input_output::TrackingDataType::downlink_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "signal_to_noise", tudat::input_output::TrackingDataType::signal_to_noise, R"doc(No documentation available.)doc" )
            .value( "spectral_max", tudat::input_output::TrackingDataType::spectral_max, R"doc(No documentation available.)doc" )
            .value( "doppler_measured_frequency",
                    tudat::input_output::TrackingDataType::doppler_measured_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_averaged_frequency",
                    tudat::input_output::TrackingDataType::doppler_averaged_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_base_frequency",
                    tudat::input_output::TrackingDataType::doppler_base_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_noise", tudat::input_output::TrackingDataType::doppler_noise, R"doc(No documentation available.)doc" )
            .value( "doppler_bandwidth", tudat::input_output::TrackingDataType::doppler_bandwidth, R"doc(No documentation available.)doc" )
            .value( "receiving_station_name",
                    tudat::input_output::TrackingDataType::receiving_station_name,
                    R"doc(No documentation available.)doc" )
            .value( "transmitting_station_name",
                    tudat::input_output::TrackingDataType::transmitting_station_name,
                    R"doc(No documentation available.)doc" )
            .value( "time_tag_delay", tudat::input_output::TrackingDataType::time_tag_delay )
            .value( "sample_number", tudat::input_output::TrackingDataType::sample_number )
            .value( "utc_day_of_year", tudat::input_output::TrackingDataType::utc_day_of_year )
            .value( "reference_body_distance", tudat::input_output::TrackingDataType::reference_body_distance )
            .value( "transmission_frequency_constant_term", tudat::input_output::TrackingDataType::transmission_frequency_constant_term )
            .value( "transmission_frequency_linear_term", tudat::input_output::TrackingDataType::transmission_frequency_linear_term )
            .value( "doppler_predicted_frequency_hz", tudat::input_output::TrackingDataType::doppler_predicted_frequency_hz )
            .value( "doppler_troposphere_correction", tudat::input_output::TrackingDataType::doppler_troposphere_correction )
            .value( "scan_nr", tudat::input_output::TrackingDataType::scan_nr )
            .export_values( );

    py::enum_< tudat::input_output::TrackingTxtFileReadFilterType >(
            m, "TrackingTxtFileReadFilterType", R"doc(No documentation available.)doc" )
            .value( "no_tracking_txt_file_filter",
                    tudat::input_output::TrackingTxtFileReadFilterType::no_tracking_txt_file_filter,
                    R"doc(No documentation available.)doc" )
            .value( "ifms_tracking_txt_file_filter",
                    tudat::input_output::TrackingTxtFileReadFilterType::ifms_tracking_txt_file_filter,
                    R"doc(No documentation available.)doc" )
            .export_values( );

    py::class_< tio::solar_activity::SolarActivityData, std::shared_ptr< tio::solar_activity::SolarActivityData > >(
            m, "SolarActivityData", R"doc(No documentation available.)doc" )
            .def_readonly( "solar_radio_flux_107_observed", &tio::solar_activity::SolarActivityData::solarRadioFlux107Observed );

    // py::class_<std::map<
    //     double,
    //     std::shared_ptr<tio::solar_activity::SolarActivityData>>>(
    //     m, "SolarActivityDataMap");

    m.def( "read_solar_activity_data",
           &tio::solar_activity::readSolarActivityData,
           py::arg( "file_path" ),
           R"doc(
 Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.

 :param file_path: Path to the space weather data file.
 )doc" );

    py::class_< tio::solar_activity::SolarActivityContainer, std::shared_ptr< tio::solar_activity::SolarActivityContainer > >(
            m, "SolarActivityContainer" )

            .def( py::init< const std::map< double, std::shared_ptr< tio::solar_activity::SolarActivityData > >& >( ),
                  py::arg( "solar_activity_data_map" ) )

            .def( "get_solar_activity_data",
                  &tio::solar_activity::SolarActivityContainer::getSolarActivityData,
                  py::arg( "time" ),
                  R"doc(
        Returns the nearest SolarActivityData (in UTC Julian days) for the given time in seconds since J2000.
        )doc" )

            .def( "get_solar_activity_data_map",
                  &tio::solar_activity::SolarActivityContainer::getSolarActivityDataMap,
                  R"doc(Returns the full map of SolarActivityData.)doc" );

    py::enum_< tio::OdfDataType >( m, "OdfDataType", R"doc(Possible data types in orbit section of ODF file)doc" )
            .value( "narrowband_spacecraft_vlbi_doppler_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_doppler_mode )
            .value( "narrowband_spacecraft_vlbi_phase_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_phase_mode )
            .value( "narrowband_quasar_vlbi_doppler_mode", tio::OdfDataType::narrowband_quasar_vlbi_doppler_mode )
            .value( "narrowband_quasar_vlbi_phase_mode", tio::OdfDataType::narrowband_quasar_vlbi_phase_mode )
            .value( "wideband_spacecraft_vlbi", tio::OdfDataType::wideband_spacecraft_vlbi )
            .value( "wideband_quasar_vlbi", tio::OdfDataType::wideband_quasar_vlbi )
            .value( "one_way_doppler", tio::OdfDataType::one_way_doppler )
            .value( "two_way_doppler", tio::OdfDataType::two_way_doppler )
            .value( "three_way_doppler", tio::OdfDataType::three_way_doppler )
            .value( "one_way_total_count_phase", tio::OdfDataType::one_way_total_count_phase )
            .value( "two_way_total_count_phase", tio::OdfDataType::two_way_total_count_phase )
            .value( "three_way_total_count_phase", tio::OdfDataType::three_way_total_count_phase )
            .value( "pra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::pra_planetary_operational_discrete_spectrum_range )
            .value( "sra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::sra_planetary_operational_discrete_spectrum_range )
            .value( "re_range", tio::OdfDataType::re_range )
            .value( "azimuth_angle", tio::OdfDataType::azimuth_angle )
            .value( "elevation_angle", tio::OdfDataType::elevation_angle )
            .value( "hour_angle", tio::OdfDataType::hour_angle )
            .value( "declination_angle", tio::OdfDataType::declination_angle )
            .value( "x_angle_east", tio::OdfDataType::x_angle_east )
            .value( "y_angle_east", tio::OdfDataType::y_angle_east )
            .value( "x_angle_south", tio::OdfDataType::x_angle_south )
            .value( "y_angle_south", tio::OdfDataType::y_angle_south );

    py::class_< tio::OdfCommonDataBlock, std::shared_ptr< tio::OdfCommonDataBlock > >(
            m, "OdfCommonDataBlock", R"doc(Base class observable-independent ODF data containers

        The data section of an ODF is split into blocks (lines), each associated with
        an observation epoch. The first elements of each block (e.g. observation epoch,
        value of the observable) are common for all the observable types, but the
        values in the remaining columns will have different meanings for different
        types of observations. This class serves as interface to the observable-independent
        part of an ODF data block. The different classes inheriting from OdfDataSpecificBlock
        provide interfaces to the observable-specific part of the blocks.
        )doc" )
            .def_property_readonly( "observable_time", &tio::OdfCommonDataBlock::getObservableTime )
            .def_property_readonly( "observable_value", &tio::OdfCommonDataBlock::getObservableValue )
            .def_property_readonly( "receiving_station_downlink_delay", &tio::OdfCommonDataBlock::getReceivingStationDownlinkDelay )
            .def_readonly( "format_id", &tio::OdfCommonDataBlock::formatId_ )
            .def_readonly( "receiving_station_id", &tio::OdfCommonDataBlock::receivingStationId_ )
            .def_readonly( "transmitting_station_id", &tio::OdfCommonDataBlock::transmittingStationId_ )
            .def_readonly( "transmitting_station_network_id", &tio::OdfCommonDataBlock::transmittingStationNetworkId_ )
            .def_readonly( "data_type", &tio::OdfCommonDataBlock::dataType_ )
            .def_readonly( "downlink_band_id", &tio::OdfCommonDataBlock::downlinkBandId_ )
            .def_readonly( "uplink_band_id", &tio::OdfCommonDataBlock::uplinkBandId_ )
            .def_readonly( "reference_band_id", &tio::OdfCommonDataBlock::referenceBandId_ )
            .def_readonly( "is_invalid", &tio::OdfCommonDataBlock::validity_ )
            .def( "print_data_block",
                  &tio::OdfCommonDataBlock::printDataBlock,
                  py::arg( "output_file" ),
                  R"doc(Write the contents of the data block to a text file

                  The file is created if it does not exist, and it can have, for example, txt extension

                  :param output_file: Contents will be written to the file defined by this path
                  )doc" );

    py::class_< tio::OdfDataSpecificBlock, std::shared_ptr< tio::OdfDataSpecificBlock > >(
            m, "OdfDataSpecificBlock", R"doc(Base class observable-dependent ODF data containers

        The data section of an ODF is split into blocks (lines), each associated with
        an observation epoch. The first elements of each block (e.g. observation epoch,
        value of the observable) are common for all the observable types, but the
        values in the remaining columns will have different meanings for different
        types of observations. This base class serves as parent for interfaces to the
        observable-specific part of an ODF data block. The interface to the common
        part of the blocks is provided by the OdfCommonDataBlock class.
        )doc" );

    py::class_< tio::OdfDopplerDataBlock, std::shared_ptr< tio::OdfDopplerDataBlock >, tio::OdfDataSpecificBlock >(
            m, "OdfDopplerDataBlock", R"doc(Container for ODF Doppler-specific data)doc" )
            .def_property_readonly( "receiver_channel", &tio::OdfDopplerDataBlock::getReceiverChannel )
            .def_property_readonly( "spacecraft_id", &tio::OdfDopplerDataBlock::getSpacecraftId )
            .def_property_readonly( "receiver_exciter_flag", &tio::OdfDopplerDataBlock::getReceiverExciterFlag )
            .def_property_readonly( "reference_frequency", &tio::OdfDopplerDataBlock::getReferenceFrequency )
            .def_property_readonly( "compression_time", &tio::OdfDopplerDataBlock::getCompressionTime )
            .def_property_readonly( "transmitting_station_uplink_delay", &tio::OdfDopplerDataBlock::getTransmittingStationUplinkDelay );

    py::class_< tio::OdfDataBlock, std::shared_ptr< tio::OdfDataBlock > >(
            m, "OdfDataBlock", R"doc(Contents of a line of the data section of an ODF)doc" )
            .def_property_readonly( "observable_specific_data_block", &tio::OdfDataBlock::getObservableSpecificDataBlock )
            .def_property_readonly( "common_data_block", &tio::OdfDataBlock::getCommonDataBlock )
            .def( "print_data_block",
                  &tio::OdfDataBlock::printDataBlock,
                  py::arg( "output_file" ),
                  R"doc(Write the contents of the data block to a text file

                  The file is created if it does not exist, and it can have, for example, txt extension

                  :param output_file: Contents will be written to the file defined by this path
                  )doc" );

    py::class_< tio::OdfRampBlock, std::shared_ptr< tio::OdfRampBlock > >(
            m, "OdfRampBlock", R"doc(Contents of a line of the ramp section of an ODF)doc" )
            .def_property_readonly( "ramp_start_frequency", &tio::OdfRampBlock::getRampStartFrequency )
            .def_property_readonly( "ramp_rate", &tio::OdfRampBlock::getRampRate )
            .def_property_readonly( "ramp_start_epoch", &tio::OdfRampBlock::getRampStartTime )
            .def_property_readonly( "ramp_end_epoch", &tio::OdfRampBlock::getRampEndTime )
            .def_property_readonly( "transmitting_station_id", &tio::OdfRampBlock::getTransmittingStationId );

    py::class_< tio::OdfRawFileContents, std::shared_ptr< tio::OdfRawFileContents > >(
            m, "OdfRawFileContents", R"doc(No documentation available.)doc" )
            .def_property_readonly( "data_blocks", &tio::OdfRawFileContents::getDataBlocks )
            .def_property_readonly( "ramp_blocks", &tio::OdfRawFileContents::getRampBlocks )
            .def_property_readonly( "clock_offset_blocks", &tio::OdfRawFileContents::getClockOffsetBlocks )
            .def_readonly( "file_reference_date", &tio::OdfRawFileContents::fileReferenceDate_ )
            .def_readonly( "file_reference_time", &tio::OdfRawFileContents::fileReferenceTime_ )
            .def( "write_to_text_file",
                  &tio::OdfRawFileContents::writeOdfToTextFile,
                  py::arg( "output_file" ),
                  R"doc(No documentation available.)doc" );

    m.def( "read_odf_file", &tio::readOdfFile, py::arg( "file_name" ), R"doc(No documentation available.)doc" );

    m.def( "set_dsn_weather_data_in_ground_stations",
           py::overload_cast< tudat::simulation_setup::SystemOfBodies&,
                              const std::vector< std::string >&,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::map< int, std::vector< std::string > >&,
                              const std::string& >( &tio::setDsnWeatherDataInGroundStations ),
           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "interpolator_settings" ) = tudat::interpolators::linearInterpolation( ),
           py::arg( "ground_stations_per_complex" ) = tudat::simulation_setup::getDefaultDsnStationNamesPerComplex( ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           R"doc(No documentation available.)doc" );

    py::class_< tudat::input_output::TrackingTxtFileContents, std::shared_ptr< tudat::input_output::TrackingTxtFileContents > >(
            m, "TrackingTxtFileContents", R"doc(No documentation available.)doc" )
            .def( py::init< const std::string, const std::vector< std::string >, const char, const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t ",
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_val",
                  py::overload_cast< tio::TrackingDataType, double >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "value" ),
                  R"doc(No documentation available.)doc" )
            .def( "get_available_datatypes",
                  &tio::TrackingTxtFileContents::getAllAvailableDataTypes,
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_str",
                  py::overload_cast< tio::TrackingDataType, const std::string& >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "str_value" ),
                  R"doc(No documentation available.)doc" )
            .def_property_readonly(
                    "column_field_types", &tio::TrackingTxtFileContents::getRawColumnTypes, R"doc(No documentation available.)doc" )
            .def_property_readonly(
                    "double_datamap", &tio::TrackingTxtFileContents::getDoubleDataMap, R"doc(No documentation available.)doc" )
            .def_property_readonly( "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap, R"doc(No documentation available.)doc" )
            .def_property_readonly( "num_rows", &tio::TrackingTxtFileContents::getNumRows, R"doc(No documentation available.)doc" );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false,
           py::arg( "data_filter_method" ) = tio::no_tracking_txt_file_filter );

    m.def( "grail_antenna_file_reader", &tio::grailAntennaFileReader, py::arg( "file_name" ), R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_0_file_reader",
           &tio::grailMassLevel0FileReader,
           py::arg( "file_name" ),
           R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_1_file_reader",
           &tio::grailMassLevel1FileReader,
           py::arg( "file_name" ),
           py::arg( "data_level" ) = "1b",
           R"doc(No documentation available.)doc" );

    m.def( "read_ifms_file",
           &tio::readIfmsFile,
           py::arg( "file_name" ),
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(Load contents of IFMS file into object

           The keys of the dictionary represent the different columns of the IFMS file, and their values are lists with all the values in the associated column as strings.

           Two of the columns of an IFMS file contain, respectively, the Doppler averaged frequency and a tropospheric correction for the station. When the `apply_tropospheric_correction` option is set to true, the content of the first column is modified by subtracting the values in the second.

           :param file_name: String representing the path to the file to be loaded
           :param apply_tropospheric_correction: Whether to modify the averaged Doppler frequency as described above (Default: True)
           :param remove_invalid_lines: Boolean defining whether a line is skipped if the transmit frequency, observed frequency, or troposphere correction is undefined (Default: True)

           :return ifms_contents: Dictionary with contents of the IFMS file as lists of strings
           )doc" );
    m.def( "set_estrack_weather_data_in_ground_stations",
           py::overload_cast< tudat::simulation_setup::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::string,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::string& >( &tio::setEstrackWeatherDataInGroundStation ),

           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "ground_station_name" ),
           py::arg( "interpolator_settings" ) = tudat::interpolators::cubicSplineInterpolation( ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           R"doc(No documentation available.)doc" );

    // Temporary fix: Asigned default to variable to avoid compilation error
    const std::vector< std::string > fdet_cols = {
        "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz"
    };

    m.def( "read_fdets_file",
           &tio::readFdetsFile,
           py::arg( "file_name" ),
           py::arg( "column_types" ) = fdet_cols,
           R"doc(Load contents of Fdets file into object

           This function loads the contents of an Fdets file into a dictionary with keys defined by the strings listed in `column_types`. If a list of columns is not specified, the following structure is assumed:

           - UTC datetime string
           - Signal to noise ratio [-]
           - Normalized spectral max (Not sure what this is)
           - Doppler measured frequency [Hz]
           - Doppler noise [Hz]

           :param file_name: String representing the path to the file to be loaded
           :param column_types: Identifiers of the columns present in the Fdets file
           :return fdets_contents: Dictionary with contents of the Fdets file as lists of strings
           )doc" );
};

}  // namespace data
}  // namespace tudatpy
