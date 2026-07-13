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
#include "expose_coma.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/comaModelInputOutput.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace data_access
{
namespace environment
{
namespace coma
{

void expose_coma( py::module& m )
{
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
           :func:`~tudatpy.data.environment.coma.coma_model_file_processor`
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
   processor = data.environment.coma.coma_model_file_processor(file_paths)

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
           :func:`~tudatpy.data.environment.coma.coma_model_file_processor`
           and its ``create_coma_stokes_dataset()`` method, or load from pre-existing CSV files.


 Examples
 --------
 Create a Stokes dataset by transforming polynomial coefficients:

 .. code-block:: python

   # Create file processor from polynomial files
   processor = data.environment.coma.coma_model_file_processor(
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
   processor = data.environment.coma.coma_model_file_processor(
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
            .def( "is_poly",
                  &tss::ComaWindDatasetCollection::isPoly,
                  R"doc(

 Check if the collection contains polynomial coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.data.environment.coma.ComaPolyDataset` objects,
     False if it contains :class:`~tudatpy.data.environment.coma.ComaStokesDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files
   processor = data.environment.coma.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_poly_coefficient_datasets()

   # Check dataset type
   if wind_datasets.is_poly():
       print("Collection contains polynomial datasets")


      )doc" )
            .def( "is_stokes",
                  &tss::ComaWindDatasetCollection::isStokes,
                  R"doc(

 Check if the collection contains Stokes coefficient datasets.

 Returns
 -------
 bool
     True if the collection contains :class:`~tudatpy.data.environment.coma.ComaStokesDataset` objects,
     False if it contains :class:`~tudatpy.data.environment.coma.ComaPolyDataset` objects.


 Examples
 --------
 .. code-block:: python

   # Create wind dataset collection from polynomial files and convert to Stokes
   processor = data.environment.coma.coma_wind_file_processor(
       x_file_paths, y_file_paths, z_file_paths)
   wind_datasets = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Check dataset type
   if wind_datasets.is_stokes():
       print("Collection contains Stokes datasets")


      )doc" );

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
           :func:`~tudatpy.data.environment.coma.coma_model_file_processor`.


 Examples
 --------
 Polynomial coefficient workflow:

 .. code-block:: python

   # Create processor from polynomial coefficient files
   poly_files = ["h2o_epoch1.txt", "h2o_epoch2.txt"]
   processor = data.environment.coma.coma_model_file_processor(poly_files)

   # Create polynomial dataset directly
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Or transform to Stokes coefficients
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

 Stokes coefficient workflow:

 .. code-block:: python

   # Create processor from existing Stokes CSV files
   processor = data.environment.coma.coma_model_file_processor(
       input_dir="coma_stokes_data",
       prefix="stokes")

   # Load Stokes dataset
   stokes_dataset = processor.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc" )
            .def( "create_poly_coefficient_dataset",
                  &tss::ComaModelFileProcessor::createPolyCoefDataset,
                  R"doc(

 Create polynomial coefficient dataset from the loaded files.

 Reads and processes polynomial coefficient files to create a :class:`~tudatpy.data.environment.coma.ComaPolyDataset`.
 This method is only available when the processor was constructed from polynomial coefficient files
 using the file path variant of :func:`~tudatpy.data.environment.coma.coma_model_file_processor`.

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
   processor = data.environment.coma.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Create polynomial dataset
   poly_dataset = processor.create_poly_coefficient_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model_from_poly_data(
       poly_data=poly_dataset,
       molecular_weight=0.018015)


      )doc" )
            .def( "create_coma_stokes_dataset",
                  py::overload_cast<>( &tss::ComaModelFileProcessor::createSHDataset, py::const_ ),
                  R"doc(

 Create Stokes coefficient dataset from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.data.environment.coma.coma_model_file_processor`.
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
   processor = data.environment.coma.coma_model_file_processor(
       input_dir="stokes_data",
       prefix="stokes")

   # Load Stokes dataset (no parameters needed)
   stokes_dataset = processor.create_coma_stokes_dataset()

   # Use in coma atmosphere model
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)

      )doc" )
            .def( "create_coma_stokes_dataset",
                  py::overload_cast< const std::vector< double >&,
                                     const std::vector< double >&,
                                     const int,
                                     const int,
                                     const bool,
                                     const bool >( &tss::ComaModelFileProcessor::createSHDataset, py::const_ ),
                  py::arg( "radii_m" ),
                  py::arg( "sol_longitudes_deg" ),
                  py::arg( "requested_max_degree" ) = -1,
                  py::arg( "requested_max_order" ) = -1,
                  py::arg( "compute_reduced_coeffs" ) = true,
                  py::arg( "is_log2" ) = true,
                  R"doc(

 Create Stokes coefficient dataset by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.data.environment.coma.coma_model_file_processor`. It transforms
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
   processor = data.environment.coma.coma_model_file_processor(
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

      )doc" )
            .def( "create_coma_spherical_harmonic_coefficient_files",
                  &tss::ComaModelFileProcessor::createSHFiles,
                  py::arg( "output_dir" ),
                  py::arg( "radii_m" ),
                  py::arg( "sol_longitudes_deg" ),
                  py::arg( "requested_max_degree" ) = -1,
                  py::arg( "requested_max_order" ) = -1,
                  py::arg( "compute_reduced_coeffs" ) = true,
                  py::arg( "is_log2" ) = true,
                  R"doc(

 Create and save Stokes coefficient CSV files from polynomial coefficients.

 Transforms polynomial coefficients to Stokes coefficients and saves them as CSV files in the specified
 output directory. This is useful for pre-computing Stokes coefficients to avoid repeated transformations
 during multiple simulation runs. The generated CSV files can later be loaded using a processor created
 with the directory-based variant of :func:`~tudatpy.data.environment.coma.coma_model_file_processor`.

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
   processor = data.environment.coma.coma_model_file_processor(
       ["h2o_poly_epoch1.txt", "h2o_poly_epoch2.txt"])

   # Pre-compute and save Stokes coefficients to CSV files
   processor.create_coma_spherical_harmonic_coefficient_files(
       output_dir="stokes_output",
       radii_m=[1000.0, 2000.0, 5000.0, 10000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0],
       requested_max_degree=10,
       requested_max_order=10)

   # Later, load from the saved CSV files
   processor_from_csv = data.environment.coma.coma_model_file_processor(
       input_dir="stokes_output",
       prefix="stokes")
   stokes_dataset = processor_from_csv.create_coma_stokes_dataset(
       radii_m=[],  # Ignored when loading from files
       sol_longitudes_deg=[])


      )doc" );

    m.def( "coma_model_file_processor",
           py::overload_cast< const std::vector< std::string >& >( &tss::comaModelFileProcessorFromPolyFiles ),
           py::arg( "file_paths" ),
           R"doc(

 Function for creating coma model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.data.environment.coma.ComaModelFileProcessor` that loads
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
   processor = data.environment.coma.coma_model_file_processor(poly_files)

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
   processor = data.environment.coma.coma_model_file_processor(poly_files)

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


      )doc" );

    m.def( "coma_model_file_processor",
           py::overload_cast< const std::string&, const std::string& >( &tss::comaModelFileProcessorFromSHFiles ),
           py::arg( "input_dir" ),
           py::arg( "prefix" ) = "stokes",
           R"doc(

 Function for creating coma model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.data.environment.coma.ComaModelFileProcessor` that loads
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
   processor = data.environment.coma.coma_model_file_processor(
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
   poly_processor = data.environment.coma.coma_model_file_processor(
       ["h2o_poly.txt"])
   poly_processor.create_coma_spherical_harmonic_coefficient_files(
       output_dir="stokes_saved",
       radii_m=[1000.0, 2000.0, 5000.0],
       sol_longitudes_deg=[0.0, 90.0, 180.0, 270.0])

   # Step 2: Later, load from saved Stokes files
   stokes_processor = data.environment.coma.coma_model_file_processor(
       input_dir="stokes_saved",
       prefix="stokes")
   stokes_dataset = stokes_processor.create_coma_stokes_dataset(
       radii_m=[],
       sol_longitudes_deg=[])

   # Step 3: Use in simulation
   coma_settings = environment_setup.atmosphere.coma_model_from_stokes_data(
       stokes_data=stokes_dataset,
       molecular_weight=0.018015)


      )doc" );

    // === Coma wind processing: file processor ===
    py::class_< tss::ComaWindModelFileProcessor, std::shared_ptr< tss::ComaWindModelFileProcessor > >( m,
                                                                                                       "ComaWindModelFileProcessor",
                                                                                                       R"doc(
        Processor for creating wind model datasets from three component file sources.

        This class manages the creation of ComaWindDatasetCollection from three sets of files
        (one for each spatial component: x, y, z). It provides a simplified interface for
        wind model setup by handling all three components together.
        )doc" )
            .def( "create_poly_coefficient_datasets",
                  &tss::ComaWindModelFileProcessor::createPolyCoefDatasets,
                  R"doc(

 Create polynomial coefficient dataset collection for all three wind components.

 Reads and processes polynomial coefficient files for x, y, and z wind velocity components to create
 a :class:`~tudatpy.data.environment.coma.ComaWindDatasetCollection`. This method is
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
       x_files, y_files, z_files)

   # Create polynomial dataset collection
   poly_datasets = wind_processor.create_poly_coefficient_datasets()


      )doc" )
            .def( "create_coma_stokes_dataset",
                  py::overload_cast<>( &tss::ComaWindModelFileProcessor::createSHDatasets, py::const_ ),
                  R"doc(

 Create Stokes coefficient dataset collection from preloaded CSV files (parameterless version).

 This method is only available when the processor was constructed from Stokes coefficient CSV files
 using :func:`~tudatpy.data.environment.coma.coma_wind_file_processor`.
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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

      )doc" )
            .def( "create_coma_stokes_dataset",
                  py::overload_cast< const std::vector< double >&, const std::vector< double >&, const int, const int >(
                          &tss::ComaWindModelFileProcessor::createSHDatasets, py::const_ ),
                  py::arg( "radii_m" ),
                  py::arg( "sol_longitudes_deg" ),
                  py::arg( "requested_max_degree" ) = -1,
                  py::arg( "requested_max_order" ) = -1,
                  R"doc(

 Create Stokes coefficient dataset collection by transforming polynomial coefficients (parameterized version).

 This method is only available when the processor was constructed from polynomial coefficient files
 using :func:`~tudatpy.data.environment.coma.coma_wind_file_processor`. It transforms
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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

      )doc" )
            .def( "create_coma_spherical_harmonic_coefficient_files",
                  &tss::ComaWindModelFileProcessor::createSHFiles,
                  py::arg( "x_output_dir" ),
                  py::arg( "y_output_dir" ),
                  py::arg( "z_output_dir" ),
                  py::arg( "radii_m" ),
                  py::arg( "sol_longitudes_deg" ),
                  py::arg( "requested_max_degree" ) = -1,
                  py::arg( "requested_max_order" ) = -1,
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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
   wind_processor_from_csv = data.environment.coma.coma_wind_file_processor(
       x_input_dir="stokes_wind/x_component",
       y_input_dir="stokes_wind/y_component",
       z_input_dir="stokes_wind/z_component",
       prefix="stokes")


      )doc" )
            .def( "is_poly_type",
                  &tss::ComaWindModelFileProcessor::isPolyType,
                  R"doc(Check if this processor works with polynomial coefficient files.)doc" )
            .def( "is_stokes_type",
                  &tss::ComaWindModelFileProcessor::isStokesType,
                  R"doc(Check if this processor works with Stokes coefficient files.)doc" );

    m.def( "coma_wind_file_processor",
           py::overload_cast< const std::vector< std::string >&, const std::vector< std::string >&, const std::vector< std::string >& >(
                   &tss::comaWindModelFileProcessorFromPolyFiles ),
           py::arg( "x_file_paths" ),
           py::arg( "y_file_paths" ),
           py::arg( "z_file_paths" ),
           R"doc(

 Function for creating coma wind model file processor from polynomial coefficient files.

 Creates a :class:`~tudatpy.data.environment.coma.ComaWindModelFileProcessor` that loads
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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


      )doc" );

    m.def( "coma_wind_file_processor",
           py::overload_cast< const std::string&, const std::string&, const std::string&, const std::string& >(
                   &tss::comaWindModelFileProcessorFromSHFiles ),
           py::arg( "x_input_dir" ),
           py::arg( "y_input_dir" ),
           py::arg( "z_input_dir" ),
           py::arg( "prefix" ) = "stokes",
           R"doc(

 Function for creating coma wind model file processor from Stokes coefficient CSV files.

 Creates a :class:`~tudatpy.data.environment.coma.ComaWindModelFileProcessor` that loads
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
   wind_processor = data.environment.coma.coma_wind_file_processor(
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
   poly_processor = data.environment.coma.coma_wind_file_processor(
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
   stokes_processor = data.environment.coma.coma_wind_file_processor(
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


      )doc" );

}  // expose_coma

}  // namespace coma
}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy
