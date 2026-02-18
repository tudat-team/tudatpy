/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_environment.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tudat/astro/aerodynamics.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/astro/gravitation.h>
#include <tudat/basics/deprecationWarnings.h>

#include "scalarTypes.h"

#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace tr = tudat::reference_frames;
namespace te = tudat::ephemerides;
namespace teo = tudat::earth_orientation;
namespace tgs = tudat::ground_stations;
namespace tr = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace trf = tudat::reference_frames;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;
namespace tsm = tudat::system_models;
namespace tom = tudat::observation_models;
namespace tem = tudat::electromagnetism;

namespace tudat
{

namespace aerodynamics
{

double getTotalSurfaceArea( const std::shared_ptr< HypersonicLocalInclinationAnalysis > coefficientGenerator )
{
    double totalSurfaceArea = 0.0;
    for( int i = 0; i < coefficientGenerator->getNumberOfVehicleParts( ); i++ )
    {
        totalSurfaceArea += std::fabs( coefficientGenerator->getVehiclePart( i )->getTotalArea( ) );
    }
    return totalSurfaceArea;
}

//! Function that saves the vehicle mesh data used for a
//! HypersonicLocalInclinationAnalysis to a file
std::pair< std::vector< Eigen::Vector3d >, std::vector< Eigen::Vector3d > > getVehicleMesh(
        const std::shared_ptr< HypersonicLocalInclinationAnalysis > localInclinationAnalysis )
{
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshPoints = localInclinationAnalysis->getMeshPoints( );
    std::vector< boost::multi_array< Eigen::Vector3d, 2 > > meshSurfaceNormals = localInclinationAnalysis->getPanelSurfaceNormals( );

    //    boost::array< int, 3 > independentVariables;
    //    independentVariables[ 0 ] = 0;
    //    independentVariables[ 1 ] = 6;
    //    independentVariables[ 2 ] = 0;

    //    std::vector< std::vector< std::vector< double > > >
    //    pressureCoefficients =
    //            localInclinationAnalysis->getPressureCoefficientList(
    //            independentVariables );

    std::vector< Eigen::Vector3d > meshPointsList;
    std::vector< Eigen::Vector3d > surfaceNormalsList;
    //    std::map< int, Eigen::Vector1d > pressureCoefficientsList;

    // int counter = 0;
    for( unsigned int i = 0; i < meshPoints.size( ); i++ )
    {
        for( unsigned int j = 0; j < meshPoints.at( i ).shape( )[ 0 ] - 1; j++ )
        {
            for( unsigned int k = 0; k < meshPoints.at( i ).shape( )[ 1 ] - 1; k++ )
            {
                meshPointsList.push_back( meshPoints[ i ][ j ][ k ] );
                surfaceNormalsList.push_back( meshSurfaceNormals[ i ][ j ][ k ] );
                //                pressureCoefficientsList[ counter ] =
                //                ( Eigen::Vector1d( ) <<
                //                pressureCoefficients[ i ][ j ][ k ]
                //                ).finished( );
                // counter++;
            }
        }
    }

    return std::make_pair( meshPointsList, surfaceNormalsList );
}

}  // namespace aerodynamics

}  // namespace tudat

namespace tudatpy
{
namespace dynamics
{
namespace environment
{

void expose_environment( py::module& m )
{
    /*!
     **************   EPHEMERIDES  ******************
     */

    py::class_< te::Ephemeris, std::shared_ptr< te::Ephemeris > >( m, "Ephemeris", R"doc(

         Object that computes the state of a body as a function of time


         Object (typically stored inside a :class:`~Body` object) that computes the state of a body as a function of time,
         both outside of a propagation, and during a propagation if the given body's translational state is not propagated.
         Note that this object computes the state w.r.t. its own origin (defined by ``frame_origin``), which need not be the same as the global frame origin
         of the environment.





      )doc" )
            .def( "cartesian_state",
                  &te::Ephemeris::getCartesianState,
                  py::arg( "current_time" ),
                  R"doc(

         This function returns the Cartesian state (position and velocity) at the given time, w.r.t. the ``frame_origin``.


         Parameters
         ----------
         current_time : float
             Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

         Returns
         -------
         numpy.ndarray
             Requested Cartesian state





     )doc" )
            .def( "cartesian_position",
                  &te::Ephemeris::getCartesianPosition,
                  py::arg( "current_time" ),
                  R"doc(

         As ``cartesian_state``, but only the three position components


         Parameters
         ----------
         current_time : float
             Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

         Returns
         -------
         numpy.ndarray
             Requested Cartesian position





     )doc" )
            .def( "cartesian_velocity",
                  &te::Ephemeris::getCartesianVelocity,
                  py::arg( "current_time" ),
                  R"doc(

         As ``cartesian_state``, but only the three velocity components


         Parameters
         ----------
         current_time : float
             Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

         Returns
         -------
         numpy.ndarray
             Requested Cartesian velocity





     )doc" )
            .def_property_readonly( "frame_origin",
                                    &te::Ephemeris::getReferenceFrameOrigin,
                                    R"doc(

         **read-only**

         Name of the reference body/point w.r.t. which this object provides its states


         :type: str
      )doc" )
            .def_property_readonly( "frame_orientation",
                                    &te::Ephemeris::getReferenceFrameOrientation,
                                    R"doc(

         **read-only**

         Name of the frame orientation w.r.t which this object provides its states



         :type: str
      )doc" );

    py::class_< te::ConstantEphemeris, std::shared_ptr< te::ConstantEphemeris >, te::Ephemeris >( m,
                                                                                                  "ConstantEphemeris",
                                                                                                  R"doc(No documentation found.)doc" )
            .def( py::init< const std::function< Eigen::Vector6d( ) >,  //<pybind11/functional.h>,<pybind11/eigen.h>
                            const std::string&,
                            const std::string& >( ),
                  py::arg( "constant_state_function" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" )
            .def( py::init< const Eigen::Vector6d,  //<pybind11/eigen.h>
                            const std::string&,
                            const std::string& >( ),
                  py::arg( "constant_state" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" )
            .def( "update_constant_state",
                  &te::ConstantEphemeris::updateConstantState,
                  py::arg( "new_state" ),
                  R"doc(No documentation found.)doc" );

    py::class_< te::KeplerEphemeris, std::shared_ptr< te::KeplerEphemeris >, te::Ephemeris >( m, "KeplerEphemeris" );

    py::class_< te::MultiArcEphemeris, std::shared_ptr< te::MultiArcEphemeris >, te::Ephemeris >( m, "MultiArcEphemeris" )
            .def( py::init< const std::map< double, std::shared_ptr< te::Ephemeris > >&, const std::string&, const std::string& >( ),
                  py::arg( "single_arc_ephemerides" ),
                  py::arg( "reference_frame_origin" ) = "SSB",
                  py::arg( "reference_frame_orientation" ) = "ECLIPJ2000" );

    py::class_< te::TabulatedCartesianEphemeris< double, double >,
                std::shared_ptr< te::TabulatedCartesianEphemeris< double, double > >,
                te::Ephemeris >( m, "TabulatedEphemeris" )
            .def_property( "interpolator",
                           &te::TabulatedCartesianEphemeris< double, double >::getDynamicVectorSizeInterpolator,
                           py::overload_cast< const std::shared_ptr< ti::OneDimensionalInterpolator< double, Eigen::VectorXd > > >(
                                   &te::TabulatedCartesianEphemeris< double, double >::resetInterpolator ) );

    py::class_< te::Tle, std::shared_ptr< te::Tle > >( m, "Tle" )
            .def( py::init<  // ctor 1
                          const std::string& >( ),
                  py::arg( "lines" ) )
            .def( py::init<  // ctor 2
                          const std::string&,
                          const std::string& >( ),
                  py::arg( "line_1" ),
                  py::arg( "line_2" ) )
            .def( "get_epoch", &te::Tle::getEpoch )
            .def( "get_b_star", &te::Tle::getBStar )
            .def( "get_epoch", &te::Tle::getEpoch )
            .def( "get_inclination", &te::Tle::getInclination )
            .def( "get_right_ascension", &te::Tle::getRightAscension )
            .def( "get_eccentricity", &te::Tle::getEccentricity )
            .def( "get_arg_of_perigee", &te::Tle::getArgOfPerigee )
            .def( "get_mean_anomaly", &te::Tle::getMeanAnomaly )
            .def( "get_mean_motion", &te::Tle::getMeanMotion );

    py::class_< te::TleEphemeris, std::shared_ptr< te::TleEphemeris >, te::Ephemeris >( m, "TleEphemeris" )
            .def( py::init< const std::string&, const std::string&, const std::shared_ptr< te::Tle >, const bool >( ),
                  py::arg( "frame_origin" ) = "Earth",
                  py::arg( "frame_orientation" ) = "J2000",
                  py::arg( "tle" ) = nullptr,
                  py::arg( "use_sdp" ) = false );

    /*!
     **************   END EPHEMERIDES  ******************
     */

    py::class_< tudat::environment::IonosphereModel, std::shared_ptr< tudat::environment::IonosphereModel > >( m,
                                                                                                               "IonosphereModel",
                                                                                                               R"doc(
Base class for ionospheric models.

Provides the vertical total electron content (VTEC) in TECU (1 TECU = 1e16 e-/m²) based on geodetic position (latitude, longitude) and time.

This is the base class from which models like TabulatedIonosphereModel or GlobalIonosphereModelVtecCalculator retrieve electron content data. The model is typically stored
inside a `Body` instance and used in observation corrections or environmental queries.
)doc" );

    py::class_< ta::AtmosphereModel, std::shared_ptr< ta::AtmosphereModel > >( m,
                                                                               "AtmosphereModel",
                                                                               R"doc(

         Object that provides the atmospheric properties of the body.

         Object that provides the atmospheric properties of the body, as a function of altitude, latitude, longitude and time. Depending on the implementation, the
         this dependence may be limited to altitude-only (e.g. standard atmosphere models). This object can be accessed directly by the user to compute atmospheric properties
         outside the loop of the propagation by calling one of its member functions. During the propagation, each body undergoing aerodynamic forces has a :class:`AtmosphericFlightConditions`
         object associated with it (accessed from a boyd through :attr:`~Body.flight_condition`) that links the atmosphere model to the aerodynamic model.

     )doc" )
            .def( "get_density",
                  &ta::AtmosphereModel::getDensity,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(

         Function to compute the atmospheric freestream density at a given location.

         Parameters
         ----------
         altitude : float
             Local altitude above the body surface at which the property is to be computed
         latitude : float
             Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         longitude : float
             Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the property is to be computed.

         Returns
         -------
         float
             Freestream density at the given time and location


     )doc" )
            .def( "get_pressure",
                  &ta::AtmosphereModel::getPressure,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(

         Function to compute the atmospheric freestream static pressure at a given location.

         Parameters
         ----------
         altitude : float
             Local altitude above the body surface at which the property is to be computed
         latitude : float
             Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         longitude : float
             Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the property is to be computed.

         Returns
         -------
         float
             Freestream static pressure at the given time and location


     )doc" )
            .def( "get_temperature",
                  &ta::AtmosphereModel::getTemperature,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(

         Function to compute the atmospheric freestream temperature at a given location.

         Parameters
         ----------
         altitude : float
             Local altitude above the body surface at which the property is to be computed
         latitude : float
             Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         longitude : float
             Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the property is to be computed.

         Returns
         -------
         float
             Freestream temperature at the given time and location


     )doc" )
            .def( "get_speed_of_sound",
                  &ta::AtmosphereModel::getSpeedOfSound,
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(

         Function to compute the atmospheric freestream speed of sound at a given location.

         Parameters
         ----------
         altitude : float
             Local altitude above the body surface at which the property is to be computed
         latitude : float
             Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         longitude : float
             Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the property is to be computed.

         Returns
         -------
         float
             Freestream speed of sound at the given time and location


     )doc" )
            .def( "get_number_density",
                  &ta::AtmosphereModel::getNumberDensity,
                  py::arg( "species" ),
                  py::arg( "altitude" ),
                  py::arg( "longitude" ),
                  py::arg( "latitude" ),
                  py::arg( "time" ),
                  R"doc(

         Function to compute the atmospheric freestream number density of a given specie at a given location.

         Parameters
         ----------
         species : AtmosphericCompositionSpecies
             Atmospheric species for which the number density is to be computed
         altitude : float
             Local altitude above the body surface at which the property is to be computed
         latitude : float
             Geographic latitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         longitude : float
             Geographic longitude (in the body-fixed frame of the body with the atmosphere) at which the property is to be computed
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the property is to be computed.

         Returns
         -------
         float
             Freestream number density of the requested specie at the given time and location


     )doc" );

    py::class_< ta::AerodynamicCoefficientInterface, std::shared_ptr< ta::AerodynamicCoefficientInterface > >(
            m,
            "AerodynamicCoefficientInterface",
            R"doc(

         Base class for computing the current aerodynamic coefficients of the body


         Base class for computing the current aerodynamic coefficients of the body. The implementation of the computation
         depends on the choice of aerodynamic coefficient model (see :ref:`aerodynamic_coefficients` for available options).
         During the propagation, this object is automatically updated to the current state by the :class:`~AtmosphericFlightConditions` object.
         The user may override the current aerodynamic coefficients when using, for instance, a custom aerodynamic guidance model
         (see `here <https://docs.tudat.space/en/latest/_src_getting_started/_src_examples/notebooks/propagation/reentry_trajectory.html>`_ for an example).
         using the member functions of this class.





      )doc" )
            .def_property_readonly( "reference_area",
                                    &ta::AerodynamicCoefficientInterface::getReferenceArea,
                                    R"doc(

         **read-only**

         The aerodynamic reference area :math:`A` of the coefficients


         :type: float
      )doc" )
            .def_property_readonly( "current_force_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                                    R"doc(

         **read-only**

         The current aerodynamic force coefficients, in the frame defined by the :attr:`~force_coefficient_frame` attribute,
         as computed by the last call to the :meth:`~update_coefficients` function.


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "current_moment_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                                    R"doc(

         **read-only**

         The current aerodynamic moment coefficients, in the frame defined by the :attr:`~moment_coefficient_frame` attribute,
         as computed by the last call to the :meth:`~update_coefficients` function.


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "current_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentAerodynamicCoefficients,
                                    R"doc(

         **read-only**

         Concatenation of :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients`


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "force_coefficient_frame",
                                    &ta::AerodynamicCoefficientInterface::getForceCoefficientsFrame,
                                    R"doc(

         **read-only**

         Reference frame in which the  :attr:`~current_force_coefficients` are defined


         :type: AerodynamicCoefficientFrames
      )doc" )
            .def_property_readonly( "moment_coefficient_frame",
                                    &ta::AerodynamicCoefficientInterface::getMomentCoefficientsFrame,
                                    R"doc(

         **read-only**

         Reference frame in which the  :attr:`~current_moment_coefficients` are defined


         :type: AerodynamicCoefficientFrames
      )doc" )
            .def_property_readonly( "independent_variable_names",
                                    &ta::AerodynamicCoefficientInterface::getIndependentVariableNames,
                                    R"doc(

         **read-only**

         List of independent variables from which the aerodynamic coefficients are computed (e.g. required input to :meth:`~update_coefficients` function).


         :type: list[AerodynamicCoefficientsIndependentVariables]
      )doc" )
            .def_property_readonly( "current_control_surface_free_force_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeForceCoefficients,
                                    R"doc(

         **read-only**

         Same as :attr:`current_force_coefficients`, but without contribution (if any) from control surfaces


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "current_control_surface_free_moment_coefficients",
                                    &ta::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeMomentCoefficients,
                                    R"doc(

         **read-only**

         Same as :attr:`current_moment_coefficients`, but without contribution (if any) from control surfaces


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "control_surface_independent_variable_names",
                                    &ta::AerodynamicCoefficientInterface::getControlSurfaceIndependentVariables,
                                    R"doc(

         **read-only**

         List of independent variables from which the aerodynamic coefficients of each control surface are computed, with dictionary key being the control surface name (e.g. required input to :meth:`~update_full_coefficients` function).


         :type: dict[str,list[AerodynamicCoefficientsIndependentVariables]]
      )doc" )
            .def( "current_control_surface_force_coefficient_increment",
                  &ta::AerodynamicCoefficientInterface::getCurrentForceCoefficientIncrement,
                  py::arg( "control_surface_name" ),
                  R"doc(

         Function to get the contribution from a single control surface to the aerodynamic force coefficient, as compute by last call to :meth:`~update_full_coefficients`



         Parameters
         ----------
         control_surface_name : str
             The name of the control surface for which the contribution is to be retrieved

         Returns
         -------
         numpy.ndarray
             Contribution from the requested control surface to the aerodynamic force coefficient





     )doc" )
            .def( "current_control_surface_moment_coefficient_increment",
                  &ta::AerodynamicCoefficientInterface::getCurrentMomentCoefficientIncrement,
                  py::arg( "control_surface_name" ),
                  R"doc(

         Function to get the contribution from a single control surface to the aerodynamic moment coefficients, as compute by last call to :meth:`~update_full_coefficients`



         Parameters
         ----------
         control_surface_name : str
             The name of the control surface for which the contribution is to be retrieved

         Returns
         -------
         numpy.ndarray
             Contribution from the requested control surface to the aerodynamic moment coefficients





     )doc" )
            .def( "set_control_surface_increments",
                  &ta::AerodynamicCoefficientInterface::setControlSurfaceIncrements,
                  py::arg( "control_surface_list" ),
                  R"doc(No documentation found.)doc" )
            .def( "update_coefficients",
                  &ta::AerodynamicCoefficientInterface::updateCurrentCoefficients,
                  py::arg( "independent_variables" ),
                  py::arg( "time" ),
                  R"doc(

         Function to update the aerodynamic coefficients of the body only


         Function to update the aerodynamic coefficients of the body only (without the control surface contribution),
         based on the current state. This function may be called by the user, but will set *only* the
         :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients` (while leaving the
         :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients` unchanged)


         Parameters
         ----------
         independent_variables : list[float]
             List of inputs from which the aerodynamic coefficients are to be computed, with each entry corresponding to the
             value of the physical variable defined by the :attr:`independent_variable_names` attribute.

         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB)

         Returns
         -------
         numpy.ndarray
             Contribution from the requested control surface to the aerodynamic moment coefficients





     )doc" )
            .def( "update_full_coefficients",
                  &ta::AerodynamicCoefficientInterface::updateFullCurrentCoefficients,
                  py::arg( "independent_variables" ),
                  py::arg( "control_surface_independent_variables" ),
                  py::arg( "time" ),
                  py::arg( "check_force_contribution" ) = true,
                  R"doc(

         Function to update the aerodynamic coefficients, from both the body and its control surfaces


         Function to update the aerodynamic coefficients of both the body and its control surfaces,
         based on the current state. This function will call the :meth:`~update_coefficients` function to update the body coefficients.
         This function may be called by the user, and will set the following attributes:
         :attr:`~current_force_coefficients`, :attr:`~current_moment_coefficients` ,
         :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients`.
         In addition, it will modify the coefficients returned by the :meth:`~current_control_surface_force_coefficient_increment` and
         :meth:`~current_control_surface_moment_coefficient_increment` functions


         Parameters
         ----------
         independent_variables : list[float]
             List of inputs from which the aerodynamic coefficients of the body are to be computed, with each entry corresponding to the
             value of the physical variable defined by the :attr:`independent_variable_names` attribute.

         control_surface_independent_variables : dict[str,list[float]]
             List of inputs from which the control surface aerodynamic coefficients are to be computed (with dictionary key the control surface name),
             with each entry corresponding to the
             value of the physical variable defined by the :attr:`control_surface_independent_variable_names` attribute.

         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB)

         check_force_contribution : bool, default = True
             Boolean that determines if the force contribution to the aerodynamic moments should be added. Note that this input is
             only used if the :attr:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings.add_force_contribution_to_moments` attribute is set to True.





     )doc" );

    py::class_< ta::AerodynamicCoefficientGenerator< 3, 6 >,
                std::shared_ptr< ta::AerodynamicCoefficientGenerator< 3, 6 > >,
                ta::AerodynamicCoefficientInterface >( m, "AerodynamicCoefficientGenerator36", "<no_doc, only_dec>" );

    py::class_< ta::HypersonicLocalInclinationAnalysis,
                std::shared_ptr< ta::HypersonicLocalInclinationAnalysis >,
                ta::AerodynamicCoefficientGenerator< 3, 6 > >( m, "HypersonicLocalInclinationAnalysis" )
            .def( py::init< const std::vector< std::vector< double > >&,
                            const std::shared_ptr< tudat::SurfaceGeometry >,
                            const std::vector< int >&,
                            const std::vector< int >&,
                            const std::vector< bool >&,
                            const std::vector< std::vector< int > >&,
                            const double,
                            const double,
                            const Eigen::Vector3d&,
                            const bool >( ),
                  py::arg( "independent_variable_points" ),
                  py::arg( "body_shape" ),
                  py::arg( "number_of_lines" ),
                  py::arg( "number_of_points" ),
                  py::arg( "invert_orders" ),
                  py::arg( "selected_methods" ),
                  py::arg( "reference_area" ),
                  py::arg( "reference_length" ),
                  py::arg( "moment_reference_point" ),
                  py::arg( "save_pressure_coefficients" ) = false,
                  R"doc(

         Class constructor, taking the shape of the vehicle, and various analysis options as input.


         Parameters
         ----------
         independent_variable_points : list[list[float]]
             List containing three lists, with each sublist containing the data points of each of the
             independent variables for the coefficient generation. The physical meaning of each of the
             three independent variables is: 0 = mach number, 1 = angle of attack, 2 = angle of sideslip.
             Each of the subvectors must be sorted in ascending order.

         body_shape : SurfaceGeometry
             Class that defines the shape of the vehicle as a continuous surface. The local inclination analysis
             discretizes the surface of the vehicle into quadrilateral panels, defined by the other inputs to
             this constructor. In case the :class:`tudat.geometry.SurfaceGeometry` object is made up of multiple
             sub-shapes, different settings may be used for each

         number_of_lines : List[ float ]
             Number of discretization points in the first independent surface variable of each of the subparts of body_shape.
             The size of this list should match the number of parts of which the body_shape is composed. The first independent
             variable of a subpart typically runs along the longitudinal vehicle direction

         number_of_points : List[ float ]
             Number of discretization points in the second independent surface variable of each of the subparts of body_shape.
             The size of this list should match the number of parts of which the body_shape is composed. The first independent
             variable of a subpart typically runs along the lateral vehicle direction

         invert_orders : List[ bool ]
             Booleans to denote whether the surface normals of the panels of each discretized body_shape subpart are to be inverted
             (i.e. inward-facing->outward facing or vice versa). The size of this list should match the number of parts of which the body_shape is composed.

         selected_methods : List[ List[ int ] ]
             Double list of selected local inclination methods, the first index (outer list) represents compression or expansion (0 and 1),
             the second index (inner list) denotes the vehicle part index. The size of this inner list should match the number of parts of which the body_shape is composed.
             The int defining the method type is interpreted as follows.
             For the compression methods, the following are available:
             *  0: Newtonian Method.
             *  1: Modified Newtonian.
             *  2 and 3: not available at this moment.
             *  4: Tangent-wedge method.
             *  5: Tangent-cone method.
             *  6: Modified Dahlem-Buck method.
             *  7: VanDyke unified pressure method.
             *  8: Smyth Delta Wing method.
             *  9: Hankey flat surface method
             The expansion method has the following options:
             *  0: Vacuum Pressure coefficient method.
             *  1: Zero Pressure function.
             *  4: High Mach base pressure method.
             *  3 or 5: Prandtl-Meyer method.
             *  6: ACM empirical pressure coefficient.

         reference_area : float
             Reference area used to non-dimensionalize aerodynamic forces and moments.

         moment_reference_point : numpy.ndarray
             Reference point wrt which aerodynamic moments are calculated.

         save_pressure_coefficients : bool
             Boolean denoting whether to save the pressure coefficients that are computed to files





     )doc" )
            .def( "clear_data", &ta::HypersonicLocalInclinationAnalysis::clearData );

    py::class_< ta::ControlSurfaceIncrementAerodynamicInterface, std::shared_ptr< ta::ControlSurfaceIncrementAerodynamicInterface > >(
            m, "ControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>" );

    py::class_< ta::CustomControlSurfaceIncrementAerodynamicInterface,
                std::shared_ptr< ta::CustomControlSurfaceIncrementAerodynamicInterface >,
                ta::ControlSurfaceIncrementAerodynamicInterface >(
            m, "CustomControlSurfaceIncrementAerodynamicInterface", "<no_doc, only_dec>" )
            .def( py::init< const std::function< Eigen::Vector6d( const std::vector< double >& ) >,
                            const std::vector< ta::AerodynamicCoefficientsIndependentVariables > >( ),
                  py::arg( "coefficient_function" ),
                  py::arg( "independent_variable_names" ) );

    m.def( "get_default_local_inclination_mach_points",
           &ta::getDefaultHypersonicLocalInclinationMachPoints,
           py::arg( "mach_regime" ) = "Full" );

    m.def( "get_default_local_inclination_angle_of_attack_points", &ta::getDefaultHypersonicLocalInclinationAngleOfAttackPoints );

    m.def( "get_default_local_inclination_sideslip_angle_points", &ta::getDefaultHypersonicLocalInclinationAngleOfSideslipPoints );

    m.def( "save_vehicle_mesh_to_file",
           &ta::saveVehicleMeshToFile,
           py::arg( "local_inclination_analysis_object" ),
           py::arg( "output_directory" ),
           py::arg( "output_file_prefix" ) = "",
           R"doc(

 Function to save the mesh used for a hypersonic local inclination analysis to a file.

 Function to save the mesh used for a hypersonic local inclination analysis to a file. This function saves
 two files to the specified directory, with filenames: "ShapeFile.dat" and "SurfaceNormalFile.dat", where these
 files names may be prefixed by an optional string (see below). The first of these files contains four columns defining
 the surface points that define mesh, with Column 0: point index; Column 1: x-position of point; Column 1: y-position of point;
 Column 2: z-position of point. The second file contains four columns with Column 0: point index; Column 1: x-component of surface normal;
 Column 1: y-position of surface normal; Column 2: z-position of surface normal.


 Parameters
 ----------
 local_inclination_analysis_object : HypersonicLocalInclinationAnalysis
     Object used to calculate the aerodynamics of the vehicle

 output_directory : str
     Directory to which the files are to be saved

 output_file_prefix : str, default=''
     Optional prefix of output file names






     )doc" );

    m.def( "get_local_inclination_total_vehicle_area", &ta::getTotalSurfaceArea, py::arg( "local_inclination_analysis_object" ) );

    m.def( "get_local_inclination_mesh", &ta::getVehicleMesh, py::arg( "local_inclination_analysis_object" ) );

    py::class_< tsm::VehicleSystems, std::shared_ptr< tsm::VehicleSystems > >( m, "VehicleSystems", R"doc(

         Object used to store physical (hardware) properties of a vehicle.






      )doc" )
            .def( py::init<>( ) )
            .def( "set_control_surface_deflection",
                  &tsm::VehicleSystems::setCurrentControlSurfaceDeflection,
                  py::arg( "control_surface_id" ),
                  py::arg( "deflection_angle" ),
                  R"doc(

         Function to set the current deflection of an aerodynamic control surface.


         Function to set the current deflection of an aerodynamic control surface,
         identified by its name. To set the control surface deflection, the control
         surface has to exist. A control surface is created whenever control surfaces are
         defined in a body's aerodynamic coefficient interface.


         Parameters
         ----------
         control_surface_id : str
             The identified (name) of the given control surface

         deflection_angle : float
             The deflection (in radians) that the control surface is to be set to. This will
             typically influence the aerodynamic coefficients of the vehicle





     )doc" )
            .def( "set_transponder_turnaround_ratio",
                  py::overload_cast< std::map< std::pair< tom::FrequencyBands, tom::FrequencyBands >, double >& >(
                          &tsm::VehicleSystems::setTransponderTurnaroundRatio ),
                  py::arg( "transponder_ratio_per_uplink_and_downlink_"
                           "frequency_band" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_default_transponder_turnaround_ratio_function",
                  &tsm::VehicleSystems::setDefaultTransponderTurnaroundRatio,
                  R"doc(Retrieve standard, DSN turnaround ratios based on the frequency bands of the link)doc" )
            .def( "set_transmitted_frequency_calculator",
                  &tsm::VehicleSystems::setTransmittedFrequencyCalculator,
                  py::arg( "transmitted_frequency_calculator" ),
                  R"doc(
                    Set the transmitted frequency calculator for the vehicle.
                    This function assigns a frequency calculator to the vehicle, which can be used
                    to determine the frequency transmitted by an onboard station or system.
                    Parameters
                    ----------
                    transmitted_frequency_calculator : StationFrequencyInterpolator
                        The frequency calculator object to be associated with the vehicle.
                    )doc" )
            .def( "get_control_surface_deflection",
                  &tsm::VehicleSystems::getCurrentControlSurfaceDeflection,
                  py::arg( "control_surface_id" ),
                  R"doc(

         Function to retrieve the current deflection of an aerodynamic control surface.


         Function to retrieve the current deflection of an aerodynamic control surface,
         identified by its name. To extract the control surface deflection, the control
         surface has to exist. A control surface is created whenever control surfaces are
         defined in a body's aerodynamic coefficient interface.


         Parameters
         ----------
         control_surface_id : str
             The identified (name) of the given control surface

         Returns
         -------
         float
             Current deflection (in radians) that the control surface





     )doc" )
            .def( "set_reference_point",
                  py::overload_cast< const std::string, const Eigen::Vector3d&, const std::string, const std::string >(
                          &tsm::VehicleSystems::setReferencePointPosition ),
                  py::arg( "reference_point" ),
                  py::arg( "location" ),
                  py::arg( "frame_origin" ) = "",
                  py::arg( "frame_orientation" ) = "",
                  R"doc(No documentation found.)doc" )
            .def( "set_reference_point",
                  py::overload_cast< const std::string, std::shared_ptr< te::Ephemeris > >(
                          &tsm::VehicleSystems::setReferencePointPosition ),
                  py::arg( "reference_point" ),
                  py::arg( "ephemeris" ),
                  R"doc(No documentation found.)doc" )

            .def( "get_engine_model",
                  &tsm::VehicleSystems::getEngineModel,
                  py::arg( "engine_name" ),
                  R"doc(

         Function to retrieve an engine model from the vehicle



         Parameters
         ----------
         engine_name : str
             The identifier for the engine model that is to be retrieved

         Returns
         -------
         EngineModel
             Model for the engine that is requested





     )doc" )
            .def( "set_timing_system", &tsm::VehicleSystems::setTimingSystem, py::arg( "timing_system" ) );

    py::class_< tss::RigidBodyProperties, std::shared_ptr< tss::RigidBodyProperties > >( m, "RigidBodyProperties", R"doc(

         Object that defines the mass, center of mass, and inertia tensor as a function of time.

         Object that defines the mass, center of mass, and inertia tensor as a function of time, typically used for evaluation of torques and non-conservative forces
         in numerical state propagation. Note that this object does *not* define properties of a gravity field (it defines the inertial mass rather than the gravitational mass)

      )doc" )
            .def( "update",
                  &tss::RigidBodyProperties::update,
                  py::arg( "time" ),
                  R"doc(

         Function to update the body properties to the current time. This function is called automatically during a propagation loop.
         In case these properties are not time-dependent (e.g. when using the :func:`~tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties` settings)
         this function does nothing (since no update is needed).

         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) to which this object is to be updated

         Returns
         -------

      )doc" )
            .def_property_readonly( "current_mass",
                                    &tss::RigidBodyProperties::getCurrentMass,
                                    R"doc(

         Mass of the object, as set by the latest call to the ``update`` function of this object.
      )doc" )
            .def_property_readonly( "current_center_of_mass",
                                    &tss::RigidBodyProperties::getCurrentCenterOfMass,
                                    R"doc(

         Position of the center of mass of the object (in the body-centered, body-fixed frame), as set by the latest call to the ``update`` function of this object.

      )doc" )
            .def_property_readonly( "current_inertia_tensor",
                                    &tss::RigidBodyProperties::getCurrentInertiaTensor,
                                    R"doc(

        Inertia tensor of the object (with axes along those of the body-fixed frame), as set by the latest call to the ``update`` function of this object.

      )doc" );

    py::class_< tsm::TimingSystem, std::shared_ptr< tsm::TimingSystem > >( m,
                                                                           "TimingSystem",
                                                                           R"doc(No documentation found.)doc" )

            .def(  // ctor 1
                    py::init< const std::vector< tudat::Time >,
                              const std::vector< double >,
                              const std::function< std::function< double( const double ) >( const double, const double, const double ) >,
                              const double >( ),
                    py::arg( "arc_times" ),
                    py::arg( "all_arcs_polynomial_drift_coefficients" ) = std::vector< double >( ),
                    py::arg( "clock_noise_generation_function" ) = nullptr,
                    py::arg( "clock_noise_time_step" ) = 1.0E-3 )
            .def(  // ctor 2
                    py::init< const std::vector< tudat::Time >,
                              const std::vector< std::vector< double > >,
                              const std::function< std::function< double( const double ) >( const double, const double, const double ) >,
                              const double >( ),
                    py::arg( "arc_times" ),
                    py::arg( "polynomial_drift_coefficients" ),
                    py::arg( "clock_noise_generation_function" ) = nullptr,
                    py::arg( "clock_noise_time_step" ) = 1.0E-3 )
            .def(  // ctor 3
                    py::init< const std::vector< std::vector< double > >,
                              const std::vector< std::function< double( const double ) > >,
                              const std::vector< tudat::Time > >( ),
                    py::arg( "polynomial_drift_coefficients" ),
                    py::arg( "stochastic_clock_noise_functions" ),
                    py::arg( "arc_times" ) );

    py::class_< tsm::EngineModel, std::shared_ptr< tsm::EngineModel > >( m, "EngineModel" )
            .def_property_readonly( "thrust_magnitude_calculator", &tsm::EngineModel::getThrustMagnitudeWrapper );

    /*!
     **************   FLIGHT CONDITIONS AND ASSOCIATED FUNCTIONALITY
     *******************
     */

    py::class_< trf::AerodynamicAngleCalculator, std::shared_ptr< trf::AerodynamicAngleCalculator > >(
            m, "AerodynamicAngleCalculator", R"doc(

         Object to calculate (aerodynamic) orientation angles, and frame transformations,
         from current vehicle state.


         Object to calculate (aerodynamic) orientation angles (list given by the :class:`~AerodynamicsReferenceFrameAngles` enum)
         and transformations between frames (list given by the :class:`~AerodynamicsReferenceFrames` enum) from current vehicle state.





      )doc" )
            .def( "get_rotation_matrix_between_frames",
                  &trf::AerodynamicAngleCalculator::getRotationMatrixBetweenFrames,
                  py::arg( "original_frame" ),
                  py::arg( "target_frame" ),
                  R"doc(

         Function to get the rotation matrix between two frames.


         Function to get the rotation matrix between two frames. This function
         is meant to be used only *during* a numerical propagation, in particular
         for the definition of a custom (e.g. guidance) model.


         Parameters
         ----------
         original_frame : AerodynamicsReferenceFrames
             The frame :math:`A` from which the rotation matrix is to be calculated

         target_frame : AerodynamicsReferenceFrames
             The frame :math:`B` to which the rotation matrix is to be calculated

         Returns
         -------
         numpy.ndarray
             Rotation matrix :math:`\mathbf{R}^{B/A}` from frame :math:`A` to frame `B`





     )doc" )
            .def( "get_angle",
                  &trf::AerodynamicAngleCalculator::getAerodynamicAngle,
                  py::arg( "angle_type" ),
                  R"doc(

         Function to get a single orientation angle


         Function to get a single orientation angle. This function
         is meant to be used only *during* a numerical propagation, in particular
         for the definition of a custom (e.g. guidance) model.


         Parameters
         ----------
         original_frame : AerodynamicsReferenceFrameAngles
             The identifier for the angle that is to be returned

         Returns
         -------
         double
             Value of requested angle





     )doc" )
            // Function removed; error is shown
            .def( "set_body_orientation_angles",
                  &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved2,
                  py::arg( "angle_of_attack" ) = TUDAT_NAN,
                  py::arg( "angle_of_sideslip" ) = TUDAT_NAN,
                  py::arg( "bank_angle" ) = TUDAT_NAN,
                  py::arg( "silence_warnings" ) = false )
            // Function removed; error is shown
            .def( "set_body_orientation_angle_functions",
                  &trf::AerodynamicAngleCalculator::setOrientationAngleFunctionsRemoved1,
                  py::arg( "angle_of_attack_function" ) = std::function< double( ) >( ),    // <pybind11/functional.h>
                  py::arg( "angle_of_sideslip_function" ) = std::function< double( ) >( ),  // <pybind11/functional.h>
                  py::arg( "bank_angle_function" ) = std::function< double( ) >( ),         // <pybind11/functional.h>
                  py::arg( "angle_update_function" ) = std::function< void( const double ) >( ),
                  py::arg( "silence_warnings" ) = false );

    py::class_< ta::FlightConditions, std::shared_ptr< ta::FlightConditions > >( m, "FlightConditions", R"doc(

         Object that calculates various state-derived quantities typically
         relevant for flight dynamics.


         Object that calculates various state-derived quantities typically
         relevant for flight dynamics, such as latitude, longitude,
         altitude, etc. It also contains an
         :py:class:`~AerodynamicAngleCalculator` that computes derived
         angles (flight path, heading angle, etc.). This object is limited
         to non-atmospheric flight. For flight through Body objects
         endowed with an atmosphere model, the derived class
         :py:class:`~AtmosphericFlightConditions` is used. This object is
         stored inside a Body object, and represents the flight conditions
         of a single body w.r.t. a single central body.





      )doc" )
            //            .def(py::init<
            //                 const
            //                 std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
            //                 const
            //                 std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>(),
            //                 py::arg("shape_model"),
            //                 py::arg("aerodynamic_angle_calculator") =
            //                 std::shared_ptr<
            //                 tr::AerodynamicAngleCalculator>())
            .def( "update_conditions", &ta::FlightConditions::updateConditions, py::arg( "current_time" ) )
            .def_property_readonly( "aerodynamic_angle_calculator",
                                    &ta::FlightConditions::getAerodynamicAngleCalculator,
                                    R"doc(

         **read-only**

         The object that is responsible for computing various relevant
         flight dynamics angles and frame rotations.


         :type: AerodynamicAngleCalculator
      )doc" )
            .def_property_readonly( "longitude",
                                    &ta::FlightConditions::getCurrentLongitude,
                                    R"doc(

         **read-only**

         The body-fixed longitude of the body w.r.t. its central body.


         :type: float
      )doc" )
            .def_property_readonly( "latitude",
                                    &ta::FlightConditions::getCurrentLatitude,
                                    R"doc(

         **read-only**

         The body-fixed geographic latitude of the body w.r.t. its
         central body.


         :type: float
      )doc" )
            .def_property_readonly( "geodetic_latitude",
                                    &ta::FlightConditions::getCurrentGeodeticLatitude,
                                    R"doc(

         **read-only**

         The body-fixed geographic latitude of the body w.r.t. its
         central body.


         :type: float
      )doc" )
            .def_property_readonly( "time", &ta::FlightConditions::getCurrentTime, R"doc(

         **read-only**

         The current time, at which this object was last updated


         :type: float
      )doc" )
            .def_property_readonly( "body_centered_body_fixed_state",
                                    &ta::FlightConditions::getCurrentBodyCenteredBodyFixedState,
                                    R"doc(

         **read-only**

         Cartesian translational state, expressed in a frame centered
         at, and fixed to, the central body. Note that, due to the
         rotation of the central body, the norm of the body-fixed,
         body-centered, velocity differs from the norm of the inertial
         body-centered velocity.


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "altitude",
                                    &ta::FlightConditions::getCurrentAltitude,
                                    R"doc(

         **read-only**

         The current time, at which this object was last updated


         :type: float
      )doc" );

    py::class_< ta::AtmosphericFlightConditions, std::shared_ptr< ta::AtmosphericFlightConditions >, ta::FlightConditions >(
            m, "AtmosphericFlightConditions", R"doc(

         Object that calculates various state-derived quantities typically
         relevant for flight dynamics, for flight in an atmosphere.


         Object that calculates various state-derived quantities typically
         relevant for flight dynamics, for flight in an atmosphere, such
         as latitude,  longitude, altitude, density, Mach number etc. It
         also contains an ``AerodynamicAngleCalculator`` that computes
         derived angles (flight path, heading angle, etc.). This object is
         derived from ``FlightConditions``, which performs computations for
         non-atmospheric flight only. This object is stored inside a Body
         object, and represents the flight conditions of a single body
         w.r.t. a single central body.





      )doc" )
            .def_property_readonly( "density",
                                    &ta::AtmosphericFlightConditions::getCurrentDensity,
                                    R"doc(

         **read-only**

         The freestream atmospheric density at the body's current
         location.


         :type: float
      )doc" )
            .def_property_readonly( "temperature",
                                    &ta::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                                    R"doc(

         **read-only**

         The freestream atmospheric temperature at the body's current
         location.


         :type: float
      )doc" )
            .def_property_readonly( "dynamic_pressure",
                                    &ta::AtmosphericFlightConditions::getCurrentDynamicPressure,
                                    R"doc(

         **read-only**

         The freestream atmospheric dynamic pressure at the body's
         current location.


         :type: float
      )doc" )
            .def_property_readonly( "pressure",
                                    &ta::AtmosphericFlightConditions::getCurrentPressure,
                                    R"doc(

         **read-only**

         The freestream atmospheric static pressure at the body's
         current location.


         :type: float
      )doc" )
            .def_property_readonly( "airspeed",
                                    &ta::AtmosphericFlightConditions::getCurrentAirspeed,
                                    R"doc(

         **read-only**

         The airspeed of the body w.r.t. the atmosphere.


         :type: float
      )doc" )
            .def_property_readonly( "mach_number",
                                    &ta::AtmosphericFlightConditions::getCurrentMachNumber,
                                    R"doc(

         **read-only**

         The freestream Mach number of the body.


         :type: float
      )doc" )
            .def_property_readonly( "airspeed_velocity",
                                    &ta::AtmosphericFlightConditions::getCurrentAirspeedBasedVelocity,
                                    R"doc(

         **read-only**

         The velocity vector of the body w.r.t. the freestream
         atmosphere (e.g. vectorial counterpart of airspeed).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "speed_of_sound",
                                    &ta::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                                    R"doc(

         **read-only**

         The freestream atmospheric speed of sound at the body's current
         location.


         :type: float
      )doc" )
            .def_property_readonly( "aero_coefficient_independent_variables",
                                    &ta::AtmosphericFlightConditions::getAerodynamicCoefficientIndependentVariables,
                                    R"doc(

         **read-only**

         List of current values of independent variables of aerodynamic
         coefficients. This list is only defined if the body has an
         :py:class:`~AerodynamicCoefficientInterface` that has
         dependencies on environmental variables (e.g. Mach number,
         angle of attack, etc.).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly(
                    "control_surface_aero_coefficient_independent_"
                    "variables",
                    &ta::AtmosphericFlightConditions::getControlSurfaceAerodynamicCoefficientIndependentVariables,
                    R"doc(

         **read-only**

         List of lists current values of independent variables of
         aerodynamic coefficients for control surfaces. The outer list
         defines the control surface, the inner list the values of the
         independent variables. This list is only defined if the body
         has an :py:class:`~AerodynamicCoefficientInterface` with
         control surfaces that have dependencies on environmental
         variables (e.g. Mach number, angle of attack, etc.).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "aerodynamic_coefficient_interface",
                                    &ta::AtmosphericFlightConditions::getAerodynamicCoefficientInterface,
                                    R"doc(

         **read-only**

         Object extracted from the same Body object as this
         :py:class:`~AtmosphericFlightConditions` object, which defines
         the aerodynamic coefficients.


         :type: AerodynamicCoefficientInterface
      )doc" );

    /*!
     **************   ROTATION MODELS  ******************
     */

    py::class_< te::RotationalEphemeris, std::shared_ptr< te::RotationalEphemeris > >( m, "RotationalEphemeris", R"doc(

         Object that stores the rotational state of the bodies.


         Object that stores the rotational state of the bodies. This object can be used to calculate rotation matrices,
         which are used to transform coordinates between reference frames.





      )doc" )
            .def( "body_fixed_to_inertial_rotation",
                  &te::RotationalEphemeris::getRotationMatrixToBaseFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get rotation matrix from body-fixed frame to inertial frame over time.


         Function to get rotation matrix from body-fixed (target) frame to inertial (base) frame over time.
         The calculation of this rotation matrix depends on the specific rotation model that has been defined,
         either from an a priori definition (see :ref:`rotation_model` submodule) or from processing
         the results of propagation of the rotational equations of motion.


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the rotation matrix is evaluated

         Returns
         -------
         numpy.ndarray
             Rotation matrix :math:`\mathbf{R}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`





     )doc" )
            .def( "time_derivative_body_fixed_to_inertial_rotation",
                  &te::RotationalEphemeris::getDerivativeOfRotationToBaseFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time.


         Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time (see ``body_fixed_to_inertial_rotation``),
         denoted :math:`\dot{\mathbf{R}}^{(I/B)}`,


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the rotation matrix derivative is evaluated

         Returns
         -------
         numpy.ndarray
             Rotation matrix :math:`\dot{\mathbf{R}}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`





     )doc" )
            .def( "inertial_to_body_fixed_rotation",
                  &te::RotationalEphemeris::getRotationMatrixToTargetFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get rotation matrix from inertial frame to body-fixed frame over time.


         Function computes the inverse (equal to transpose) rotation of the ``body_fixed_to_inertial_rotation`` function.


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the rotation matrix is evaluated





     )doc" )
            .def( "time_derivative_inertial_to_body_fixed_rotation",
                  &te::RotationalEphemeris::getDerivativeOfRotationToTargetFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time.


         Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time (see ``inertial_to_body_fixed_rotation``),
         denoted :math:`\dot{\mathbf{R}}^{(B/I)}`,


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the rotation matrix derivative is evaluated

         Returns
         -------
         numpy.ndarray
             Rotation matrix :math:`\dot{\mathbf{R}}^{B/I}` from inertial frame `I` to body-fixed frame :math:`B`





     )doc" )
            .def( "angular_velocity_in_body_fixed_frame",
                  &te::RotationalEphemeris::getRotationalVelocityVectorInTargetFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get the body's angular velocity vector, expressed in the body-fixed frame.


         Function to get the body's angular velocity vector :math:`\boldsymbol{\omega}^{(B)}`, expressed in the body-fixed frame :math:`B`.
         The calculation of the angular velocity depends on the specific rotation model that has been defined,
         either from an a priori definition (see :ref:`rotation_model` submodule) or from processing
         the results of propagation of the rotational equations of motion.
         Note that when numerically propagating rotational dynamics, this angular velocity vector is typically directly defined
         in the last three entries of the state vector.


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the angular velocity vector is evaluated

         Returns
         -------
         numpy.ndarray
             Angular velocity vector of the body  :math:`\boldsymbol{\omega}^{(B)}` expressed in the body-fixed frame :math:`B`





     )doc" )
            .def( "angular_velocity_in_inertial_frame",
                  &te::RotationalEphemeris::getRotationalVelocityVectorInBaseFrame,
                  py::arg( "time" ),
                  R"doc(

         Function to get the body's angular velocity vector, expressed in the inertial frame.


         Function to get the body's angular velocity vector :math:`\boldsymbol{\omega}^{(I)}`, expressed in the body-fixed frame :math:`I`.
         This quantity is computed from :math:`\mathbf{R}^{I/B}\boldsymbol{\omega}^{(B)}`, see the ``angular_velocity_in_body_fixed_frame`` and
         ``body_fixed_to_inertial_rotation`` functions.


         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the angular velocity vector is evaluated

         Returns
         -------
         numpy.ndarray
             Angular velocity vector of the body  :math:`\boldsymbol{\omega}^{(B)}` expressed in the body-fixed frame :math:`B`





     )doc" )
            .def_property_readonly( "body_fixed_frame_name",
                                    &te::RotationalEphemeris::getTargetFrameOrientation,
                                    R"doc(

         **read-only**

         The identifier of the body-fixed frame, used in other parts of the simulation to identify it.


         :type: str
      )doc" )
            .def_property_readonly( "inertial_frame_name",
                                    &te::RotationalEphemeris::getBaseFrameOrientation,
                                    R"doc(

         **read-only**

         The identifier of the inertial frame, used in other parts of the simulation to identify it.


         :type: str
      )doc" );

    m.def( "transform_to_inertial_orientation",
           &te::transformStateToInertialOrientation< double, double >,
           py::arg( "state_in_body_fixed_frame" ),
           py::arg( "current_time" ),
           py::arg( "rotational_ephemeris" ),
           R"doc(

 Function to convert a Cartesian state vector from a body-fixed to an inertial frame

 Function to convert a Cartesian state vector from a body-fixed to an inertial frame, using a :class:`~tudatpy.dynamics.environment.RotationalEphemeris`
 object as a model for the rotation. The body-fixed frame from which the conversion takes place is the :attr:`~tudatpy.dynamics.environment.RotationalEphemeris.body_fixed_frame_name` frame,
 the (assumedly) inertial frame to which the conversion is done is :attr:`~tudatpy.dynamics.environment.RotationalEphemeris.inertial_frame_name`.

 This function calls :func:`~tudatpy.astro.element_conversion.rotate_state_to_frame` (with frame :math:`A` the inertial frame, and frame :math:`B` the body-fixed frame). The present function
 computes the required rotation matrix and its time derivative from the ``rotational_ephemeris`` input given here.

 Parameters
 ----------
 state_in_body_fixed_frame : numpy.ndarray[numpy.float64[6, 1]]
     Cartesian state (position and velocity) in the body-fixed frame

 current_time : astro.time_representation.Time
     Time object representing seconds since J2000 (TDB) at which the transformation is to be computed

 rotational_ephemeris : RotationalEphemeris
     Boy rotation model that is to be used to convert the body-fixed state to inertial state

 Returns
 -------
 numpy.ndarray[numpy.float64[6, 1]]
     Cartesian state transformed to inertial frame, using ``rotational_ephemeris`` model, from body-fixed ``state_in_body_fixed_frame``






     )doc" );

    py::class_< te::LongitudeLibrationCalculator, std::shared_ptr< te::LongitudeLibrationCalculator > >( m,
                                                                                                         "LongitudeLibrationCalculator" );

    py::class_< te::DirectLongitudeLibrationCalculator,
                std::shared_ptr< te::DirectLongitudeLibrationCalculator >,
                te::LongitudeLibrationCalculator >( m, "DirectLongitudeLibrationCalculator" )
            .def( py::init< const double >( ), py::arg( "scaled_libration_amplitude" ) );

    py::class_< te::SynchronousRotationalEphemeris, std::shared_ptr< te::SynchronousRotationalEphemeris >, te::RotationalEphemeris >(
            m, "SynchronousRotationalEphemeris" )
            .def_property( "libration_calculator",
                           &te::SynchronousRotationalEphemeris::getLongitudeLibrationCalculator,
                           &te::SynchronousRotationalEphemeris::setLibrationCalculation );

    py::class_< te::AerodynamicAngleRotationalEphemeris,
                std::shared_ptr< te::AerodynamicAngleRotationalEphemeris >,
                te::RotationalEphemeris >( m, "AerodynamicAngleRotationalEphemeris" )
            .def( "reset_aerodynamic_angle_function", &te::AerodynamicAngleRotationalEphemeris::setAerodynamicAngleFunction );

    py::class_< teo::EarthOrientationAnglesCalculator, std::shared_ptr< teo::EarthOrientationAnglesCalculator > >(
            m, "EarthOrientationAnglesCalculator", R"doc(

         Object for computing high-accuracy Earth orientation angles

     )doc" )
            .def( "get_gcrs_to_itrs_rotation_angles",
                  &teo::EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs< TIME_TYPE >,
                  py::arg( "epoch" ),
                  py::arg( "time_scale" ) = tba::tdb_scale,
                  R"doc(

         Function to compute high-accuracy Earth orientation angles

         Function to compute high-accuracy Earth orientation angle quantities :math:`X,Y,s,x_{p},y_{p}` and UT1 (from which :math:`\theta_{E}` is computed)
         as described in :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`

         Parameters
         ----------
         epoch : astro.time_representation.Time
             Time object representing seconds since J2000 at which the Earth orientation angles are to be compute
         time_scale : TimeScales
             Time scale in which the input epoch is given

         Returns
         -------
         tuple[list[float],float]
             Pair (tuple of size two) with the first entry a list of orientation angles :math:`X,Y,s,x_{p},y_{p}` (in that order) and the second entry the current UT1.

     )doc" );

    py::class_< te::GcrsToItrsRotationModel, std::shared_ptr< te::GcrsToItrsRotationModel >, te::RotationalEphemeris >(
            m, "GcrsToItrsRotationModel", R"doc(

         Object for high-accuracy GCRS<->ITRS rotation.

         Object derived from :class:`~RotationalEphemeris` that implements the high-accuracy GCRS<->ITRS rotation as per the IERS 2010 Conventions. The details of the model are described in
         :func:`~tudatpy.dynamics.environment_setup.rotation_model.gcrs_to_itrs`
         With the exception of :math:`s'`, the list of angles used to compute the full rotation are computed by an object of type :class:`~EarthOrientationAnglesCalculator` (which can be retrieved from this rotation model
         through :attr:`~GcrsToItrsRotationModel.angles_calculator`.

     )doc" )
            .def_property_readonly( "angles_calculator",
                                    &te::GcrsToItrsRotationModel::getAnglesCalculator,
                                    R"doc(

         **read-only**

         Object that computes the Earth rotation angles :math:`X,Y,s,\theta_{E},x_{p},y_{p}`


         :type: EarthOrientationAnglesCalculator

     )doc" );

    py::class_< te::DirectionBasedRotationalEphemeris, std::shared_ptr< te::DirectionBasedRotationalEphemeris >, te::RotationalEphemeris >(
            m, "CustomInertialDirectionBasedRotationalEphemeris" )
            .def_property_readonly( "inertial_body_axis_calculator",
                                    &te::DirectionBasedRotationalEphemeris::getInertialBodyAxisDirectionCalculator );

    py::class_< te::InertialBodyFixedDirectionCalculator, std::shared_ptr< te::InertialBodyFixedDirectionCalculator > >(
            m, "InertialBodyFixedDirectionCalculator" );

    py::class_< te::CustomBodyFixedDirectionCalculator,
                std::shared_ptr< te::CustomBodyFixedDirectionCalculator >,
                te::InertialBodyFixedDirectionCalculator >( m, "CustomBodyFixedDirectionCalculator" )
            .def_property( "inertial_body_axis_direction_function",
                           &te::CustomBodyFixedDirectionCalculator::getInertialBodyAxisDirectionFunction,
                           &te::CustomBodyFixedDirectionCalculator::resetInertialBodyAxisDirectionFunction );

    /*!
     **************   GRAVITY FIELD  ******************
     */

    py::class_< tg::GravityFieldModel, std::shared_ptr< tg::GravityFieldModel > >( m,
                                                                                   "GravityFieldModel",
                                                                                   R"doc(

         Object that provides the gravity field of a body


         Object (typically stored inside a :class:`~Body` object) that provides the gravity field of a body, typically (but not exclusively) for
         use in gravitational acceleration and torque models. This base class allows access to the gravitational parameter of the body.
         Specific derived classes are implemented to provide models for more detailed gravity field models (e.g. spherical harmonics, polyhedron).




     )doc" )
            .def( py::init< const double, const std::function< void( ) > >( ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "update_inertia_tensor" ) = std::function< void( ) >( )  // <pybind11/functional.h>
                  )
            .def( "get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter )
            .def_property( "gravitational_parameter",
                           &tg::GravityFieldModel::getGravitationalParameter,
                           &tg::GravityFieldModel::resetGravitationalParameter,
                           R"doc(

         Value of the gravity field's gravitational parameters :math:`\mu`


         :type: float





      )doc" );

    py::class_< tg::SphericalHarmonicsGravityField, std::shared_ptr< tg::SphericalHarmonicsGravityField >, tg::GravityFieldModel >(
            m, "SphericalHarmonicsGravityField", R"doc(

         Object that provides a spherical harmonic gravity field of a body.

         Object (typically stored inside a :class:`~Body` object) that provides a spherical harmonic gravity field of a body, typically (but not exclusively) for
         use in gravitational acceleration and torque models. This class is derived from :class:`~GravityFieldModel`.  This object is typically created using the :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`
         settings function. If any time variations of the gravity field are provided, an object of the derived class :class:`~TimeVariableSphericalHarmonicsGravityField` is created.

     )doc" )
            .def_property_readonly( "reference_radius",
                                    &tg::SphericalHarmonicsGravityField::getReferenceRadius,
                                    R"doc(

         **read-only**

         Reference radius :math:`R` of the gravity field

         :type: float
      )doc" )
            .def_property_readonly( "maximum_degree",
                                    &tg::SphericalHarmonicsGravityField::getDegreeOfExpansion,
                                    R"doc(

         **read-only**

         Maximum spherical harmonic degree :math:`l_{max}` for which coefficients are defined

         :type: int
      )doc" )
            .def_property_readonly( "maximum_order",
                                    &tg::SphericalHarmonicsGravityField::getOrderOfExpansion,
                                    R"doc(

         **read-only**

         Maximum spherical harmonic order :math:`m_{max}` for which coefficients are defined

         :type: int
      )doc" )
            .def_property( "cosine_coefficients",
                           &tg::SphericalHarmonicsGravityField::getCosineCoefficients,
                           &tg::SphericalHarmonicsGravityField::setCosineCoefficients,
                           R"doc(

         Matrix with cosine spherical harmonic coefficients :math:`\bar{C}_{lm}` (geodesy normalized). Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.

         :type: numpy.ndarray[numpy.float64[l, m]]
      )doc" )
            .def_property( "sine_coefficients",
                           &tg::SphericalHarmonicsGravityField::getSineCoefficients,
                           &tg::SphericalHarmonicsGravityField::setSineCoefficients,
                           R"doc(

         Matrix with sine spherical harmonic coefficients :math:`\bar{S}_{lm}` (geodesy normalized). Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.

         :type: numpy.ndarray[numpy.float64[l, m]]
      )doc" );

    py::class_< tg::TimeDependentSphericalHarmonicsGravityField,
                std::shared_ptr< tg::TimeDependentSphericalHarmonicsGravityField >,
                tg::SphericalHarmonicsGravityField >( m, "TimeDependentSphericalHarmonicsGravityField", R"doc(

            Derived class of :class:`~SphericalHarmonicsGravityField` that is created when any gravity field variations are detected.

            Derived class of :class:`~SphericalHarmonicsGravityField` that is created instead when any gravity field variations are detected during object creation
            (typically during a call of :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`)
            This object computes the time-variability of spherical harmonic coefficients from a list of :class:`~GravityFieldVariationModel` objects.
            The ``cosine_coefficients`` and ``sine_coefficients`` attributes provide the instantaneous coefficients (including the time-variability)
            The ``nominal_cosine_coefficients`` and ``nominal_sine_coefficients`` provide the static (e.g. without time-variations) coefficients.

            )doc" )

            .def_property( "nominal_cosine_coefficients",
                           &tg::TimeDependentSphericalHarmonicsGravityField::getNominalCosineCoefficients,
                           &tg::TimeDependentSphericalHarmonicsGravityField::setNominalCosineCoefficients,
                           R"doc(

         Matrix with cosine spherical harmonic coefficients :math:`\bar{C}_{lm}` (geodesy normalized) *excluding* time-variations. Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.

         :type: numpy.ndarray[numpy.float64[l, m]]
      )doc" )
            .def_property( "nominal_sine_coefficients",
                           &tg::TimeDependentSphericalHarmonicsGravityField::getNominalSineCoefficients,
                           &tg::TimeDependentSphericalHarmonicsGravityField::setNominalSineCoefficients,
                           R"doc(

         Matrix with sine spherical harmonic coefficients :math:`\bar{S}_{lm}` (geodesy normalized) *excluding* time-variations. Entry :math:`(i,j)` denotes coefficient at degree :math:`i` and order :math:`j`.

         :type: numpy.ndarray[numpy.float64[l, m]]
      )doc" )

            .def_property_readonly( "gravity_field_variation_models",
                                    &tg::TimeDependentSphericalHarmonicsGravityField::getGravityFieldVariations,
                                    R"doc(

         **read-only**

         List of gravity field variation models that the object uses to update the spherical harmonic coefficients at every time step

         :type: list[GravityFieldVariationModel]
        )doc" );

    py::class_< tg::PolyhedronGravityField, std::shared_ptr< tg::PolyhedronGravityField >, tg::GravityFieldModel >(
            m, "PolyhedronGravityField" )
            .def_property_readonly( "volume", &tg::PolyhedronGravityField::getVolume )
            .def_property_readonly( "vertices_coordinates", &tg::PolyhedronGravityField::getVerticesCoordinates )
            .def_property_readonly( "vertices_defining_each_facet", &tg::PolyhedronGravityField::getVerticesDefiningEachFacet );

    py::class_< tg::GravityFieldVariations, std::shared_ptr< tg::GravityFieldVariations > >( m, "GravityFieldVariationModel", R"doc(

        Object that computes a single type of gravity field variation.

        Object that computes a single type of gravity field variation. This object is typically not used directly, but internally by the :class:`~TimeDependentSphericalHarmonicsGravityField` class.
    )doc" );

    /*!
     **************   RADIATION MODELS  ******************
     */
    py::class_< tem::RadiationPressureTargetModel, std::shared_ptr< tem::RadiationPressureTargetModel > >( m,
                                                                                                           "RadiationPressureTargetModel" );

    py::class_< tem::CannonballRadiationPressureTargetModel,
                std::shared_ptr< tem::CannonballRadiationPressureTargetModel >,
                tem::RadiationPressureTargetModel >( m, "CannonballRadiationPressureTargetModel" )
            .def_property( "radiation_pressure_coefficient",
                           &tem::CannonballRadiationPressureTargetModel::getCoefficient,
                           &tem::CannonballRadiationPressureTargetModel::resetCoefficient );

    py::class_< tem::RadiationSourceModel, std::shared_ptr< tem::RadiationSourceModel > >( m, "RadiationSourceModel" );
    /*!
     **************   SHAPE MODELS  ******************
     */

    py::class_< tba::BodyShapeModel, std::shared_ptr< tba::BodyShapeModel > >( m,
                                                                               "BodyShapeModel",
                                                                               R"doc(

         Object that provides a shape model for a natural body.

         Object (typically stored inside a :class:`~Body` object) that provides a shape model for a body, for instance to compute the altitude from a body-centered state, or w.r.t. which
         to place ground stations. This shape model is typically only associated with natural bodies. Shape models for spacecraft (for non-conservative force models) use properties stored inside the
         :class:`~VehicleSystems` object.

     )doc" )
            .def_property_readonly( "average_radius", &tba::BodyShapeModel::getAverageRadius, R"doc(

         **read-only**

         Average radius of the body, for use in computations that assume a spherical body shape.

         :type: float


     )doc" );

    /*!
     **************   GROUND STATION FUNCTIONALITY
     *******************
     */

    py::class_< tgs::GroundStationState, std::shared_ptr< tgs::GroundStationState > >( m, "GroundStationState", R"doc(

         Object that performs computations of the current (body-fixed) position and frame conversions of the ground station.

         Object that performs computations of the current (body-fixed) position and frame conversions of the ground station. In the simplest situation,
         only a Cartesian position is provided, which is then assumed constant. If time variations (for instance due to tides or plate motion) are present,
         their impact on station position is computed in this object.

     )doc" )
            .def( "get_cartesian_state",
                  &tgs::GroundStationState::getCartesianStateInTime,
                  py::arg( "current_time" ),
                  py::arg( "target_frame_origin" ) = "" )
            .def( "get_cartesian_position",
                  &tgs::GroundStationState::getCartesianPositionInTime,
                  py::arg( "current_time" ),
                  py::arg( "target_frame_origin" ) = "",
                  R"doc(

         This function computes the position of the station as a function of time.

         This function computes the position of the station as a function of time, in a frame with body-fixed orientation.
         Some time-variations of the station position depend on the *origin* of the frame in which the computation is to be
         used. For instance, relativistic correction to the Earth-fixed position is different in a geocentric or barycentric frame.
         However, the output of this function is always given in the body-fixed, body-centered frame.

         Parameters
         ----------
         current_time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the position is to be computed.

         target_frame_origin: str, default = ""
             Identifier for the frame origin w.r.t. which the computed position is to be used.

         Returns
         -------
         numpy.ndarray
             Cartesian position of the station at the current epoch, in a body-centered, body-fixed frame

     )doc" )
            .def_property_readonly( "cartesian_position_at_reference_epoch",
                                    &tgs::GroundStationState::getNominalCartesianPosition,
                                    R"doc(


         **read-only**

         Cartesian position of the ground station, at the reference epoch, in a body-fixed, body-centered frame.

         :type: numpy.ndarray[numpy.float64[3, 1]]

     )doc" )
            .def_property_readonly( "spherical_position_at_reference_epoch",
                                    &tgs::GroundStationState::getNominalSphericalPosition,
                                    R"doc(

         **read-only**

         Spherical position of the ground station (distance w.r.t. body center, latitude, longitude), at the reference epoch, in a body-fixed, body-centered frame.

         :type: numpy.ndarray[numpy.float64[3, 1]]

     )doc" )
            .def_property_readonly( "geodetic_position_at_reference_epoch",
                                    &tgs::GroundStationState::getNominalGeodeticPosition,
                                    R"doc(

         **read-only**

         Geodetic position of the ground station (altitude w.r.t. body shape model, geodetic latitude, longitude), at the reference epoch, in a body-fixed, body-centered frame.

         :type: numpy.ndarray[numpy.float64[3, 1]]

     )doc" )
            .def_property_readonly( "rotation_matrix_body_fixed_to_topocentric",
                                    &tgs::GroundStationState::getConstantRotationMatrixFromBodyFixedToTopocentricFrame,
                                    R"doc(

         **read-only**

         Rotation matrix from the body-fixed frame (of the station's central body) to the topocentric frame
         of the ground station. The body-fixed frame is defined by the rotation model of the body object (:attr:`~tudatpy.dynamics.environment.Body.rotation_model`).
         The axes of the topocentric frame are defined such that the x-axis is in East direction, the z-direction is upwards, perpendicular to the body's surface sphere
         (with properties defined by the central body's shape model :attr:`~tudatpy.dynamics.environment.Body.shape_model`). The y-axis completes the frame, and is in northern direction.
         For time-varying ground station positions, this function uses the station position at reference epoch for the computation of the axes.

         :type: numpy.ndarray[numpy.float64[3, 3]]


     )doc" );

    py::class_< tgs::GroundStation, std::shared_ptr< tgs::GroundStation > >( m, "GroundStation", R"doc(

         Object used to define and store properties of a ground station.

         Object (typically stored inside a :class:`~Body` object) used to define and store properties of a ground station, typically used in modelling tracking observations to/from a ground station.

     )doc" )
            .def( "set_transmitting_frequency_calculator",
                  &tgs::GroundStation::setTransmittingFrequencyCalculator,
                  py::arg( "transmitting_frequency_calculator" ) )
            //            .def( "set_water_vapor_partial_pressure_function",
            //                  &tgs::GroundStation::setWaterVaporPartialPressureFunction,
            //                  py::arg( "water_vapor_partial_pressure_function" ) )
            //            .def( "set_temperature_function",
            //            &tgs::GroundStation::setTemperatureFunction, py::arg(
            //            "temperature_function" ) ) .def( "set_pressure_function",
            //            &tgs::GroundStation::setPressureFunction, py::arg( "pressure_function" ) )
            //            .def( "set_relative_humidity_function",
            //                  &tgs::GroundStation::setRelativeHumidityFunction,
            //                  py::arg( "relative_humidity_function" ) )
            .def_property( "transmitting_frequency_calculator",
                           &tgs::GroundStation::getTransmittingFrequencyCalculator,
                           &tgs::GroundStation::setTransmittingFrequencyCalculator,
                           R"doc(

         Object that provides the transmission frequency as a function of time for (radio) tracking stations. This attribute is typically set automatically when loading tracking data files (e.g. ODF, IFMS, TNF, etc.)

         :type: TransmittingFrequencyCalculator

     )doc" )
            .def_property_readonly( "temperature_function", &tgs::GroundStation::getTemperatureFunction, R"doc(

         Function that provides the local temperature at the ground station (typically use for media corrections) as a function of time

         :type: :type: callable[[float], float]

     )doc" )
            .def_property_readonly( "pressure_function", &tgs::GroundStation::getPressureFunction, R"doc(

         Function that provides the local pressure at the ground station (typically use for media corrections) as a function of time

         :type: :type: callable[[float], float]

     )doc" )
            .def_property_readonly( "relative_humidity_function",
                                    &tgs::GroundStation::getRelativeHumidityFunction,
                                    R"doc(

         Function that provides the local relative humidity at the ground station (typically use for media corrections) as a function of time

         :type: :type: callable[[float], float]

     )doc" )
            .def_property_readonly( "pointing_angles_calculator",
                                    &tgs::GroundStation::getPointingAnglesCalculator,
                                    R"doc(

         **read-only**

         Object that performs computations of the azimuth and elevation of an arbitrary target as observed by the ground station

         :type: PointingAnglesCalculator

     )doc" )
            .def_property_readonly( "station_state", &tgs::GroundStation::getNominalStationState, R"doc(

         **read-only**

         Object that performs computations of the current (body-fixed) position and frame conversions of the ground station.

         :type: GroundStationState

     )doc" )
            .def( "set_timing_system", &tgs::GroundStation::setTimingSystem, py::arg( "timing_system" ) )

            .def( "set_station_meteo_data", &tudat::ground_stations::GroundStation::setMeteoData, py::arg( "meteo_data" ) );

    py::class_< tgs::StationFrequencyInterpolator, std::shared_ptr< tgs::StationFrequencyInterpolator > >(
            m, "TransmittingFrequencyCalculator", R"doc(No documentation found.)doc" );

    py::class_< tgs::ConstantFrequencyInterpolator,
                std::shared_ptr< tgs::ConstantFrequencyInterpolator >,
                tgs::StationFrequencyInterpolator >( m, "ConstantTransmittingFrequencyCalculator" )
            .def( py::init< double >( ), py::arg( "frequency" ) );

    py::class_< tgs::PiecewiseLinearFrequencyInterpolator,
                std::shared_ptr< tgs::PiecewiseLinearFrequencyInterpolator >,
                tgs::StationFrequencyInterpolator >( m, "PiecewiseLinearFrequencyInterpolator" )
            .def( py::init< const std::vector< tudat::Time >&,
                            const std::vector< tudat::Time >&,
                            const std::vector< double >&,
                            const std::vector< double >& >( ),
                  py::arg( "start_times" ),
                  py::arg( "end_times" ),
                  py::arg( "ramp_rates" ),
                  py::arg( "start_frequency" ) )
            .def_property_readonly( "start_times", &tgs::PiecewiseLinearFrequencyInterpolator::getStartTimes )
            .def_property_readonly( "end_times", &tgs::PiecewiseLinearFrequencyInterpolator::getEndTimes )
            .def_property_readonly( "ramp_rates", &tgs::PiecewiseLinearFrequencyInterpolator::getRampRates )
            .def_property_readonly( "start_frequencies", &tgs::PiecewiseLinearFrequencyInterpolator::getStartFrequencies )
            .def( "compute_current_frequency",
                  &tgs::PiecewiseLinearFrequencyInterpolator::computeCurrentFrequency< double, tudat::Time >,
                  py::arg( "lookup_time_original" ) );

    py::class_< tgs::PointingAnglesCalculator, std::shared_ptr< tgs::PointingAnglesCalculator > >( m, "PointingAnglesCalculator" )
            .def( "calculate_elevation_angle",
                  py::overload_cast< const Eigen::Vector3d&, const double >(
                          &tgs::PointingAnglesCalculator::calculateElevationAngleFromInertialVector ),
                  py::arg( "inertial_vector_to_target" ),
                  py::arg( "time" ),
                  R"doc(

         Calculate the elevation angle of a target object.

         Parameters
         ----------
         inertial_vector_to_target : numpy.ndarray
             Vector from ground station to target in inertial frame
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which to calculate the angle

         Returns
         -------
         float
             Elevation angle in radians

         )doc" )
            .def( "calculate_azimuth_angle",
                  py::overload_cast< const Eigen::Vector3d&, const double >(
                          &tgs::PointingAnglesCalculator::calculateAzimuthAngleFromInertialVector ),
                  py::arg( "inertial_vector_to_target" ),
                  py::arg( "time" ),
                  R"doc(

         Calculate the azimuth angle of a target object.

         Parameters
         ----------
         inertial_vector_to_target : numpy.ndarray
             Vector from ground station to target in inertial frame
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which to calculate the angle

         Returns
         -------
         float
             Azimuth angle in radians

         )doc" )
            .def( "convert_inertial_vector_to_topocentric",
                  &tgs::PointingAnglesCalculator::convertVectorFromInertialToTopocentricFrame,
                  py::arg( "inertial_vector" ),
                  py::arg( "time" ) );

    py::enum_< tudat::ground_stations::MeteoDataEntries >( m, "MeteoDataEntries" )
            .value( "temperature_meteo_data", tudat::ground_stations::temperature_meteo_data )
            .value( "pressure_meteo_data", tudat::ground_stations::pressure_meteo_data )
            .value( "water_vapor_pressure_meteo_data", tudat::ground_stations::water_vapor_pressure_meteo_data )
            .value( "relative_humidity_meteo_data", tudat::ground_stations::relative_humidity_meteo_data )
            .value( "dew_point_meteo_data", tudat::ground_stations::dew_point_meteo_data )
            .export_values( );

    py::class_< tudat::ground_stations::StationMeteoData, std::shared_ptr< tudat::ground_stations::StationMeteoData > >(
            m, "StationMeteoData" );

    py::class_< tudat::ground_stations::ContinuousInterpolatedMeteoData,
                std::shared_ptr< tudat::ground_stations::ContinuousInterpolatedMeteoData >,
                tudat::ground_stations::StationMeteoData >( m, "ContinuousInterpolatedMeteoData" )
            .def( py::init< std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > >,
                            std::map< tudat::ground_stations::MeteoDataEntries, int > >( ),
                  py::arg( "interpolator" ),
                  py::arg( "vector_entries" ) );

    /*!
     **************   BODY OBJECTS AND ASSOCIATED FUNCTIONALITY
     *******************
     */

    py::class_< tss::Body, std::shared_ptr< tss::Body > >( m,
                                                           "Body",
                                                           R"doc(

         Object that stores the environment properties and current state of
         a single body.


         Object that stores the environment properties and current state
         of a single celestial body (natural or artificial). Each separate
         environment model (gravity field, ephemeris, etc.) is stored as a
         member object in this class. During each time step, the Body gets
         updated to the current time/propagated state, and the current
         properties, in as much as they are time-dependent, can be
         extracted from this object





      )doc" )
            //            .def_property(
            //                    "ephemeris_frame_to_base_frame",
            //                    &tss::Body::getEphemerisFrameToBaseFrame,
            //                    &tss::Body::setEphemerisFrameToBaseFrame )
            .def_property_readonly( "state", &tss::Body::getState, R"doc(

         **read-only**

         The translational state of the Body, as set during the current
         step of the numerical propagation. The translational state
         stored here is always in Cartesian elements, w.r.t. the global
         frame origin, with axes along the global frame orientation. If
         the body's translational state is numerically propagated, this
         property gets extracted from the propagated state vector. If it
         is not propagated, the state is extracted from this body's
         ephemeris. In both cases, any required state transformations
         are automatically applied. Note that this function  is *only*
         valid during the numerical propagation if any aspects of the
         dynamics or dependent variables require the body's state.


         :type: numpy.ndarray
      )doc" )
            .def( "get_ionosphere_model", &tudat::simulation_setup::Body::getIonosphereModel )
            .def_property_readonly( "position",
                                    &tss::Body::getPosition,
                                    R"doc(

         **read-only**

         The translational position of the Body, as set during the
         current step of the numerical propagation
         (see :py:attr:`~state`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "velocity",
                                    &tss::Body::getVelocity,
                                    R"doc(

         **read-only**

         The translational velocity of the Body, as set during the
         current step of the numerical propagation
         (see :py:attr:`~state`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "inertial_to_body_fixed_frame",
                                    &tss::Body::getCurrentRotationMatrixToLocalFrame,
                                    R"doc(

         **read-only**

         The rotation from inertial frame (with global frame
         orientation) to this Body's body-fixed frame. The rotation is
         always returned here as a rotation matrix.  If the body's
         rotational state is numerically propagated, this property gets
         extracted from the propagated state vector. If it is not
         propagated, the state is extracted from this body's rotational
         ephemeris.

         .. note:: This function is **only** valid during the
                   numerical propagation if any aspects of the dynamics
                   or dependent variables require the body's rotational
                   state.


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "body_fixed_to_inertial_frame",
                                    &tss::Body::getCurrentRotationMatrixToGlobalFrame,
                                    R"doc(

         **read-only**

         The rotation from this Body's body-fixed frame to inertial
         frame (see :py:attr:`~inertial_to_body_fixed_frame`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "inertial_to_body_fixed_frame_derivative",
                                    &tss::Body::getCurrentRotationMatrixDerivativeToLocalFrame,
                                    R"doc(

         **read-only**

         Time derivative of rotation matrix from inertial frame to this
         Body's body-fixed frame
         (see :py:attr:`~inertial_to_body_fixed_frame`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "body_fixed_to_inertial_frame_derivative",
                                    &tss::Body::getCurrentRotationMatrixDerivativeToGlobalFrame,
                                    R"doc(

         **read-only**

         Time derivative of rotation matrix from this Body's body-fixed
         frame to inertial frame
         (see :py:attr:`~inertial_to_body_fixed_frame`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "inertial_angular_velocity",
                                    &tss::Body::getCurrentAngularVelocityVectorInGlobalFrame,
                                    R"doc(

         **read-only**

         Angular velocity vector of the body, expressed in inertial
         frame (see :py:attr:`~inertial_to_body_fixed_frame`).


         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "body_fixed_angular_velocity",
                                    &tss::Body::getCurrentAngularVelocityVectorInLocalFrame,
                                    R"doc(

         **read-only**

         Angular velocity vector of the body, expressed in body-fixed
         frame (see :py:attr:`~inertial_to_body_fixed_frame`).


         :type: numpy.ndarray
      )doc" )
            .def_property( "mass", &tss::Body::getBodyMass, &tss::Body::setConstantBodyMass, R"doc(

         The current mass :math:`m` of the vehicle, as used in the calculation of
         non-conservative acceleration. This attribute is a shorthand for accessing the
         mass as computed/stored in the :attr:`Body.rigid_body_properties` attribute. For certain
         types of rigid-body properties, this attribute cannot be used to (re)set the current
         mass. If the body has no    :attr:`Body.rigid_body_properties`, and this function is used to
         set a mass, a new object is automatically created, with settings analogous to the
         the :func:`~tudatpy.dynamics.environment_setup.rigid_body.constant_rigid_body_properties` setting.

         Unlike the attributes containing the state, orientation, angular velocity
         of the Body, this attribute may be used to retrieve the state during the
         propagation *and* to define the mass of a vehicle.


         :type: float
      )doc" )
            .def_property( "inertia_tensor",
                           &tss::Body::getBodyInertiaTensor,
                           py::overload_cast< const Eigen::Matrix3d& >( &tss::Body::setBodyInertiaTensor ),
                           R"doc(

         The current inertia tensor :math:`\mathbf{I}` of the vehicle, as used in the calculation of
         (for instance) the response to torques. This attribute is a shorthand for accessing the
         inertia tensor as computed/stored in the :attr:`~Body.rigid_body_properties` attribute. For certain
         types of rigid-body properties, this attribute cannot be used to (re)set the current
         mass.

         Unlike the attributes containing the state, orientation, angular velocity
         of the Body, this attribute may be used to retrieve the state during the
         propagation *and* to define the mass of a vehicle.


         :type: numpy.ndarray
      )doc" )
            .def( "state_in_base_frame_from_ephemeris",
                  &tss::Body::getStateInBaseFrameFromEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "time" ),
                  R"doc(

         This function returns the body's state, as computed from its ephemeris model (extracted from :attr:`~Body.ephemeris`) at the current time, and (if needed)
         translates this state to the global frame origin. For the case where the origin of the body's ephemeris (extracted from :attr:`~Ephemeris.frame_origin`) is equal to the
         global frame origin of the system of bodies it is in (extracted from :attr:`SystemOfBodies.global_frame_origin`), this function is equal to ``Body.ephemeris.cartesian_state( time )``.
         Where the global frame origin and ephemeris origin is not equal, other bodies' ephemerides are queried as needed to provide this body's state w.r.t. the global frame origin


         Parameters
         ----------
         time : astro.time_representation.Time
             Time object representing seconds since J2000 (TDB) at which the state is to be computed
         Returns
         -------
         numpy.ndarray
             Cartesian state (position and velocity) of the body w.r.t. the global frame origin at the requested time.





     )doc" )
            .def_property( "ephemeris", &tss::Body::getEphemeris, &tss::Body::setEphemeris, R"doc(

         Object defining the ephemeris model of this body, used to calculate its current state as a function of time.
         Depending on the selected type of model, the type of this attribute
         is of type :class:`~Ephemeris`, or a derived class thereof.


         :type: Ephemeris
      )doc" )
            .def_property( "atmosphere_model",
                           &tss::Body::getAtmosphereModel,
                           &tss::Body::setAtmosphereModel,
                           R"doc(

         Object defining the atmosphere model of this body, used to calculate density, temperature, etc. at a given
         state/time. Depending on the selected type of model, the type of this attribute
         is of type :class:`~AtmosphereModel`, or a derived class thereof.


         :type: AtmosphereModel
      )doc" )
            .def_property( "shape_model", &tss::Body::getShapeModel, &tss::Body::setShapeModel, R"doc(

         Object defining the a shape model of this body, used to define the exterior shape of the body, for instance for
         the calculation of vehicle's altitude. Depending on the selected type of model, the type of this attribute
         is of type BodyShapeModel, or a derived class thereof.


         :type: BodyShapeModel
      )doc" )
            .def_property( "gravity_field_model",
                           &tss::Body::getGravityFieldModel,
                           &tss::Body::setGravityFieldModel,
                           R"doc(

         Object defining the a gravity field model of this body, used to define the exterior gravitational potential, and
         its gradient(s). Depending on the selected type of model, the type of this attribute
         is of type GravityFieldModel, or a derived class thereof.


         :type: GravityFieldModel
      )doc" )
            .def_property( "aerodynamic_coefficient_interface",
                           &tss::Body::getAerodynamicCoefficientInterface,
                           &tss::Body::setAerodynamicCoefficientInterface,
                           R"doc(

         Object defining the aerodynamic coefficients of the body (force-only, or force and moment)
         as a function of any number of independent variables. Depending on the selected type of model, the type of this attribute
         is of type AerodynamicCoefficientInterface, or a derived class thereof.


         :type: AerodynamicCoefficientInterface
      )doc" )
            .def_property( "flight_conditions",
                           &tss::Body::getFlightConditions,
                           &tss::Body::setFlightConditions,
                           R"doc(

         Object used to calculated and store the current flight conditions of a vehicle (altitude, latitude, longitude,
         flight-path angle, etc.) w.r.t. a central body. In case the central body contains an atmosphere, this object
         also stores current local density, Mach number, etc. This object is typically used for aerodynamic accelerations,
         guidance models or other central-body-related custom models.


         :type: FlightConditions
      )doc" )
            .def_property( "rotation_model",
                           &tss::Body::getRotationalEphemeris,
                           &tss::Body::setRotationalEphemeris,
                           R"doc(

         Object defining the orientation of the body, used to calculate the rotation to/from a body-fixed
         frame (and its derivate). Depending on the selected type of model, the type of this attribute
         is of type RotationalEphemeris, or a derived class thereof.


         :type: RotationalEphemeris
      )doc" )
            .def_property( "system_models",
                           &tss::Body::getVehicleSystems,
                           &tss::Body::setVehicleSystems,
                           R"doc(

         Object used to store physical (hardware) properties of a vehicle, such as engines, control surfaces, etc. This
         object is typically created automatically whenever such a hardware model needs to be assigned to a vehicle.


         :type: VehicleSystems
      )doc" )
            .def_property( "rigid_body_properties",
                           &tss::Body::getMassProperties,
                           &tss::Body::setMassProperties,
                           R"doc(

        Object defining the mass, center of mass and inertia tensor of the body. This object is distinct from
        the gravity field of a body (defined by the :attr:`Body.gravity_field` object). A body endowed with this property does *not*
        automatically have a gravity field created for it. However, the whenever a body is endowed with a gravity field,
        a rigid body properties attribute is created to be consistent with this gravity field (e.g. for a spherical harmonic gravity field
        the mass, center of mass and inertia tensor are created from the gravitational parameter, degree-1 coefficients, and degree-2 coefficients plus mean moment of inertia, respectively).

        :type: RigidBodyProperties
     )doc" )
            .def_property( "radiation_pressure_target_models",
                           &tss::Body::getRadiationPressureTargetModels,
                           &tss::Body::setRadiationPressureTargetModels,
                           R"doc(

        List of radiation pressure target models that exist in the body. These objects define how incoming radiation interacts with the body to produce
        a force/torque. A single body may be endowed with multiple target models, which may be selected
        for an acceleration depending on the application. For instance, a body may have a cannonball target model and a panelled target model available,
        and use one for solar radiation pressure acceleration, and the other for planetary radiation pressure acceleration (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure`).
        This attribute is a list of :class:`~RadiationPressureTargetModel`, or a derived class thereof.


        :type: list[RadiationPressureTargetModel]

     )doc" )
            .def_property( "radiation_pressure_source_model",
                           &tss::Body::getRadiationSourceModel,
                           &tss::Body::setRadiationSourceModel,
                           R"doc(

        Object that defines the radiation that a body emits, primarily for the calculation of radiation pressure acceleration.
        It computes the irradiance at a given target location.

        This attribute is a list of :class:`~RadiationSourceModel`, or a derived class thereof.


        :type: RadiationSourceModel

     )doc" )
            .def_property_readonly( "gravitational_parameter", &tss::Body::getGravitationalParameter, R"doc(
         **read-only**

         Attribute of convenience, equivalent to ``.gravity_field_model.gravitational_parameter``


         :type: float
      )doc" )
            .def( "get_ground_station",
                  &tss::Body::getGroundStation,
                  py::arg( "station_name" ),
                  R"doc(

         This function extracts a ground station object from the body.

         This function extracts a ground station object, for a station of a given name, from the body.
         If no station of this name exists, an exception is thrown


         Parameters
         ----------
         station_name : str
             Name of the ground station that is to be retrieved.

         Returns
         -------
         GroundStation
             Ground station object of the station of requested name





     )doc" )
            .def_property_readonly( "ground_station_list",
                                    &tss::Body::getGroundStationMap,
                                    R"doc(

         Dictionary of all ground stations that exist in the body, with dictionary key being the name of the station,
         and the ground station object the key of the dictionary.


         :type: dict[str,GroundStation]
      )doc" );

    py::class_< tss::SystemOfBodies, std::shared_ptr< tss::SystemOfBodies > >( m, "SystemOfBodies", R"doc(

         Object that contains a set of Body objects and associated frame
         information.


         Object that contains a set of Body objects and associated frame
         information. This object stored the entire environment for a
         typical Tudat numerical simulation, and is fundamental for the
         overall Tudat architecture.





      )doc" )
            .def( "get",
                  &tss::SystemOfBodies::getBody,
                  py::arg( "body_name" ),
                  R"doc(

         This function extracts a single Body object from the SystemOfBodies.


         Parameters
         ----------
         body_name : str
             Name of the Body that is to be retrieved.

         Returns
         -------
         Body
             Body object of the requested name





     )doc" )
            .def( "get_body",
                  &tss::SystemOfBodies::getBody,
                  py::arg( "body_name" ),
                  R"doc(

         Deprecated version of :py:func:`~get`





     )doc" )
            .def( "create_empty_body",
                  &tss::SystemOfBodies::createEmptyBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = 1,
                  R"doc(

         This function creates a new empty body.

         This function creates a new empty body, and adds it to the
         :py:class:`~SystemOfBodies`. Since the body is empty, it will
         not have any environment models defined. These must all be
         added manually by a user.


         Parameters
         ----------
         body_name : string
             Name of the Body that is to be added

         process_body : bool, default=True
             Variable that defines whether this new Body will have its
             global frame origin/orientation set to conform to rest of
             the environment.

             .. warning:: Only in very rare cases should
                          this variable be anything other than ``True``.
                          Users are recommended to keep this default value
                          intact.





         Examples
         --------

         This function is often used early on in the environment
         creation segment of a simulation, following the creation of
         a :py:class:`~SystemOfBodies` from the default settings
         for celestial bodies.

         .. code-block:: python
            :emphasize-lines: 18

            # Define string names for bodies to be created from default.
            bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

            # Use "Earth"/"J2000" as global frame origin and orientation.
            global_frame_origin = "Earth"
            global_frame_orientation = "J2000"

            # Create default body settings, usually from `spice`.
            body_settings = environment_setup.get_default_body_settings(
                bodies_to_create,
                global_frame_origin,
                global_frame_orientation)

            # Create system of selected celestial bodies
            bodies = environment_setup.create_system_of_bodies(body_settings)

            # Create vehicle objects.
            bodies.create_empty_body("Delfi-C3")

     )doc" )
            .def( "does_body_exist",
                  &tss::SystemOfBodies::doesBodyExist,
                  py::arg( "body_name" ),
                  R"doc(

         Function to check if a body with a given name exists in the SystemOfBodies

         Parameters
         ----------
         body_name : string
             Name of the Body whose existence is to be checked

         Returns
         -------
         bool
             True if the body exists in this object, false if not




     )doc" )
            .def( "list_of_bodies",
                  &tss::SystemOfBodies::getListOfBodies,
                  R"doc(

         List of names of bodies that are stored in this SystemOfBodies
     )doc" )
            //            .def("get_body_dict",
            //            &tss::SystemOfBodies::getMap,
            //                 get_docstring("SystemOfBodies.get_body_dict").c_str())
            .def( "add_body",
                  &tss::SystemOfBodies::addBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_to_add" ),
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = 1,
                  R"doc(

         This function adds an existing body, which the user has
         separately created, to the :py:class:`~SystemOfBodies`.



         Parameters
         ----------
         body_to_add : Body
             Body object that is to be added.

         body_name : numpy.ndarray
             Name of the Body that is to be added.

         process_body : bool, default=True
             Variable that defines whether this new Body will have its
             global frame origin/orientation set to conform to rest of
             the environment.

             .. warning:: Only in very rare cases should this variable be
                          anything other than ``True``. Users are
                          recommended to keep this default value intact.





     )doc" )
            .def( "remove_body",
                  &tss::SystemOfBodies::deleteBody,
                  py::arg( "body_name" ),
                  R"doc(

         This function removes an existing body from the
         :py:class:`~SystemOfBodies`.



         .. warning:: This function does *not* necessarily delete the
                      Body object, it only removes it from this object.
                      If any existing models in the simulation refer to
                      this Body, it will persist in memory.


         Parameters
         ----------
         body_name : numpy.ndarray
             Name of the Body that is to be removed.





     )doc" )
            .def( "global_frame_orientation",
                  &tss::SystemOfBodies::getFrameOrientation,
                  R"doc(

         Common global frame orientation for all bodies in this SystemOfBodies, described in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#frame-orientation>`__.

     )doc" )
            .def( "global_frame_origin",
                  &tss::SystemOfBodies::getFrameOrigin,
                  R"doc(

         Common global frame origin for all bodies in this SystemOfBodies, described in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#global-origin>`__.

     )doc" );

    //            .def_property_readonly("number_of_bodies",
    //            &tss::SystemOfBodies::getNumberOfBodies,
    //                                   get_docstring("number_of_bodies").c_str()
    //                                   );

    /*!
     **************   SUPPORTING FUNCTIONS USED ENVIRONMENT MODELS
     *******************
     */
}

}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy
