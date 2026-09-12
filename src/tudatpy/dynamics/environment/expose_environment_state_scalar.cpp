/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_environment_state_scalar.h"

#include "scalarTypes.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding )
{
    bodyBinding.def( "state_in_base_frame_from_ephemeris",
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





     )doc" );

    systemOfBodiesBinding
            .def( "create_empty_body",
                  &tss::SystemOfBodies::createEmptyBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = true,
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
            .def( "add_body",
                  &tss::SystemOfBodies::addBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_to_add" ),
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = true,
                  R"doc(

         This function adds an existing body, which the user has
         separately created, to the :py:class:`~SystemOfBodies`.



         Parameters
         ----------
         body_to_add : Body
             Body object that is to be added.

         body_name : str
             Name of the Body that is to be added.

         process_body : bool, default=True
             Variable that defines whether this new Body will have its
             global frame origin/orientation set to conform to rest of
             the environment.

             .. warning:: Only in very rare cases should this variable be
                          anything other than ``True``. Users are
                          recommended to keep this default value intact.





     )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy
