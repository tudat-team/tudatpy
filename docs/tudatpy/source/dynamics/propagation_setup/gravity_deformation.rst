Gravity-coefficient propagation
===============================

Tudat propagates degree-two gravity-field *variations*, rather than replacing
the nominal gravity field.  For each body, the five-element numerical state is

.. math::

   [\Delta C_{20},\Delta C_{21},\Delta C_{22},\Delta S_{21},\Delta S_{22}].

All five propagated entries are unnormalised.  Spherical-harmonic environment
settings, the coefficient matrices exposed by a body, and the
``static_coefficients`` accepted by
:func:`~tudatpy.dynamics.propagation_setup.propagator.maxwell_deformation` use
geodesy normalisation.  Conversion occurs at the deformation-model and
integrated-variation boundaries.  ``static_coefficients`` identify a static
contribution already present in the nominal environment field; this
contribution is excluded from the deformable baseline and is not propagated.
Consequently, the instantaneous total field is the nominal field plus this
integrated variation plus any other configured gravity-field variations.

The propagation setup automatically creates an internal integrated
gravity-field variation object.  A static spherical-harmonic field is promoted
to a time-dependent field without changing its gravitational parameter,
reference radius, nominal coefficients, body-fixed frame, or scaled mean
moment of inertia.  No direct mutation of body coefficients is required.

During propagation, environment updates read the current five-element slice
directly from the numerical state.  If
:attr:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorProcessingSettings.set_integrated_result`
is enabled, post-processing installs the complete state history and its
configured interpolator in the same variation object and refreshes the cached
total field at the greatest saved epoch.  Outside propagation, later
environment updates evaluate that interpolator; during a subsequent
propagation they switch back to the live state.  This is the same lifecycle
distinction used for propagated synchronous rotations.

Inertia coupling
----------------

The body's rigid-body properties object owns both the current inertia tensor
and its derivative.  A gravity field can exist without providing inertia.  A
spherical-harmonic field provides inertia only when it has complete degree-two
coefficients and a finite scaled mean moment of inertia.  Check
``body.rigid_body_properties.inertia_tensor_available`` and
``inertia_tensor_derivative_available`` before accessing the corresponding
values.  When available, the tensor is refreshed automatically from the
current total degree-two field, including ordinary gravity variations computed
during environment updates.  Its derivative is refreshed from a live
integrated coefficient rate during gravity propagation, or from an installed
coefficient history after that propagation.  Existing instantaneous
gravity-variation models do not yet expose coefficient derivatives; their
contribution to the inertia derivative remains unsupported.  If an inertia
tensor is available without a supported rate source, rotational propagation
retains the legacy zero-inertia-derivative behavior while the public derivative
remains unavailable.

If gravity-coefficient rates and angular acceleration depend on one another,
Tudat solves the resulting algebraic derivative system.  The default is a
scaled direct affine solve with a checked fixed-point fallback.  Configure
tolerances, iteration limits, and non-convergence handling through
``propagator_settings.coupled_state_derivative_solver_settings``.

Supported modes
---------------

The current Maxwell implementation supports degree and order 2, with exactly
the five coefficients listed above.  It may be propagated alone or as one
component of a single-arc multi-type propagation, including coupled rotational
dynamics.  Multiple bodies may be included in one gravity propagation; their
five-element states are concatenated in ``bodies_to_integrate`` order.
Single-arc result interpolation and environment installation are supported.
Multi-arc and hybrid-arc installation of a single continuous gravity history
is not currently part of the public contract.

Python setup
------------

Create deformation settings with
:func:`~tudatpy.dynamics.propagation_setup.propagator.maxwell_deformation`,
turn them into models with
:func:`~tudatpy.dynamics.propagation_setup.create_gravity_deformation_models`,
and create propagation settings with
:func:`~tudatpy.dynamics.propagation_setup.propagator.gravity_deformation`.
