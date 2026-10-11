Planetary barycenter propagation and estimation
==============================================

The dedicated barycenter registry implements prescribed internal ephemerides
with an estimated common translational state. Existing propagation, physical
force partials and hybrid-arc solvers are reused through small adapters.

Scope and invariants
--------------------

* A barycenter is not a Body and cannot own physical environment models or be
  selected as an acceleration source. SystemOfBodies keeps a separate registry.
* Members are physical bodies with prescribed internal ephemerides. Fixed,
  normalized, positive weights define the center; default weights are the member
  gravitational parameters sampled at registration. Estimating a member GM must
  not silently change the coordinate definition.
* Registration snapshots nominal member ephemerides in the common inertial frame.
  The nominal weighted center and each centered member offset are evaluated from
  these snapshots, never from the subsequently reconstructed member ephemerides.
  Thus x_i = x_B + d_i(t) and sum(w_i*d_i) = 0 without a dependency cycle.
* For the first implementation, the environment origin and propagated barycenter
  origin are SSB and the barycenter uses Cartesian Cowell propagation.
  Barycenter dynamics use selected external point-mass gravity
  evaluated at the barycenter. Unsupported configurations must fail explicitly.
* Overlapping memberships, duplicate names/members, invalid weights, and
  independently propagated prescribed members must be rejected.
* Scalar and epoch precision must be preserved through the ephemeris interface;
  conversions cannot silently replace an integrated long-double state by double.

State access and lifecycle
--------------------------

Barycenter stores the state, epoch/validity of the current integration stage,
ephemeris, and nominal ephemeris snapshots. SystemOfBodies retains the configured
external force map separately, preventing ownership cycles when models capture
their barycenter. All ownership uses shared pointers. SystemOfBodies adds
explicit translational-state operations (existence,
ephemeris get/set, current-state set, ephemeris-state evaluation). Existing body
lookup and body enumeration remain physical-body-only.

Member ephemerides compose the prescribed offset with the barycenter trajectory.
During environment updates they use a current integrated barycenter state when
available. Otherwise they use the barycenter ephemeris. State validity is reset
at every stage, including repeated epochs; leaving propagation invalidates the
stage state. Initial-state and observation queries use the ephemeris explicitly.

CentralBodyData resolves a prescribed member origin to its parent barycenter
plus an offset. If the parent is in the integrated vector, origin reconstruction
uses that vector directly, before environment updates. This avoids consulting an
old ephemeris or the previous Runge--Kutta stage. Dependency ordering uses the
parent while the user-facing origin remains Jupiter.

Frame managers and integrated-state processors accept the separate registry via
translational-state access. The barycenter's nominal ephemeris is replaced by the
appropriate integrated ephemeris when propagation results are installed. The
normal hybrid-arc solver retains one continuous barycenter trajectory, including
between spacecraft arcs.

Relative dynamics
-----------------

For Juno relative to Jupiter, r_J = r_B + rho_J(t) and
r_Juno = r_B + rho_J(t) + q. Consequently

    q'' = a_Juno - a_B - rho_J''(t).

The acceleration factory creates direct physical accelerations for a spacecraft
whose origin is a prescribed barycenter member. A distinct origin-acceleration
model supplies -a_B-rho_J''. This replaces the usual per-source third-body
subtractions; mixing both would subtract the origin acceleration twice.

The external barycenter force list must be configured before constructing its
derived-origin dynamics (the factory performs a barycenter pass before physical
bodies in a combined request). A missing force definition is a setup error.
The origin model also works when only the spacecraft is currently integrated,
as in the original multi-arc component of a hybrid solver.

The prescribed offset acceleration is obtained from the centered ephemeris
velocity using a documented finite-difference interval. Its state derivative is
zero because the internal ephemerides and weights are fixed. Its numerical
accuracy is checked separately from integrated-state sensitivities.

Variational equations and observations
--------------------------------------

SystemOfBodies supplies a single member-to-parent dependency map. The
variational assembler requests each physical member's partial when assembling
the parent's column, as well as partials directly with respect to that parent.
Central-origin scaling resolves a member origin to the parent and includes
position and velocity dependencies while leaving the kinematic q'=v rows intact.

The new origin acceleration has analytic state and parameter partials formed by
negating the configured barycenter acceleration partials. Existing physical
acceleration partial implementations remain authoritative for their own forces.

The state-transition interface uses the same resolved central dependency when
reconstructing absolute spacecraft sensitivities. Direct observation partials of
a member resolve to its parent's state. Thus Jupiter, moons, and Juno can be
observed without falsely requiring independently estimated member states.

A local Newtonian Jupiter or moon acceleration depends only on fixed internal
geometry and q. Its barycenter derivative must cancel between the accelerated
and accelerating absolute-state contributions. External tides and inertial
tracking geometry retain their barycenter sensitivity.

Public interface
----------------

After creating the physical environment and its nominal ephemerides::

    center = bodies.add_barycenter(
        "Jupiter_barycenter", ["Jupiter", "Io", "Europa", "Ganymede", "Callisto"])

The optional ``weights`` argument overrides the default member GMs; registration
normalizes and freezes these weights. ``center.members``, ``center.weights``,
``center.ephemeris`` and ``center.get_member_offset(member, time)`` expose the
definition and trajectory. ``bodies.get_barycenter(name)`` retrieves the point;
ordinary body lookup and enumeration do not include it.

Use ``Jupiter_barycenter`` in the existing translational propagation and initial
state estimation settings. Configure its external point-mass accelerations
before constructing Juno's accelerations if using separate factory calls. Keep
``Jupiter`` as Juno's central body and as its physical gravity source. For hybrid
arcs, put the barycenter in the single-arc settings and Juno in the independent
multi-arc settings. The usual hybrid initial-state parameter factory then creates
one barycenter state and one spacecraft state per encounter.

The optional ``offset_acceleration_step`` is the interval in seconds for a
five-point derivative of prescribed offset velocity (default 10 s). Choose it
consistent with the smoothness and time resolution of the nominal ephemerides.
Prescribed members cannot also be propagated independently. Physical models
continue to belong to their member bodies; a barycenter cannot exert a force.

Tests and acceptance
--------------------

The focused C++ barycenter-estimation test target uses a controlled Jovian
system with Jupiter and all four Galilean moons, realistic masses and orbital
scales, solar-system external gravity, and a massless Juno-like spacecraft with
close Jupiter passages. Synthetic ephemerides keep the test reproducible without
downloading spacecraft or satellite kernels.

1. Registry/centering tests: normalized weights, correct nominal state, membership
   validation, physical-body-only enumeration, no acceleration from a barycenter,
   and weighted-center identity before/after propagation and same-epoch updates.
2. Single arc: jointly estimate six barycenter and six Jupiter-relative Juno
   initial-state components from ideal observations around a close passage.
   Validate propagated and observation partials against finite differences of
   complete repropagations, including individual member observation mappings.
3. Hybrid arc: one continuous barycenter state over three flybys at 0, 90, and
   180 Julian days (the deterministic three-month convention for this test).
   Each flyby has an independent six-component Juno initial state. Verify 24
   initial-state parameters, recovery from injected errors, shared barycenter
   sensitivity in all three arcs, and zero dependence on other Juno arcs.
4. Check the local-gravity common-translation cancellation and agreement between
   relative propagation and reconstructed absolute-state equations.
5. Run relevant existing translational, variational, single-/multi-/hybrid-arc
   estimation regressions, and compile/smoke-test the Python exposure.

Ideal Cartesian spacecraft observations may be used for a well-conditioned
end-to-end state-recovery test; Earth--Juno tracking observables must additionally
exercise the light-time/observation partial path. Numerical thresholds must be
justified by integration/interpolation and finite-difference convergence, not
relaxed simply to hide a failed derivative or rank-deficient estimation.
