# PR 913 split: decisions and open issues

The architecture branch is based on `develop`. The gravity-deformation branch is
stacked on the architecture branch. These branches are local until the split is
reviewed.

1. **Compatibility boundary — decided.** Preserve the existing Python interface,
   including the legacy scaled-mean-moment setting and deprecation path. Direct
   C++ construction and interfaces may change where the new ownership model
   requires it. The architecture branch owns this migration.
2. **Rates of existing gravity-field variations — required.** Add coefficient-rate
   support for the existing periodic, polynomial, tabulated, tidal, and custom
   variation models in the gravity-deformation branch. Their contributions to
   the inertia-tensor derivative must be combined with the integrated variation
   rate. The precise contract for custom and tabulated data still needs a choice
   between explicit derivative input and a numerical fallback.
3. **Maxwell variation input — open.** The current implementation subtracts
   nominal coefficients from the total field. With other variations present,
   this includes their corrections. Decide whether the Maxwell constitutive
   equation should evolve only its integrated variation or the combined
   correction, and test that decision.
4. **Coupled derivatives — open.** The same-epoch loop between gravity-coefficient
   rates, inertia derivatives, and angular acceleration still needs a focused
   solution. The resulting state derivatives and variational equations must use
   the same coupling.

The architecture branch contains gravity-derived rigid-body settings and
properties, their gravity-field link and synchronization, inertia availability,
Python compatibility, and related estimation and torque-partial changes. The
gravity-deformation branch contains integrated variations, the Maxwell model,
gravity-deformation propagation, coefficient-rate implementations, and their
tests and Python exposure.
