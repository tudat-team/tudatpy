.. _observation_corrections:

``observation_corrections``
===========================

This module contains functionality to compute and apply several observation corrections.
These corrections are computed once and applied on an existing set of observations before
handing them to the estimation. As such, a reasonably accurate a priori ephemeris of the
object should be known to accurately compute the corrections. These can be retrieved from
e.g. JPL Horizons.


.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations/observation_corrections/photocenter_correction
   /estimation/observations/observation_corrections/light_deflection_correction
