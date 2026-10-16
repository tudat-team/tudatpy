INPOP17a time kernel and prism convention
=========================================

The published lunar reflector reference point remains **unresolved**.
`Manche (2011) <https://www.imcce.fr/content/medias/publications/publications-recherche/theses-habilitations/docs/Manche_These.pdf>`_
chapter 9, pp.92–94, Eq.9.1 lists geometric, relativistic and atmospheric
light time. Section 9.3 specifies principal-axis reflector coordinates,
but neither defines a front-face centre, vertex or effective optical point.
The `published 2018 article <https://doi.org/10.1093/mnras/sty096>`_ section
2.3 refers to that thesis and Moyer (2003); it does not resolve the prism
reference. The INPOP17a/Table9 and INPOP13c/Table13 coordinate sources
also provide no recovered definition. Thermal refractive-index discussion
concerns return detectability. Missing an explicit glass term is not proof
that optical delay is absorbed in the coordinates. Predeclared ``n*h/c``
per-leg and zero-explicit-glass controls must retain this ambiguity; residual
improvement cannot determine the convention.

The official ``_time.bsp`` is usable as an independent **geocentric** clock
comparison. The segment is named ``TT-TDB``, target 1000000001 relative to
1000000000, SPK type2, frame1. The release README specifies seconds. Its
first pseudo-position component is signed TT−TDB in seconds; no kilometre
conversion applies. It covers ±3234816000 TDB seconds since J2000, with
18720 four-day records, 38 words each, 12 coefficients per component.
The original TDB header/raw unit flags and primary source hashes are retained
in ``summary.json``.

Direct long-double Chebyshev evaluation and SPICE evaluation agree within
4.6e-18 s. INPOP minus SOFA TT−TDB at the three saved launches is respectively
+0.42146 ns (APOLLO March2), +0.32837 ns (Grasse March9) and +1.39241 ns
(APOLLO April21). At a declared 100 m/s radial speed the largest absolute
launch sensitivity is 0.140 micrometres. This compares values, not complete
RTT clock-rate differences. The topocentric station term is separate;
loading the file does not replace the native converter. No native module,
observation run, build or model change was used for this coefficient audit.

``audit_time_spk.py`` extracts the SPK directory and each selected coefficient
record, using the standard `NAIF type2 layout <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%202:%20Chebyshev%20(position%20only)>`_.
``summary.json`` preserves the three actual X-coefficient vectors as well
as source identities and signed results. Reproduce this static check::

    python analysis/laser_ranging/results/inpop17a_time_and_prism_audit/audit_time_spk.py \
      --time-spk /path/to/inpop17a_TDB_m100_p100_spice_time.bsp \
      --output /new/path/time-spk-audit.json

The default launch evidence is the adjacent committed absolute-clock audit.
The portable script only changes path handling/formatting and adds recorded
coefficients; packaging performed no additional numerical runs.
