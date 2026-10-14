# IERS ocean pole-tide regression data

Sources and SHA256 hashes are pinned in `sources.json`. The published test is
unchanged. The grid extract retains four original coefficient rows at longitudes
232.25/232.75 and latitudes -43.75/-43.25 degrees; only the axis counts and extent
headers are adjusted for the extract. The full grid remains external.

The published series explicitly uses the IERS 2003 mean pole and the older
G=6.673e-11, K=5339.4043696 metres/radian. Tests use those values solely to
reproduce its published reference. Production retains the 2018 secular pole and
IERS 2010 chapter 1 constants (K=5340.428562229677 metres/radian).

The 731 published rows print both wobble and displacement with seven significant
digits. Regression tolerance is 1 nm; maximum residual using the printed wobble
is 0.7894 nm (including the rounded displacement).

To regenerate or verify from a local directory containing the two pinned downloads:

```sh
python3 extract_reference_grid.py /path/to/reference-directory
python3 extract_reference_grid.py /path/to/reference-directory --check
```

The script checks both SHA256 hashes before extracting original rows; it performs
no network requests and evaluates no production formula.

The geocentric coordinate interpretation is supported by Desai (2002),
*Observing the pole tide with satellite altimetry*, DOI10.1029/2001JC001224,
Eqs.(3)-(6): Earth-centred Cartesian coordinates are rewritten as spherical
colatitude/longitude, and Eq.(14) uses the same angular coordinates for loading.
Thus the model uses geocentric latitude (pi/2-colatitude) and a spherical
radial/north/east basis. The raw provider header remains neutral.
Primary paper mirror: https://igsac-cnes.cls.fr/documents/gins/Doc_Itrf/Desai2002_poletide.pdf
Paper SHA256: 6f72013f5abc12accf926ed5442d78cde4a85993ed8bc22684b52c92ef6724c9

The paper's Eq.(7) distinguishes reported and instantaneous rotation poles using
a pole derivative. This implementation follows the conventional IERS Eq.(25)
wobble from reported xp/yp minus the 2018 secular pole; no derivative correction
from the research paper is applied.
