# Station-loading source fixtures

`imls_merra2_extract.eph` is a reproducible extract of the original IMLS file
`atm_merra2_20200101_0000.eph`, retrieved 14 October 2026 from the URL recorded
in `extract_manifest.json`. The original has 724403 bytes, SHA256
`0a9a2e507fbc6eb66901c2653fbb16ef06fe1acaaafa6c73466168c1c82e8ca8`.
The extract retains provider labels `00NA`, `weuc`, `yela`, their unchanged
Cartesian S rows and radial/east/north D rows, the unchanged A/T records,
signatures and the first 170 source lines of descriptive comments. Only the
P site's/displacement counts change from 4341 to 3. The manifest records that
exact change and a SHA256/source line for every unchanged extract line.

The full original provider file violates the specification's uniqueness rule:
it defines eight site labels twice (AGGO, AIRA, APEX, Arel, HOFN, IITK, NOTO,
Sisl). The strict reader rejects it at the first duplicated AGGO definition,
source line 563. `imls_duplicate_identifier.eph` retains the unchanged AGGO
site/sample rows as a reproducible negative fixture; only declared site/sample
counts change to two. The manifest records all duplicate-label source lines
and every unchanged negative-fixture row. The positive extract selects three
unambiguous labels; no original label or source value is relabeled or repaired.

Regenerate offline:

```sh
python extract.py /path/to/atm_merra2_20200101_0000.eph
```

Sidecars identify the original 2005.06.30 EPHEDISP specification and the IMLS
atmosphere/tips pages, including their retrieval timestamps and hashes. The
provider's actual file has A before the three T records, whereas the format
specification lists T before A. The reader permits those two complete layouts,
and rejects interleaving or duplicate sections. The file also demonstrates a
valid singleton with sample interval zero and matching begin/end epochs.

The native tests check fixed-field values against the unchanged source rows.
Additional synthetic records independently exercise partial/empty per-site
coverage, arbitrary internal identifiers, ignored informational fields,
fractional source epochs, both section layouts, CR/LF record delimiters, counts,
sorting, duplicates, invalid intervals, nonfinite fields and internal gaps.
Synthetic epoch/index modifications are explicitly separate from the real
provider extract; they do not claim additional observed provider data.

The named MERRA2 atmospheric profile is explicit, never selected by a filename.
The pinned description removes mean and eight diurnal/semidiurnal harmonics;
it retains annual, semiannual and linear drift. The source's Love-number origin
is whole-Earth CM including atmosphere. IMLS's loose CE/CF terminology is not
used to equate those physical origins. Product release remains caller supplied.

The standalone CMake project builds the reader with only the standard library;
Boost belongs only to the native test executable.
