The annual observation runs stopped before final CSV output at the FCULa domain
guard. Source inspection then found two normal-point records whose printed
round-trip times cannot represent lunar ranging. This is a source flight-time
problem established without a computed residual or change to that guard.

| Monthly source line | UTC on 2017-02-20 | Wavelength | Printed RTT | Half-range |
| --- | --- | --- | --- | --- |
| `apollo15_201702.npt:168` | 13:23:08.387557778100 | 1064.2 nm | 24.976 ns | 3.743808215504 m |
| `apollo15_201702.npt:211` | 13:23:20.476420038600 | 532.1 nm | 17.271 ns | 2.588857771059 m |

The [official daily product](https://edc.dgfi.tum.de/pub/slr/data/npt_crd/apollo15/2017/apollo15_20170220.npt)
contains byte-identical NP lines at daily lines 12 and 26. Both monthly sessions
label the target Apollo15 and contain CRD record 11, range type 2, transmit
epoch event 2, station-system-delay-applied flag 1, and two returns. These normal
metadata checks therefore did not establish that their flight times were lunar.
Their short ranges suggest local or calibration contamination, but the source
does not establish the cause. No replacement flight time can be inferred.

An exhaustive printed-field audit of all 1,304 CRD record-11 rows in the 60
annual products, including the out-of-scope Matera records, finds only these
two outside a broad 1–4 s lunar flight-time screen. Every other Grasse RTT is
2.338143316423–2.687189638736 s. All 311 APOLLO Group-e RTTs are
2.3492204584470–2.6770120814302 s. This audit does not establish that all remaining
source records are accurate; it identifies the two physically impossible paths.

The [evidence JSON](nonlunar_flight_time_audit.json) preserves exact source fields,
daily/monthly hashes and the distinct selection proposal. Original provider
bytes, annual inventory, and frozen raw manifest counts remain unchanged:
311 APOLLO + 1,271 Grasse source rows, with the prior single requested Grasse
exclusion. Separately classifying these two rows as invalid lunar flight times
would leave 311 + 1,268 = 1,579 solver-eligible rows. This proposal was subsequently approved as an exact source-validity exclusion
and is now reflected in the annual manifest and final eligibility markers;
it is separate from the user-requested March 9 exclusion and any residual
selection. The evidence JSON preserves the original pre-approval proposal. The repeated August source records
remain distinct source-line entries. No new observation solve was performed.
