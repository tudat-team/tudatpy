# MRO TNF short-arc fixture

`mro_short_arc_2012_077.tnf` contains complete SFDU records extracted from the
NASA PDS file `mromagr2012_076_0840xmmmv1.tnf`. It retains the 759 earlier ramp
records needed to evaluate the uplink frequency at light-time-shifted transmit
epochs, followed by source-file bytes 62,155,862 through 67,155,899, inclusive.
The bounded Doppler observations span 2012-03-17 11:43:25.5 through 12:50:42.5
UTC.

Source:
https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/tnf/mromagr2012_076_0840xmmmv1.tnf

The fixture is 5,109,334 bytes. Its SHA-256 digest is
`a48ddb1a9c7e1ab44213b9f34021e8ed330f1ef7cd49630da5b8f3449efd7720`.
The bounded file keeps the calibrated link-delay and residual checks unchanged
while avoiding the full 78,995,114-byte daily TNF file and thousands of repeated
single-observation deletions during the test.
