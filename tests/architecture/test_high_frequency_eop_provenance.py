"""Check the built-in source bytes without importing the compiled kernel."""

import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]


class HighFrequencyEopProvenance(unittest.TestCase):
    def test_embedded_table_matches_pinned_primary_source(self):
        manifest = json.loads(
            (
                ROOT / "tests/test_tudat/src/astro/earth_orientation/"
                "desai_sibois_reference/sources.json"
            ).read_text()
        )
        source = manifest["files"]["desai_model_jgrb51665-sup-0002-ds01.txt"]
        embedded = (ROOT / "src/tudat/io/desaiSibois2016EopData.inc").read_bytes()
        table = embedded.split(b'R"DS_EOP(', 1)[1].split(b')DS_EOP"', 1)[0]
        self.assertEqual(len(table), source["bytes"])
        self.assertEqual(hashlib.sha256(table).hexdigest(), source["sha256"])
        self.assertEqual(
            sum(bool(line.strip()) and not line.startswith(b"#") for line in table.splitlines()),
            159,
        )


if __name__ == "__main__":
    unittest.main()
