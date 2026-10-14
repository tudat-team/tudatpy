"""Verify pinned extract hashes, source ranges and edits without a kernel."""

import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests/test_tudat/src/io/station_loading_fixtures"


class StationLoadingProvenance(unittest.TestCase):
    def check_ranges_and_edits(self, table, ranges, edits):
        lines = table.splitlines(keepends=True)
        source_lines = [number for first, last in ranges for number in range(first, last + 1)]
        self.assertEqual(source_lines, sorted(set(source_lines)))
        self.assertEqual(len(lines), len(source_lines) + len(edits))
        for edit in edits:
            self.assertNotIn(edit["source_line"], source_lines)
            self.assertEqual(lines[edit["extract_line"] - 1].decode(), edit["replacement"])
            self.assertEqual(edit["original"], "P T 3 S       4341 E      1 D       4341\n")
            self.assertTrue(edit["replacement"].startswith("P T 3 S "))

    def test_exact_extract_hashes_source_ranges_and_count_edits(self):
        manifest = json.loads((FIXTURES / "extract_manifest.json").read_text())
        table = (FIXTURES / "imls_merra2_extract.eph").read_bytes()
        self.assertEqual(hashlib.sha256(table).hexdigest(), manifest["extract_sha256"])
        self.assertEqual(len(table), manifest["extract_bytes"])
        self.check_ranges_and_edits(table, manifest["unchanged_source_ranges"], manifest["edits"])
        self.assertEqual(manifest["unchanged_source_ranges"][0], [1, 170])
        self.assertEqual(manifest["source_bytes"], 724403)
        self.assertEqual(
            manifest["source_sha256"],
            "0a9a2e507fbc6eb66901c2653fbb16ef06fe1acaaafa6c73466168c1c82e8ca8",
        )
        negative_manifest = manifest["negative_duplicate_fixture"]
        negative = (FIXTURES / "imls_duplicate_identifier.eph").read_bytes()
        self.assertEqual(hashlib.sha256(negative).hexdigest(), negative_manifest["sha256"])
        self.assertEqual(len(negative), negative_manifest["bytes"])
        self.check_ranges_and_edits(
            negative, negative_manifest["unchanged_source_ranges"], negative_manifest["edits"]
        )
        self.assertIn([562, 563], negative_manifest["unchanged_source_ranges"])
        self.assertEqual(manifest["source_duplicate_site_labels"]["AGGO"], [562, 563])


if __name__ == "__main__":
    unittest.main()
