"""Verify real-provider extract provenance without loading the compiled kernel."""

import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests/test_tudat/src/io/station_loading_fixtures"


class StationLoadingProvenance(unittest.TestCase):
    def test_every_unchanged_extract_record_matches_its_source_digest(self):
        manifest = json.loads((FIXTURES / "extract_manifest.json").read_text())
        table = (FIXTURES / "imls_merra2_extract.eph").read_bytes()
        self.assertEqual(hashlib.sha256(table).hexdigest(), manifest["extract_sha256"])
        self.assertEqual(len(table), manifest["extract_bytes"])
        lines = table.splitlines(keepends=True)
        for record in manifest["unchanged_records"]:
            with self.subTest(source_line=record["source_line"]):
                line = lines[record["extract_line"] - 1]
                self.assertEqual(hashlib.sha256(line).hexdigest(), record["sha256"])
        changes = manifest["changed_counts"]
        self.assertEqual(
            [line.decode().strip() for line in lines if line.startswith(b"P ")],
            [changes["extract"]],
        )
        self.assertEqual(changes["original"], "P T 3 S       4341 E      1 D       4341")
        self.assertEqual(
            manifest["source_sha256"],
            "0a9a2e507fbc6eb66901c2653fbb16ef06fe1acaaafa6c73466168c1c82e8ca8",
        )
        negative_manifest = manifest["negative_duplicate_fixture"]
        negative = (FIXTURES / "imls_duplicate_identifier.eph").read_bytes()
        self.assertEqual(hashlib.sha256(negative).hexdigest(), negative_manifest["sha256"])
        negative_lines = negative.splitlines(keepends=True)
        for record in negative_manifest["unchanged_records"]:
            self.assertEqual(
                hashlib.sha256(negative_lines[record["extract_line"] - 1]).hexdigest(),
                record["sha256"],
            )
        self.assertEqual(manifest["source_duplicate_site_labels"]["AGGO"], [562, 563])


if __name__ == "__main__":
    unittest.main()
