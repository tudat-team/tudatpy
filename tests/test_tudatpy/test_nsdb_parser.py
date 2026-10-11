"""
Unit tests for the NSDB metadata parser (``_parse_contents_section``).

Tests cover:
- Single-line fields
- Multi-line fields (continuation lines)
- Values containing colons
- Blank lines within sections
- Section boundaries (continuation stops at next section header)
- Duplicate key detection
- Empty input
- Mixed single and multi-line fields
"""
import importlib.util
import os

import pytest

# Import NSDBQuery directly from source to avoid triggering the
# tudatpy.data package init (which requires a compiled C extension).
_nsdb_path = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "src", "tudatpy", "data", "nsdb", "nsdb_data.py",
)
_spec = importlib.util.spec_from_file_location("nsdb_data", _nsdb_path)
_nsdb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nsdb)
NSDBQuery = _nsdb.NSDBQuery


# ---------------------------------------------------------------------------
# Fixture: an instance of NSDBQuery (the method is static, but we call it
# through the class for clarity).
# ---------------------------------------------------------------------------
@pytest.fixture
def parser():
    return NSDBQuery._parse_contents_section


# ---------------------------------------------------------------------------
# Single-line fields
# ---------------------------------------------------------------------------
class TestSingleLineFields:
    def test_single_key_value(self, parser):
        text = "observatory: 083\n"
        assert parser(text) == {"observatory": "083"}

    def test_multiple_single_line(self, parser):
        text = "observatory: 083\ntype: absolute\n"
        assert parser(text) == {"observatory": "083", "type": "absolute"}

    def test_value_with_spaces(self, parser):
        text = "observatory: 083 - Golosseevo-Kiev\n"
        assert parser(text) == {"observatory": "083 - Golosseevo-Kiev"}

    def test_key_with_dots(self, parser):
        text = "number of observations: 42\n"
        assert parser(text) == {"number of observations": "42"}


# ---------------------------------------------------------------------------
# Multi-line fields
# ---------------------------------------------------------------------------
class TestMultiLineFields:
    def test_two_line_value(self, parser):
        text = (
            "observatory: 083 - Golosseevo-Kiev\n"
            "             188 - Majdanak\n" \
            "             095 - Simeiz\n"
        )
        assert parser(text) == {
            "observatory": ["083 - Golosseevo-Kiev", "188 - Majdanak", "095 - Simeiz"]
        }

    def test_three_line_value(self, parser):
        text = (
            "observatory: 083 - Golosseevo-Kiev\n"
            "             188 - Majdanak\n"
            "             095 - Simeiz\n"
        )
        assert parser(text) == {
            "observatory": [
                "083 - Golosseevo-Kiev",
                "188 - Majdanak",
                "095 - Simeiz",
            ]
        }

    def test_mixed_single_and_multi(self, parser):
        text = (
            "observatory: 083 - Golosseevo-Kiev\n"
            "             188 - Majdanak\n"
            "type: absolute\n"
        )
        assert parser(text) == {
            "observatory": ["083 - Golosseevo-Kiev", "188 - Majdanak"],
            "type": "absolute",
        }

    def test_continuation_with_leading_whitespace(self, parser):
        """Continuation lines may have varying indentation."""
        text = (
            "satellites: Triton\n"
            "  Nereid\n"
            "    Proteus\n"
        )
        assert parser(text) == {
            "satellites": ["Triton", "Nereid", "Proteus"]
        }


# ---------------------------------------------------------------------------
# Values containing colons
# ---------------------------------------------------------------------------
class TestValuesWithColons:
    def test_colon_in_value_single_line(self, parser):
        text = "reference: Smith et al. (2020): A study\n"
        assert parser(text) == {"reference": "Smith et al. (2020): A study"}

    def test_colon_in_value_multiline(self, parser):
        text = (
            "reference: Smith et al. (2020): A study\n"
            "          Jones (2019): Another work\n"
        )
        assert parser(text) == {
            "reference": [
                "Smith et al. (2020): A study",
                "Jones (2019): Another work",
            ]
        }

    def test_url_with_colons(self, parser):
        text = "source: https://nsdb.imcce.fr/obspos\n"
        assert parser(text) == {"source": "https://nsdb.imcce.fr/obspos"}


# ---------------------------------------------------------------------------
# Blank lines
# ---------------------------------------------------------------------------
class TestBlankLines:
    def test_blank_lines_between_fields(self, parser):
        text = (
            "observatory: 083\n"
            "\n"
            "type: absolute\n"
        )
        assert parser(text) == {"observatory": "083", "type": "absolute"}

    def test_blank_lines_within_multiline(self, parser):
        """Blank lines inside a multi-line block are skipped."""
        text = (
            "observatory: 083 - Golosseevo-Kiev\n"
            "\n"
            "             188 - Majdanak\n"
        )
        assert parser(text) == {
            "observatory": ["083 - Golosseevo-Kiev", "188 - Majdanak"]
        }

    def test_leading_blank_lines(self, parser):
        text = "\n\nobservatory: 083\n"
        assert parser(text) == {"observatory": "083"}

    def test_trailing_blank_lines(self, parser):
        text = "observatory: 083\n\n\n"
        assert parser(text) == {"observatory": "083"}


# ---------------------------------------------------------------------------
# Section boundaries
# ---------------------------------------------------------------------------
class TestSectionBoundaries:
    def test_stops_at_section_header(self, parser):
        """Continuation lines must not cross a section header."""
        text = (
            "observatory: 083 - Golosseevo-Kiev\n"
            "Reference.\n"
            "type: absolute\n"
        )
        # "Reference." is a section header, so "083 - Golosseevo-Kiev" is
        # the only value for "observatory".
        assert parser(text) == {"observatory": "083 - Golosseevo-Kiev"}

    def test_multiple_sections(self, parser):
        text = (
            "observatory: 083\n"
            "Reference.\n"
            "author: Smith\n"
            "Format.\n"
            "1. column one\n"
        )
        assert parser(text) == {"observatory": "083"}


# ---------------------------------------------------------------------------
# Duplicate keys
# ---------------------------------------------------------------------------
class TestDuplicateKeys:
    def test_duplicate_key_raises(self, parser):
        text = (
            "observatory: 083\n"
            "observatory: 188\n"
        )
        with pytest.raises(ValueError, match="Duplicate key.*observatory"):
            parser(text)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------
class TestEdgeCases:
    def test_empty_text(self, parser):
        assert parser("") == {}

    def test_only_whitespace(self, parser):
        assert parser("   \n  \n  ") == {}

    def test_only_section_header(self, parser):
        """A section header alone should not produce any entries."""
        text = "Contents.\n"
        assert parser(text) == {}

    def test_key_without_value(self, parser):
        """A key followed by colon and nothing else."""
        text = "observatory:\n"
        result = parser(text)
        assert result == {"observatory": ""}

    def test_key_without_value_then_continuation(self, parser):
        """A key with empty first line followed by continuation lines."""
        text = (
            "observatory:\n"
            "  083 - Golosseevo-Kiev\n"
            "  188 - Majdanak\n"
        )
        assert parser(text) == {
            "observatory": ["", "083 - Golosseevo-Kiev", "188 - Majdanak"]
        }