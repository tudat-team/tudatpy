import pytest
import pandas as pd
from tudatpy.data.mpc import BatchMPC
from tudatpy.data.mpc.parser_80col.parsers import (
    parse_80cols_data,
    parse_80cols_file,
    identify_object
)
from tudatpy.data.mpc.parser_80col.unpackers import (
    unpack_permanent_minor_planet,
    unpack_provisional_minor_planet,
    unpack_provisional_comet_or_satellite,
    unpack_permanent_natural_satellite
)

# ==============================================================================
# SECTION A: UNIT TESTS (Logic Only, Offline)
# ==============================================================================
# These tests verify the unpacking logic for various object types in isolation.

# --- 1. Permanent Minor Planet Tests (Asteroids) ---
@pytest.mark.parametrize("packed_input, expected_output", [
    # Case 3: Basic Range (0 - 99,999) - Straight numeric conversion
    ("00001", "1"),
    ("00433", "433"),
    ("99999", "99999"),

    # Case 2: Standard Range (100,000 - 619,999) - Single char prefix + 4 digits
    # A0000 = 10 * 10000 + 0 = 100,000
    ("A0000", "100000"),
    ("A0001", "100001"),
    ("Z9999", "359999"), # Z=35
    ("a0000", "360000"), # a=36
    ("z9999", "619999"), # z=61

    # Case 1: Extended Range (> 620,000) - Tilde prefix + Base62 chars
    ("~0000", "620000"),
    ("~0001", "620001"),
    ("~AZaz", "3140113"), # Complex Base62 check
])
def test_unpack_permanent_minor_planet_success(packed_input, expected_output):
    """Verifies unpacking of permanent minor planet IDs (Asteroid Numbers)."""
    assert unpack_permanent_minor_planet(packed_input) == expected_output

@pytest.mark.parametrize("invalid_input, error_msg", [
    ("123", "Invalid packed length"),       # Too short
    ("123456", "Invalid packed length"),    # Too long
    ("A123x", "Suffix must be numeric"),    # Non-digit suffix in Case 2
    ("?0000", "Invalid start character"),   # Invalid Prefix
])
def test_unpack_permanent_minor_planet_errors(invalid_input, error_msg):
    """Ensures invalid formats raise appropriate ValueError."""
    with pytest.raises(ValueError, match=error_msg):
        unpack_permanent_minor_planet(invalid_input)


# --- 2. Provisional Minor Planet Tests ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("PLS2040", "2040 P-L"),    # Survey: Palomar-Leiden
    ("J98M00M", "1998 MM"),     # Standard: 1998 + M + 0 + M
    ("K07T08A", "2007 TA8"),    # Numeric Cycle: 2007 + T + 8 + A
    ("K07Tf8A", "2007 TA418"),  # Base62 Cycle: 'f8' -> 418
])
def test_unpack_provisional_minor_planet_success(packed_input, expected_output):
    """Verifies unpacking of provisional minor planet designations."""
    assert unpack_provisional_minor_planet(packed_input) == expected_output

@pytest.mark.parametrize("invalid_input", ["J98", "J98M0012", "X99A00A"])
def test_unpack_provisional_minor_planet_errors(invalid_input):
    """Ensures invalid provisional formats raise ValueError."""
    with pytest.raises(ValueError):
        unpack_provisional_minor_planet(invalid_input)


# --- 3. Provisional Comet & Satellite Tests ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("K07C01A", "2007 C1-A"), # Comet: 2007 C1 Fragment A
    ("J95X020", "1995 X2"),   # Satellite: 1995 X2
    ("K21A05b", "2021 A5-B"), # Fragment casing check (b -> B)
])
def test_unpack_provisional_comet_sat_success(packed_input, expected_output):
    """Verifies unpacking of provisional comet and satellite designations."""
    assert unpack_provisional_comet_or_satellite(packed_input) == expected_output


# --- 4. Permanent Natural Satellite Tests ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("J001S", "Jupiter I"),   # J = Jupiter
    ("S010S", "Saturn X"),    # S = Saturn
    ("N002S", "Neptune II"),  # N = Neptune
])
def test_unpack_permanent_natural_satellite_success(packed_input, expected_output):
    """Verifies unpacking of permanent natural satellite IDs (Roman Numerals)."""
    assert unpack_permanent_natural_satellite(packed_input) == expected_output


# ==============================================================================
# SECTION B: PARSER LOGIC TESTS (Data Handling, Offline)
# ==============================================================================
# These tests ensure the parser correctly handles lists of strings, reads files,
# filters out metadata lines, and correctly extracts the object ID.

def test_80cols_line_parser_logic():
    """
    Tests filtering, column extraction, and format handling on raw string lists.
    No file I/O or internet required.
    """
    # 1. Setup Data: Real observation lines
    line_atlas = '0003I         C2025 07 03.94821918 01 44.983-18 39 35.17         18.0 GW#0JSN598'

    # Eros: Observation ('S') line followed by Parallax ('s') line
    line_eros_valid = '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500'
    line_eros_skip = '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'

    line_charon = 'D4341J79M00A*4A1979 06 25.66181 20 27 06.64 -15 37 11.5          19.0   M4986413'
    line_2025FA22 = '     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'

    lines_to_parse = [line_atlas, line_eros_valid, line_eros_skip, line_charon, line_2025FA22]

    # 2. Parse
    parsed_table = parse_80cols_data(lines_to_parse)

    # 3. Assertions
    # A. Check filtering (satellite parallax 's' line should be dropped)
    # Input has 5 lines, but 1 is a parallax ('s') line, so we expect 4 valid obs.
    assert len(parsed_table) == 4, "Failed to filter 's' line (satellite parallax)"

    # B. Check ID extraction (Strings)
    ids = parsed_table['number'].astype(str).tolist()
    assert ids[0] == '3I'          # Comet Unpacking
    assert ids[1] == '433'         # Asteroid Unpacking
    assert ids[2] == '134341'      # 'D4341' -> 134341 Unpacking
    assert ids[3] == '2025 FA22'   # Provisional Unpacking

def test_80cols_malformed_lines():
    """Tests that invalid line lengths or formats raise ValueError."""
    bad_len = '    K25F22A  C2025 10 13.24277' # Too short
    # Valid length (80) but missing separator in RA (cols 34 or 37)
    bad_fmt = '     K25F22A  C2025 10 13.24277 0020 45.76 +25 53 06.1          18.3 RrET147718 '

    with pytest.raises(ValueError, match="Line Length Error"):
        parse_80cols_data([bad_len])

    with pytest.raises(ValueError, match="Invalid Separator"):
        parse_80cols_data([bad_fmt])

def test_parse_80cols_file_io(tmp_path):
    """Tests file reading capability using a temporary file (clean cleanup)."""
    # 1. Create temp file
    p = tmp_path / "test_obs.txt"
    p.write_text(
        '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500\n'
        '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'
    ) # this is a satellite observation from eros

    # 2. Parse file path
    parsed_table = parse_80cols_file(str(p))

    # 3. Assert
    assert len(parsed_table) == 1 # must be one single observation line, since the second line is a satellite parallax line
    assert str(parsed_table['number'][0]) == '433'


# ==============================================================================
# SECTION C: VALIDATION LOGIC TESTS (Detailed Failure Modes)
# ==============================================================================
# These tests verify that the parser's internal validation logic (range checks,
# separator checks) correctly identifies and reports bad data.

def test_validation_separator_error():
    """Ensures parser catches invalid spaces between RA/DEC components."""
    # Modified from a real line to introduce error: RA Sep '00 20' -> '00x20'
    bad_line = '     K25F22A  C2025 10 13.24277 00x20 45.76 +25 53 06.1          18.3 RrET147718'
    with pytest.raises(ValueError, match="Invalid Separator"):
        parse_80cols_data([bad_line])

def test_validation_range_error_seconds():
    """Ensures parser catches invalid seconds (>= 60)."""
    # Modified from a real line: Seconds '45.76' -> '60.00'
    bad_line = '     K25F22A  C2025 10 13.24277 00 20 60.00 +25 53 06.1          18.3 RrET147718'
    with pytest.raises(ValueError, match="RA Seconds '60.0' >= 60"):
        parse_80cols_data([bad_line])

def test_validation_range_error_month():
    """Ensures parser catches invalid months."""
    # Modified from a real line: Month '10' -> '13'
    bad_line = '     K25F22A  C2025 13 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'
    with pytest.raises(ValueError, match="Month '13' out of range"):
        parse_80cols_data([bad_line])


# ==============================================================================
# SECTION D: SATELLITE STRUCTURE TESTS
# ==============================================================================
# The MPC format for natural satellites requires pairs of lines:
# 1. Observation line (Flag 'S' in col 15)
# 2. Parallax/Offset line (Flag 's' in col 15)
# These tests ensure the parser enforces this structure.

def test_satellite_mismatch_counts():
    """Ensures parser catches count mismatch between 'S' and 's' lines."""
    line_S = '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500'
    line_s = '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'

    # Case: 2 Observations ('S'), 1 Parallax ('s')
    with pytest.raises(ValueError, match="Mismatch in satellite lines"):
        parse_80cols_data([line_S, line_S, line_s])

def test_satellite_orphan_S():
    """Ensures parser catches an 'S' line not immediately followed by an 's' line."""
    line_S = '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500'
    line_s = '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'
    line_std = '0003I         C2025 07 03.94821918 01 44.983-18 39 35.17         18.0 GW#0JSN598'

    # We provide valid counts (1 S, 1 s) but broken sequence: S -> Std -> s
    # This triggers "Observation 'S' not followed by Parallax"
    with pytest.raises(ValueError, match="Observation 'S' not followed by Parallax"):
        parse_80cols_data([line_S, line_std, line_s])


# ==============================================================================
# SECTION E: IDENTIFICATION LOGIC TESTS (Object Types)
# ==============================================================================
# Tests the `identify_object` internal function directly to ensure correct
# categorization based on regex patterns (Comet vs Asteroid vs Satellite).

@pytest.mark.parametrize("row_data, expected_type, expected_name", [
    ({"number": "0001I", "provisional_designation": ""}, "Interstellar", "1I"),
    ({"number": "0009P", "provisional_designation": ""}, "Comet", "9P"),
    ({"number": "J013S", "provisional_designation": ""}, "Natural Satellite", "Jupiter XIII"),
    ({"number": "00433", "provisional_designation": ""}, "Minor Planet", "(433)"),
])
def test_identify_object_types(row_data, expected_type, expected_name):
    """Tests the identification logic directly."""
    row = pd.Series(row_data)
    result = identify_object(row)
    assert result['obj_type'] == expected_type
    assert result['unpacked_name'] == expected_name

def test_identify_object_missing_ids():
    """Ensures error is raised if both ID columns are empty."""
    row = pd.Series({"number": "     ", "provisional_designation": "       "})
    # Updated match string to reflect actual error message in parser
    with pytest.raises(ValueError, match="Missing both permanent ID and provisional ID"):
        identify_object(row)


# ==============================================================================
# SECTION F: INTEGRATION TEST (BatchMPC Consistency, Online)
# ==============================================================================
# This test requires an internet connection. It queries the actual MPC via BatchMPC
# and compares the results against our parser's output for the same objects.
# This ensures that our unpacking logic matches Tudat/BatchMPC's expectations.

def test_batch_mpc_vs_parser_consistency():
    """Checks consistency between Parser output and BatchMPC expectations."""
    batch = BatchMPC()

    # We use a mix of object types to stress-test the consistency:
    # - '3I': Interstellar Comet (Short format)
    # - 433: Numbered Asteroid (Integer input)
    # - 134341: High-numbered Asteroid (Pluto/Charon) requiring packed conversion
    # - '2025 FA22': Provisional Designation (String input)
    target_objects = ['3I', 134341, '2025 FA22', 433]
    target_types = ['comet_number', 'asteroid_number', 'asteroid_designation', 'asteroid_number']

    try:
        # Fetch expected names from MPC (The Source of Truth)
        batch.get_observations(target_objects, id_types=target_types)
    except Exception as e:
        pytest.skip(f"Skipping online integration test: {e}")

    # Raw lines matching the query above
    lines_main = [
        '0003I         C2025 07 03.94821918 01 44.983-18 39 35.17         18.0 GW#0JSN598', # 3I
        'D4341J79M00A*4A1979 06 25.66181 20 27 06.64 -15 37 11.5          19.0   M4986413', # 134341
        '     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'  # 2025 FA22
    ]
    lines_sat = [
        '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500', # 433 (S)
        '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'  # 433 (s)
    ]

    # Run Parser
    table_main = parse_80cols_data(lines_main)
    table_sat = parse_80cols_data(lines_sat)

    # Combine results (Order: 3I, 134341, 2025 FA22, 433)
    parsed_ids = table_main['number'].astype(str).tolist() + table_sat['number'].astype(str).tolist()

    # Normalize BatchMPC objects to strings for comparison
    batch_ids = [str(x) for x in batch.MPC_objects]

    # Assertions
    print(f"\nBatchMPC: {batch_ids}")
    print(f"Parser:   {parsed_ids}")

    assert len(parsed_ids) == len(batch_ids)
    for pid, bid in zip(parsed_ids, batch_ids):
        assert pid == bid, f"Mismatch! Parser: '{pid}' vs BatchMPC: '{bid}'"