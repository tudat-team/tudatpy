import pytest
from tudatpy.data.mpc import BatchMPC
from tudatpy.data.mpc.parser_80col.parsers import parse_80cols_data, parse_80cols_file
from tudatpy.data.mpc.parser_80col.unpackers import (
    unpack_permanent_minor_planet,
    unpack_provisional_minor_planet,
    unpack_provisional_comet_or_satellite,
    unpack_permanent_natural_satellite
)

# ==============================================================================
# SECTION A: UNIT TESTS (Logic Only, Offline)
# ==============================================================================

# --- 1. Permanent Minor Planet Tests (Asteroids) ---
@pytest.mark.parametrize("packed_input, expected_output", [
    # Case 3: Basic Range (0 - 99,999)
    ("00001", "1"),
    ("00433", "433"),
    ("99999", "99999"),
    # Case 2: Standard Range (100,000 - 619,999)
    ("A0000", "100000"),
    ("A0001", "100001"),
    ("Z9999", "359999"), # Z=35
    ("a0000", "360000"), # a=36
    ("z9999", "619999"), # z=61
    # Case 1: Extended Range (> 620,000)
    ("~0000", "620000"),
    ("~0001", "620001"),
    ("~AZaz", "3140113"), # Complex Base62 check
])

def test_unpack_permanent_minor_planet_success(packed_input, expected_output):
    assert unpack_permanent_minor_planet(packed_input) == expected_output

@pytest.mark.parametrize("invalid_input, error_msg", [
    ("123", "Invalid packed length"),
    ("123456", "Invalid packed length"),
    ("A123x", "Suffix must be numeric"),
    ("?0000", "Invalid start character"),
])
def test_unpack_permanent_minor_planet_errors(invalid_input, error_msg):
    with pytest.raises(ValueError, match=error_msg):
        unpack_permanent_minor_planet(invalid_input)


# --- 2. Provisional Minor Planet Tests ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("PLS2040", "2040 P-L"),    # Survey
    ("J98M00M", "1998 MM"),     # Standard
    ("K07T08A", "2007 TA8"),    # Numeric Cycle
    ("K07Tf8A", "2007 TA418"),  # Base62 Cycle
])
def test_unpack_provisional_minor_planet_success(packed_input, expected_output):
    assert unpack_provisional_minor_planet(packed_input) == expected_output

@pytest.mark.parametrize("invalid_input", ["J98", "J98M0012", "X99A00A"])
def test_unpack_provisional_minor_planet_errors(invalid_input):
    with pytest.raises(ValueError):
        unpack_provisional_minor_planet(invalid_input)


# --- 3. Provisional Comet & Satellite Tests ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("K07C01A", "2007 C1-A"), # Comet
    ("J95X020", "1995 X2"),   # Satellite
    ("K21A05b", "2021 A5-B"), # Fragment casing
])
def test_unpack_provisional_comet_sat_success(packed_input, expected_output):
    assert unpack_provisional_comet_or_satellite(packed_input) == expected_output


# --- 4. Permanent Natural Satellite Tests (Uses TudatPy Logic) ---
@pytest.mark.parametrize("packed_input, expected_output", [
    ("J001S", "Jupiter I"),
    ("S010S", "Saturn X"),
    ("N002S", "Neptune II"),
])
def test_unpack_permanent_natural_satellite_success(packed_input, expected_output):
    assert unpack_permanent_natural_satellite(packed_input) == expected_output


# ==============================================================================
# SECTION B: PARSER LOGIC TESTS (Data Handling, Offline)
# ==============================================================================

def test_80cols_line_parser_logic():
    """
    Tests filtering, column extraction, and format handling on raw string lists.
    No file I/O or internet required.
    """
    # 1. Setup Data
    line_atlas = '0003I         C2025 07 03.94821918 01 44.983-18 39 35.17         18.0 GW#0JSN598'
    line_eros_valid = '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500'
    line_eros_skip = '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57' # 's' line
    line_charon = 'D4341J79M00A*4A1979 06 25.66181 20 27 06.64 -15 37 11.5          19.0   M4986413'
    line_2025FA22 = '     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'

    lines_to_parse = [line_atlas, line_eros_valid, line_eros_skip, line_charon, line_2025FA22]

    # 2. Parse
    parsed_table = parse_80cols_data(lines_to_parse)

    # 3. Assertions
    # A. Check filtering (satellite parallax 's' line should be dropped)
    assert len(parsed_table) == 4, "Failed to filter 's' line"

    # B. Check ID extraction (Strings)
    # Adjust strict equality based on your specific parser return types
    ids = parsed_table['number'].astype(str).tolist()
    assert ids[0] == '3I'          # Comet Unpacking
    assert ids[1] == '433'         # Asteroid Unpacking
    assert ids[2] == '134341'      # 'D4341' -> 134341 Unpacking
    assert ids[3] == '2025 FA22'   # Provisional Unpacking

def test_80cols_malformed_lines():
    """Tests that invalid line lengths or formats raise ValueError."""
    bad_len = '    K25F22A  C2025 10 13.24277'
    bad_fmt = '      K25F22A  C2025 10 13.24277 00 20 4576 +25 53 06.1          18.3 RrET147718' # No RA space

    with pytest.raises(ValueError):
        parse_80cols_data([bad_len])
    with pytest.raises(ValueError):
        parse_80cols_data([bad_fmt])

def test_parse_80cols_file_io(tmp_path):
    """Tests file reading using a temporary file (clean cleanup)."""
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
# SECTION C: INTEGRATION TEST (BatchMPC Consistency, Online (Needs BatchMPC))
# ==============================================================================

def test_batch_mpc_vs_parser_consistency():
    """
    CRITICAL: Ensures the Parser naming convention matches BatchMPC/Tudat naming convention.
    If this fails, downstream filters will mismatch (e.g. str vs int, or packed vs unpacked).
    """
    # 1. Setup "Ground Truth" (BatchMPC)
    batch = BatchMPC()
    target_objects = ['3I', 134341, '2025 FA22', 433]
    target_types = ['comet_number', 'asteroid_number', 'asteroid_designation', 'asteroid_number']

    try:
        # Fetch expected names from MPC
        batch.get_observations(target_objects, id_types=target_types)
    except Exception as e:
        pytest.skip(f"Skipping online integration test: {e}")

    # 2. Setup "Subject" (Raw Lines matching the query above)
    ground_station_lines = [
        '0003I         C2025 07 03.94821918 01 44.983-18 39 35.17         18.0 GW#0JSN598', # 3I
        'D4341J79M00A*4A1979 06 25.66181 20 27 06.64 -15 37 11.5          19.0   M4986413', # 134341
        '     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'  # 2025 FA22
    ]

    satellite_lines = [
        '00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500',
        '00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57'
    ] # this is a satellite observation

    # 3. Run Parser
    parsed_ground_station_table = parse_80cols_data(ground_station_lines)
    parsed_satellite_table = parse_80cols_data(satellite_lines)
    parsed_ids = parsed_ground_station_table['number'].astype(str).tolist() + parsed_satellite_table['number'].astype(str).tolist()
    batch_ids = [str(x) for x in batch.MPC_objects]

    # 4. Compare
    print(f"\nBatchMPC: {batch_ids}")
    print(f"Parser:   {parsed_ids}")

    assert len(parsed_ids) == len(batch_ids)
    for pid, bid in zip(parsed_ids, batch_ids):
        assert pid == bid, f"Mismatch! Parser: '{pid}' vs BatchMPC: '{bid}'"