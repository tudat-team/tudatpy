"""
Internal helper functions for unpacking various MPC 80-column designation formats
Users are referred to the Minor Planet Center website:
https://minorplanetcenter.net/iau/info/OpticalObs.html
for info about the designation format for Optical Astrometric Observations Of Comets,
Minor Planets and Natural Satellites.

Packing:
Converting a long, human-readable designation into a short, coded version that fits the 80-column space.
Unpacking:
Converting the short, packed form back into a human-readable designation.
"""

from tudatpy.util._support import transform_integer_to_roman_number

_BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
BASE62_MAP = {char: i for i, char in enumerate(_BASE62)}
""" Base-62 characters (0-9, A-Z, a-z) for unpacking.
  Letter,     Dates        |        Letter ,      Dates
    A   ,   Jan. 1-15      |       B       ,  Jan. 16-31
    C   ,   Feb. 1-15      |       D       ,  Feb. 16-29
    E   ,   Mar. 1-15      |       F       ,  Mar. 16-31
    G   ,   Apr. 1-15      |       H       ,  Apr. 16-30
    J   ,   May  1-15      |       K       ,  May  16-31
    L   ,   June 1-15      |       M       ,  June 16-30
    N   ,   July 1-15      |       O       ,  July 16-31
    P   ,   Aug. 1-15      |       Q       ,  Aug. 16-31
    R   ,   Sept.1-15      |       S       ,  Sept.16-30
    T   ,   Oct. 1-15      |       U       ,  Oct. 16-31
    V   ,   Nov. 1-15      |       W       ,  Nov. 16-30
    X   ,   Dec. 1-15      |       Y       , Dec. 16-31
             
             (Note: I is omitted and Z is unused)
             
The order within the month is indicated using letters as follows:
    A = 1st     B = 2nd     C = 3rd     D = 4th     E = 5th
    F = 6th     G = 7th     H = 8th     J = 9th     K = 10th
    L = 11th    M = 12th    N = 13th    O = 14th    P = 15th
    Q = 16th    R = 17th    S = 18th    T = 19th    U = 20th
    V = 21st    W = 22nd    X = 23rd    Y = 24th    Z = 25th

                        (Note: I is omitted)
"""

CENTURY_MAP = {'I': '18', 'J': '19', 'K': '20'}
""" Mapping for packed designations: century mapping, e.g., 'I' -> '18th century, etc...'] """

SURVEY_MAP = {'PLS': 'P-L', 'T1S': 'T-1', 'T2S': 'T-2', 'T3S': 'T-3'}
""" Mapping for packed designations: survey map.    
   Survey,                         Identifier
   Palomar-Leiden (1960),             P-L
   First Trojan Survey (1971),        T-1
   Second Trojan Survey (1973),       T-2
   Third Trojan Survey (1977),        T-3 
"""

PLANET_MAP = {'J': 'Jupiter', 'S': 'Saturn', 'U': 'Uranus', 'N': 'Neptune'}
""" Mapping for packed designations: planet map. 
Permanent designations for natural satellites are the Roman numerals. 
A permanent natural-satellite packed designation starts with a letter indicating to which
planet the satellite belongs ("J" = Jupiter, "S" = Saturn, "U" = Uranus, "N" = Neptune),
followed by a three-character zero-padded right-justified string containing
the numerical representation of the Roman numeral, followed by "S".
E.g., Jupiter XIII is represented as "J013S", Saturn X as "S010S" and Neptune II as "N002S".
"""

OBS_TYPES_TO_DROP = ["x", "X", "V", "v", "W", "w", "R", "r", "Q", "q", "O"]
""" Observation types to be ignored during parsing.
    See: https://minorplanetcenter.net/iau/info/OpticalObs.html (Note 2)
"""

map_reason_to_drop = {
    "x": "Replaced discovery observation (suppressed in residual blocks).",
    "X": "Replaced discovery observation (suppressed in residual blocks).",
    "V": "Roving Observer observation.",
    "v": "Roving Observer observation.",
    "W": "Roving Observer observation (converted from XML).",
    "w": "Roving Observer observation (converted from XML).",
    "R": "Radar observation.",
    "r": "Radar observation.",
    "Q": "Radar observation (converted from XML).",
    "q": "Radar observation (converted from XML).",
    "O": "Offset observation (natural satellites)."
}

# --- Unpacking Functions ---
def unpack_permanent_minor_planet(packed: str) -> str:
    """
    Unpacks a 5-character packed permanent minor planet designation.
    Ref: https://minorplanetcenter.net/iau/info/PackedDes.html#perm

    Format Rules:
    1. 00001 - 99999 : Straight digits (0-9)
    2. A0000 - z9999 : Base62 char + 4 digits (100,000 - 619,999)
    3. ~0000 - ~zzzz : Tilde + 4 Base62 chars (620,000 +)
    """

    # 1. Validation: Length Check
    if len(packed) != 5:
        raise ValueError(f"Invalid packed length: '{packed}'. Must be exactly 5 characters.")

    first_char = packed[0]

    # --- Case 1: Extended Range (> 620,000) ---
    # Format: ~ + 4 Base62 characters
    if first_char == '~':
        value = 0
        for char in packed[1:]:
            if char not in BASE62_MAP:
                raise ValueError(f"Invalid Base62 character '{char}' in extended packed designation '{packed}'.")
            value = value * 62 + BASE62_MAP[char]

        return str(value + 620000)

    # --- Case 2: Standard Range (100,000 - 619,999) ---
    # Format: Base62 char (A-z) + 4 Digits
    elif first_char.isalpha():
        # Validate suffix is numeric
        if not packed[1:].isdigit():
            raise ValueError(f"Invalid format for range 100k+: '{packed}'. Suffix must be numeric digits.")

        prefix_val = BASE62_MAP[first_char]
        suffix_val = int(packed[1:])

        return str(prefix_val * 10000 + suffix_val)

    # --- Case 3: Basic Range (0 - 99,999) ---
    # Format: 5 Digits
    elif first_char.isdigit():
        # Validate entire string is numeric
        if not packed.isdigit():
            raise ValueError(f"Invalid format for basic range: '{packed}'. Must be all digits.")

        return str(int(packed))

    # --- Case 4: Invalid Start Character ---
    else:
        raise ValueError(f"Invalid start character '{first_char}' in packed designation '{packed}'. Expected digit, letter, or '~'.")


def unpack_provisional_minor_planet(packed: str) -> str:
    """
    Unpacks a 7-character packed provisional minor planet designation.

    Parameters
    ----------
    packed : str
        The 7-character packed provisional minor planet designation.

    Returns
    -------
    str
        The unpacked provisional minor planet designation.
    """
    # 1. Validation: Length Check
    if len(packed) != 7:
        raise ValueError(f"Invalid packed length: '{packed}'. Must be exactly 7 characters.")

    # Case 1: Survey designations (e.g., PLS2040 for 2040 P-L)
    if packed[0:3] in SURVEY_MAP:
        number = int(packed[3:])
        survey_type = SURVEY_MAP[packed[0:3]]
        return f"{number} {survey_type}"

    # Case 2: Standard provisional designations (e.g., K07Tf8A for 2007 TA418)
    if CENTURY_MAP.get(packed[0], '??') == '??':
        raise ValueError(f'Invalid letter {packed[0]} does not exist in CENTURY MAP: {CENTURY_MAP}.')
    year = CENTURY_MAP.get(packed[0], '??') + packed[1:3]
    half_month_1 = packed[3]
    half_month_2 = packed[6]
    cycle_packed = packed[4:6]

    cycle = 0
    if cycle_packed[0].isdigit():
        cycle = int(cycle_packed)
    else:  # Letter-digit format (e.g., 'f8' -> 418)
        val1 = BASE62_MAP[cycle_packed[0]]
        val2 = int(cycle_packed[1])
        cycle = val1 * 10 + val2

    if cycle == 0:
        return f"{year} {half_month_1}{half_month_2}"
    else:
        return f"{year} {half_month_1}{half_month_2}{cycle}"


def unpack_provisional_comet_or_satellite(packed: str) -> str: 
    """
    Unpacks a 7-char provisional comet/satellite designation.
    
    Parameters
    ---------- 
    packed : str
        The 7-character packed provisional comet or satellite designation.
        Examples: 'K07C01A' for 2007 C1 A, 'J95X020' for 1995 X2.

    Returns
    -------
    str
        The unpacked provisional comet or satellite designation.
    """

    if CENTURY_MAP.get(packed[0], '??') == '??':
        raise ValueError(f'Invalid letter {packed[0]} does not exist in CENTURY MAP: {CENTURY_MAP}.')
    year = CENTURY_MAP.get(packed[0], '??') + packed[1:3]
    half_month = packed[3]
    order = int(packed[4:6])
    fragment = packed[6]

    desig = f"{year} {half_month}{order}"
    if fragment != '0':  # Append fragment letter if it exists
        desig += f"-{fragment.upper()}"
    return desig

def unpack_permanent_natural_satellite(packed: str) -> str:
    """
    Unpacks a 5-char permanent natural satellite designation.

    Parameters
    ----------
    packed : str
        The 5-character packed permanent natural satellite designation.
        Examples: 'J013S' for Jupiter XIII, 'S010S' for Saturn X, 'N002S' for Neptune II.

    Returns
    -------
    str
        The unpacked permanent natural satellite designation.
    """
    if len(packed) != 5:
        raise ValueError(f"Invalid packed length: '{packed}'. Must be exactly 5 characters.")

    if PLANET_MAP.get(packed[0], 'Unknown Planet') == 'Unknown Planet':
        raise ValueError(f'Invalid planet code {packed[0]} does not exist in PLANET MAP: {PLANET_MAP}.')
    planet_name = PLANET_MAP.get(packed[0], "Unknown Planet")


    if planet_name == "Unknown Planet":
        raise ValueError("Unknown planet with code " + packed[0] + " found in PLANET_MAP. "
                         "Make sure your lines are compliant with the MPC format "
                         "(https://minorplanetcenter.net/iau/info/OpticalObs.html), or update the PLANET_MAP."
                        )
    number = int(packed[1:4])
    return f"{planet_name} {transform_integer_to_roman_number(number)}"