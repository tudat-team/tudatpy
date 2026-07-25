import pytest
from unittest.mock import MagicMock, patch
from tudatpy.data.discos import DiscosQuery

# ------------------------------------------------
# Mock Responses Taken from discos API web request
MOCK_RESPONSE_DISCOS_ID = {
    "data": {
        "type": "object",
        "id": "25544",
        "attributes": {
            "cosparId": "1998-068B",
            "vimpelId": None,
            "satno": 25547,
            "name": "Delta II second stage (Delta 7925)",
            "objectClass": "Rocket Body",
            "mass": 942.52,
            "shape": "Cyl + Cone",
            "width": None,
            "height": 5.97,
            "depth": None,
            "diameter": 2.44,
            "span": 5.97,
            "xSectMax": 15.2988933574707,
            "xSectMin": 4.67594650560305,
            "xSectAvg": 13.7787112193795,
            "firstEpoch": "1998-11-22",
            "mission": None,
            "predDecayDate": None,
            "active": None,
            "cataloguedFragments": 0,
            "onOrbitCataloguedFragments": 0,
        },
        "relationships": {},
        "links": {},
    }
}

MOCK_RESPONSE_NORAD_ID = {
    "data": [
        {
            "type": "object",
            "id": "44371",
            "attributes": {
                "cosparId": "1998-067A",
                "vimpelId": None,
                "satno": 25544,
                "name": "International Space Station",
                "objectClass": "Payload",
                "mass": 450000,
            },
            "relationships": {},
            "links": {},
        }
    ]
}
# ------------------------------------------------


@pytest.fixture
def discos_query():
    return DiscosQuery(token="fake_token_123")


@patch("tudatpy.data.discos.discos.requests.get")
def test_query_norad_success(mock_get, discos_query):
    """
    Tests Query by NORAD ID (Default).
    API Behavior: Returns a LIST of objects inside 'data'.
    """
    mock_response = MagicMock()
    mock_response.ok = True
    mock_response.json.return_value = MOCK_RESPONSE_NORAD_ID
    mock_get.return_value = mock_response

    result = discos_query.query_object("25544", is_discos_id=False, verbose=False)

    # Verify request
    mock_get.assert_called_once()
    args, kwargs = mock_get.call_args
    assert "filter=eq(satno,25544)" in args[0]
    assert "Authorization" in kwargs.get("headers", {})

    # Verify result
    assert isinstance(result, dict)
    assert result.get("mass") == 450000
    assert result.get("name") == "International Space Station"


@patch("tudatpy.data.discos.discos.requests.get")
def test_query_discos_id_success(mock_get, discos_query):
    """
    Tests query by DISCOS ID using is_discos_id = True.
    API Behavior: Returns a SINGLE DICT inside 'data'.
    """
    mock_response = MagicMock()
    mock_response.ok = True
    mock_response.json.return_value = MOCK_RESPONSE_DISCOS_ID
    mock_get.return_value = mock_response

    result = discos_query.query_object("25544", is_discos_id=True, verbose=False)

    # Verify request
    mock_get.assert_called_once()
    args, kwargs = mock_get.call_args
    assert "/api/objects/25544" in args[0]
    assert "Authorization" in kwargs.get("headers", {})

    # Verify result
    assert isinstance(result, dict)
    assert result.get("mass") == 942.52
    assert result.get("name") == "Delta II second stage (Delta 7925)"


@patch("tudatpy.data.discos.discos.requests.get")
def test_query_not_found(mock_get, discos_query):
    """
    Tests invalid satellite norad and discos id.
    API Behavior: Returns an empty list inside 'data'.
    """
    mock_response = MagicMock()
    mock_response.ok = True
    mock_response.json.return_value = {"data": []}
    mock_get.return_value = mock_response

    result_norad = discos_query.query_object("00000", verbose=False)
    result_discos = discos_query.query_object("00000", is_discos_id=True, verbose=False)

    assert result_norad is None
    assert result_discos is None


@patch("tudatpy.data.discos.discos.requests.get")
def test_api_error(mock_get, discos_query):
    """
    Scenario: API server error or Auth failure.
    API Behavior: Returns ok=False and an error message.
    """
    mock_response = MagicMock()
    mock_response.ok = False
    mock_response.json.return_value = {"errors": "This is an Error"}
    mock_get.return_value = mock_response

    result = discos_query.query_object("12345", verbose=False)

    mock_get.assert_called_once()
    assert result is None
