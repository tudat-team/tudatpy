def test_import_from():
    from tudatpy import constants

    _ = constants.JULIAN_DAY


def test_import_as():
    import tudatpy.constants as constants

    _ = constants.JULIAN_DAY
