include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/tests/data)
add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)
add_subdirectory(src/support)

# Unzip compressed data and place in build tree.
add_custom_target( unpack_test_data ALL)
add_custom_command(
    TARGET unpack_test_data PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/tests/data/
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/data.zip
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    COMMENT "Unpacking data.zip"
    VERBATIM
)

# Extract the archived reference headers and data before building their test consumers.
foreach(test_target IN ITEMS
        test_earth_orientation_EarthOrientationDataEvaluator
        test_earth_orientation_EopReader
        test_earth_orientation_HighFrequencyEarthOrientation
        test_earth_orientation_ShortPeriodEopCorrections
        test_ground_stations_AtmosphericTideDeformation
        test_ground_stations_AtmosphericTideLoading
        test_ground_stations_EarthDeformationTimeRouting
        test_ground_stations_Iers2010DeformationModel
        test_ground_stations_OceanPoleTideDisplacement
        test_ground_stations_OceanTideDisplacement
        test_ground_stations_StationCoordinateHistory
        test_ground_stations_StationHistoryFrameTransformation
        test_ground_stations_StationReferencePoint
        test_reference_frames_TerrestrialFrameTransformation
        test_io_CrdAndSinexReaders
        test_io_SinexStationHistory
        test_io_StationLoadingData)
    add_dependencies(${test_target} unpack_test_data)
endforeach()
