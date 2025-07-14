

function(TUDAT_ADD_DATA_DIR arg1)
    get_filename_component(dir_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${arg1} DESTINATION ${INSTALL_DATA_DIR}/${dir_name})
endfunction()


function(TUDAT_ADD_TEST_DATA arg1)
    # We move the test files into the testing directory
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${arg1}/" DESTINATION "${CMAKE_BINARY_DIR}/${arg1}/")
endfunction()
