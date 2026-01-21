cmake_minimum_required(VERSION 3.20)

# Display basic cmake configuration
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

# Remove caches if present [Trick to recompile changes for BUILD-TREE]
# From tudat cmake, might not be needed.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()


# Policies
cmake_policy(SET CMP0069 NEW)  # Missing documentation
cmake_policy(SET CMP0048 NEW)  # Missing documentation

# Get version from file
file(STRINGS "version" VERSION_LABEL)
list(GET VERSION_LABEL 0 VERSION_LABEL)
string(
    REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)"
    TUDATPY_VERSION "${VERSION_LABEL}"
)

# Project declaration.
project(tudatpy VERSION "${TUDATPY_VERSION}" LANGUAGES CXX C)
set(BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Set position independent code [From tudat cmake]
set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Allow access to cmake modules
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudat")
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudat/yolo")
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudatpy")

# Additional statements from original cmake files (might not be needed)
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
enable_testing()
if (MINGW)
    # https://github.com/tudat-team/tudatpy-feedstock/issues/3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_hypot=hypot")
    # https://stackoverflow.com/questions/16596876/object-file-has-too-many-sections
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    # https://sourceforge.net/p/mingw-w64/mailman/message/33154210/
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--large-address-aware")
endif ()

if (MSVC)
    # https://docs.microsoft.com/en-us/cpp/build/reference/bigobj-increase-number-of-sections-in-dot-obj-file?view=vs-2019
    # This is equivalent to the big-obj for MinGW, this should stop error C1128.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()

# +============================================================================
# BUILD OPTIONS
#  Offer the user the choice of defining the build variation.
# +============================================================================
option(
    TUDAT_BUILD_TESTS
    "Build C++ tests for tudat."
    ON
)
option(
    TUDAT_DOWNLOAD_AND_BUILD_BOOST
    "Downloads and builds boost"
    OFF
)
option(
    TUDAT_BUILD_WITH_PROPAGATION_TESTS
    "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)"
    OFF
)
option(
    TUDAT_BUILD_WITH_ESTIMATION_TOOLS
    "Build tudat with estimation tools."
    ON
)
option(
    TUDAT_BUILD_TUDAT_TUTORIALS
    "Build C++ tutorials."
    OFF
)
option(
    TUDAT_BUILD_STATIC_LIBRARY
    "Build tudat as a static library, instead of dynamic."
    ON
)
option(
    TUDAT_BUILD_WITH_FILTERS
    "Build Tudat with preliminary Kalman Filter functionality."
    OFF
)
option(
    TUDAT_BUILD_WITH_SOFA_INTERFACE
    "Build Tudat with sofa interface."
    ON
)
option(
    TUDAT_BUILD_WITH_JSON_INTERFACE
    "Build Tudat with json interface."
    OFF
)
option(
    TUDAT_BUILD_WITH_FFTW3
    "Build Tudat with FFTW3."
    OFF
)
option(
    TUDAT_BUILD_WITH_PAGMO
    "Build Tudat with pagmo."
    OFF
)
if(CMAKE_CXX_SIMULATE_ID MATCHES "MSVC")
    # Build extended precision propagation tools.
    option(
        TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS
        "Build tudat with extended precision propagation tools."
        ON
    )
else()
    option(
        TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS
        "Build tudat with extended precision propagation tools."
        ON
    )
endif()

# Build as part of a GitHub Actions workflow
# Option enables the use of ccache by MSVC
# see https://github.com/ccache/ccache/wiki/MS-Visual-Studio
option(
    TUDAT_BUILD_GITHUB_ACTIONS
    "Build as part of GitHub Actions workflow."
    OFF
)

# +============================================================================
# EMSCRIPTEN / WEBASSEMBLY CONFIGURATION
#  Detect Emscripten and configure appropriate build settings
# +============================================================================
if (CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR EMSCRIPTEN)
    message(STATUS "")
    message(STATUS "*** EMSCRIPTEN/WEBASSEMBLY BUILD DETECTED ***")
    message(STATUS "")

    # Mark this as a WASM build
    set(TUDAT_BUILD_WASM ON)

    # Force settings that are incompatible or unnecessary for WASM
    set(TUDAT_BUILD_TESTS OFF)
    set(TUDAT_BUILD_TUDAT_TUTORIALS OFF)
    set(TUDAT_BUILD_STATIC_LIBRARY ON)

    # SOFA is required for WASM builds due to unconditional includes in source
    # It will be fetched automatically
    set(TUDAT_BUILD_WITH_SOFA_INTERFACE ON)
    set(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS OFF)

    # CSPICE's f2c-generated code has function signature mismatches (e.g., s_cat, s_copy
    # return void but are called expecting int). WASM is strict about this and will crash
    # without emulation. This flag must be set globally to apply to all linked libraries.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sEMULATE_FUNCTION_POINTER_CASTS=1")

    message(STATUS "WASM build: Tests disabled")
    message(STATUS "WASM build: Tutorials disabled")
    message(STATUS "WASM build: Static library forced")
    message(STATUS "")

    # Create update-emscripten target
    # This allows users to update the Emscripten SDK version
    set(EMSDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.emsdk")
    set(EMSDK_INSTALL_DIR "${EMSDK_DIR}/emsdk")

    if(WIN32)
        set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk.bat")
    else()
        set(EMSDK_EXE "${EMSDK_INSTALL_DIR}/emsdk")
    endif()

    add_custom_target(update-emscripten
        COMMAND git -C "${EMSDK_INSTALL_DIR}" pull
        COMMAND "${EMSDK_EXE}" install ${EMSDK_VERSION}
        COMMAND "${EMSDK_EXE}" activate ${EMSDK_VERSION}
        COMMAND ${CMAKE_COMMAND} -E echo "${EMSDK_VERSION}" > "${EMSDK_DIR}/.installed_version"
        WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
        COMMENT "Updating Emscripten SDK to version ${EMSDK_VERSION}"
        VERBATIM
    )

    # Target to list available Emscripten versions
    add_custom_target(list-emscripten-versions
        COMMAND "${EMSDK_EXE}" list
        WORKING_DIRECTORY "${EMSDK_INSTALL_DIR}"
        COMMENT "Listing available Emscripten SDK versions"
        VERBATIM
    )
else()
    set(TUDAT_BUILD_WASM OFF)
endif()

if (TUDAT_BUILD_GITHUB_ACTIONS AND MSVC)
    find_program(ccache_exe ccache)
    if(ccache_exe)
        message(STATUS "CONFIGURING MSVC FOR CCACHE")
        file(COPY_FILE ${ccache_exe} ${CMAKE_BINARY_DIR}/cl.exe ONLY_IF_DIFFERENT)

        # By default Visual Studio generators will use /Zi which is not compatible
        # with ccache, so tell Visual Studio to use /Z7 instead.
        message(STATUS "Setting MSVC debug information format to 'Embedded'")
        set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

        set(CMAKE_VS_GLOBALS
        "CLToolExe=cl.exe"
        "CLToolPath=${CMAKE_BINARY_DIR}"
        "TrackFileAccess=false"
        "UseMultiToolTask=true"
        "DebugInformationFormat=OldStyle"
        )
    endif()
endif()

include(show_build_configuration)

# +============================================================================
# TUDAT DEFINITIONS (Based on user configuration)
# +============================================================================
set(Tudat_DEFINITIONS "")
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FILTERS=${TUDAT_BUILD_WITH_FILTERS}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_SOFA_INTERFACE=${TUDAT_BUILD_WITH_SOFA_INTERFACE}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FFTW3=${TUDAT_BUILD_WITH_FFTW3}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_JSON_INTERFACE=${TUDAT_BUILD_WITH_JSON_INTERFACE}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}"
)

# TODO: NOT SURE IF THIS BELONGS HERE (Comes from tudatpy cmake)
if (_ENABLE_EXTENDED_ALIGNED_STORAGE)
    add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()
if (NOT TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=1)
endif ()
if (NOT TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=1)
endif ()

# +============================================================================
# INSTALL TREE CONFIGURATION (Project name independent)
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
# TODO: Consider if this can work with <tudat/paths.hpp> going forward.

include(install_tree_configuration)

# +============================================================================
# DEPENDENCIES SPECIFIC TO BUILD OPTIONS
# +============================================================================
if (NOT TUDAT_BUILD_STATIC_LIBRARY)
    # https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

# Set required boost components.
set(
    _TUDAT_REQUIRED_BOOST_LIBS
    filesystem system regex date_time thread chrono atomic unit_test_framework
)
# if (TUDAT_BUILD_TESTS)
#     # Add the unit test framework, if needed.
#     list(APPEND _TUDAT_REQUIRED_BOOST_LIBS unit_test_framework)
# endif ()

# Download boost if option is used.
if (
    TUDAT_DOWNLOAD_AND_BUILD_BOOST
    OR (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1)
)
    set(TUDAT_DOWNLOAD_AND_BUILD_BOOST ON)
    # Set the required boost components (for tudat)
    set(BoostComponents ${_TUDAT_REQUIRED_BOOST_LIBS})
    # Pick a boost version
    if (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1)
        set(BoostVersion ${TUDAT_BOOST_VERSION})
        set(BoostSHA1 ${TUDAT_BOOST_SHA1})
    else ()
        set(BoostVersion 1.67.0)
        set(BoostSHA1 694ae3f4f899d1a80eb7a3b31b33be73c423c1ae)
    endif ()
    message(STATUS "Downloading Boost Version ${BoostVersion}")
    include(cmake_modules/add_boost.cmake)
    message(STATUS "Downloading Boost Version ${BoostVersion} - done")
    #    set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/boost/stage/lib)
    #    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/boost)
endif ()

# Set cmake build-type.
# If not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file), the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

# Set compiler based on user input
include(compiler)

# Find Boost and check GLBICXX_USE_CXX11_ABI (TODO: Better description)
include(TudatFindBoost)
message(STATUS "Checking for _GLIBCXX_USE_CXX11_ABI definition...")
if (DEFINED _GLIBCXX_USE_CXX11_ABI)
    if (_GLIBCXX_USE_CXX11_ABI)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=1")
        set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=1")
        message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=1 added to CMAKE_CXX_FLAGS!")
    elseif (NOT _GLIBCXX_USE_CXX11_ABI)
        message(Boost_VERSION_MAJOR: ${Boost_VERSION_MAJOR})
        message(Boost_VERSION_MINOR: ${Boost_VERSION_MINOR})
        message(TUDAT_BUILD_GNU: ${TUDAT_BUILD_GNU})
        if (
            (${Boost_VERSION_MAJOR} LESS_EQUAL "1")
            AND (${Boost_VERSION_MINOR} LESS_EQUAL "67") AND ${TUDAT_BUILD_GNU}
        )
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(
                Tudat_DEFINITIONS
                "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0"
            )
            message(
                STATUS
                "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!"
            )
        elseif (
            (${Boost_VERSION_MAJOR} GREATER_EQUAL "1")
            AND (${Boost_VERSION_MINOR} GREATER_EQUAL "70")
            AND ${TUDAT_BUILD_GNU}
        )
            #            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!")
        elseif (${TUDAT_BUILD_MSVC})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(
                Tudat_DEFINITIONS
                "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0"
            )
            message(
                STATUS
                "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!"
            )
        else ()
            message(STATUS "_GLIBCXX_USE_CXX11_ABI was set as 0, but since the boost version is >=1.68|>1.70 [${Boost_VERSION}] and compiler_id is GNU, it's being ignored!")
            add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
        endif ()
    endif ()
else ()
    message(STATUS "_GLIBCXX_USE_CXX11_ABI was not found.")
endif ()

if ((UNIX OR (WIN32 AND CLANG)) AND NOT TUDAT_BUILD_WASM)
    # Skip -isystem for WASM - Emscripten's sysroot already has correct include paths
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# Get tudat resource manager dependency.
if(TUDAT_BUILD_WASM)
    # For WASM builds, TudatResources is optional - data can be loaded via virtual filesystem
    find_package(TudatResources QUIET)
    if(TudatResources_FOUND)
        message(STATUS "TudatResources found for WASM build")
        include_directories(TudatResources_INCLUDE_DIRS)
    else()
        message(STATUS "TudatResources not found - WASM build will proceed without it")
        set(TudatResources_INCLUDE_DIRS "")
    endif()
else()
    find_package(TudatResources REQUIRED)
    message(STATUS "TudatResources_INCLUDE_DIRS: ${TudatResources_INCLUDE_DIRS}")
    include_directories(TudatResources_INCLUDE_DIRS)
endif()

# Find Eigen3 library on local system and include directories.
# For WASM builds, download Eigen if not found (it's header-only)
if(TUDAT_BUILD_WASM)
    find_package(Eigen3 QUIET)
    if(NOT Eigen3_FOUND)
        message(STATUS "Eigen3 not found, downloading for WASM build...")
        include(FetchContent)
        FetchContent_Declare(
            eigen
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
            URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
        )
        # Use FetchContent_Populate instead of FetchContent_MakeAvailable
        # This downloads Eigen without adding it as a subdirectory, avoiding
        # Eigen's build system which requires Fortran for BLAS/LAPACK tests
        FetchContent_GetProperties(eigen)
        if(NOT eigen_POPULATED)
            FetchContent_Populate(eigen)
        endif()
        set(EIGEN3_INCLUDE_DIR "${eigen_SOURCE_DIR}")
        # Add Eigen includes immediately after setting the variable
        include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}")
        message(STATUS "Eigen3 downloaded to ${EIGEN3_INCLUDE_DIR}")
    else()
        include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}")
    endif()
else()
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}")
endif()
if (APPLE OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# PAGMO dependency
if (TUDAT_BUILD_WITH_PAGMO)
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO=1)
    find_package(Pagmo REQUIRED)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO=0)
endif ()

# CSPICE dependency
add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
# For WASM builds, automatically download and build CSpice with Emscripten
if(TUDAT_BUILD_WASM)
    find_package(CSpice QUIET 1.0.0)
    if(NOT CSpice_FOUND)
        message(STATUS "CSpice not found, downloading and building for WASM...")
        include(FetchContent)
        FetchContent_Declare(
            cspice
            GIT_REPOSITORY https://github.com/tudat-team/cspice-cmake.git
            GIT_TAG master
        )
        # Use FetchContent_Populate to download first, then apply WASM patches
        FetchContent_GetProperties(cspice)
        if(NOT cspice_POPULATED)
            # Silence CMP0169 warning - we need Populate to patch before add_subdirectory
            if(POLICY CMP0169)
                cmake_policy(PUSH)
                cmake_policy(SET CMP0169 OLD)
            endif()
            FetchContent_Populate(cspice)
            if(POLICY CMP0169)
                cmake_policy(POP)
            endif()

        # Apply WASM compatibility patches to CSPICE
        # These fix f2c function signature mismatches that cause WASM crashes
        set(CSPICE_PATCH_FILE "${CMAKE_SOURCE_DIR}/cmake/patches/cspice-wasm.patch")
        if(EXISTS "${CSPICE_PATCH_FILE}")
            message(STATUS "Applying WASM patches to CSpice...")
            # Patch sig_die.c - change abort() to exit(1) for WASM
            file(READ "${cspice_SOURCE_DIR}/src/cspice/sig_die.c" SIG_DIE_CONTENT)
            if(NOT SIG_DIE_CONTENT MATCHES "__EMSCRIPTEN__")
                string(REPLACE
                    "fflush(stderr);\n      /* now get a core */\n#ifdef SIGIOT"
                    "fflush(stderr);\n#ifdef __EMSCRIPTEN__\n      /* In WASM, abort() causes a trap instruction that results in\n         \"RuntimeError: unreachable\". Use exit(1) instead for clean shutdown. */\n      exit(1);\n#else\n      /* now get a core */\n#ifdef SIGIOT"
                    SIG_DIE_CONTENT "${SIG_DIE_CONTENT}")
                string(REPLACE
                    "abort();\n      }"
                    "abort();\n#endif\n      }"
                    SIG_DIE_CONTENT "${SIG_DIE_CONTENT}")
                file(WRITE "${cspice_SOURCE_DIR}/src/cspice/sig_die.c" "${SIG_DIE_CONTENT}")
                message(STATUS "  Patched sig_die.c")
            endif()

            # Patch s_copy.c - change void to int return type
            file(READ "${cspice_SOURCE_DIR}/src/cspice/s_copy.c" S_COPY_CONTENT)
            if(NOT S_COPY_CONTENT MATCHES "return 0;")
                string(REPLACE
                    "/* assign strings:  a = b */\n\n#ifdef KR_headers\nVOID s_copy"
                    "/* assign strings:  a = b */\n/* Modified for WASM: Changed return type from void to int to match f2c-generated caller expectations */\n\n#ifdef KR_headers\nint s_copy"
                    S_COPY_CONTENT "${S_COPY_CONTENT}")
                string(REPLACE
                    "void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)"
                    "int s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)"
                    S_COPY_CONTENT "${S_COPY_CONTENT}")
                string(REPLACE
                    "*a++ = ' ';\n\t\t}\n\t}"
                    "*a++ = ' ';\n\t\t}\n\treturn 0;\n\t}"
                    S_COPY_CONTENT "${S_COPY_CONTENT}")
                file(WRITE "${cspice_SOURCE_DIR}/src/cspice/s_copy.c" "${S_COPY_CONTENT}")
                message(STATUS "  Patched s_copy.c")
            endif()

            # Patch s_cat.c - change VOID to int return type
            file(READ "${cspice_SOURCE_DIR}/src/cspice/s_cat.c" S_CAT_CONTENT)
            if(NOT S_CAT_CONTENT MATCHES "return 0;")
                string(REPLACE
                    "#endif /* NO_OVERWRITE */\n\n VOID\n#ifdef KR_headers"
                    "#endif /* NO_OVERWRITE */\n\n/* Modified for WASM: Changed return type from VOID to int to match f2c-generated caller expectations */\n int\n#ifdef KR_headers"
                    S_CAT_CONTENT "${S_CAT_CONTENT}")
                string(REPLACE
                    "free(lp1);\n\t\t}\n#endif\n\t}"
                    "free(lp1);\n\t\t}\n#endif\n\treturn 0;\n\t}"
                    S_CAT_CONTENT "${S_CAT_CONTENT}")
                file(WRITE "${cspice_SOURCE_DIR}/src/cspice/s_cat.c" "${S_CAT_CONTENT}")
                message(STATUS "  Patched s_cat.c")
            endif()

            # Patch rsfe.c - change zzsetnnread_ void to int
            file(READ "${cspice_SOURCE_DIR}/src/cspice/rsfe.c" RSFE_CONTENT)
            if(RSFE_CONTENT MATCHES "void.*zzsetnnread_" AND NOT RSFE_CONTENT MATCHES "int.*zzsetnnread_")
                string(REPLACE
                    "void    zzsetnnread_( logical * on );"
                    "int     zzsetnnread_( logical * on );"
                    RSFE_CONTENT "${RSFE_CONTENT}")
                string(REPLACE
                    "void zzsetnnread_( logical * on )\n   {\n   read_non_native = *on;\n   }"
                    "/* Modified for WASM: Changed return type from void to int to match f2c-generated caller expectations */\nint zzsetnnread_( logical * on )\n   {\n   read_non_native = *on;\n   return 0;\n   }"
                    RSFE_CONTENT "${RSFE_CONTENT}")
                file(WRITE "${cspice_SOURCE_DIR}/src/cspice/rsfe.c" "${RSFE_CONTENT}")
                message(STATUS "  Patched rsfe.c")
            endif()
        endif()

        # Now add the subdirectory to build
        set(CSPICE_BUILD_SHARED_LIBS OFF CACHE BOOL "Build static CSpice" FORCE)
        add_subdirectory(${cspice_SOURCE_DIR} ${cspice_BINARY_DIR})

        # Suppress additional Clang/Emscripten warnings for CSPICE legacy f2c-generated C code
        if(TARGET cspice)
            target_compile_options(cspice PRIVATE
                -Wno-unused-but-set-variable
                -Wno-deprecated-non-prototype
            )
        endif()
        endif() # if(NOT cspice_POPULATED)
        set(CSpice_LIBRARIES cspice)
        # Headers are in include/cspice/, so include dir is just include/
        set(CSpice_INCLUDE_DIRS "${cspice_SOURCE_DIR}/include")
        # Add include directory for the build
        include_directories(SYSTEM "${CSpice_INCLUDE_DIRS}")
        message(STATUS "CSpice downloaded to ${cspice_SOURCE_DIR}")
        message(STATUS "CSpice include dir: ${CSpice_INCLUDE_DIRS}")
    endif()
else()
    find_package(CSpice REQUIRED 1.0.0)
endif()

# SOFA dependency
if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    if(TUDAT_BUILD_WASM)
        # For WASM, always download SOFA to avoid using incompatible native builds
        message(STATUS "SOFA: Downloading and building for WASM...")
        include(FetchContent)
        FetchContent_Declare(
            sofa
            GIT_REPOSITORY https://github.com/tudat-team/sofa-cmake.git
            GIT_TAG master
        )
        # Use FetchContent_Populate to download without adding subdirectory
        # Then manually create a static library target
        FetchContent_GetProperties(sofa)
        if(NOT sofa_POPULATED)
            FetchContent_Populate(sofa)
        endif()
        # Collect all SOFA source files
        file(GLOB SOFA_SOURCES "${sofa_SOURCE_DIR}/src/*.c")
        # Create static library manually
        add_library(sofa STATIC ${SOFA_SOURCES})
        target_include_directories(sofa PUBLIC "${sofa_SOURCE_DIR}/include")
        set(Sofa_LIBRARIES sofa)
        set(Sofa_INCLUDE_DIRS "${sofa_SOURCE_DIR}/include")
        include_directories(SYSTEM "${Sofa_INCLUDE_DIRS}")
        message(STATUS "SOFA downloaded to ${sofa_SOURCE_DIR}")
        message(STATUS "SOFA include dir: ${Sofa_INCLUDE_DIRS}")
    else()
        find_package(Sofa REQUIRED 1.0.0)
    endif()
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
endif ()

# NRLMSISE-00 dependency
# For WASM builds, automatically download and build nrlmsise00 with Emscripten
if(TUDAT_BUILD_WASM)
    find_package(nrlmsise00 QUIET 0.1)
    if(NOT nrlmsise00_FOUND)
        message(STATUS "nrlmsise00 not found, downloading and building for WASM...")
        include(FetchContent)
        FetchContent_Declare(
            nrlmsise00_fetch
            GIT_REPOSITORY https://github.com/tudat-team/nrlmsise-00-cmake.git
            GIT_TAG master
        )
        FetchContent_GetProperties(nrlmsise00_fetch)
        if(NOT nrlmsise00_fetch_POPULATED)
            FetchContent_Populate(nrlmsise00_fetch)
        endif()
        # Create include directory structure expected by tudat (nrlmsise00/nrlmsise-00.h)
        set(NRLMSISE00_INCLUDE_BASE "${CMAKE_CURRENT_BINARY_DIR}/nrlmsise00_include")
        file(MAKE_DIRECTORY "${NRLMSISE00_INCLUDE_BASE}/nrlmsise00")
        file(COPY "${nrlmsise00_fetch_SOURCE_DIR}/nrlmsise-00.h"
             DESTINATION "${NRLMSISE00_INCLUDE_BASE}/nrlmsise00")
        # Build nrlmsise00 as a static library (repo doesn't have CMakeLists.txt)
        add_library(nrlmsise00 STATIC
            ${nrlmsise00_fetch_SOURCE_DIR}/nrlmsise-00.c
            ${nrlmsise00_fetch_SOURCE_DIR}/nrlmsise-00_data.c
        )
        target_include_directories(nrlmsise00 PUBLIC
            ${nrlmsise00_fetch_SOURCE_DIR}
            ${NRLMSISE00_INCLUDE_BASE}
        )
        set(NRLMSISE00_LIBRARIES nrlmsise00)
        set(NRLMSISE00_INCLUDE_DIRS ${NRLMSISE00_INCLUDE_BASE})
        include_directories(SYSTEM "${NRLMSISE00_INCLUDE_DIRS}")
        message(STATUS "nrlmsise00 downloaded to ${nrlmsise00_fetch_SOURCE_DIR}")
        message(STATUS "nrlmsise00 include dir: ${NRLMSISE00_INCLUDE_DIRS}")
    else()
        message(STATUS ${NRLMSISE00_LIBRARIES})
    endif()
    add_definitions(-DTUDAT_BUILD_WITH_NRLMSISE=1)
else()
    find_package(nrlmsise00 REQUIRED 0.1)
    message(STATUS ${NRLMSISE00_LIBRARIES})
    add_definitions(-DTUDAT_BUILD_WITH_NRLMSISE=1)
endif()

# CALCEPH dependency (for direct binary SPK reading in WASM)
# CALCEPH is a pure C library that can read SPK files without CSPICE's furnsh_c()
# which crashes in WASM due to f2c FORTRAN I/O issues.
if(TUDAT_BUILD_WASM)
    message(STATUS "Downloading CALCEPH for WASM binary SPK reading...")
    include(FetchContent)
    FetchContent_Declare(
        calceph
        URL https://www.imcce.fr/content/medias/recherche/equipes/asd/calceph/calceph-4.0.5.tar.gz
        URL_HASH SHA256=3460d8a3e10a86e7fe0228d5d9abcda589713b8ed3ee007ce061ae01f8c2e1ea
    )
    FetchContent_GetProperties(calceph)
    if(NOT calceph_POPULATED)
        FetchContent_Populate(calceph)
    endif()

    # Collect all CALCEPH C source files (exclude interfaces_C_F90.c which requires Fortran)
    file(GLOB CALCEPH_SOURCES "${calceph_SOURCE_DIR}/src/*.c")
    # Remove the Fortran interface file if present
    list(FILTER CALCEPH_SOURCES EXCLUDE REGEX "interfaces_C_F90\\.c$")

    # Generate calcephconfig.h with WASM-compatible settings
    set(CALCEPH_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}/calceph_config")
    file(MAKE_DIRECTORY "${CALCEPH_CONFIG_DIR}")
    file(WRITE "${CALCEPH_CONFIG_DIR}/calcephconfig.h" [[
/* CALCEPH configuration for WASM builds */
#ifndef __CALCEPH_CONFIG_H__
#define __CALCEPH_CONFIG_H__

/* ANSI C standard headers (for stdarg.h, stdlib.h, etc.) */
#define STDC_HEADERS 1

/* Standard headers available in WASM/Emscripten */
#define HAVE_STDIO_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STDARG_H 1
#define HAVE_LIMITS_H 1
#define HAVE_STRING_H 1
#define HAVE_ERRNO_H 1
#define HAVE_MATH_H 1
#define HAVE_STDBOOL_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_UNISTD_H 1
#define HAVE_TIME_H 1

/* WASM doesn't have mmap */
/* #undef HAVE_SYS_MMAN_H */
/* #undef HAVE_MMAP */

/* Available functions in Emscripten */
#define HAVE_RINT 1
#define HAVE_SNPRINTF 1
#define HAVE_VSNPRINTF 1
#define HAVE_STRNCMP 1
#define HAVE_STRCASECMP 1
#define HAVE_FSEEKO 1

/* Not available in WASM */
/* #undef HAVE_VASPRINTF */
/* #undef HAVE_STRTOD_L */
/* #undef HAVE_NEWLOCALE */
/* #undef HAVE_FREELOCALE */
/* #undef HAVE_SETLOCALE */
/* #undef HAVE_XLOCALE_H */
/* #undef HAVE_LOCALE_H */

/* C99 bool type */
#define HAVE_BOOL 1

/* GCC attribute for unused variables */
#define HAVE_VAR_ATTRIBUTE_UNUSED 1

/* Disable thread support for WASM single-threaded mode */
/* #undef HAVE_PTHREAD */

#endif /* __CALCEPH_CONFIG_H__ */
]])

    # Create static library
    add_library(calceph STATIC ${CALCEPH_SOURCES})
    target_include_directories(calceph PUBLIC
        "${calceph_SOURCE_DIR}/src"
        "${CALCEPH_CONFIG_DIR}"
    )
    target_compile_definitions(calceph PRIVATE
        HAVE_CONFIG_H=1
    )
    # Suppress warnings in third-party code
    target_compile_options(calceph PRIVATE -w)

    set(CALCEPH_LIBRARIES calceph)
    set(CALCEPH_INCLUDE_DIRS "${calceph_SOURCE_DIR}/src" "${CALCEPH_CONFIG_DIR}")
    include_directories(SYSTEM "${CALCEPH_INCLUDE_DIRS}")
    add_definitions(-DTUDAT_BUILD_WITH_CALCEPH=1)
    set(TUDAT_BUILD_WITH_CALCEPH TRUE)  # CMake variable for conditional source inclusion
    message(STATUS "CALCEPH downloaded to ${calceph_SOURCE_DIR}")
    message(STATUS "CALCEPH include dirs: ${CALCEPH_INCLUDE_DIRS}")
else()
    # For non-WASM builds, CALCEPH is optional
    set(CALCEPH_LIBRARIES "")
    set(CALCEPH_INCLUDE_DIRS "")
endif()

# FFTW3 dependency
if (TUDAT_BUILD_WITH_FFTW3)
    if(NOT APPLE)
        find_package(FFTW3)
        message(STATUS "FFTW " ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        if(WIN32)
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/fftw3.lib")
        else()
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.so")
        endif()
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    else( )
        find_package(FFTW3)
        message(STATUS ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        set(
            CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -isystem \"${FFTW3_INCLUDE_DIR}\""
        )
        set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.dylib")
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    endif( )
else( )
    add_definitions(-DTUDAT_BUILD_WITH_FFTW3=0)
endif( )

# JSON dependency
if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# TODO: Not sure if this is needed. Comes from CMake of tudatpy
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)

# Include YOLO functionality
# TODO: Check what this does and if we need it
include(YOLOProjectAddTestCase)
include(YOLOProjectAddLibrary)
include(YOLOProjectAddExecutable)
include(YOLOProjectAddExternalData)

# Define list of libraries in tudat and run setup for submodules
include(TudatLinkLibraries)
add_subdirectory(src/tudat/astro)
add_subdirectory(src/tudat/basics)
add_subdirectory(src/tudat/interface)
add_subdirectory(src/tudat/math)
add_subdirectory(src/tudat/simulation)
add_subdirectory(src/tudat/io)
add_subdirectory(src/tudat/utils)
if (TUDAT_BUILD_TESTS)
    add_subdirectory(tests/test_tudat)
endif ()

# WASM-specific tests (don't require Boost.Test)
if (TUDAT_BUILD_WASM)
    add_subdirectory(tests/wasm)
endif ()

if (TUDAT_BUILD_TUDAT_TUTORIALS)
    add_subdirectory(examples/tudat)
endif()

# Additional contents from original tudat cmake
# ========================================================================

file(
    RELATIVE_PATH TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX
    "${CMAKE_INSTALL_PREFIX}"
    "${INSTALL_DATA_DIR}/tudat"
)

# Set up config.hpp file of tudat for current version
configure_file(
        "${BASE_DIRECTORY}/configuration/tudat/config.hpp.in"
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp" @ONLY
)

# Installation of the header files.
install(
    DIRECTORY include/
    DESTINATION "${INSTALL_INCLUDE_DIR}"
)
install(
    FILES "${PROJECT_BINARY_DIR}/include/tudat/config.hpp"
    DESTINATION ${INSTALL_INCLUDE_DIR}/tudat
)

# Add all targets to the build-tree export set
# For WASM builds, skip export as dependencies are built locally via FetchContent
if(NOT TUDAT_BUILD_WASM)
    export(
        TARGETS ${Tudat_LIBRARIES}
        NAMESPACE Tudat::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/tudat_export.cmake"
    )

    # Installation of the export modules.
    install(EXPORT tudat_export NAMESPACE Tudat:: DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}")
endif()

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE Tudat)

# +============================================================================
# PROJECT CONFIG
#  For external use of project.
# +============================================================================
include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${BASE_DIRECTORY}/configuration/tudat/tudat-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake"
        INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}"
)

configure_file(
        "${BASE_DIRECTORY}/configuration/tudat/tudat-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake" @ONLY
)

# Configure file path for tudat data loading.
#set(Tudat_DATA_DIRS ${INSTALL_CMAKE_DIR}/${RELATIVE_INSTALL_DATA_DIR})


file(
        COPY "${BASE_DIRECTORY}/cmake/tudat/TudatFindBoost.cmake"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake"
        "${BASE_DIRECTORY}/cmake/tudat/TudatFindBoost.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}"
)

# Take care of versioning.

# NOTE: SameMinorVersion available only
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}"
)

# ========
message(
    STATUS
    "TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX: ${TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX}"
)

# +============================================================================
# Tudatpy specific configuration
# +============================================================================

message(STATUS Tudat: [Tudat_PROPAGATION_LIBRARIES] ${Tudat_PROPAGATION_LIBRARIES})
message(STATUS Tudat: [Tudat_INCLUDE_DIRS] ${Tudat_INCLUDE_DIRS})

# Pybind11 dependency - skip for WASM builds (no Python bindings in browser)
if(NOT TUDAT_BUILD_WASM)
    include(YACMAPythonSetup)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)
    set(pybind11_INCLUDE_DIRS ${pybind11_BUILD_INCLUDE_DIRS}) # TODO: Check
    set(pybind11_INCLUDE_DIR ${pybind11_BUILD_INCLUDE_DIRS})  # TODO: Check

    # Source directory for tudatpy
    set(TUDATPY_SOURCE_DIR ${BASE_DIRECTORY}/src/tudatpy)

    # Add docs to build directory
    file(
        COPY ${BASE_DIRECTORY}/docs/tudatpy
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/docs
    )

    # Include tudatpy CMake configuration
    include_directories(include/tudatpy)
    add_subdirectory(src/tudatpy)
endif()

# WASM bindings subdirectory (only for Emscripten builds)
if(TUDAT_BUILD_WASM)
    message(STATUS "Adding WASM bindings subdirectory")
    add_subdirectory(src/tudatpy_wasm)
endif()

# +============================================================================
# CLEAN UP (Project name independent)
#  Cleanup in case project is not top level.
# +============================================================================
unset(INSTALL_LIB_DIR)
unset(INSTALL_BIN_DIR)
unset(INSTALL_INCLUDE_DIR)
unset(INSTALL_DATA_DIR)
unset(INSTALL_MAN_DIR)
unset(INSTALL_TESTS_DIR)
unset(INSTALL_CMAKE_DIR)

unset(RELATIVE_INSTALL_LIB_DIR)
unset(RELATIVE_INSTALL_BIN_DIR)
unset(RELATIVE_INSTALL_INCLUDE_DIR)
unset(RELATIVE_INSTALL_DATA_DIR)
unset(RELATIVE_INSTALL_MAN_DIR)
unset(RELATIVE_INSTALL_TESTS_DIR)
unset(RELATIVE_INSTALL_CMAKE_DIR)
