 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #
 #    References
 #       compiler.cmake GPLv3
 #         Software Competence Center Hagenberg GmbH (SCCH)
 #         <thomas.natschlaeger@scch.at>, <office@scch.at>

# Provide options to force building with GNU or Clang, if the standard compiler is not desired
option(USE_CLANG "Force build with clang (if gcc is standard)" OFF) # OFF is the default
option(USE_GNU   "Force build with gcc (if clang is standard)" OFF) # OFF is the default
set( CLANG_C_COMPILER   "/usr/bin/clang" CACHE FILEPATH "Path to clang C compiler" )
set( CLANG_CXX_COMPILER "/usr/bin/clang++" CACHE FILEPATH "Path to clang C++ compiler" )
set( GNU_C_COMPILER     "/usr/bin/gcc" CACHE FILEPATH "Path to C compiler" )
set( GNU_CXX_COMPILER   "/usr/bin/g++" CACHE FILEPATH "Path to C++ compiler" )

# Set the platform type and override compiler if necessary
if( USE_CLANG OR USE_GNU )
    message(STATUS "Guessing compiler executables:")
    if( USE_GNU )
        message(STATUS "  GNU C     : ${GNU_C_COMPILER}")
        message(STATUS "  GNU C++   : ${GNU_CXX_COMPILER}")
        set( CMAKE_C_COMPILER        "${GNU_C_COMPILER}" )
        set( CMAKE_CXX_COMPILER      "${GNU_CXX_COMPILER}" )
        set( TUDAT_BUILD_GNU ON)
    elseif( USE_CLANG )
        message(STATUS "  Clang C   : ${CLANG_C_COMPILER}")
        message(STATUS "  Clang C++ : ${CLANG_CXX_COMPILER}")
        set( CMAKE_C_COMPILER        "${CLANG_C_COMPILER}" )
        set( CMAKE_CXX_COMPILER      "${CLANG_CXX_COMPILER}" )
        set( TUDAT_BUILD_CLANG ON)

        # Gets unset by CMake, because "since llvm doesn't have its own binutils
        # but uses the regular ar, objcopy, etc. (instead of llvm-objcopy etc.)".
        # Note: -D_CMAKE_TOOLCHAIN_PREFIX=llvm- will persist.
        set ( _CMAKE_TOOLCHAIN_PREFIX "llvm-" )
    endif()
    message(STATUS "  (if incorrect, please set these variables manually)")
else()
    if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$" )
        set( TUDAT_BUILD_CLANG ON)
    elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
        set( TUDAT_BUILD_GNU   ON)
    elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
        set( TUDAT_BUILD_MSVC  ON)
    endif()
endif()

# Set the compile flags
if( TUDAT_BUILD_CLANG )
    message(STATUS "Using clang compiler.")
    set ( CMAKE_C_FLAGS                "-Wall -std=c11" )
    set ( CMAKE_C_FLAGS_DEBUG          "-g" )
    set ( CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG" )
    set ( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
    set ( CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" )

    set ( CMAKE_CXX_FLAGS                "-Wall -Wextra -Wno-unused-parameter -std=c++11" )
    set ( CMAKE_CXX_FLAGS_DEBUG          "-g" )
    set ( CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" )
    set ( CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" )
    set ( CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" )

elseif( TUDAT_BUILD_GNU )
    message(STATUS "Using gnucxx compiler.")
    include( CheckCXXCompilerFlag )
    check_cxx_compiler_flag( "-std=c++11" CXX_SUPPORTS_CXX11 )
    if ( CXX_SUPPORTS_CXX11 )
        set ( CMAKE_CXX_FLAGS     "-Wall -std=c++11" )
    else()
        check_cxx_compiler_flag( "-std=c++0x" CXX_SUPPORTS_CXX0x )
        if ( CXX_SUPPORTS_CXX0x )
            set ( CMAKE_CXX_FLAGS "-Wall -std=c++0x" )
        endif()
    endif()

    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG          "-g")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unused-parameter -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")

    # MinGW fixes
    if( MINGW AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
      # MinGW fails to build with O2 or O3 optimization on several math.h function
      # http://ehc.ac/p/mingw/bugs/2250/
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__NO_INLINE__")
      # MinGW gives some c11plus.xe out of memory messages:
      # http://sourceforge.net/p/mingw-w64/mailman/message/33182613/
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0")
      # MinGW32 4.8.1 has no defenitions for _aligned_malloc/realloc/free
      #
      add_definitions(-DEIGEN_MALLOC_ALREADY_ALIGNED=1)
      add_definitions(-DEIGEN_DONT_ALIGN=1)
    endif()

elseif( TUDAT_BUILD_MSVC )
    message(STATUS "Using msvc compiler.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach(flag_var
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        # Find all dynamic runtime (MD) references and replace with static (MT)
        if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
    if( MSVC_VERSION GREATER  1500 )
        # Multiprocessor support during compilation
        add_definitions( "/MP" )
    endif()
else()
    message(STATUS "Compiler not identified: ${CMAKE_CXX_COMPILER_ID}" )
    message(STATUS "  Path: ${CMAKE_CXX_COMPILER}")
endif()

message(STATUS "Building with flags: ${CMAKE_CXX_FLAGS}.")
