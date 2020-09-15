 #    Copyright (c) 2010-2019, Delft University of Technology
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

 add_compile_definitions(CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
 # Provide options to force building with GNU or Clang, if the standard compiler is not desired
 option(USE_CLANG "Force build with clang (if gcc is standard)" OFF) # OFF is the default
 option(USE_GNU "Force build with gcc (if clang is standard)" OFF) # OFF is the default
 set(CLANG_C_COMPILER "/usr/bin/clang" CACHE FILEPATH "Path to clang C compiler")
 set(CLANG_CXX_COMPILER "/usr/bin/clang++" CACHE FILEPATH "Path to clang C++ compiler")
 set(GNU_C_COMPILER "/usr/bin/gcc" CACHE FILEPATH "Path to C compiler")
 set(GNU_CXX_COMPILER "/usr/bin/g++" CACHE FILEPATH "Path to C++ compiler")

 # Set the platform type and override compiler if necessary
 if (USE_CLANG OR USE_GNU)
     message(STATUS "Guessing compiler executables:")
     if (USE_GNU)
         message(STATUS "  GNU C     : ${GNU_C_COMPILER}")
         message(STATUS "  GNU C++   : ${GNU_CXX_COMPILER}")
         set(CMAKE_C_COMPILER "${GNU_C_COMPILER}")
         set(CMAKE_CXX_COMPILER "${GNU_CXX_COMPILER}")
         set(TUDAT_BUILD_GNU ON)
     elseif (USE_CLANG)
         message(STATUS "  Clang C   : ${CLANG_C_COMPILER}")
         message(STATUS "  Clang C++ : ${CLANG_CXX_COMPILER}")
         set(CMAKE_C_COMPILER "${CLANG_C_COMPILER}")
         set(CMAKE_CXX_COMPILER "${CLANG_CXX_COMPILER}")
         set(TUDAT_BUILD_CLANG ON)

         # Gets unset by cmake, because "since llvm doesn't have its own binutils
         # but uses the regular ar, objcopy, etc. (instead of llvm-objcopy etc.)".
         # Note: -D_CMAKE_TOOLCHAIN_PREFIX=llvm- will persist.
         set(_CMAKE_TOOLCHAIN_PREFIX "llvm-")
     endif ()
     message(STATUS "  (if incorrect, please set these variables manually)")
 else ()
     if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
         set(TUDAT_BUILD_CLANG ON)
         set(TUDAT_BUILD_GNU OFF)
         set(TUDAT_BUILD_MSVC OFF)
     elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
         set(TUDAT_BUILD_CLANG OFF)
         set(TUDAT_BUILD_GNU ON)
         set(TUDAT_BUILD_MSVC OFF)
     elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
         set(TUDAT_BUILD_CLANG OFF)
         set(TUDAT_BUILD_GNU OFF)
         set(TUDAT_BUILD_MSVC ON)
     endif ()
 endif ()

 # Set the compile flags
 if (TUDAT_BUILD_CLANG)
     add_compile_definitions(TUDAT_BUILD_CLANG)
     if(WIN32)
      add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
     endif()
     message(STATUS "Using clang compiler.")
     set(CMAKE_C_FLAGS "-Wall -std=c11")
     set(CMAKE_C_FLAGS_DEBUG "-g")
     set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")
     set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
     set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
     if (APPLE)
         set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -std=c++11 -stdlib=libc++")
     elseif (NOT WIN32)
         set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -std=c++11 -stdlib=libstdc++")
     else()
         set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -std=c++11")
     endif ()
     set(CMAKE_CXX_FLAGS_DEBUG "-g")
     set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
     set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

 elseif (TUDAT_BUILD_GNU)
     add_compile_definitions(TUDAT_BUILD_GNU)
     message(STATUS "Using gnucxx compiler.")
     include(CheckCXXCompilerFlag)
     check_cxx_compiler_flag("-std=c++11" CXX_SUPPORTS_CXX11)
     if (CXX_SUPPORTS_CXX11)
         set(CMAKE_CXX_FLAGS "-Wall -std=c++11")
     else ()
         check_cxx_compiler_flag("-std=c++0x" CXX_SUPPORTS_CXX0x)
         if (CXX_SUPPORTS_CXX0x)
             set(CMAKE_CXX_FLAGS "-Wall -std=c++0x")
         endif ()
     endif ()

     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
     set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
     set(CMAKE_CXX_FLAGS_DEBUG "-Og")

     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unused-parameter -Wno-unused-variable -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")

     # MinGW fixes
     if (MINGW AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
         # MinGW fails to build with O2 or O3 optimization on several math.h function
         # http://ehc.ac/p/mingw/bugs/2250/
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__NO_INLINE__")

         # MinGW32 4.8.1 has no defenitions for _aligned_malloc/realloc/free
         #
         add_definitions(-DEIGEN_MALLOC_ALREADY_ALIGNED=1)
         add_definitions(-DEIGEN_DONT_ALIGN=1)
     endif ()

     # Fix exceptions not being caught
     if (MINGW AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.4))
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
     endif ()

     if (MINGW AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
         # MinGW gives some c11plus.xe out of memory messages:
         # http://sourceforge.net/p/mingw-w64/mailman/message/33182613/
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0")
         #         set(CMAKE_EXE_LINKER_FLAGS "-Wl,--large-address-aware")
     endif ()

 elseif (TUDAT_BUILD_MSVC)
     add_compile_definitions(TUDAT_BUILD_MSVC)
     add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
     message(STATUS "Using MSVC compiler.")
     # problem: https://dev.azure.com/tudat-team/feedstock-builds/_build/results?buildId=95&view=logs&j=00f5923e-fdef-5026-5091-0d5a0b3d5a2c&t=3cc4a9ed-60e1-5810-6eb3-5f9cd4a26dba
     # solution: https://stackoverflow.com/questions/1091662/vc-internal-compiler-error
     #     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
     if (TUDAT_FORCE_DYNAMIC_RUNTIME)
         # This is needed for conda builds, as the prebuilt libraries are MD.
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")
         message(STATUS "Forcing MD_DynamicRuntime on MSVC")
     else ()
         # The following applied to legacy Tudat.
         # Because we are using static boost libraries, with static runtime, we need to force MSVC to
         # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
         foreach (flag_var
                 CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
                 CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
             # Find all dynamic runtime (MD) references and replace with static (MT)
             if (${flag_var} MATCHES "/MD")
                 string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
             endif (${flag_var} MATCHES "/MD")
         endforeach (flag_var)
     endif ()
     if (MSVC_VERSION GREATER 1500)
         # Multiprocessor support during compilation
         add_definitions("/MP")
     endif ()
 else ()
     message(STATUS "Compiler not identified: ${CMAKE_CXX_COMPILER_ID}")
     message(STATUS "  Path: ${CMAKE_CXX_COMPILER}")
 endif ()

 set(CMAKE_POSITION_INDEPENDENT_CODE ON)
 message(STATUS "Building with flags: ${CMAKE_CXX_FLAGS}.")


 if (MSVC)
     #    add_definitions(-Dinline=__inline)
     message(STATUS "Using [${CMAKE_C_COMPILER_ID}] compiler")
     if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
         set(MSVC_DISABLED_WARNINGS_LIST
                 "C4305" # 'initializing': truncation from 'int' to 'bool'
                 "C4101" # : unreferenced local variable
                 "C4018" # 'expression' : signed/unsigned mismatch
                 "C4057" # 'operator' : 'identifier1' indirection to
                 # slightly different base types from 'identifier2'
                 "C4068" # : unknown pragma directive
                 "C4100" # 'identifier' : unreferenced formal parameter
                 "C4127" # conditional expression is constant
                 "C4146" # unary minus operator applied to unsigned type,
                 # result still unsigned
                 "C4244" # 'argument' : conversion from 'type1' to 'type2',
                 # possible loss of data
                 "C4245" # 'conversion' : conversion from 'type1' to 'type2',
                 # signed/unsigned mismatch
                 "C4267" # 'var' : conversion from 'size_t' to 'type',
                 # possible loss of data
                 "C4389" # 'operator' : signed/unsigned mismatch
                 "C4706" # assignment within conditional expression
                 "C4996" # The POSIX name for this item is deprecated.
                 # Instead, use the ISO C and C++ conformant name
                 )
     elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
         add_definitions(-D_CRT_SUPPRESS_RESTRICT)
         set(MSVC_DISABLED_WARNINGS_LIST
                 "C111"  # Unreachable statement
                 "C128"  # Unreachable loop
                 "C167"  # Unexplict casting unsigned to signed
                 "C186"  # Pointless comparison of unsigned int with zero
                 "C188"  # Enumerated type mixed with another type
                 "C344"  # Redeclared type
                 "C556"  # Unexplict casting signed to unsigned
                 "C869"  # Unreferenced parameters
                 "C1786" # Deprecated functions
                 "C2545" # Empty else statement
                 "C2557" # Comparing signed to unsigned
                 "C2722" # List init syntax is c++11 feature
                 "C3280" # Declaration hides variable
                 )
     endif ()
     string(REPLACE "C" " -wd" MSVC_DISABLED_WARNINGS_STR
             ${MSVC_DISABLED_WARNINGS_LIST})
     string(REGEX REPLACE "[/-]W[1234][ ]?" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -MP -W4 ${MSVC_DISABLED_WARNINGS_STR}")

     message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
     add_definitions(${MSVC_DISABLED_WARNINGS_STR})
 endif ()
