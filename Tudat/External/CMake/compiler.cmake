#
# ./cmake/compiler.cmake
# ------------------------------------------------------------------------------
# Copyright (C) 2012-2014 Software Competence Center Hagenberg GmbH (SCCH)
# <thomas.natschlaeger@scch.at>, <office@scch.at>
# -----------------------------------------------------------------------------
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
# This code is subject to dual-licensing. Please contact office@scch.at
# if you are interested in obtaining a differently licensed version.
#


option(USE_CLANG "build application with clang" OFF) # OFF is the default
SET (CLANG_C_COMPILER   "/usr/bin/clang" CACHE FILEPATH "Path to clang C compiler" )
SET (CLANG_CXX_COMPILER "/usr/bin/clang++" CACHE FILEPATH "Path to clang C++ compiler" )

if(USE_CLANG)
    message(STATUS "Using clang compiler.")
    set ( CMAKE_C_COMPILER             "${CLANG_C_COMPILER}" )
    set ( CMAKE_C_FLAGS                "-Wall -std=c11" )
    set ( CMAKE_C_FLAGS_DEBUG          "-g" )
    set ( CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG" )
    set ( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
    set ( CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" )

    set ( CMAKE_CXX_COMPILER             "${CLANG_CXX_COMPILER}" )
    set ( CMAKE_CXX_FLAGS                "-Wall -std=c++11" )
    set ( CMAKE_CXX_FLAGS_DEBUG          "-g" )
    set ( CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" )
    set ( CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" )
    set ( CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" )

    # Gets unset by CMake, because "since llvm doesn't have its own binutils
    # but uses the regular ar, objcopy, etc. (instead of llvm-objcopy etc.)".
    # Note: -D_CMAKE_TOOLCHAIN_PREFIX=llvm- will persist.
    set ( _CMAKE_TOOLCHAIN_PREFIX "llvm-" )

elseif( CMAKE_COMPILER_IS_GNUCXX )
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

    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g")
    set(CMAKE_CXX_FLAGS_RELEASE        "-g -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG          "-g")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
	
	# MinGW fixes
	if( MINGW )
	  # MinGW fails to build with O2 or O3 optimization on several math.h function
	  # http://ehc.ac/p/mingw/bugs/2250/
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__NO_INLINE__ -ftrack-macro-expansion=0")
	  # MinGW gives some c11plus.xe out of memory messages:
	  # http://sourceforge.net/p/mingw-w64/mailman/message/33182613/
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0")
	  # MinGW32 4.8.1 has no defenitions for _aligned_malloc/realloc/free
	  # 
	  add_definitions(-DEIGEN_MALLOC_ALREADY_ALIGNED=1) 
	  add_definitions(-DEIGEN_DONT_ALIGN=1)
	endif()

elseif( MSVC )
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
endif()
