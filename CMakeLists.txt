cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Standards C++11 requis.
set_property(GLOBAL PROPERTY CXX_STANDARD 11)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED ON)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
endif()

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo")

project (Formal_Maths)

# The version number.
set (FORMAL_MATHS_VERSION_MAJOR 1)
set (FORMAL_MATHS_VERSION_MINOR 0)
 
configure_file (
  "${PROJECT_SOURCE_DIR}/config/Formal_Maths_config.h.in"
  "${PROJECT_BINARY_DIR}/config/Formal_Maths_config.h"
  )

# Activer dossiers dans Visual
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Niveau d'avertissements 4 sous Windows.
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long -pedantic")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")

file(
	GLOB_RECURSE
	source_files
	${PROJECT_BINARY_DIR}/config/*
	inc/*
	src/*
)

source_group(
	"Headers"
	FILES
	inc/Formal_Maths.h
	${PROJECT_BINARY_DIR}/config/Formal_Maths_config.h
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
include_directories(${PROJECT_BINARY_DIR}/config)

add_library(Formal_Maths STATIC ${source_files})
set_target_properties(Formal_Maths PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(tests)

# add the install targets
install (TARGETS Formal_Maths DESTINATION bin)
install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/inc/Formal_Maths.h
         DESTINATION include)
install (FILES ${PROJECT_BINARY_DIR}/config/Formal_Maths_config.h
         DESTINATION include)
 
 set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing that CMake has configured at least once.")