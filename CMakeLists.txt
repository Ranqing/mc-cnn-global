cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)

set (PROJECT_NAME mc-cnn-global)
project( ${PROJECT_NAME} )

option(WITH_OPENCV "" ON)
option(WITH_PCL "" OFF)

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(use_cxx11)

use_cxx11()

file(GLOB PROJECT_FILES
  "*.h"
  "*.cpp"
)

if(WITH_OPENCV)
	find_package(OpenCV REQUIRED)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
	link_directories( ${OpenCV_LIB_DIR} )
endif(WITH_OPENCV)

if(WITH_PCL)
	find_package(PCL 1.8 REQUIRED)
        include_directories(${PCL_INCLUDE_DIRS})
        link_directories(${PCL_LIBRARY_DIRS})
        add_definitions(${PCL_DEFINITIONS})
endif(WITH_PCL)

add_executable ( ${PROJECT_NAME} ${PROJECT_FILES})

if(WITH_OPENCV)
	target_link_libraries ( ${PROJECT_NAME} ${OpenCV_LIBS} )
endif(WITH_OPENCV)

if(WITH_PCL)
	target_link_libraries ( ${PROJECT_NAME} ${PCL_LIBRARIES} )
endif(WITH_PCL)

