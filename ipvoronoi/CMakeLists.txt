# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project( ipvoronoi_Demo )

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL QUIET)


if ( CGAL_FOUND )
  include(${CGAL_USE_FILE})

  # ipvoronoi target
  add_executable  (ipvoronoi ipvoronoi.cpp)
  target_link_libraries(ipvoronoi ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS ipvoronoi )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(ipvoronoi)

else()

  message(STATUS "NOTICE: This requires CGAL and will not be compiled.")

endif()
