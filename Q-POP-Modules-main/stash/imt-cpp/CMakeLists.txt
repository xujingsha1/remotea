find_package(DOLFIN REQUIRED)
find_package(Boost REQUIRED COMPONENTS timer)

include(${DOLFIN_USE_FILE})

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Do not throw error for 'multi-line comments' (these are typical in
# rst which includes LaTeX)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wno-comment" HAVE_NO_MULTLINE)
if (HAVE_NO_MULTLINE)
  set(CMAKE_CXX_FLAGS "-Wno-comment ${CMAKE_CXX_FLAGS}")
endif()

add_executable(qpop-imt)
target_sources(
  qpop-imt PRIVATE main.cpp ${PROJECT_SOURCE_DIR}/external/pugixml/pugixml.cpp
                   imt.cpp)
target_include_directories(
  qpop-imt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                  ${PROJECT_SOURCE_DIR}/external/pugixml ${DOLFIN_INCLUDE_DIRS})
target_link_libraries(qpop-imt ${DOLFIN_LIBRARIES} Boost::timer Boost::boost)

add_custom_target(
  compile_ufl ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/imt.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/imt.cpp
  COMMENT "Compile imt.ufl file")

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/imt.h ${CMAKE_CURRENT_SOURCE_DIR}/imt.cpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/imt.ufl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND ffc -l dolfin -O -f split -f quadrature_degree=6 imt.ufl
  COMMENT "Create imt.h and imt.cpp by: ffc -l dolfin -O -f split -f quadrature_degree=6 ${CMAKE_CURRENT_SOURCE_DIR}/imt.ufl")

add_dependencies(qpop-imt compile_ufl)
