cmake_minimum_required(VERSION 3.18)
project(DiagPointTests LANGUAGES CXX)
enable_testing()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(PHYSICS_ROOT ${CMAKE_CURRENT_LIST_DIR}/../../.. ABSOLUTE)
set(AMReX_FORTRAN OFF CACHE BOOL "" FORCE)
set(AMReX_PARTICLES OFF CACHE BOOL "" FORCE)
set(AMReX_LINEAR_SOLVERS OFF CACHE BOOL "" FORCE)
set(AMReX_AMRLEVEL OFF CACHE BOOL "" FORCE)
set(AMReX_EB OFF CACHE BOOL "" FORCE)
set(AMReX_PLOTFILE_TOOLS OFF CACHE BOOL "" FORCE)
add_subdirectory(${PHYSICS_ROOT}/Submodules/amrex amrex)
add_executable(test_diagpoint
  ${CMAKE_CURRENT_LIST_DIR}/test_diagpoint.cpp
  ${PHYSICS_ROOT}/Source/Utility/Diagnostics/DiagPoint.cpp
  ${PHYSICS_ROOT}/Source/Utility/Diagnostics/DiagBase.cpp
  ${PHYSICS_ROOT}/Source/Utility/Diagnostics/DiagFilter.cpp)
target_include_directories(test_diagpoint PRIVATE
  ${PHYSICS_ROOT}/Source
  ${PHYSICS_ROOT}/Source/Utility/Diagnostics)
target_link_libraries(test_diagpoint PRIVATE AMReX::amrex)

function(add_point_test scenario)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${scenario})
  if(AMReX_MPI)
    add_test(NAME ${scenario} COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_diagpoint> ${MPIEXEC_POSTFLAGS} scenario=${scenario}
      amrex.throw_exception=1 amrex.signal_handling=0)
  else()
    add_test(NAME ${scenario} COMMAND test_diagpoint scenario=${scenario}
      amrex.throw_exception=1 amrex.signal_handling=0)
  endif()
  set_tests_properties(${scenario} PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${scenario} TIMEOUT 45)
endfunction()

add_point_test(history)
add_point_test(no_volume_fraction)
foreach(scenario upper_boundary lower_boundary internal_boundary outside
    bad_dimension nonfinite bad_interval period filters empty_fields
    missing_field covered changed_fields changed_position older_restart
    inconsistent_step broken_history unwritable)
  add_point_test(${scenario})
endforeach()
