# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests (GTest-based)
file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(REMOVE_ITEM all_tests
  atom/test_mcmc_atom.cpp
  atom/test_mcmc_correlator_atom.cpp
  atom/test_correlator_atom.cpp
  dimer/test_mcmc_dimer.cpp
  dimer/test_mcmc_dimer_correlator.cpp
)

foreach(test ${all_tests})
  get_filename_component(test_stem ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  if(test_dir STREQUAL "")
    set(test_name ${test_stem})
  else()
    string(REPLACE "/" "_" test_prefix ${test_dir})
    set(test_name ${test_prefix}_${test_stem})
  endif()
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest gtest_main)
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TARGET ${test_name} PROPERTY OUTPUT_NAME ${test_stem})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++20
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
  endif()
endforeach()

# MPI + GTest tests (custom main for MPI init)
add_executable(test_mcmc_atom atom/test_mcmc_atom.cpp)
target_link_libraries(test_mcmc_atom ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest)
set_property(TARGET test_mcmc_atom PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/atom)
add_test(NAME test_mcmc_atom COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/atom/test_mcmc_atom)

add_executable(test_mcmc_dimer dimer/test_mcmc_dimer.cpp)
target_link_libraries(test_mcmc_dimer ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest)
set_property(TARGET test_mcmc_dimer PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dimer)
add_test(NAME test_mcmc_dimer COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/dimer/test_mcmc_dimer)

add_executable(test_mcmc_dimer_correlator dimer/test_mcmc_dimer_correlator.cpp)
target_link_libraries(test_mcmc_dimer_correlator ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest)
set_property(TARGET test_mcmc_dimer_correlator PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dimer)
add_test(NAME test_mcmc_dimer_correlator COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/dimer/test_mcmc_dimer_correlator)

add_executable(test_correlator_atom atom/test_correlator_atom.cpp)
target_link_libraries(test_correlator_atom ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest)
set_property(TARGET test_correlator_atom PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/atom)
add_test(NAME test_correlator_atom COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/atom/test_correlator_atom)

add_executable(test_mcmc_correlator_atom atom/test_mcmc_correlator_atom.cpp)
target_link_libraries(test_mcmc_correlator_atom ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest)
set_property(TARGET test_mcmc_correlator_atom PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/atom)
add_test(NAME test_mcmc_correlator_atom COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/atom/test_mcmc_correlator_atom)

