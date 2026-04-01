# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests (GTest-based)
file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(REMOVE_ITEM all_tests
  test_diagram.cpp
  test_diagram_dimer.cpp
  test_diagram2.cpp
  dimer_order_4.cpp
  test_hubbard_atom.cpp
  test_atomic_cumulant.cpp
)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest gtest_main)
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
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

# Standalone tests (no GTest dependency)
add_executable(test_diagram2 test_diagram2.cpp)
target_link_libraries(test_diagram2 ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)
add_test(NAME test_diagram2 COMMAND test_diagram2)
