# MC executables (not tests — these are standalone applications)

add_executable(mcmc_atom mcmc_atom.cpp)
target_link_libraries(mcmc_atom ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)

add_executable(mcmc_dimer mcmc_dimer.cpp)
target_link_libraries(mcmc_dimer ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)
