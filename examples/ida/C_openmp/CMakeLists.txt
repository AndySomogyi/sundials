# ---------------------------------------------------------------
# $Revision:  $
# $Date:  $
# ---------------------------------------------------------------
# Programmer: Ting Yan @ SMU
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the IDA OpenMP examples


# Only include tests if OpenMP is enabled
IF(OPENMP_FOUND)
  # Add variable ida_examples with the names of the openmp IDA examples
  SET(IDA_examples
    idaFoodWeb_bnd_omp
    idaFoodWeb_kry_omp
    )

  SET(IDA_examples_BL
  )
ENDIF(OPENMP_FOUND)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_static)
  SET(NVECS_LIB sundials_nvecopenmp_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(IDA_LIB sundials_ida_shared)
  SET(NVECS_LIB sundials_nvecopenmp_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")


# Set-up linker flags and link libraries
SET(SUNDIALS_LIBS ${IDA_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

IF(KLU_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_FOUND)

IF(PETSC_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${PETSC_LIBRARIES})
ENDIF(PETSC_FOUND)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

# Add the build and install targets for each IDA example

FOREACH(example ${IDA_examples})
  ADD_EXECUTABLE(${example} ${example}.c)
  SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
  SUNDIALS_ADD_TEST(${example} ${example} TEST_ARGS 4)
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
  IF(EXAMPLES_INSTALL)
    INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/C_openmp)
  ENDIF(EXAMPLES_INSTALL)
ENDFOREACH(example ${IDA_examples})

# If Lapack support is enabled, add the build and install targets for
# the examples using Lapack

IF(LAPACK_FOUND)
  FOREACH(example ${IDA_examples_BL})
    ADD_EXECUTABLE(${example} ${example}.c)
    SET_TARGET_PROPERTIES(${example} PROPERTIES FOLDER "Examples")
    SUNDIALS_ADD_TEST(${example} ${example} TEST_ARGS 4)
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    IF(EXAMPLES_INSTALL)
      INSTALL(FILES ${example}.c ${example}.out DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/C_openmp)
    ENDIF(EXAMPLES_INSTALL)
  ENDFOREACH(example ${IDA_examples_BL})
ENDIF(LAPACK_FOUND)

IF(EXAMPLES_INSTALL)

  # Install the README file
  INSTALL(FILES README DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/C_openmp)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  SET(SOLVER "IDA")
  SET(SOLVER_LIB "sundials_ida")
  LIST2STRING(IDA_examples EXAMPLES)
  IF(LAPACK_FOUND)
    LIST2STRING(IDA_examples_BL EXAMPLES_BL)
  ELSE(LAPACK_FOUND)
    SET(EXAMPLES_BL "")
  ENDIF(LAPACK_FOUND)
  
  STRING (REPLACE ";" " " TMP_STR ${EXAMPLES})
  SET(EXAMPLES ${TMP_STR})
  
  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmp_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/C_openmp/CMakeLists.txt
      @ONLY
      )

  # install CMakelists.txt
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/examples/ida/C_openmp/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/C_openmp 
    )

  IF(UNIX)
    # generate Makefile and place it in the binary dir
    CONFIGURE_FILE(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmp_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/C_openmp/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    INSTALL(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/C_openmp/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/C_openmp
      RENAME Makefile
      )
  ENDIF(UNIX)

ENDIF(EXAMPLES_INSTALL)
